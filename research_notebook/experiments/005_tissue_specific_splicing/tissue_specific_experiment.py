#!/usr/bin/env python3
"""
Experiment 005 (v2): Tissue-specific vs averaged SPLICE_JUNCTIONS for rare disease genes.

Hypothesis: For brain-expressed rare disease genes (SCN1A, SCN2A, MECP2),
brain-tissue-specific AlphaGenome splicing scores will have different
performance characteristics than the tissue-averaged scores we used in
the original benchmark.

Key fix: SPLICE_JUNCTIONS returns AnnData with multiple rows per variant
(one per affected junction). We must aggregate per-variant by taking
the max-abs across (a) junctions × (b) tracks.
"""
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon, ttest_1samp
from sklearn.metrics import roc_auc_score, average_precision_score

from alphagenome.data import genome, ontology
from alphagenome.atlas import atlas

OUTPUT_DIR = Path("research_notebook/experiments/005_tissue_specific_splicing")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

key = open('/Users/hermes/projects/alphagenome-work/.alphagenome_key').read().strip()
os.environ['ALPHAGENOME_API_KEY'] = key
atlas_client = atlas.create(os.environ['ALPHAGENOME_API_KEY'])

# Identify brain-related tracks in SPLICE_JUNCTIONS
sm = atlas_client.scorer_metadata()
sj_meta = sm['SPLICE_JUNCTIONS'].track_metadata
brain_terms = ['brain', 'neuron', 'cortex', 'hippocamp', 'cerebell',
               'neural', 'glia', 'spinal', 'forebrain', 'midbrain',
               'dorsolateral prefrontal', 'anterior cingulate']
brain_mask = sj_meta['biosample_name'].str.contains(
    '|'.join(brain_terms), case=False, na=False)
brain_tracks = sj_meta[brain_mask]
print(f"SPLICE_JUNCTIONS: {len(sj_meta)} total tracks, {len(brain_tracks)} brain-related")

# Get unique ontology terms
brain_ats = []
for _, row in brain_tracks.iterrows():
    curie = str(row['ontology_curie'])
    try:
        brain_ats.append(ontology.from_curie(curie))
    except Exception:
        pass
unique_brain_ats = list(set(brain_ats))
print(f"Using {len(unique_brain_ats)} unique brain ontology terms")


def load_variants(gene):
    if gene == 'scn1a':
        df = pd.read_csv("outputs/benchmark_scn1a_live_api_raw.csv")
    else:
        df = pd.read_csv(f"outputs/cross_disease_{gene}_raw.csv")
    df = df[(df['ref'].str.len() == 1) & (df['alt'].str.len() == 1)].copy()
    return df


def make_variant(row):
    chrom = str(row['chrom']).strip()
    if not chrom.startswith('chr'):
        chrom = 'chr' + chrom
    return genome.Variant(
        chromosome=chrom,
        position=int(row['pos']),
        reference_bases=str(row['ref']).strip(),
        alternate_bases=str(row['alt']).strip(),
    )


def score_variants(df, brain_only=False):
    """Score variants. Returns per-variant max-abs score."""
    variants = [make_variant(r) for _, r in df.iterrows()]
    
    kwargs = dict(
        variants=variants,
        requested_scorers=['SPLICE_JUNCTIONS'],
        progress_bar=False,
        max_workers=4,
    )
    if brain_only:
        kwargs['ontology_terms'] = unique_brain_ats
    
    cond = 'brain' if brain_only else 'avg'
    print(f"  ({cond}) scoring {len(variants)} variants...")
    
    t0 = time.time()
    result = atlas_client.query_variants(**kwargs)
    elapsed = time.time() - t0
    
    # result['SPLICE_JUNCTIONS'] is AnnData with shape (n_junctions_total, n_tracks)
    # obs['variant'] is a Variant object per row indicating which input variant
    # Group by variant and take max-abs score per variant across (junctions × tracks)
    
    ad = result['SPLICE_JUNCTIONS']
    obs = ad.obs
    X = ad.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    # Build per-variant max-abs
    variant_scores = {}
    for i in range(len(obs)):
        var = obs.iloc[i]['variant']
        # Convert Variant to canonical string for dict key
        var_key = str(var) if not hasattr(var, 'chromosome') else (
            f"{var.chromosome}:{var.position}:{var.reference_bases}>{var.alternate_bases}"
        )
        if var_key not in variant_scores:
            variant_scores[var_key] = -np.inf
        row_max = float(np.abs(X[i]).max())
        if row_max > variant_scores[var_key]:
            variant_scores[var_key] = row_max
    
    # Map back to df
    scores = []
    for _, row in df.iterrows():
        chrom = str(row['chrom']).strip()
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        key_str = f"{chrom}:{int(row['pos'])}:{row['ref']}>{row['alt']}"
        if key_str in variant_scores:
            scores.append(variant_scores[key_str])
        else:
            scores.append(np.nan)
    
    return scores, elapsed


def auroc_auprc(scores, labels):
    auc = roc_auc_score(labels, scores)
    ap = average_precision_score(labels, scores)
    return auc, ap


def main():
    genes = ['scn1a', 'scn2a', 'mecp2', 'cftr', 'dmd']
    
    all_results = []
    
    for brain_only in [False, True]:
        cond = 'brain' if brain_only else 'avg'
        for gene in genes:
            try:
                df = load_variants(gene)
                scores, elapsed = score_variants(df, brain_only=brain_only)
                df[f'SJ_{cond}_score'] = scores
                
                labels = (df['clnsig_category'] == 'pathogenic').astype(int).values
                valid = ~np.isnan(scores)
                if valid.sum() > 0 and labels[valid].sum() > 0 and (1 - labels[valid]).sum() > 0:
                    auc, ap = auroc_auprc(
                        np.array(scores)[valid], labels[valid])
                    row = {
                        'gene': gene,
                        'condition': cond,
                        'n_total': int(valid.sum()),
                        'n_path': int(labels[valid].sum()),
                        'n_ben': int((1 - labels[valid]).sum()),
                        'auroc': auc,
                        'auprc': ap,
                        'runtime_s': elapsed,
                    }
                    all_results.append(row)
                    print(f"  {gene.upper()} ({cond}): AUROC={auc:.4f}, AUPRC={ap:.4f} "
                          f"(n_path={row['n_path']}, n_ben={row['n_ben']}, "
                          f"elapsed={elapsed:.1f}s)")
                
                # Save per-variant scores
                df[['chrom', 'pos', 'ref', 'alt', 'clnsig_category', f'SJ_{cond}_score']].to_csv(
                    OUTPUT_DIR / f'per_variant_{gene}_{cond}.csv', index=False)
            except Exception as e:
                import traceback
                print(f"Error on {gene} ({cond}): {e}")
                traceback.print_exc()
    
    # Save summary
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUTPUT_DIR / 'tissue_vs_averaged.csv', index=False)
    print(f"\n=== Final comparison ===")
    if len(results_df) == 0:
        print("No results to compare!")
        return
    
    pivot = results_df.pivot_table(
        index='gene', columns='condition', values='auprc')
    print("\nAUPRC pivot:")
    print(pivot.to_string())
    
    if 'brain' in pivot.columns and 'avg' in pivot.columns:
        pivot['delta'] = pivot['brain'] - pivot['avg']
        print("\nDelta (brain - avg):")
        print(pivot[['delta']].to_string())
        
        deltas = pivot['delta'].dropna()
        if len(deltas) >= 2:
            w_stat, w_p = wilcoxon(deltas)
            t_stat, t_p = ttest_1samp(deltas, 0)
            print(f"\nWilcoxon on deltas: W={w_stat}, p={w_p:.4f}")
            print(f"One-sample t-test: t={t_stat:.3f}, p={t_p:.4f}")
            print(f"Mean delta: {deltas.mean():.4f}")
            if deltas.mean() > 0:
                print("→ Brain-filtered is BETTER on average")
            elif deltas.mean() < 0:
                print("→ Brain-filtered is WORSE on average")
            else:
                print("→ No difference")


if __name__ == '__main__':
    main()