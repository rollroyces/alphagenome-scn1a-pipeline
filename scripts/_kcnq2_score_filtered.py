#!/usr/bin/env python3
"""
Exp 008 follow-up: re-score KCNQ2 with the SAME consequence filter as the
other 5 cross-disease genes (pathogenic + splicing-related vs benign + intronic).
This is the apples-to-apples comparison to answer "does KCNQ2 fit the pattern?"
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


SCORER_NAMES = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
NEG_CAP = 200  # cap benign+intronic to match other genes
POS_CAP = 46   # all pathogenic+splicing available
SEED = 42


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set")
        return 1

    # Re-extract from VCF with MC info
    import pysam
    tb = pysam.TabixFile("data/clinvar_grch38.vcf.gz",
                         index="data/clinvar_grch38.vcf.gz.tbi")
    rows = []
    for line in tb.fetch("20:63400679-63472909"):
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        info = {e.split("=", 1)[0]: e.split("=", 1)[1]
                for e in fields[7].split(";") if "=" in e}
        if "KCNQ2" not in info.get("GENEINFO", ""):
            continue
        if len(fields[3]) != 1:
            continue
        alts = fields[4].split(",") if fields[4] != "." else []
        snv_alt = next((a for a in alts if len(a) == 1), None)
        if snv_alt is None:
            continue
        clnsig = info.get("CLNSIG", "")
        s = clnsig.lower()
        if "conflicting" in s:
            cat = "conflicting"
        elif "pathogenic" in s and "benign" not in s:
            cat = "pathogenic"
        elif "benign" in s and "pathogenic" not in s:
            cat = "benign"
        elif "uncertain" in s or "vus" in s:
            cat = "uncertain"
        else:
            cat = "other"
        mc = info.get("MC", "").lower()
        rows.append({
            "chrom": 20, "pos": int(fields[1]),
            "ref": fields[3], "alt": snv_alt,
            "clnsig_category": cat, "mc": mc,
        })
    tb.close()

    df = pd.DataFrame(rows)
    pos = df[(df["clnsig_category"] == "pathogenic")
             & (df["mc"].str.contains("splice"))].copy()
    neg = df[(df["clnsig_category"] == "benign")
             & (df["mc"].str.contains("intron"))].copy()

    np.random.seed(SEED)
    if len(pos) > POS_CAP:
        pos = pos.sample(n=POS_CAP, random_state=SEED)
    if len(neg) > NEG_CAP:
        neg = neg.sample(n=NEG_CAP, random_state=SEED)

    pos["label"] = "positive"
    neg["label"] = "negative"
    bench = pd.concat([pos, neg], ignore_index=True)
    bench = bench.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    print(f"Apples-to-apples benchmark: {len(bench)} variants "
          f"({(bench['label']=='positive').sum()} pos, "
          f"{(bench['label']=='negative').sum()} neg)")

    dna_model = dna_client.create(api_key)
    out_rows = []
    t0 = time.time()
    n_succ = n_fail = 0
    for i, r in bench.iterrows():
        if i % 25 == 0 or i == len(bench) - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            remaining = (len(bench) - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{len(bench)}] {elapsed:.0f}s, "
                  f"success={n_succ}, fail={n_fail}, ETA {remaining:.0f}s",
                  flush=True)
        try:
            v = genome.Variant(
                chromosome="chr" + str(r["chrom"]),
                position=int(r["pos"]),
                reference_bases=str(r["ref"]),
                alternate_bases=str(r["alt"]),
            )
            interval = v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
            scorers = [
                variant_scorers.RECOMMENDED_VARIANT_SCORERS[n] for n in SCORER_NAMES
            ]
            scores = dna_model.score_variant(
                interval=interval, variant=v, variant_scorers=scorers,
            )
            row = {
                "chrom": r["chrom"], "pos": r["pos"],
                "ref": r["ref"], "alt": r["alt"],
                "label": r["label"],
                "score_success": True,
            }
            for ann, name in zip(scores, SCORER_NAMES):
                row[f"{name}_score"] = float(ann.X.sum())
            n_succ += 1
        except Exception as e:
            row = {
                "chrom": r["chrom"], "pos": r["pos"],
                "ref": r["ref"], "alt": r["alt"],
                "label": r["label"],
                "score_success": False,
                "score_error": f"{type(e).__name__}: {str(e)[:80]}",
            }
            n_fail += 1
        out_rows.append(row)

    elapsed = time.time() - t0
    print(f"\nDone. {n_succ}/{len(bench)} OK, {n_fail} fail, {elapsed:.0f}s")
    raw = pd.DataFrame(out_rows)
    raw.to_csv("outputs/cross_disease_kcnq2_filtered_raw.csv", index=False)

    # Metrics
    for name in SCORER_NAMES:
        ok = raw[raw["score_success"]].copy()
        y_true = (ok["label"] == "positive").astype(int).values
        y_score = ok[f"{name}_score"].values
        auroc = float(roc_auc_score(y_true, y_score))
        auprc = float(average_precision_score(y_true, y_score))
        # bootstrap CI
        rng = np.random.default_rng(SEED)
        aucs = []
        for _ in range(1000):
            idx = rng.integers(0, len(y_true), size=len(y_true))
            if len(np.unique(y_true[idx])) < 2:
                continue
            aucs.append(average_precision_score(y_true[idx], y_score[idx]))
        ci = np.percentile(aucs, [2.5, 97.5]) if aucs else (np.nan, np.nan)
        n_top = max(1, int(np.ceil(len(y_score) * 0.05)))
        top5 = float(y_true[np.argsort(y_score)[-n_top:]].mean())
        print(f"  {name}: AUPRC={auprc:.4f} [{ci[0]:.3f},{ci[1]:.3f}], "
              f"AUROC={auroc:.4f}, top-5%={top5:.3f}, n={len(ok)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
