# AlphaGenome SDK cheatsheet

Things I learned reading the source code that aren't obvious from the tutorials.

## Two clients

```python
from alphagenome.models import dna_client     # Live prediction API
from alphagenome.atlas import atlas          # Pre-computed scores (Atlas)
```

`dna_client.create(API_KEY)` → use when you have a custom sequence or a variant not in Atlas.
`atlas.create(API_KEY)` → use when scoring a variant in the human reference genome.

## Hard limits (all in `dna_client`)

```python
SEQUENCE_LENGTH_16KB    # 16_384
SEQUENCE_LENGTH_100KB   # 102_400
SEQUENCE_LENGTH_500KB   # 524_288
SEQUENCE_LENGTH_1MB     # 1_048_576
MAX_VARIANT_SCORERS_PER_REQUEST = 20
MAX_ISM_INTERVAL_WIDTH = 10   # base pairs per ISM chunk (internal batching)
```

The 4 supported sequence lengths are the only ones the model accepts. Any
input must be padded/centered to one of these.

## Output types

11 modalities. The `_ACTIVE` variants are tuned for variant scoring (filter
for active chromatin / expressed genes); they often work better than the
base versions when the question is "does this variant disrupt something?":

```
ATAC, ATAC_ACTIVE
CAGE, CAGE_ACTIVE
CHIP_HISTONE, CHIP_HISTONE_ACTIVE
CHIP_TF, CHIP_TF_ACTIVE
CONTACT_MAPS
DNASE, DNASE_ACTIVE
PROCAP, PROCAP_ACTIVE
RNA_SEQ, RNA_SEQ_ACTIVE
SPLICE_JUNCTIONS, SPLICE_SITE_USAGE, SPLICE_SITES
POLYADENYLATION
```

## Organism

```python
Organism.HOMO_SAPIENS   # default
Organism.MUS_MUSCULUS   # also supported, but most scorers may not
```

Mouse is "supported" but limited — many variant scorers don't have mouse
track data. Don't assume parity.

## Model versions

```python
ModelVersion.ALL_FOLDS  # = 1, default — uses all 5 folds
ModelVersion.FOLD_0 ... FOLD_3  # individual folds for ensembling
```

5-fold cross-validation. The default uses the ensemble (better predictions,
no calibration breakdown). For benchmarking, you want individual folds so
you can measure variance.

## Predict vs score vs ISM

| Function | Returns | Use case |
|---|---|---|
| `predict_sequence(seq, ...)` | `Output` (all modalities) | Custom sequence, raw track values |
| `predict_interval(interval, ...)` | `Output` | Real genomic region, raw track values |
| `predict_variant(interval, variant, ...)` | `VariantOutput` (ref + alt) | Compare ref vs alt predictions |
| `score_variant(interval, variant, [scorers], ...)` | `list[AnnData]` | Aggregate ref/alt delta into a per-gene score |
| `score_ism_variants(interval, ism_interval, ...)` | `list[list[AnnData]]` | Mutate every base in a sub-region |
| `AtlasClient.get_interval(interval, scorers, ...)` | `Mapping[str, AnnData]` | Pre-computed scores for all variants in a region |

## Variant scorers — pick wisely

`variant_scorers.RECOMMENDED_VARIANT_SCORERS` gives you 18 sensible defaults,
which is under the per-request limit of 20. To submit them all in one call:

```python
scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
# or pick a subset for specific modalities:
scorers = [
    variant_scorers.RECOMMENDED_VARIANT_SCORERS[k]
    for k in ['SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS']
]
```

For bulk processing, ALWAYS use all 18 — it's free, gives you per-modality
predictions, and lets you filter post-hoc.

## Atlas internals

```python
_INTERVAL_CHUNK_SIZE = 32   # bp per chunk when querying an interval
DEFAULT_MAX_WORKERS = 10    # parallel queries
```

`AtlasClient.get_interval()` chunks your interval into 32 bp pieces and
queries each in parallel with 10 workers. Total work ≈ interval_width × 3
variants (ref + 3 alts per position).

If you want to score a 1 Mb region, that's ~3M variants split into ~31,250
chunks. Set `max_workers` higher (50–100) if you have network headroom.

## Output shapes

| Output | Shape | Notes |
|---|---|---|
| `TrackData.values` | `(sequence_length, n_tracks)` | For 1 Mb + 1 tissue = (1_048_576, 1) |
| `VariantOutput.ref.X` | `(n_tracks,)` or `(n_tracks, n_bins)` | Per-tissue/per-bin ref prediction |
| `AnnData` from `score_variant` | `(1, n_tracks)` or `(n_genes, n_tracks)` | The `.var` is track metadata |
| `AnnData` from `get_interval` | `(n_variants × 3_alts, n_tracks)` | `.obs` has variant info |

## Common gotchas

1. **Don't pass `int` args to subprocess.run in Python 3.14** — pass `'60'`
   as a string. This is a Hermes quirk in our env.

2. **Sequence padding must center the real sequence** in the 1 Mb window:
   `'ACGT'.center(SEQUENCE_LENGTH_1MB, 'N')` puts ACGT in the middle,
   not the start. The model treats the whole window as the input.

3. **`interval.resize(N)` centers on the interval's midpoint**. To pad a
   variant to 1 Mb, use `variant.reference_interval.resize(...)`.

4. **`ontology_terms` is REQUIRED** for predict calls (no default). Pass
   at least one UBERON / EFO term or you'll get an error.

5. **API key in shell env, never in code**. Use `os.environ.get(...)` and
   fail loudly if missing.

6. **Rate limits are undocumented but real**. Watch for `RESOURCE_EXHAUSTED`
   gRPC errors. The free tier is for "1000s of predictions, not millions".

## File locations that matter

- GENCODE GTF for hg38: `https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather` (cached locally at `data/gencode.v46.annotation.gtf.gz.feather`)
- Model version info: `CHANGELOG.md` in the upstream repo
- Variant scorer spec: `src/alphagenome/models/variant_scorers.py`
- Atlas query protos: `src/alphagenome/protos/atlas_service.proto`
