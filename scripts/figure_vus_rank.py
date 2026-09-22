#!/usr/bin/env python
"""Figure 3: VUS re-scoring — rank vs score plot showing the 61 high-impact cluster."""
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/Users/hermes/projects/alphagenome-work/outputs/vus_rescored.csv')
df = df.sort_values('SPLICE_SITES_score', ascending=False).reset_index(drop=True)
df['rank'] = range(1, len(df) + 1)

fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(df['rank'], df['SPLICE_SITES_score'],
           s=4, alpha=0.4, color='steelblue', label='All VUS')
high = df[df['SPLICE_SITES_score'] >= 0.5]
ax.scatter(high['rank'], high['SPLICE_SITES_score'],
           s=10, alpha=0.7, color='crimson', label=f'High-impact (score ≥ 0.5, n={len(high)})')
ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, label='Threshold')

ax.set_xlabel('Variant rank (descending AlphaGenome score)', fontsize=11)
ax.set_ylabel('SPLICE_SITES_score (AlphaGenome)', fontsize=11)
ax.set_title(f'SCN1A VUS re-scoring — {len(high)} of {len(df)} VUS flagged as high-impact', fontsize=13)
ax.legend(loc='upper right', fontsize=10)
ax.set_ylim(-0.05, 2.0)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('/Users/hermes/projects/alphagenome-work/figures/vus_rank_score.png', dpi=150)
print(f"Saved → figures/vus_rank_score.png")
print(f"  {len(high)}/{len(df)} variants above 0.5 threshold")
