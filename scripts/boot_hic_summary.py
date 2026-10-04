"""Summarise the donor x genomic-block bootstrap of boot_hic.py into estimates and 95% percentile intervals (Fig. 4d,f; ST8D)."""
import numpy as np, pandas as pd
from poslayers.config import OUTDIR
rows = []
for cell in ('GM12878', 'IMR90'):
    d = pd.read_csv(OUTDIR + f'boot_hic_{cell}.csv'); pt = d[d.rep == -1].iloc[0]; b = d[d.rep >= 0]
    for k in ['contact_coef', 'k_distance', 'k_contact', 'k_t1', 'k_t2', 'k_t3']:
        lo, hi = np.nanpercentile(b[k], [2.5, 97.5])
        rows.append({'cell': cell, 'stat': k, 'estimate': pt[k], 'ci_low': lo, 'ci_high': hi, 'boot_se': b[k].std(), 'B': b[k].notna().sum()})
    g = b.k_t1 - b.k_t3
    rows.append({'cell': cell, 'stat': 'k_t1_minus_k_t3', 'estimate': pt.k_t1 - pt.k_t3, 'ci_low': np.nanpercentile(g, 2.5), 'ci_high': np.nanpercentile(g, 97.5), 'boot_se': g.std(), 'B': g.notna().sum()})
    rows.append({'cell': cell, 'stat': 'P(k_t1>k_t2>k_t3) across replicates', 'estimate': np.mean((b.k_t1 > b.k_t2) & (b.k_t2 > b.k_t3)), 'ci_low': np.nan, 'ci_high': np.nan, 'boot_se': np.nan, 'B': len(b)})
R = pd.DataFrame(rows); R.to_csv(OUTDIR + 'boot_hic_summary.csv', index=False); print(R.round(3).to_string(index=False))
