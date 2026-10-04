"""Saturating-response prediction: the contact exponent k should be lower for highly expressed (closer to saturation) pairs."""
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
pw = lambda c, r0, A, k: r0 + A * c ** k
for cell, tissue in [('GM12878', 'cells_ebv-transformed_lymphocytes'), ('IMR90', 'cells_cultured_fibroblasts')]:
    C = pd.read_csv(f'/home/claude/repo/data/external/gtex/gene_reads_adult_gtex_v11_{tissue}_gct.gz', sep='\t', skiprows=2, index_col=0).drop(columns='Description')
    C.index = C.index.str.split('.').str[0]; C = C[~C.index.duplicated()]; lv = np.log2(C / C.sum() * 1e6 + 1).mean(1)
    H = pd.read_csv(f'/home/claude/atlas/hic_coupling_{cell}.csv.gz'); H = H[H.contact_KR > 0].copy()
    H['expr'] = np.minimum(lv.reindex(H.g1).values, lv.reindex(H.g2).values); H = H.dropna(subset=['expr'])
    H['et'] = pd.qcut(H.expr, 3, labels=['low', 'mid', 'high'])
    for t, d in H.groupby('et', observed=True):
        d = d.copy(); d['cb'] = pd.qcut(np.log10(d.contact_KR), 15, labels=False)
        b = d.groupby('cb').agg(c=('contact_KR', 'median'), r=('r', 'mean'), se=('r', lambda x: x.std() / np.sqrt(len(x))))
        p, cov = curve_fit(pw, b.c, b.r, p0=[0.0, 0.002, 0.5], sigma=b.se, maxfev=50000, bounds=([-0.05, 0, 0.01], [0.05, 5, 3]))
        print(f'{cell} expression tertile {t:>4}: k = {p[2]:.2f} ± {np.sqrt(cov[2, 2]):.2f} | mean coupling {d.r.mean():.3f} | pairs {len(d):,}')
