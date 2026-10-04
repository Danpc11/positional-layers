import os
DATA = os.environ.get('POSLAYERS_DATA', 'data').rstrip('/') + '/'
"""Test of the hub model: r = r0 + A*(1 + (c*/c)^(2/3))^(-3/2) versus a pure power law r = r0 + A*c^k.
Pairs binned by Hi-C contact (all distances 25 kb-2 Mb); local slope of the coupling-contact relation across contact."""
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
hub = lambda c, r0, A, cs: r0 + A * (1 + (cs / c) ** (2 / 3)) ** (-1.5)
pw = lambda c, r0, A, k: r0 + A * c ** k
out = []
for cell in ['GM12878', 'IMR90']:
    H = pd.read_csv(fDATA + 'hic_coupling_{cell}.csv.gz'); H = H[H.contact_KR > 0]
    H['cb'] = pd.qcut(np.log10(H.contact_KR), 20, labels=False)
    b = H.groupby('cb').agg(c=('contact_KR', 'median'), r=('r', 'mean'), se=('r', lambda x: x.std() / np.sqrt(len(x))), d=('tss_distance', 'median'), n=('r', 'size'))
    w = 1 / b.se ** 2
    p1, _ = curve_fit(hub, b.c, b.r, p0=[0.005, 0.2, 500], sigma=b.se, maxfev=50000, bounds=([-0.05, 0, 1], [0.05, 5, 1e6]))
    p2, _ = curve_fit(pw, b.c, b.r, p0=[0.0, 0.002, 0.5], sigma=b.se, maxfev=50000, bounds=([-0.05, 0, 0.01], [0.05, 5, 3]))
    chi1 = np.sum(w * (b.r - hub(b.c, *p1)) ** 2); chi2 = np.sum(w * (b.r - pw(b.c, *p2)) ** 2)
    x = (p1[2] / b.c) ** (2 / 3); b['k_model'] = x / (1 + x)
    # empirical local slope of (r - r0) vs contact, from adjacent contact bins
    y = np.log(np.clip(b.r - p1[0], 1e-5, None)); lc = np.log(b.c); b['k_empirical'] = np.gradient(y, lc)
    print(f'== {cell}: hub model r0={p1[0]:.4f} A={p1[1]:.3f} c*={p1[2]:.0f}  chi2={chi1:.1f} | power law k={p2[2]:.2f} chi2={chi2:.1f}  (20 bins)')
    dstar = np.interp(np.log(p1[2]), np.log(b.c.values[::-1]), b.d.values[::-1])
    print(f'   distance at which contact = c* (separation = hub size): ~{dstar / 1e3:.0f} kb')
    print(b[['n', 'd', 'c', 'r', 'k_model', 'k_empirical']].round(4).to_string())
    b['cell'] = cell; out.append(b)
pd.concat(out).to_csv('derivation_hub_model.csv')
