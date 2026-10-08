import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import A, OI, OUTDIR, W, lab, save, placeholder
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
R = pd.read_csv(A + 'atlas_results.csv'); P1 = pd.read_csv(OUTDIR + 'p1_spectra.csv')
R = pd.read_csv(A + 'atlas_results.csv'); L = pd.read_csv(OUTDIR + 'isochore_v2.csv'); L = L.rename(columns={**{f'held_pred_L{k}': f'pred_L{k}' for k in (1, 2, 5, 10, 20, 30)}, **{f'held_obs_L{k}': f'obs_gc_component_L{k}' for k in (1, 2, 5, 10, 20, 30)}}); P4 = pd.read_csv(OUTDIR + 'p4_domain_scale_tests.csv'); P5 = pd.read_csv(OUTDIR + 'p5a_liver_gc.csv'); S = pd.read_csv(OUTDIR + 'sim_decay_recovery.csv'); P1 = pd.read_csv(OUTDIR + 'p1_spectra.csv')
fig = plt.figure(figsize=(W, W * 0.66)); gs = fig.add_gridspec(2, 4, hspace=0.32, wspace=0.6, height_ratios=[1.55, 1])
# a: what has to be removed before co-expression can be read (schematic)
ax = fig.add_subplot(gs[0, :]); placeholder(ax, 'Schematic: separating sources of expression variation', 'Tissue average, GC content, copy number, remaining variation.', key='fig1a'); lab(ax, 'a', -0.02)
# b: reproducible peaks do not depend on real gene order
ax = fig.add_subplot(gs[1, 0]); ax.scatter(R.universal_peaks_random_order, R.universal_peaks_real_order, s=10, color=OI['blue'], lw=0)
m = max(R.universal_peaks_real_order.max(), R.universal_peaks_random_order.max()) * 1.05; ax.plot([0, m], [0, m], 'k:', lw=0.7)
ax.set_xlabel('Peaks, randomised gene order'); ax.set_ylabel('Peaks, real gene order'); ax.text(0.05, 0.9, f'median ratio {np.median(R.universal_peaks_real_order / R.universal_peaks_random_order):.2f}', transform=ax.transAxes, fontsize=6); lab(ax, 'b')
# c: correlation left after each correction (liver)
ax = fig.add_subplot(gs[1, 1]); lbl = ['no correction', '− cell composition', '− GC trend', '− 5 PCs']
cis = [P4.cis_L1.iloc[0], P5.cis_L1.iloc[0], P5.cis_L1.iloc[1], P4.cis_L1.iloc[3]]; dom = [P4.domain_L10_30.iloc[0], P5.domain_L10_30.iloc[0], P5.domain_L10_30.iloc[1], P4.domain_L10_30.iloc[3]]
x = np.arange(4); ax.bar(x - 0.2, cis, 0.4, color=OI['blue'], label='adjacent genes'); ax.bar(x + 0.2, dom, 0.4, color=OI['red'], label='10–30 genes apart')
ax.set_xticks(x); ax.set_xticklabels(lbl, fontsize=5.6, rotation=35, ha='right'); ax.set_ylabel('Correlation above\nrandomised gene order'); ax.set_ylim(0, 0.34); ax.legend(loc='upper right', fontsize=5.5); lab(ax, 'c')
# d: GC correction removes most domain-scale correlation in every tissue
ax = fig.add_subplot(gs[1, 2])
for _, r in R.iterrows(): ax.plot([0, 1], [r.domain_L10_30_raw, r.domain_minusGC], color='0.75', lw=0.5)
ax.plot([0, 1], [R.domain_L10_30_raw.median(), R.domain_minusGC.median()], 'o-', color=OI['red'], lw=1.4, ms=4, label='median, 36 tissues')
ax.set_xticks([0, 1]); ax.set_xticklabels(['Before GC\ncorrection', 'After GC\ncorrection']); ax.set_xlim(-0.3, 1.3); ax.set_ylabel('Correlation, genes 10–30 apart'); ax.legend(loc='upper right', fontsize=5.5)
lab(ax, 'd')
# e: the GC contribution is predicted on held-out chromosomes
ax = fig.add_subplot(gs[1, 3])
for lg, col in [(1, OI['blue']), (5, OI['green']), (20, OI['red'])]:
    ax.scatter(L[f'pred_L{lg}'], L[f'obs_gc_component_L{lg}'], s=9, color=col, lw=0, label=f'{lg} gene' + ('s' if lg > 1 else '') + f' apart (r = {np.corrcoef(L[f"pred_L{lg}"], L[f"obs_gc_component_L{lg}"])[0, 1]:.2f})')
m = max(L[[c for c in L.columns if c.startswith('obs')]].max().max(), L[[c for c in L.columns if c.startswith('pred')]].max().max()) * 1.05
ax.plot([0, m], [0, m], 'k:', lw=0.7); ax.set_xlabel('Predicted (odd chr.)'); ax.set_ylabel('Observed (even chr.)'); ax.legend(loc='lower right', fontsize=5.5, handletextpad=0.1)
lab(ax, 'e')
save(fig, 'Fig1_nonpositional_layers'); print('Fig1_nonpositional_layers ok')
