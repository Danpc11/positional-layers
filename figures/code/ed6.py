"""Extended Data Fig. 4: coupling is a reproducible, local property of gene pairs (coupling_robustness.py)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import OI, OUTDIR, W, lab, save, tissue_label
import matplotlib.pyplot as plt

S = pd.read_csv(OUTDIR + 'coupling_robustness.csv'); P = pd.read_csv(OUTDIR + 'coupling_robustness_pairs_thyroid.csv.gz')
TIS = ['thyroid', 'cells_ebv-transformed_lymphocytes', 'muscle_skeletal', 'lung']
fig = plt.figure(figsize=(W, W * 0.6)); gs = fig.add_gridspec(2, 6, hspace=0.75, wspace=1.6)

def hexpanel(ax, x, y, xl, yl, title, cmap):
    lim = (-0.4, 0.8)
    ax.hexbin(x, y, gridsize=45, extent=lim + lim, cmap=cmap, bins='log', mincnt=1, linewidths=0)
    ax.plot(lim, lim, color='0.5', lw=0.6, ls=':'); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title, fontsize=6.5)
    ax.text(0.04, 0.96, f'r = {np.corrcoef(x, y)[0, 1]:.2f}\n{len(x):,} pairs', transform=ax.transAxes, va='top', fontsize=5.5)
adj, dis = P[P.type == 'adjacent'], P[P.type == 'matched distant']
ax = fig.add_subplot(gs[0, 0:2]); hexpanel(ax, adj.pearson, adj.spearman, 'Pearson coupling', 'Spearman coupling', 'Adjacent pairs, thyroid', 'Blues'); lab(ax, 'a')
ax = fig.add_subplot(gs[0, 2:4]); hexpanel(ax, adj.half_a, adj.half_b, 'Coupling, donor half A', 'Coupling, donor half B', 'Adjacent pairs, thyroid', 'Blues'); lab(ax, 'b')
ax = fig.add_subplot(gs[0, 4:6]); hexpanel(ax, dis.half_a, dis.half_b, 'Coupling, donor half A', 'Coupling, donor half B', 'Matched distant pairs, thyroid', 'Greys'); lab(ax, 'c')
x = np.arange(len(TIS)); labels = [tissue_label(t) for t in TIS]
ax = fig.add_subplot(gs[1, 0:3])
for k, (col, name, off) in enumerate([('split_half_reproducibility_adjacent', 'adjacent', -0.12), ('split_half_reproducibility_distant', 'matched distant', 0.12)]):
    c = OI['blue'] if k == 0 else '0.45'
    for K, mk, fill in ((0, 'o', 'white'), (15, 'o', c)):
        v = [S[(S.tissue == t) & (S.pcs_removed == K)][col].iloc[0] for t in TIS]
        ax.scatter(x + off, v, s=22, marker=mk, facecolor=fill, edgecolor=c, lw=1, zorder=3, label=f'{name}, {"15 PCs removed" if K else "no PCs removed"}')
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=5.5); ax.set_ylim(0, 1.05); ax.set_ylabel('Split-half reproducibility\n(r across pairs)')
lg = ax.legend(loc='lower left', fontsize=5.5, ncol=2, handletextpad=0.2, columnspacing=0.8); lg.keep_position = True; lab(ax, 'd')
ax = fig.add_subplot(gs[1, 3:6]); wbar = 0.35
S15 = S[S.pcs_removed == 15].set_index('tissue').loc[TIS]
ax.bar(x - wbar / 2, 100 * S15.frac_adjacent_q05_positive, wbar, color=OI['blue'], label='adjacent, positive')
ax.bar(x + wbar / 2, 100 * S15.frac_distant_q05_positive, wbar, color='0.6', label='matched distant, positive')
ax.bar(x - wbar / 2, -100 * S15.frac_adjacent_q05_negative, wbar, color=OI['blue'], alpha=0.45, label='adjacent, negative')
ax.bar(x + wbar / 2, -100 * S15.frac_distant_q05_negative, wbar, color='0.6', alpha=0.45, label='matched distant, negative')
ax.axhline(0, color='k', lw=0.6); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=5.5); ax.set_ylabel('Pairs with q < 0.05 (%)\n(negative values: negative coupling)')
ax.set_ylim(-25, 60); lg = ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), fontsize=5.5, ncol=2, handletextpad=0.3, columnspacing=0.8, borderaxespad=0); lg.keep_position = True; lab(ax, 'e')
save(fig, 'ExtData_Fig6_coupling_robustness'); print('ed6 ok')
