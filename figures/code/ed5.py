"""Extended Data Fig. 5: replication and robustness of the two architectural barriers."""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import OI, OUTDIR, W, lab, save
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
fig = plt.figure(figsize=(W, W * 0.62)); gs = fig.add_gridspec(2, 2, hspace=0.75, wspace=0.85)
TIS = [('cells_ebv-transformed_lymphocytes', 'lymphoblastoid', OI['blue']), ('thyroid', 'thyroid', OI['green']), ('cells_cultured_fibroblasts', 'fibroblasts', OI['red'])]
# a: active intervening gene, absolute coupling by quartile
ax = fig.add_subplot(gs[0, 0])
for t, nm, col in TIS:
    d = pd.read_csv(OUTDIR + f'active_gene_barrier_means_{t}.csv'); ax.errorbar(d.quartile, d.coupling, yerr=[d.coupling - d.ci_low, d.ci_high - d.coupling], fmt='o-', color=col, ms=3.5, lw=1, capsize=1.5, label=nm)
ax.set_xticks([1, 2, 3, 4]); ax.set_xlabel('Expression of the intervening gene (quartile)'); ax.set_ylabel('Coupling of flanking genes'); ax.legend(fontsize=5.5, loc='lower left'); lab(ax, 'a')
# b: coefficients, primary and robust models
ax = fig.add_subplot(gs[0, 1]); rows = []
for t, nm, col in TIS:
    d = pd.read_csv(OUTDIR + f'active_gene_barrier_{t}.csv').set_index('term')
    for term, lbl, mk in [('expression of the intervening gene (primary)', 'primary model', 'o'), ('robust: expression of j', '+ distance, gene length, orientation', 's'), ('robust: log length of j', 'gene length of j', 'D')]:
        q = d.loc[term]; rows.append((f'{nm}: {lbl}', q.coef_per_SD, q.ci_low, q.ci_high, col, mk))
for i, (n, e, l, h, col, mk) in enumerate(rows[::-1]): ax.errorbar(e, i, xerr=[[e - l], [h - e]], fmt=mk, color=col, ms=3, capsize=1.5, lw=0.8)
ax.axvline(0, color='k', lw=0.6); ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=5.5); ax.set_xlabel('Effect on coupling of flanking genes (per s.d.)'); lab(ax, 'b', -0.9)
# c: CTCF loops crossed, absolute coupling
ax = fig.add_subplot(gs[1, 0])
for cell, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    d = pd.read_csv(OUTDIR + f'ctcf_barrier_means_{cell}.csv'); ax.errorbar(d.crossed_loops, d.coupling, yerr=[d.coupling - d.ci_low, d.ci_high - d.coupling], fmt='o-', color=col, ms=3.5, lw=1, capsize=1.5, label=cell)
ax.set_xticks(range(4)); ax.set_xticklabels(['0', '1', '2', '≥3']); ax.set_xlabel('Convergent CTCF loops crossed'); ax.set_ylabel('Coupling (adjusted for\ndistance and contact)'); ax.legend(fontsize=5.5); lab(ax, 'c')
# d: per-loop coefficients with and without contact adjustment
ax = fig.add_subplot(gs[1, 1]); rows = []
for cell, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    d = pd.read_csv(OUTDIR + f'ctcf_barrier_{cell}.csv')
    for adj, lbl, mk in [(False, 'equal distance', 'o'), (True, 'equal distance and contact', 's')]:
        q = d[(d.adjusted_for_contact == adj) & (d.term == 'per crossed loop')].iloc[0]; rows.append((f'{cell}: {lbl}', q.coef, q.ci_low, q.ci_high, col, mk))
for i, (n, e, l, h, col, mk) in enumerate(rows[::-1]): ax.errorbar(e, i, xerr=[[e - l], [h - e]], fmt=mk, color=col, ms=3, capsize=1.5, lw=0.8)
ax.axvline(0, color='k', lw=0.6); ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=5.5); ax.set_xlabel('Change in coupling per crossed loop'); lab(ax, 'd', -0.9)
save(fig, 'ExtData_Fig5_barriers'); print('ed5 ok')
