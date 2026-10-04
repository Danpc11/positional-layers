import sys, glob; sys.path.insert(0, '/home/claude/figs'); from style import *
E = pd.read_csv(A + 'eqtl_cis_test.csv'); Dz = pd.read_csv(A + 'coloc_dose_response.csv'); SP = pd.read_csv(A + 'eqtl_tissue_specificity.csv'); Cc = pd.read_csv(A + 'coloc_cis_test.csv')
PC = pd.concat([pd.read_csv(f) for f in glob.glob(A + 'predcal_*.csv.gz')])
import schem
fig = plt.figure(figsize=(W, W * 0.55)); gs = fig.add_gridspec(2, 4, hspace=0.6, wspace=0.65)

ax = fig.add_subplot(gs[0, 0:2]); schem.shared_variant(ax); lab(ax, 'a', -0.04)
ax = fig.add_subplot(gs[0, 2]); cats = [('r_not_both_eGenes', 'not both\neGenes', '0.6'), ('r_eGenes_no_share', 'no shared\nvariant', OI['sky']), ('r_share_same', 'shared,\nsame sign', OI['blue'])]
for i, (c, n, col) in enumerate(cats): ax.scatter(i + rng.uniform(-0.15, 0.15, len(E)), E[c], s=6, color=col, lw=0); ax.hlines(E[c].median(), i - 0.25, i + 0.25, color='k', lw=1.2)
ax.set_xticks(range(3)); ax.set_xticklabels([n for _, n, _ in cats], fontsize=5.4); ax.set_ylabel('Adjacent-pair correlation'); lab(ax, 'b')
ax = fig.add_subplot(gs[0, 3]); s = E.sort_values('delta_same_adj'); y = np.arange(len(s))
ax.scatter(s.delta_same_adj, y, s=8, color=OI['blue'], label='same direction'); ax.scatter(s.delta_opposite_only_adj, y, s=8, color=OI['red'], label='opposite direction only')
ax.axvline(0, color='k', lw=0.6); ax.set_yticks([]); ax.set_ylabel('36 tissues'); ax.set_xlabel('Effect on coupling (adjusted)')
lab(ax, 'c', -0.08)
ax = fig.add_subplot(gs[1, 0]); order = ['0', '0-0.1', '0.1-0.5', '0.5-0.8', '>0.8']; g = Dz.groupby('p_coloc_bin').mean_r
med = [g.get_group(b).median() for b in order]; q1 = [g.get_group(b).quantile(0.25) for b in order]; q3 = [g.get_group(b).quantile(0.75) for b in order]
ax.errorbar(range(5), med, yerr=[np.subtract(med, q1), np.subtract(q3, med)], fmt='o-', color=OI['purple'], ms=4, capsize=2, lw=1.2)
ax.set_xticks(range(5)); ax.set_xticklabels(['none', '<0.1', '0.1–0.5', '0.5–0.8', '>0.8'], fontsize=5.4, rotation=30); ax.set_xlabel('P(same causal variant)'); ax.set_ylabel('Adjacent-pair correlation')
lab(ax, 'd')
ax = fig.add_subplot(gs[1, 1]); PC['bin'] = pd.qcut(PC.pred_genetic_r, 8, duplicates='drop'); b = PC.groupby('bin', observed=True).agg(p=('pred_genetic_r', 'mean'), o=('obs_r', 'mean'), se=('obs_r', lambda x: x.std() / np.sqrt(len(x))))
ax.errorbar(b.p, b.o, yerr=1.96 * b.se, fmt='o', color=OI['blue'], ms=4, capsize=2); sl = np.polyfit(PC.pred_genetic_r, PC.obs_r, 1); xx = np.linspace(b.p.min(), b.p.max(), 10)
ax.set_xscale('symlog', linthresh=0.01); ax.set_xlabel('Predicted genetic correlation 2p(1−p)β₁β₂'); ax.set_ylabel('Observed coupling')
ax.text(0.05, 0.9, f'slope {sl[0]:.2f}; r = {np.corrcoef(PC.pred_genetic_r, PC.obs_r)[0, 1]:.2f}', transform=ax.transAxes, fontsize=6); lab(ax, 'e')
ax = fig.add_subplot(gs[1, 2:4]); ax.scatter(SP.b_other, SP.b_own, s=4, color=OI['green'], alpha=0.5, lw=0); m = max(SP.b_own.max(), SP.b_other.max()); mn = min(SP.b_own.min(), SP.b_other.min())
ax.plot([mn, m], [mn, m], 'k:', lw=0.7); ax.set_xlabel('Effect of sharing in another tissue'); ax.set_ylabel('Effect of sharing in the same tissue')
lab(ax, 'f')
save(fig, 'Fig3_genetic_component'); print('Fig3_genetic_component ok')
