import sys; sys.path.insert(0, '/home/claude/figs'); from style import *
R = pd.read_csv(A + 'atlas_results.csv'); O = pd.read_csv(A + 'orientation_by_tissue.csv'); F = pd.read_csv(A + 'robust_families.csv'); Lz = pd.read_csv(A + 'lorentz_results.csv')
import schem
fig = plt.figure(figsize=(W, W * 0.55)); gs = fig.add_gridspec(2, 4, hspace=0.55, wspace=0.65)

ax = fig.add_subplot(gs[0, 0:2]); schem.two_scales(ax); lab(ax, 'a', -0.04)
ax = fig.add_subplot(gs[0, 2]); s = R.sort_values('cis_L1_minusGC_tech'); y = np.arange(len(s)); cl = s.tissue.str.startswith('cells')
ax.barh(y, s.cis_L1_minusGC_tech, color=np.where(cl, OI['orange'], OI['blue']), height=0.75)
ax.set_yticks(y[::3]); ax.set_yticklabels([tissue_label(t) for t in s.tissue][::3], fontsize=5.3); ax.set_xlabel('Correlation of adjacent genes\n(GC- and batch-corrected)')
ax.text(1.0, 1.01, 'orange: cell lines', transform=ax.transAxes, ha='right', va='bottom', fontsize=5.8, color=OI['orange']); lab(ax, 'b', -0.62)
ax = fig.add_subplot(gs[0, 3]); O['excess'] = O.mean_r - O.random_pair_floor
order = ['0-1kb', '1-5kb', '5-20kb', '20-100kb', '100-500kb', '>500kb']; mids = [0.5, 3, 12, 50, 250, 1000]
for o, col in [('divergent', OI['green']), ('tandem', OI['blue']), ('convergent', OI['red'])]:
    m = O[O.orientation == o].groupby('dist_bin').excess.median().reindex(order); ax.plot(mids, m.values, 'o-', color=col, ms=3, lw=1, label=o)
ax.set_xscale('log'); ax.set_xlabel('Intergenic distance (kb)'); ax.set_ylabel('Correlation above random pairs'); ax.legend(loc='upper right'); lab(ax, 'c')
ax = fig.add_subplot(gs[1, 0]); cats = [('r_same_stem', 'same gene family'), ('r_readthrough', 'read-through'), ('r_clusters', 'classic clusters'), ('r_all', 'all pairs'), ('r_clean', 'clean pairs')]
med = [F[c].median() for c, _ in cats]; lo = [F[c].quantile(0.25) for c, _ in cats]; hi = [F[c].quantile(0.75) for c, _ in cats]
cols = ['0.6', '0.6', '0.6', OI['blue'], OI['green']]; x = np.arange(len(cats)); ax.bar(x, med, color=cols, width=0.65); ax.errorbar(x, med, yerr=[np.subtract(med, lo), np.subtract(hi, med)], fmt='none', ecolor='k', lw=0.7, capsize=2)
ax.set_xticks(x); ax.set_xticklabels([n for _, n in cats], rotation=30, ha='right', fontsize=6); ax.set_ylabel('Adjacent-pair correlation'); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 1]); sp = np.load(A + 'spec_thyroid.npy'); ax.loglog(sp[0], sp[1], 'o', ms=2.5, color='0.25', label='thyroid, 684 samples'); ax.loglog(sp[0], sp[2], color=OI['blue'], lw=1.2, label='two Lorentzians')
ax.loglog(sp[0], sp[3], color=OI['red'], lw=0.9, ls='--', label='power law'); ax.set_xlabel('Spatial frequency (cycles per gene)'); ax.set_ylabel('Covariance spectrum S(f)'); ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3, fontsize=5.2, columnspacing=0.8)
lab(ax, 'e')
ax = fig.add_subplot(gs[1, 2:4]); xs = np.arange(len(Lz))
ax.scatter(xs - 0.12, Lz.lam_acf_short, s=14, color=OI['blue'], label='short, autocorrelation'); ax.scatter(xs + 0.12, Lz.lam_spec_short, s=14, marker='s', color=OI['sky'], label='short, spectrum')
ax.scatter(xs - 0.12, Lz.lam_acf_long, s=14, color=OI['red'], label='long, autocorrelation'); ax.set_yscale('log')
ax.set_xticks(xs); ax.set_xticklabels([tissue_label(t) for t in Lz.tissue], rotation=35, ha='right', fontsize=5.8); ax.set_ylabel('Decay length (genes)'); ax.set_ylim(0.6, 300); ax.legend(loc='upper right', fontsize=5.4, ncol=3)
lab(ax, 'f')
save(fig, 'Fig2_cis_layer'); print('Fig2_cis_layer ok')
