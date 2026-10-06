import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import A, OI, OUTDIR, W, lab, save, tissue_label, placeholder
import matplotlib.pyplot as plt
R = pd.read_csv(A + 'atlas_results.csv'); O = pd.read_csv(A + 'orientation_by_tissue.csv'); F = pd.read_csv(A + 'robust_families.csv'); Lz = pd.read_csv(OUTDIR + 'lorentz_v2.csv')
fig = plt.figure(figsize=(W, W * 1.1)); gs = fig.add_gridspec(3, 4, hspace=0.75, wspace=0.7, height_ratios=[1, 1.25, 1.05])

ax = fig.add_subplot(gs[0, 0:2]); placeholder(ax, 'Schematic: coupling between neighbouring genes', 'Coupling = correlation of corrected expression between adjacent genes across individuals.\nTwo genomic scales: short range (~1 gene) and domain scale (7–24 genes).\nPanel b ranks tissues by coupling.', key='fig2a'); p_ = ax.get_position(); ax.set_position([p_.x0, p_.y0, p_.width * 0.8, p_.height]); lab(ax, 'a', -0.04)   # leave room for the tissue names of b
# b: ranked dot plot -- coupling of every tissue on one common scale; tumours, a related measure, on their own scale
CV = pd.read_csv(OUTDIR + 'cohesin_v2_results.csv').drop_duplicates('cohort').set_index('cohort').mean_cis_excess_wt
BA = pd.read_csv(OUTDIR + 'beataml_freedman_lane.csv').iloc[0].cis_excess_wt
gt = R.sort_values('cis_L1_minusGC_tech').reset_index(drop=True)
tum = [('BLCA', 'Bladder'), ('UCEC', 'Endometrium'), ('STAD', 'Stomach'), ('COAD', 'Colon'), ('AML', 'AML')]
tumv = pd.Series({nm: (BA if c == 'AML' else CV.get(c, np.nan)) for c, nm in tum}).sort_values()
sub_b = gs[0:2, 2:4].subgridspec(2, 1, height_ratios=[len(gt), len(tumv) + 1.2], hspace=0.28)
ax = fig.add_subplot(sub_b[0]); ax.set_zorder(3); y = np.arange(len(gt)); cult = gt.tissue.str.startswith('cells').values   # zorder: tissue names drawn above panel c
ax.grid(axis='y', color='0.92', lw=0.5); ax.set_axisbelow(True)
ax.scatter(gt.cis_L1_minusGC_tech[~cult], y[~cult], s=11, color=OI['blue'], zorder=3, label='GTEx tissue')
ax.scatter(gt.cis_L1_minusGC_tech[cult], y[cult], s=13, color=OI['orange'], marker='D', zorder=3, label='cultured cells')
ax.axvline(gt.cis_L1_minusGC_tech.median(), color='0.6', lw=0.6, ls=':')
ax.set_yticks(y); ax.set_yticklabels([tissue_label(t) for t in gt.tissue], fontsize=5.5); ax.set_ylim(-0.7, len(gt) - 0.3); ax.set_xlim(0.09, 0.19)
ax.set_xlabel('Adjacent-gene coupling (after mean, GC and technical correction)', fontsize=5.6); ax.tick_params(axis='y', length=0, pad=1.5)
ax.legend(loc='lower right', fontsize=5.5, handletextpad=0.2, borderaxespad=0.2); lab(ax, 'b')
ax = fig.add_subplot(sub_b[1]); yt = np.arange(len(tumv))
ax.grid(axis='y', color='0.92', lw=0.5); ax.set_axisbelow(True); ax.scatter(tumv.values, yt, s=12, color=OI['red'], zorder=3)
ax.set_yticks(yt); ax.set_yticklabels(tumv.index, fontsize=5.5); ax.tick_params(axis='y', length=0, pad=1.5); ax.set_xlim(0.04, 0.1); ax.set_ylim(-0.7, len(tumv) - 0.3)
ax.set_xlabel('Cis excess, wild-type tumours (own scale)', fontsize=5.6)
ax = fig.add_subplot(gs[1, 0:2]); O['excess'] = O.mean_r - O.random_pair_floor
order = ['0-1kb', '1-5kb', '5-20kb', '20-100kb', '100-500kb', '>500kb']; mids = [0.5, 3, 12, 50, 250, 1000]
for o, col in [('divergent', OI['green']), ('tandem', OI['blue']), ('convergent', OI['red'])]:
    m = O[O.orientation == o].groupby('dist_bin').excess.median().reindex(order); ax.plot(mids, m.values, 'o-', color=col, ms=3, lw=1, label=o)
ax.set_xscale('log'); ax.set_xlabel('Intergenic distance (kb)'); ax.set_ylabel('Excess correlation'); ax.legend(loc='lower left'); lab(ax, 'c')
# d: UpSet of pair attributes (family/cluster/read-through, shared eQTL, same TAD) with distance-adjusted coupling
U = pd.read_csv(OUTDIR + 'upset_pair_attributes.csv'); ks = ['family', 'shared_eqtl', 'same_tad']
M = U.groupby(ks).agg(frac=('frac', 'median'), r=('r_adj', 'median'), q1=('r_adj', lambda x: x.quantile(.25)), q3=('r_adj', lambda x: x.quantile(.75))).reset_index().sort_values('r').reset_index(drop=True)
sub = gs[2, 0:2].subgridspec(3, 1, height_ratios=[1.5, 0.75, 0.75], hspace=0.08)
ax = fig.add_subplot(sub[0]); x = np.arange(len(M))
ax.errorbar(x, M.r, yerr=[M.r - M.q1, M.q3 - M.r], fmt='o', color=OI['blue'], ms=3.5, capsize=1.5, lw=0.8)
ax.axhline(M.r.iloc[0], color='0.7', lw=0.5, ls=':'); ax.set_xticks([]); ax.set_xlim(-0.6, len(M) - 0.4); ax.set_ylabel('Coupling', fontsize=6)
lab(ax, 'd')
ax2 = fig.add_subplot(sub[1], sharex=ax); ax2.bar(x, 100 * M.frac, color='0.55', width=0.6); ax2.set_yscale('log'); ax2.set_ylabel('% pairs', fontsize=5.5, rotation=0, ha='right', va='center'); ax2.set_xticks([])
ax3 = fig.add_subplot(sub[2], sharex=ax)
for k, key in enumerate(ks):
    for i in range(len(M)):
        on = bool(M[key].iloc[i]); ax3.plot(i, k, 'o', ms=3.5, color='k' if on else '0.85')
    ax3.plot(x[M[key].astype(bool)], [k] * int(M[key].sum()), 'k-', lw=0) 
for i in range(len(M)):
    ons = [k for k, key in enumerate(ks) if M[key].iloc[i]]
    if len(ons) > 1: ax3.plot([i, i], [min(ons), max(ons)], 'k-', lw=0.8)
ax3.set_yticks(range(3)); ax3.set_yticklabels(['gene family / cluster', 'shared eQTL', 'same TAD'], fontsize=5.5); ax3.set_ylim(-0.6, 2.6)
names = {(0, 0, 0): 'none', (0, 0, 1): 'TAD', (0, 1, 0): 'eQTL', (0, 1, 1): 'eQTL + TAD', (1, 0, 0): 'family', (1, 0, 1): 'family + TAD', (1, 1, 0): 'family + eQTL', (1, 1, 1): 'all three'}
ax3.set_xticks(x); ax3.set_xticklabels([names[tuple(int(bool(M[k].iloc[i])) for k in ks)] for i in range(len(M))], rotation=40, ha='right', fontsize=5.5); ax3.tick_params(axis='x', length=0)
ax3.set_xlabel('Attributes of the adjacent pair (filled dot = present)', fontsize=5.6)
for sp_ in ('top', 'right', 'bottom'): ax3.spines[sp_].set_visible(False)
ax3.invert_yaxis()
ax = fig.add_subplot(gs[2, 2]); sp = np.load(OUTDIR + 'spec_v2_thyroid.npy'); ax.loglog(sp[0], sp[1], 'o', ms=2.5, color='0.25', label='thyroid, 684 samples'); ax.loglog(sp[0], sp[2], color=OI['blue'], lw=1.2, label='two Lorentzians')
ax.loglog(sp[0], sp[3], color=OI['red'], lw=0.9, ls='--', label='power law'); ax.set_xlabel('Spatial frequency (cycles per gene)'); ax.set_ylabel('Covariance spectrum S(f)'); ax.legend(loc='lower left', fontsize=5.5)
lab(ax, 'e')
ax = fig.add_subplot(gs[2, 3]); xs = np.arange(len(Lz))
ax.scatter(xs - 0.12, Lz.lam_acf_short, s=14, color=OI['blue'], label='short, autocorrelation'); ax.scatter(xs + 0.12, Lz.lam_spec_short, s=14, marker='s', color=OI['sky'], label='short, spectrum')
ax.scatter(xs - 0.12, Lz.lam_acf_long, s=14, color=OI['red'], label='long, autocorrelation'); ax.set_yscale('log')
ax.set_xticks(xs); ax.set_xticklabels([tissue_label(t) for t in Lz.tissue], rotation=45, ha='right', fontsize=5.5); ax.set_ylabel('Decay length (genes)'); ax.set_ylim(0.6, 60); ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), fontsize=5.5, ncol=1, borderaxespad=0)
lab(ax, 'f')
save(fig, 'Fig2_cis_layer'); print('Fig2_cis_layer ok')
