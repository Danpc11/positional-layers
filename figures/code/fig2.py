import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import A, OI, OUTDIR, W, lab, save, tissue_label, placeholder
import matplotlib.pyplot as plt
R = pd.read_csv(A + 'atlas_results.csv'); O = pd.read_csv(A + 'orientation_by_tissue.csv'); F = pd.read_csv(A + 'robust_families.csv'); Lz = pd.read_csv(OUTDIR + 'lorentz_v2.csv')
fig = plt.figure(figsize=(W, W * 0.95)); gs = fig.add_gridspec(3, 4, hspace=0.75, wspace=0.7, height_ratios=[1, 1, 1.05])

ax = fig.add_subplot(gs[0, 0:2]); placeholder(ax, 'Schematic: coupling between neighbouring genes', 'Coupling = correlation of corrected expression between adjacent genes across individuals.\nTwo genomic scales: short range (~1 gene) and domain scale (7–24 genes).\nRings in b summarise four measures per tissue.'); lab(ax, 'a', -0.04)
# b: circular overview of every tissue and tumour cohort (replaces the per-tissue bar chart)
EQ = pd.read_csv(A + 'eqtl_cis_test.csv').set_index('tissue'); TB = pd.read_csv(OUTDIR + 'tad_block_bootstrap.csv').set_index('tissue')
CV = pd.read_csv(OUTDIR + 'cohesin_v2_results.csv').drop_duplicates('cohort').set_index('cohort').mean_cis_excess_wt
BA = pd.read_csv(OUTDIR + 'beataml_freedman_lane.csv').iloc[0].cis_excess_wt
TT = pd.read_csv(OUTDIR + 'tad_tumour_results.csv').drop_duplicates('cohort').set_index('cohort')
nrm = lambda t: t.lower().replace('-', '_')
gt = R.assign(k=R.tissue.map(nrm)).sort_values('cis_L1_minusGC_tech', ascending=False).reset_index(drop=True)
tum = ['BLCA', 'UCEC', 'STAD', 'COAD', 'AML']; tumv = [CV.get(c, np.nan) for c in tum[:-1]] + [BA]
n = len(gt) + len(tum) + 2; th = lambda i: np.pi / 2 - 2 * np.pi * (i + 0.5) / n; wdt = 2 * np.pi / n * 0.82
ax = fig.add_subplot(gs[0:2, 2:4], projection='polar'); ax.set_theta_zero_location('E'); ax.set_ylim(0, 1.55); ax.axis('off')
RINGS = [(1.10, 0.36, 'adjacent coupling', OI['blue']), (0.80, 0.26, 'domain scale, raw → GC-corrected', OI['orange']),
         (0.56, 0.20, 'shared-eQTL effect', OI['green']), (0.36, 0.17, 'same-TAD effect', OI['purple'])]
for r0, h, _, _c in RINGS: ax.bar(np.linspace(0, 2 * np.pi, 200), h, width=2 * np.pi / 200, bottom=r0, color='0.96', lw=0)
v1 = gt.cis_L1_minusGC_tech.values; dr, dc = gt.domain_L10_30_raw.values, gt.domain_minusGC_tech.values
eq = EQ.reindex(gt.tissue).delta_same_adj.values; td = gt.k.map(lambda k: TB.b_same_tad.get(k, np.nan)).values
mx = [v1.max(), dr.max(), np.nanmax(eq), np.nanmax(td)]
for i in range(len(gt)):
    a = th(i)
    ax.bar(a, RINGS[0][1] * v1[i] / mx[0], width=wdt, bottom=RINGS[0][0], color=OI['orange'] if gt.tissue[i].startswith('cells') else OI['blue'], lw=0)
    ax.bar(a, RINGS[1][1] * dr[i] / mx[1], width=wdt, bottom=RINGS[1][0], color='#F3C98B', lw=0)
    ax.bar(a, RINGS[1][1] * max(dc[i], 0) / mx[1], width=wdt * 0.55, bottom=RINGS[1][0], color=OI['orange'], lw=0)
    if np.isfinite(eq[i]): ax.bar(a, RINGS[2][1] * eq[i] / mx[2], width=wdt, bottom=RINGS[2][0], color=OI['green'], lw=0)
    if np.isfinite(td[i]): ax.bar(a, RINGS[3][1] * td[i] / mx[3], width=wdt, bottom=RINGS[3][0], color=OI['purple'], lw=0)
    deg = np.degrees(a) % 360; flip = 90 < deg < 270
    ax.text(a, 1.50, tissue_label(gt.tissue[i]), rotation=deg + (180 if flip else 0), rotation_mode='anchor', ha='right' if flip else 'left', va='center', fontsize=3.7)
tmx = np.nanmax(tumv)
for j, (c, v) in enumerate(zip(tum, tumv)):
    a = th(len(gt) + 2 + j); ax.bar(a, RINGS[0][1] * v / tmx, width=wdt, bottom=RINGS[0][0], color=OI['red'], lw=0)
    if c in TT.index: ax.bar(a, RINGS[3][1] * min((TT.within_wt[c] - TT.stable_wt[c]) / 0.12, 1), width=wdt, bottom=RINGS[3][0], color=OI['red'], lw=0, alpha=0.7)
    deg = np.degrees(a) % 360; flip = 90 < deg < 270
    ax.text(a, 1.50, c + ' (tumour)', rotation=deg + (180 if flip else 0), rotation_mode='anchor', ha='right' if flip else 'left', va='center', fontsize=3.9, color=OI['red'])
ax.text(0, 0, '36 GTEx tissues\n+ 5 tumour cohorts', ha='center', va='center', fontsize=4.6)
hs = [plt.Rectangle((0, 0), 1, 1, color=c) for *_, c in RINGS] + [plt.Rectangle((0, 0), 1, 1, color=OI['red'])]
fig.legend(hs, [f'outer: adjacent coupling (bar = 0 to {mx[0]:.2f})', f'2nd: domain scale, raw and GC-corrected (0 to {mx[1]:.3f})', f'3rd: shared-eQTL effect (0 to {mx[2]:.2f})', f'inner: same-TAD effect (0 to {mx[3]:.3f})',
                f'red: tumours, cis excess (0 to {tmx:.2f});\ninner ring, TAD contrast (0 to 0.12)'], loc='center left', bbox_to_anchor=(ax.get_position().x0 - 0.33, ax.get_position().y0 + 0.62 * ax.get_position().height),
           bbox_transform=fig.transFigure, fontsize=4.8, frameon=False, handlelength=0.8)
lab(ax, 'b')
ax = fig.add_subplot(gs[1, 0:2]); O['excess'] = O.mean_r - O.random_pair_floor
order = ['0-1kb', '1-5kb', '5-20kb', '20-100kb', '100-500kb', '>500kb']; mids = [0.5, 3, 12, 50, 250, 1000]
for o, col in [('divergent', OI['green']), ('tandem', OI['blue']), ('convergent', OI['red'])]:
    m = O[O.orientation == o].groupby('dist_bin').excess.median().reindex(order); ax.plot(mids, m.values, 'o-', color=col, ms=3, lw=1, label=o)
ax.set_xscale('log'); ax.set_xlabel('Intergenic distance (kb)'); ax.set_ylabel('Excess correlation'); ax.legend(loc='upper right'); lab(ax, 'c')
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
ax3.set_yticks(range(3)); ax3.set_yticklabels(['gene family / cluster', 'shared eQTL', 'same TAD'], fontsize=5.3); ax3.set_xticks([]); ax3.set_ylim(-0.6, 2.6)
for sp_ in ('top', 'right', 'bottom'): ax3.spines[sp_].set_visible(False)
ax3.invert_yaxis()
ax = fig.add_subplot(gs[2, 2]); sp = np.load(OUTDIR + 'spec_v2_thyroid.npy'); ax.loglog(sp[0], sp[1], 'o', ms=2.5, color='0.25', label='thyroid, 684 samples'); ax.loglog(sp[0], sp[2], color=OI['blue'], lw=1.2, label='two Lorentzians')
ax.loglog(sp[0], sp[3], color=OI['red'], lw=0.9, ls='--', label='power law'); ax.set_xlabel('Spatial frequency (cycles per gene)'); ax.set_ylabel('Covariance spectrum S(f)'); ax.legend(loc='lower left', fontsize=4.8)
lab(ax, 'e')
ax = fig.add_subplot(gs[2, 3]); xs = np.arange(len(Lz))
ax.scatter(xs - 0.12, Lz.lam_acf_short, s=14, color=OI['blue'], label='short, autocorrelation'); ax.scatter(xs + 0.12, Lz.lam_spec_short, s=14, marker='s', color=OI['sky'], label='short, spectrum')
ax.scatter(xs - 0.12, Lz.lam_acf_long, s=14, color=OI['red'], label='long, autocorrelation'); ax.set_yscale('log')
ax.set_xticks(xs); ax.set_xticklabels([tissue_label(t) for t in Lz.tissue], rotation=45, ha='right', fontsize=5.2); ax.set_ylabel('Decay length (genes)'); ax.set_ylim(0.6, 60); ax.legend(loc='upper left', fontsize=4.8, ncol=1)
lab(ax, 'f')
save(fig, 'Fig2_cis_layer'); print('Fig2_cis_layer ok')
