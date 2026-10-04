import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from style import A, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys, tissue_label
T = pd.read_csv(A + 'tad_cis_test.csv'); TT = pd.read_csv(OUTDIR + 'tad_tumour_results.csv'); CO = pd.read_csv(OUTDIR + 'cohesin_freedman_lane.csv'); MT = pd.read_csv(OUTDIR + 'stag2_meta.csv').iloc[0]; TB = pd.read_csv(OUTDIR + 'tad_block_bootstrap.csv').set_index('tissue')
BA = pd.read_csv(OUTDIR + 'beataml_freedman_lane.csv'); AN = pd.read_csv(OUTDIR + 'aneuploidy_continuousCN.csv')
H = {c: pd.read_csv(A + f'hic_coupling_{c}.csv.gz') for c in ['GM12878', 'IMR90']}
import schem
fig = plt.figure(figsize=(W, W * 0.62)); gs = fig.add_gridspec(2, 4, hspace=0.75, wspace=0.9)

ax = fig.add_subplot(gs[0, 0]); schem.saturation(ax); lab(ax, 'a', -0.25)
ax = fig.add_subplot(gs[0, 1]); s = T.sort_values('b_same_tad'); y = np.arange(len(s))
tb = TB.reindex(s.tissue); ax.errorbar(s.b_same_tad, y + 0.15, xerr=[s.b_same_tad.values - tb.ci_low.values, tb.ci_high.values - s.b_same_tad.values], fmt='o', ms=3, color=OI['blue'], lw=0.6, capsize=0, label='all pairs (95% block CI)'); ax.scatter(s.b_same_tad_no_shared_eQTL, y - 0.15, s=12, color=OI['sky'], marker='s', label='no shared eQTL')
ax.axvline(0, color='k', lw=0.6); ax.set_yticks(y); ax.set_yticklabels([tissue_label(t) for t in s.tissue], fontsize=5.3); ax.set_xlabel('Same-TAD increase in coupling\n(at equal distance)'); ax.set_xlim(-0.01, 0.11); ax.legend(loc='lower left', bbox_to_anchor=(-0.05, 1.0), fontsize=5.0, ncol=1, handletextpad=0.2, borderaxespad=0)
lab(ax, 'b', -0.6)
bins = [25e3, 50e3, 1e5, 2e5, 5e5, 1e6, 2e6]; mids = [37.5, 75, 150, 350, 750, 1500]
ax = fig.add_subplot(gs[0, 2]); kk = {}
for c, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    d = H[c].copy(); d['db'] = pd.cut(d.tss_distance, bins); lo, hi, mr, mc = [], [], [], []
    for b, g in d.groupby('db', observed=True):
        q = pd.qcut(g.oe.rank(method='first'), 5, labels=False); lo.append(g.r[q == 0].mean()); hi.append(g.r[q == 4].mean()); mr.append(g.r.mean()); mc.append(g.contact_KR.mean())
    ax.plot(mids, hi, 'o-', color=col, ms=3, lw=1.1, label=f'{c}, most contact'); ax.plot(mids, lo, 'o--', color=col, ms=3, lw=0.9, mfc='white', label=f'{c}, least contact'); kk[c] = (mc, mr)
ax.set_xscale('log'); ax.set_xlabel('Distance between promoters (kb)'); ax.set_ylabel('Coupling'); ax.set_ylim(-0.005, 0.16); ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=1, fontsize=5); lab(ax, 'c')
ax = fig.add_subplot(gs[0, 3])
for c, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    mc, mr = kk[c]; k = np.polyfit(np.log10(mc), np.log10(np.clip(mr, 1e-4, None)), 1)[0]; ax.loglog(mc, mr, 'o-', color=col, ms=3.5, lw=1, label=f'{c}: k = {k:.2f}')
ax.set_xlabel('Mean Hi-C contact (KR)'); ax.set_ylabel('Mean coupling'); ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), fontsize=5.0, borderaxespad=0); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 0]); r = TT.drop_duplicates('cohort')
x = np.arange(len(r)); ax.bar(x - 0.2, r.within_wt, 0.4, color=OI['green'], label='same TAD'); ax.bar(x + 0.2, r.stable_wt, 0.4, color='0.6', label='across stable boundary')
ax.set_xticks(x); ax.set_xticklabels(['Bladder', 'Endometrium'], fontsize=6.3); ax.set_ylabel('Cis excess (adjacent genes)'); ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), fontsize=5.2, ncol=1); ax.set_ylim(0, 0.22)
lab(ax, 'e', -0.3)
ax = fig.add_subplot(gs[1, 1]); BS = pd.read_csv(OUTDIR + 'boot_hic_summary.csv').set_index(['cell', 'stat']); xk = np.arange(3)
for cell, col, off in [('GM12878', OI['blue'], -0.15), ('IMR90', OI['red'], 0.15)]:
    d = BS.loc[cell].loc[['k_t1', 'k_t2', 'k_t3']]; ax.errorbar(xk + off, d.estimate, yerr=[d.estimate - d.ci_low, d.ci_high - d.estimate], fmt='o', color=col, ms=4, capsize=2, label=cell)
ax.set_xticks(xk); ax.set_xticklabels(['low', 'mid', 'high']); ax.set_xlabel('Expression of the pair (tertile)'); ax.set_ylabel('Contact exponent k'); ax.set_ylim(0, 1.05); ax.legend(loc='lower left', fontsize=5.4)
lab(ax, 'f')
ax = fig.add_subplot(gs[1, 2:4]); rows = []
for _, q in CO.iterrows(): rows.append((f"{'STAG2' if q.gene == 'STAG2' else 'CTCF'} {'BLCA' if q.cohort == 'BLCA' else 'UCEC'} {'trunc' if q['class'] == 'truncating' else 'any'}{', strict' if q.tmb_q == 0.7 else ''}", q.adj_pct, q.ci_low, q.ci_high, OI['red'] if q.gene == 'STAG2' else '0.5'))
for _, q in BA.iterrows(): rows.append((f"{q.group.replace(', ', ' ').replace('any coding', 'any').replace('truncating', 'trunc')} AML", q.adj_pct, q.ci_low, q.ci_high, OI['orange']))
rows.append(('STAG2 BLCA + AML, meta', MT.meta_pct, MT.ci_low, MT.ci_high, 'k'))
for i, (n, e, l, h, col) in enumerate(rows[::-1]): ax.errorbar(e, i, xerr=[[e - l], [h - e]], fmt='o', color=col, ms=3, capsize=1.5, lw=0.8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([n for n, *_ in rows[::-1]], fontsize=5.2); ax.yaxis.tick_right(); ax.spines['right'].set_visible(True); ax.spines['left'].set_visible(False); ax.axvline(0, color='k', lw=0.6); ax.set_xlabel('Change in cis excess (%), adjusted')
lab(ax, 'g', -0.12)
save(fig, 'Fig4_architecture_cohesin'); print('Fig4_architecture_cohesin ok')
