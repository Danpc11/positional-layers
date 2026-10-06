import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import A, OI, OUTDIR, W, lab, save, tissue_label, placeholder
import matplotlib.pyplot as plt
T = pd.read_csv(A + 'tad_cis_test.csv'); TT = pd.read_csv(OUTDIR + 'tad_tumour_results.csv'); CO = pd.read_csv(OUTDIR + 'cohesin_freedman_lane.csv'); MT = pd.read_csv(OUTDIR + 'stag2_meta.csv').iloc[0]; TB = pd.read_csv(OUTDIR + 'tad_block_bootstrap.csv').set_index('tissue')
BA = pd.read_csv(OUTDIR + 'beataml_freedman_lane.csv'); AN = pd.read_csv(OUTDIR + 'aneuploidy_continuousCN.csv')
H = {c: pd.read_csv(A + f'hic_coupling_{c}.csv.gz') for c in ['GM12878', 'IMR90']}
fig = plt.figure(figsize=(W, W * 0.62)); gs = fig.add_gridspec(2, 4, hspace=0.75, wspace=0.9)

ax = fig.add_subplot(gs[0, 0]); placeholder(ax, 'Schematic: architecture', 'Domains, convergent CTCF loop anchors and an active intervening gene limit which neighbouring genes share a regulatory source.', key='fig4a'); p_ = ax.get_position(); ax.set_position([0.005, p_.y0 - 0.12 * p_.height, p_.x1 - 0.005 + 0.01, p_.height * 1.12]); lab(ax, 'a', -0.25)   # use the empty left margin
ax = fig.add_subplot(gs[0, 1]); s = T.sort_values('b_same_tad'); y = np.arange(len(s))
tb = TB.reindex(s.tissue); ax.errorbar(s.b_same_tad, y + 0.15, xerr=[s.b_same_tad.values - tb.ci_low.values, tb.ci_high.values - s.b_same_tad.values], fmt='o', ms=3, color=OI['blue'], lw=0.6, capsize=0, label='all pairs (95% block CI)'); ax.scatter(s.b_same_tad_no_shared_eQTL, y - 0.15, s=12, color=OI['sky'], marker='s', label='no shared eQTL')
ax.axvline(0, color='k', lw=0.6); ax.set_yticks(y); ax.set_yticklabels([tissue_label(t) for t in s.tissue], fontsize=5.5); ax.set_xlabel('Same-TAD increase in coupling\n(at equal distance)'); ax.set_xlim(-0.01, 0.11); ax.legend(loc='lower left', bbox_to_anchor=(-0.05, 1.0), fontsize=5.5, ncol=1, handletextpad=0.2, borderaxespad=0)
lab(ax, 'b', -0.6)
bins = [25e3, 50e3, 1e5, 2e5, 5e5, 1e6, 2e6]; mids = [37.5, 75, 150, 350, 750, 1500]
ax = fig.add_subplot(gs[0, 2]); kk = {}
for c, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    d = H[c].copy(); d['db'] = pd.cut(d.tss_distance, bins); lo, hi, mr, mc = [], [], [], []
    for b, g in d.groupby('db', observed=True):
        q = pd.qcut(g.oe.rank(method='first'), 5, labels=False); lo.append(g.r[q == 0].mean()); hi.append(g.r[q == 4].mean())
        gp = g[g.contact_KR > 0]; mr.append(gp.r.mean()); mc.append(gp.contact_KR.mean())   # same pairs (contact > 0) as boot_hic.py
    ax.plot(mids, hi, 'o-', color=col, ms=3, lw=1.1, label=f'{c}, most contact'); ax.plot(mids, lo, 'o--', color=col, ms=3, lw=0.9, mfc='white', label=f'{c}, least contact'); kk[c] = (mc, mr)
ax.set_xscale('log'); ax.set_xlabel('Distance between promoters (kb)'); ax.set_ylabel('Coupling'); ax.set_ylim(-0.005, 0.16); ax.legend(loc='upper right', fontsize=5.5, handlelength=1.6); lab(ax, 'c')
ax = fig.add_subplot(gs[0, 3])
for c, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    mc, mr = kk[c]; k = np.polyfit(np.log10(mc), np.log10(np.clip(mr, 1e-4, None)), 1)[0]; ax.loglog(mc, mr, 'o-', color=col, ms=3.5, lw=1, label=f'{c}: k = {k:.2f}')
ax.set_xlabel('Mean Hi-C contact (KR)'); ax.set_ylabel('Mean coupling'); ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), fontsize=5.5, borderaxespad=0); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 0])                          # ctcf_barrier.py
for cell, col, off in [('GM12878', OI['blue'], -0.08), ('IMR90', OI['red'], 0.08)]:
    d = pd.read_csv(OUTDIR + f'ctcf_barrier_means_{cell}.csv')   # ratio intervals resample numerator and denominator together
    ax.errorbar(d.crossed_loops + off, d.relative, yerr=[d.relative - d.relative_ci_low, d.relative_ci_high - d.relative], fmt='o-', color=col, ms=3.5, lw=1, capsize=1.5, label=cell)
ax.axhline(1, color='0.7', lw=0.5, ls=':'); ax.set_xticks(range(4)); ax.set_xticklabels(['0', '1', '2', '≥3']); ax.set_xlabel('Convergent CTCF loops crossed'); ax.set_ylabel('Coupling relative to\nno loop crossed'); ax.set_ylim(0.3, 1.3); ax.legend(loc='lower left', fontsize=5.5)
lab(ax, 'e', -0.3)
ax = fig.add_subplot(gs[1, 1])                          # active_gene_barrier.py
for tis, col, off, nm in [('cells_ebv-transformed_lymphocytes', OI['blue'], -0.1, 'lymphoblastoid'), ('thyroid', OI['green'], 0, 'thyroid'), ('cells_cultured_fibroblasts', OI['red'], 0.1, 'fibroblasts')]:
    d = pd.read_csv(OUTDIR + f'active_gene_barrier_means_{tis}.csv')
    ax.errorbar(d.quartile + off, d.relative, yerr=[d.relative - d.relative_ci_low, d.relative_ci_high - d.relative], fmt='o-', color=col, ms=3.5, lw=1, capsize=1.5, label=nm)
ax.axhline(1, color='0.7', lw=0.5, ls=':'); ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(['1\nlowest', '2', '3', '4\nhighest']); ax.set_xlabel('Expression of the intervening gene (quartile)'); ax.set_ylabel('Coupling of flanking genes\nrelative to quartile 1'); ax.set_ylim(0.3, 1.3); ax.legend(loc='lower left', fontsize=5.5)
lab(ax, 'f')
ax = fig.add_subplot(gs[1, 2:4]); rows = []
for _, q in CO.iterrows(): rows.append((f"{'STAG2' if q.gene == 'STAG2' else 'CTCF'} {'BLCA' if q.cohort == 'BLCA' else 'UCEC'} {'trunc' if q['class'] == 'truncating' else 'any'}{', strict' if q.tmb_q == 0.7 else ''}", q.adj_pct, q.ci_low, q.ci_high, OI['red'] if q.gene == 'STAG2' else '0.5'))
for _, q in BA.iterrows(): rows.append((f"{q.group.replace(', ', ' ').replace('any coding', 'any').replace('truncating', 'trunc')} AML", q.adj_pct, q.ci_low, q.ci_high, OI['orange']))
rows.append(('STAG2 BLCA + AML, meta', MT.meta_pct, MT.ci_low, MT.ci_high, 'k'))
for i, (n, e, l, h, col) in enumerate(rows[::-1]): ax.errorbar(e, i, xerr=[[e - l], [h - e]], fmt='o', color=col, ms=3, capsize=1.5, lw=0.8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([n for n, *_ in rows[::-1]], fontsize=5.5); ax.yaxis.tick_right(); ax.spines['right'].set_visible(True); ax.spines['left'].set_visible(False); ax.axvline(0, color='k', lw=0.6); ax.set_xlabel('Change in cis excess (%), adjusted')
lab(ax, 'g', -0.12)
save(fig, 'Fig4_architecture_cohesin'); print('Fig4_architecture_cohesin ok')
