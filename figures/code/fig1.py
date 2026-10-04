import os
DATA = os.environ.get('POSLAYERS_DATA', 'data').rstrip('/') + '/'
import sys; sys.path.insert(0, os.path.dirname(__file__)); from style import *
import pyannotables as pa
exec(open(DATA + 'theory_sim.py').read().split('out = {}')[0])
R = pd.read_csv(A + 'atlas_results.csv'); P1 = pd.read_csv(DATA + 'p1_spectra.csv')
BM = pd.read_csv('/mnt/user-data/uploads/mart_export__1_.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene type': 'type'}).drop_duplicates('gid').set_index('gid')
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]; G = G[G.Chromosome.astype(str).isin([str(i) for i in range(1, 23)])].join(BM[['gc', 'type']], how='inner'); G = G[G.type == 'protein_coding'].sort_values(['Chromosome', 'Start'])
lags = np.arange(1, 61); ac = []
for L_ in lags:
    v = []
    for c, g in G.groupby('Chromosome'):
        x = (g.gc.values - g.gc.mean()) / g.gc.std()
        if len(x) > L_ + 5: v.append(np.mean(x[:-L_] * x[L_:]))
    ac.append(np.mean(v))
R = pd.read_csv(A + 'atlas_results.csv'); L = pd.read_csv(A + 'isochore_law.csv'); P4 = pd.read_csv(DATA + 'p4_domain_scale_tests.csv'); P5 = pd.read_csv(DATA + 'p5a_liver_gc.csv'); S = pd.read_csv(DATA + 'sim_decay_recovery.csv'); P1 = pd.read_csv(DATA + 'p1_spectra.csv')
import schem
fig = plt.figure(figsize=(W, W * 0.85)); gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.5)
# B identity in simulation
Xs, chrs, GC, mu_, z = simulate(seed=1, land_sd=4); c0 = chrs == 0; Xc = Xs[:, c0]; nn = Xc.shape[1]; mbar = Xc.mean(0); Dv = Xc - mbar
lhs = np.mean([np.abs(np.fft.fft(x)) ** 2 for x in Xc], 0); Cl = np.array([np.mean(np.sum(Dv * np.roll(Dv, -L, axis=1), 1)) for L in range(nn)])
land = np.abs(np.fft.fft(mbar)) ** 2; cov = np.real(np.fft.fft(Cl)); k = slice(1, nn // 2)
ax = fig.add_subplot(gs[0, 0:2]); schem.decomposition(ax); lab(ax, 'a', -0.04)
ax = fig.add_subplot(gs[0, 2]); ax.loglog(lhs[k], (land + cov)[k], '.', ms=1.4, color='0.35', rasterized=True); lim = [lhs[k].min(), lhs[k].max()]; ax.plot(lim, lim, color=OI['red'], lw=0.7)
ax.set_xlabel('Mean periodogram'); ax.set_ylabel('|M(f)|² + S(f)'); ax.text(0.95, 0.06, f'landscape {np.sum(land) / np.sum(lhs):.0%}', transform=ax.transAxes, fontsize=6, ha='right'); lab(ax, 'b')
# D universal peaks real vs random gene order
ax = fig.add_subplot(gs[1, 0]); ax.scatter(R.universal_peaks_random_order, R.universal_peaks_real_order, s=10, color=OI['blue'], lw=0)
m = max(R.universal_peaks_real_order.max(), R.universal_peaks_random_order.max()) * 1.05; ax.plot([0, m], [0, m], 'k:', lw=0.7)
ax.set_xlabel('Peaks, random gene order'); ax.set_ylabel('Peaks, real gene order'); ax.text(0.05, 0.9, f'median ratio {np.median(R.universal_peaks_real_order / R.universal_peaks_random_order):.2f}', transform=ax.transAxes, fontsize=6); lab(ax, 'c')
# E liver: consensus spectrum vs landscape spectrum
ax = fig.add_subplot(gs[1, 1]); ax.loglog(P1.w_static, P1.cons_full, '.', ms=1, color='0.45', rasterized=True); u = P1.universal.fillna(False).astype(bool)
ax.loglog(P1.w_static[u], P1.cons_full[u], '.', ms=3, color=OI['red'], label='186 universal peaks')
ax.set_xlabel('Landscape spectrum'); ax.set_ylabel('Consensus spectrum'); ax.legend(loc='upper left', bbox_to_anchor=(0, 0.92), markerscale=2)
ax.text(0.05, 0.9, f"r = {np.corrcoef(np.log(P1.w_static), np.log(P1.cons_full))[0, 1]:.3f}", transform=ax.transAxes, fontsize=6); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 2]); lbl = ['none', '−comp.', '−GC', '−5 PCs']
cis = [P4.cis_L1.iloc[0], P5.cis_L1.iloc[0], P5.cis_L1.iloc[1], P4.cis_L1.iloc[3]]; dom = [P4.domain_L10_30.iloc[0], P5.domain_L10_30.iloc[0], P5.domain_L10_30.iloc[1], P4.domain_L10_30.iloc[3]]
x = np.arange(4); ax.bar(x - 0.2, cis, 0.4, color=OI['blue'], label='adjacent genes (cis)'); ax.bar(x + 0.2, dom, 0.4, color=OI['red'], label='10–30 genes (domain)')
ax.set_xticks(x); ax.set_xticklabels(lbl, fontsize=5.6, rotation=25); ax.set_ylabel('Correlation vs permuted order'); ax.set_ylim(0, 0.34); ax.legend(loc='upper right', fontsize=5.4); lab(ax, 'e')
ax = fig.add_subplot(gs[2, 0])
for _, r in R.iterrows(): ax.plot([0, 1], [r.domain_L10_30_raw, r.domain_minusGC], color='0.75', lw=0.5)
ax.plot([0, 1], [R.domain_L10_30_raw.median(), R.domain_minusGC.median()], 'o-', color=OI['red'], lw=1.4, ms=4, label='median, 36 tissues')
ax.set_xticks([0, 1]); ax.set_xticklabels(['Naive', 'GC-corrected']); ax.set_xlim(-0.3, 1.3); ax.set_ylabel('Domain correlation (10–30 genes)'); ax.legend(loc='lower left', fontsize=5.4)
lab(ax, 'f')
ax = fig.add_subplot(gs[2, 1])
for lg, col in [(1, OI['blue']), (5, OI['green']), (20, OI['red'])]:
    ax.scatter(L[f'pred_L{lg}'], L[f'obs_gc_component_L{lg}'], s=9, color=col, lw=0, label=f'{lg} gene' + ('s' if lg > 1 else '') + f'  (r = {np.corrcoef(L[f"pred_L{lg}"], L[f"obs_gc_component_L{lg}"])[0, 1]:.2f})')
m = max(L[[c for c in L.columns if c.startswith('obs')]].max().max(), L[[c for c in L.columns if c.startswith('pred')]].max().max()) * 1.05
ax.plot([0, m], [0, m], 'k:', lw=0.7); ax.set_xlabel('Predicted GC component'); ax.set_ylabel('Observed GC component'); ax.legend(loc='upper left', fontsize=5.5)
lab(ax, 'g')
ax = fig.add_subplot(gs[2, 2]); g = S.groupby('true_lambda')[['naive', 'corrected']].median()
ax.plot(g.index, g.naive, 'o-', color=OI['red'], ms=3, lw=1, label='naive'); ax.plot(g.index, g.corrected, 's-', color=OI['blue'], ms=3, lw=1, label='GC-corrected')
ax.plot([0.8, 9], [0.8, 9], 'k:', lw=0.7, label='identity'); ax.set_xlabel('True cis length (genes)'); ax.set_ylabel('Estimated length (genes)'); ax.legend(loc='center left', bbox_to_anchor=(0.0, 0.52))
lab(ax, 'h')
save(fig, 'Fig1_nonpositional_layers'); print('Fig1_nonpositional_layers ok')
