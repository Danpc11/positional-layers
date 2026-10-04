import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from style import A, DATA, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys, tissue_label
# ---- Extended Data Fig. 1: robustness of the cis layer
P4 = pd.read_csv(OUTDIR + 'p4_domain_scale_tests.csv'); Lz = pd.read_csv(OUTDIR + 'lorentz_v2.csv'); E = pd.read_csv(A + 'eqtl_cis_test.csv'); F = pd.read_csv(A + 'robust_families.csv')
import pyannotables as pa
BMg = pd.read_csv(DATA + 'biomart_GRCh38_gene_gc.txt', sep='\t', low_memory=False).rename(columns={'Gene stable ID': 'gid', 'Gene % GC content': 'gc', 'Gene type': 'type'}).drop_duplicates('gid').set_index('gid')
Gg = pa.tables()['homo_sapiens-GRCh38-ensembl100']; Gg = Gg[~Gg.index.duplicated()]; Gg = Gg[Gg.Chromosome.astype(str).isin([str(i) for i in range(1, 23)])].join(BMg[['gc', 'type']], how='inner'); Gg = Gg[Gg.type == 'protein_coding'].sort_values(['Chromosome', 'Start'])
lags_ = np.arange(1, 61); ac = [np.mean([np.mean(((g.gc.values - g.gc.mean()) / g.gc.std())[:-L] * ((g.gc.values - g.gc.mean()) / g.gc.std())[L:]) for c, g in Gg.groupby('Chromosome') if len(g) > L + 5]) for L in lags_]
fig, axs = plt.subplots(1, 4, figsize=(W, W * 0.33)); plt.subplots_adjust(wspace=0.75)
ax = axs[0]; x = np.arange(len(P4)); ax.bar(x - 0.2, P4.cis_L1, 0.4, color=OI['blue'], label='adjacent (cis)'); ax.bar(x + 0.2, P4.domain_L10_30, 0.4, color=OI['red'], label='10–30 genes (domain)')
ax.set_xticks(x); ax.set_xticklabels(['raw', '−comp.', '−props', '−5 PCs', '−10 PCs', '−20 PCs'], fontsize=5.4, rotation=45, ha='right'); ax.set_ylabel('Correlation vs permuted order'); ax.set_ylim(0, 0.3); ax.legend(fontsize=5.2, loc='upper right'); lab(ax, 'a')
ax = axs[1]; dA = Lz.AIC_powerlaw - Lz.AIC_lor2; dB = Lz.AIC_lor1 - Lz.AIC_lor2; ax.scatter(dB, dA, s=16, color=OI['purple'])
for t, b, a in zip(Lz.tissue, dB, dA): ax.annotate(tissue_label(t), (b, a), fontsize=5.0, xytext=(3, 2), textcoords='offset points', ha='left' if b < 40 else 'right')
ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5); ax.set_xlabel('ΔAIC, 1 vs 2 Lorentzians'); ax.set_ylabel('ΔAIC, power law vs 2 Lorentzians'); lab(ax, 'b')
ax = axs[2]; ax.scatter(E.both_eGenes, 100 * E.share_of_cis_from_shared_eQTL, s=12, color=OI['orange']); ax.set_xlabel('Pairs with both genes eGenes'); ax.set_ylabel('Cis explained by shared eQTLs (%)')
lab(ax, 'c')
ax = axs[3]; ax.plot(lags_, ac, color=OI['red'], lw=1.1); ax.axhline(0, color='0.8', lw=0.5); ax.set_xlabel('Distance (genes)'); ax.set_ylabel('GC autocorrelation'); ax.text(0.95, 0.9, f'{ac[0]:.2f} at 1 gene\n{ac[9]:.2f} at 10\n{ac[29]:.2f} at 30', transform=ax.transAxes, ha='right', va='top', fontsize=5.6); lab(ax, 'd')
save(fig, 'ExtData_Fig1_robustness')
# ---- Extended Data Fig. 2: tumours and gene dosage
AN = pd.read_csv(OUTDIR + 'aneuploidy_continuousCN.csv'); FR2 = pd.read_csv(OUTDIR + 'friction2_tumour_tad_clustering.csv')
fig, axs = plt.subplots(1, 3, figsize=(W, W * 0.3)); plt.subplots_adjust(wspace=0.5)
p = AN.pivot(index='cohort', columns='expression', values='mean_far'); x = np.arange(len(p)); ax = axs[0]
ax.bar(x - 0.2, p['raw'], 0.4, color=OI['red'], label='uncorrected'); ax.bar(x + 0.2, p['continuous CN corrected'], 0.4, color=OI['blue'], label='copy-number corrected')
ax.set_xticks(x); ax.set_xticklabels(p.index, fontsize=6); ax.set_ylabel('Long-range floor (20–30 genes)'); ax.set_ylim(0, 0.14); ax.legend(fontsize=5.3, loc='upper center', ncol=2); lab(ax, 'a')
q = AN.pivot(index='cohort', columns='expression', values='rho_far_cna'); ax = axs[1]
ax.bar(x - 0.2, q['raw'], 0.4, color=OI['red']); ax.bar(x + 0.2, q['continuous CN corrected'], 0.4, color=OI['blue']); ax.axhline(0, color='k', lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(q.index, fontsize=6); ax.set_ylabel('Spearman ρ, floor vs CNA burden'); lab(ax, 'b')
ax = axs[2]; x = np.arange(len(FR2)); w = 0.2
ax.bar(x - 1.5 * w, FR2.z_normal, w, color='0.7', label='normal'); ax.bar(x - 0.5 * w, FR2.z_tumour, w, color=OI['red'], label='tumour')
ax.bar(x + 0.5 * w, FR2.z_normal_same_exclusion, w, color='0.7', hatch='///', label='normal, CNA genes removed'); ax.bar(x + 1.5 * w, FR2.z_tumour_CNA_genes_removed, w, color=OI['red'], hatch='///', label='tumour, CNA genes removed')
ax.set_xticks(x); ax.set_xticklabels(FR2.cohort, fontsize=6); ax.set_ylabel('TAD clustering of active genes (z)'); ax.legend(fontsize=5, loc='upper right'); lab(ax, 'c')
save(fig, 'ExtData_Fig2_dosage')
# ---- Extended Data Fig. 3: TAD effect across distance, with and without contact
FR1 = pd.read_csv(OUTDIR + 'friction1_tad_vs_contact.csv'); order = ['(25000.0, 50000.0]', '(50000.0, 100000.0]', '(100000.0, 200000.0]', '(200000.0, 500000.0]', '(500000.0, 1000000.0]', '(1000000.0, 2000000.0]']
labels = ['25–50 kb', '50–100 kb', '100–200 kb', '0.2–0.5 Mb', '0.5–1 Mb', '1–2 Mb']
fig, axs = plt.subplots(1, 2, figsize=(W * 0.7, W * 0.3)); plt.subplots_adjust(wspace=0.45)
for cell, col in [('GM12878', OI['blue']), ('IMR90', OI['red'])]:
    d = FR1[FR1.cell == cell].set_index('distance').reindex(order)
    axs[0].plot(range(6), d.r_same / d.r_diff, 'o-', color=col, ms=3.5, lw=1, label=cell)
    axs[1].plot(range(6), d.TAD_effect, 'o-', color=col, ms=3.5, lw=1, label=f'{cell}, distance only'); axs[1].plot(range(6), d.TAD_effect_given_contact, 's--', color=col, ms=3.5, lw=1, mfc='white', label=f'{cell}, + contact')
for ax in axs: ax.set_xticks(range(6)); ax.set_xticklabels(labels, rotation=30, fontsize=5.8)
axs[0].axhline(1.2, color='0.6', lw=0.8, ls=':'); axs[0].text(5.0, 1.0, 'constant ×1.2 (ref. 8)', fontsize=5.3, color='0.4', ha='right', va='top'); axs[0].set_ylim(0.8, None); axs[0].set_ylabel('Coupling, same TAD / different TAD'); axs[0].legend(fontsize=5.6); lab(axs[0], 'a')
axs[1].set_ylabel('Same-TAD effect on coupling'); axs[1].legend(fontsize=5.2); lab(axs[1], 'b')
save(fig, 'ExtData_Fig3_TAD_vs_contact'); print('ED ok')
