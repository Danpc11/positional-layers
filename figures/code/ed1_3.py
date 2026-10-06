import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import A, OI, OUTDIR, W, lab, save, tissue_label
import matplotlib.pyplot as plt
# ---- Extended Data Fig. 1: robustness of the cis layer
P4 = pd.read_csv(OUTDIR + 'p4_domain_scale_tests.csv'); Lz = pd.read_csv(OUTDIR + 'lorentz_v2.csv'); E = pd.read_csv(A + 'eqtl_cis_test.csv'); F = pd.read_csv(A + 'robust_families.csv')
GA = pd.read_csv(OUTDIR + 'gc_autocorrelation.csv'); lags_ = GA.lag.values; ac = GA.autocorrelation.tolist()   # scripts/gc_autocorrelation.py
fig = plt.figure(figsize=(W, W * 0.7)); gs = fig.add_gridspec(2, 4, wspace=0.75, hspace=0.6, height_ratios=[1, 1.05])
axs = [fig.add_subplot(gs[0, k]) for k in range(4)]
ax = axs[0]; x = np.arange(len(P4)); ax.bar(x - 0.2, P4.cis_L1, 0.4, color=OI['blue'], label='adjacent (cis)'); ax.bar(x + 0.2, P4.domain_L10_30, 0.4, color=OI['red'], label='10–30 genes (domain)')
ax.set_xticks(x); ax.set_xticklabels(['raw', '−comp.', '−props', '−5 PCs', '−10 PCs', '−20 PCs'], fontsize=5.5, rotation=45, ha='right'); ax.set_ylabel('Correlation vs permuted order'); ax.set_ylim(0, 0.3); lg = ax.legend(fontsize=5.5, loc='upper right', bbox_to_anchor=(1.16, 1.03), handlelength=1.0, handletextpad=0.4, borderaxespad=0); lg.keep_position = True; lab(ax, 'a')
ax = axs[1]; dA = Lz.AIC_powerlaw - Lz.AIC_lor2; dB = Lz.AIC_lor1 - Lz.AIC_lor2; ax.scatter(dB, dA, s=16, color=OI['purple'])
pts = sorted(zip(Lz.tissue, dB, dA), key=lambda z: z[2]); dy = {}
for i, (t_, b_, a_) in enumerate(pts):                    # separate labels of points that lie close together
    close = [p for p in pts[:i] if abs(p[1] - b_) < 12 and abs(p[2] - a_) < 8]
    dy[t_] = (8 if close else 0) if close or any(abs(p[1] - b_) < 12 and abs(p[2] - a_) < 8 for p in pts[i + 1:]) else 2
for t_, b_, a_ in zip(Lz.tissue, dB, dA): ax.annotate(tissue_label(t_), (b_, a_), fontsize=5.5, xytext=(3 if b_ < 40 else -3, dy[t_]), textcoords='offset points', ha='left' if b_ < 40 else 'right', va='center')
ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5); ax.set_xlabel('ΔAIC, 1 vs 2 Lorentzians'); ax.set_ylabel('ΔAIC, power law vs 2 Lorentzians'); lab(ax, 'b')
ax = axs[2]; ax.scatter(E.both_eGenes, 100 * E.share_of_cis_from_shared_eQTL, s=12, color=OI['orange']); ax.set_xlabel('Pairs with both genes eGenes'); ax.set_ylabel('Cis explained by shared eQTLs (%)')
lab(ax, 'c')
ax = axs[3]; ax.plot(lags_, ac, color=OI['red'], lw=1.1); ax.axhline(0, color='0.8', lw=0.5); ax.set_xlabel('Distance (genes)'); ax.set_ylabel('GC autocorrelation'); ax.text(0.95, 0.9, f'{ac[0]:.2f} at 1 gene\n{ac[9]:.2f} at 10\n{ac[29]:.2f} at 30', transform=ax.transAxes, ha='right', va='top', fontsize=5.6); lab(ax, 'd')
# e: coupling against physical distance, with the two-exponential fit (scripts/spectral_rigour.py)
SR = pd.read_csv(OUTDIR + 'spectral_rigour.csv').set_index('tissue')
T6 = ['thyroid', 'nerve_tibial', 'skin_sun_exposed_lower_leg', 'cells_cultured_fibroblasts', 'cells_ebv-transformed_lymphocytes', 'muscle_skeletal']
C6 = [OI['blue'], OI['green'], OI['orange'], OI['red'], OI['purple'], OI['sky']]
ax = fig.add_subplot(gs[1, 0:2]); xx = np.logspace(np.log10(5e3), np.log10(5e6), 200)
for t, col in zip(T6, C6):
    d = pd.read_csv(OUTDIR + f'coupling_by_bp_{t}.csv'); r = SR.loc[t]
    ax.scatter(d.distance_bp / 1e6, d.coupling_minus_far, s=7, color=col, zorder=3, label=tissue_label(t))
    ax.plot(xx / 1e6, r.bp_short_amplitude * np.exp(-xx / (r.bp_short_length_kb * 1e3)) + r.bp_long_amplitude * np.exp(-xx / (r.bp_long_length_Mb * 1e6)), color=col, lw=0.8)
ax.set_xscale('log'); ax.axhline(0, color='0.7', lw=0.5); ax.set_xlabel('Distance between genes (Mb)'); ax.set_ylabel('Coupling minus that of pairs\n20–40 Mb apart')
ax.text(0.98, 0.62, f"short scale {SR.bp_short_length_kb.min():.0f}–{SR.bp_short_length_kb.max():.0f} kb\nlong scale {SR.bp_long_length_Mb.min():.1f}–{SR.bp_long_length_Mb.max():.1f} Mb", transform=ax.transAxes, ha='right', va='top', fontsize=5.5)
lg = ax.legend(loc='upper right', fontsize=5.5, ncol=2, handletextpad=0.1, columnspacing=0.6, markerscale=1.2); lg.keep_position = True; lab(ax, 'e')
# f: spectrum of the mean profile against autocorrelation-preserving nulls (scripts/spectral_rigour.py)
L = pd.read_csv(OUTDIR + 'landscape_null_thyroid_chr1.csv')
ax = fig.add_subplot(gs[1, 2:4])
ax.plot(L.frequency, L.landscape_power, color='0.35', lw=0.5, label='mean profile, thyroid chr. 1')
ax.plot(L.frequency, L.red_noise_threshold, color=OI['red'], lw=1.0, label='red-noise threshold')
ax.axhline(L.block30_threshold.iloc[0], color=OI['blue'], lw=1.0, ls='--', label='block-permutation threshold')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Spatial frequency (cycles per gene)'); ax.set_ylabel('Power')
ax.set_ylim(top=max(L.red_noise_threshold.max(), L.block30_threshold.iloc[0]) * 60)   # room for the legend above the thresholds
ax.text(0.03, 0.04, f'frequencies above either threshold: {int(((L.landscape_power > L.red_noise_threshold) | (L.landscape_power > L.block30_threshold)).sum())} of {len(L):,}', transform=ax.transAxes, fontsize=5.5)
lg = ax.legend(loc='upper right', fontsize=5.5); lg.keep_position = True; lab(ax, 'f')
save(fig, 'ExtData_Fig1_robustness')
# ---- Extended Data Fig. 2: tumours and gene dosage
AN = pd.read_csv(OUTDIR + 'aneuploidy_continuousCN.csv'); FR2 = pd.read_csv(OUTDIR + 'friction2_tumour_tad_clustering.csv')
fig, axs = plt.subplots(2, 2, figsize=(W, W * 0.62)); plt.subplots_adjust(wspace=0.35, hspace=0.75); axs = axs.ravel()
p = AN.pivot(index='cohort', columns='expression', values='mean_far'); x = np.arange(len(p)); ax = axs[0]
ax.bar(x - 0.2, p['raw'], 0.4, color=OI['red'], label='uncorrected'); ax.bar(x + 0.2, p['continuous CN corrected'], 0.4, color=OI['blue'], label='copy-number corrected')
ax.set_xticks(x); ax.set_xticklabels(p.index, fontsize=6); ax.set_ylabel('Long-range floor (20–30 genes)'); ax.set_ylim(0, 0.14); lg = ax.legend(fontsize=5.5, loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=1, borderaxespad=0); lg.keep_position = True; lab(ax, 'a')
q = AN.pivot(index='cohort', columns='expression', values='rho_far_cna'); ax = axs[1]
ax.bar(x - 0.2, q['raw'], 0.4, color=OI['red']); ax.bar(x + 0.2, q['continuous CN corrected'], 0.4, color=OI['blue']); ax.axhline(0, color='k', lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(q.index, fontsize=6); ax.set_ylabel('Spearman ρ, floor vs CNA burden'); lab(ax, 'b')
ax = axs[2]; x = np.arange(len(FR2)); w = 0.2
ax.bar(x - 1.5 * w, FR2.z_normal, w, color='0.7', label='normal'); ax.bar(x - 0.5 * w, FR2.z_tumour, w, color=OI['red'], label='tumour')
ax.bar(x + 0.5 * w, FR2.z_normal_same_exclusion, w, color='0.7', hatch='////', edgecolor='white', linewidth=0, label='normal, CNA genes removed'); ax.bar(x + 1.5 * w, FR2.z_tumour_CNA_genes_removed, w, color=OI['red'], hatch='////', edgecolor='white', linewidth=0, label='tumour, CNA genes removed')   # hatch drawn in the edge colour
ax.set_xticks(x); ax.set_xticklabels(FR2.cohort, fontsize=6); ax.set_ylabel('TAD clustering of active genes (z)'); ax.legend(fontsize=5.5, loc='upper right'); lab(ax, 'c')
# d: controls for the copy-number layer (scripts/dosage_controls.py)
DC = pd.read_csv(OUTDIR + 'dosage_controls.csv'); order = ['true', 'permuted across tumours', 'segments shifted', 'purity-adjusted']
names = ['true copy\nnumber', 'permuted\nacross tumours', 'segments\nshifted', 'purity-\nadjusted']; ax = axs[3]; x = np.arange(len(order)); w = 0.36
for k, (coh, col, lbl) in enumerate([('BLCA', OI['red'], 'bladder'), ('UCEC', OI['orange'], 'endometrium')]):
    v = DC[DC.cohort == coh].set_index('copy_number').loc[order].observed_relative_to_true * 100
    ax.bar(x + (k - 0.5) * w, v.values, w, color=col, label=lbl)
ax.axhline(100, color='0.6', lw=0.6, ls=':'); ax.axhline(0, color='k', lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(names, fontsize=5.5); ax.set_ylabel('Dosage covariance\n(% of the value with true copy number)'); ax.set_ylim(-10, 125)
lg = ax.legend(loc='center', bbox_to_anchor=(0.5, 0.62), fontsize=5.5); lg.keep_position = True; lab(ax, 'd')
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
axs[0].axhline(1.2, color='0.6', lw=0.8, ls=':'); axs[0].text(5.0, 1.0, 'constant ×1.2 (ref. 8)', fontsize=5.5, color='0.4', ha='right', va='top'); axs[0].set_ylim(0.8, None); axs[0].set_ylabel('Coupling, same TAD / different TAD'); axs[0].legend(fontsize=5.6); lab(axs[0], 'a')
axs[1].set_ylabel('Same-TAD effect on coupling'); axs[1].legend(fontsize=5.5); lab(axs[1], 'b')
save(fig, 'ExtData_Fig3_TAD_vs_contact'); print('ED ok')
