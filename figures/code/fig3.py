import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import A, OI, OUTDIR, W, lab, save, tissue_label, placeholder
import matplotlib.pyplot as plt
E = pd.read_csv(A + 'eqtl_cis_test.csv'); Dz = pd.read_csv(A + 'coloc_dose_response.csv'); SP = pd.read_csv(A + 'eqtl_tissue_specificity.csv'); Cc = pd.read_csv(A + 'coloc_cis_test.csv')
PC = pd.read_csv(OUTDIR + 'eqtl_law_v2_pairs.csv').rename(columns={'pred_r': 'pred_genetic_r', 'obs_r_int_pc15': 'obs_r'})
rng = np.random.default_rng(0)
fig = plt.figure(figsize=(W, W * 1.12)); gs = fig.add_gridspec(3, 4, hspace=0.62, wspace=0.65, height_ratios=[1, 1, 1.55])

ax = fig.add_subplot(gs[0, 0:2]); placeholder(ax, 'Schematic: coupling as the sum of shared sources', 'Each source (regulatory element, genotype, copy-number segment) varies across\nindividuals and reaches nearby genes with a footprint. Covariance between two genes =\nsum over shared sources of source variance × footprint on gene i × footprint on gene j.'); lab(ax, 'a', -0.04)
DS = pd.read_csv(OUTDIR + 'dosage_prediction.csv')                     # dosage_prediction.py: out-of-sample prediction of the copy-number layer
ax = fig.add_subplot(gs[0, 2:4])
for (c, d), col in zip(DS.groupby('cohort'), [OI['blue'], OI['red']]):
    ax.plot(d.lag, d.predicted_dosage_cov, '-', color=col, lw=1.2, label=f'{dict(BLCA="bladder", UCEC="endometrium")[c]}, predicted')
    ax.plot(d.lag, d.observed_dosage_cov, 'o', color=col, ms=3.5, mfc='white', mew=1, label=f'{dict(BLCA="bladder", UCEC="endometrium")[c]}, observed')
ax.set_xscale('log'); ax.set_xlabel('Distance between genes (genes)'); ax.set_ylabel('Dosage covariance')
ax.text(0.97, 0.95, 'no fitted parameter;\nobserved/predicted 0.99–1.07', transform=ax.transAxes, ha='right', va='top', fontsize=5.6); ax.legend(loc='lower left', fontsize=5.2); lab(ax, 'b')
ax = fig.add_subplot(gs[1, 0]); cats = [('r_not_both_eGenes', 'not both\neGenes', '0.6'), ('r_eGenes_no_share', 'no shared\nvariant', OI['sky']), ('r_share_same', 'shared,\nsame sign', OI['blue'])]
for i, (c, n, col) in enumerate(cats): ax.scatter(i + rng.uniform(-0.15, 0.15, len(E)), E[c], s=6, color=col, lw=0); ax.hlines(E[c].median(), i - 0.25, i + 0.25, color='k', lw=1.2)
ax.set_xticks(range(3)); ax.set_xticklabels([n for _, n, _ in cats], fontsize=5.4); ax.set_ylabel('Adjacent-pair correlation'); lab(ax, 'c')
ax = fig.add_subplot(gs[1, 1]); s = E.sort_values('delta_same_adj'); y = np.arange(len(s))
ax.scatter(s.delta_same_adj, y, s=8, color=OI['blue'], label='same direction'); ax.scatter(s.delta_opposite_only_adj, y, s=8, color=OI['red'], label='opposite direction only')
ax.axvline(0, color='k', lw=0.6); ax.set_yticks([]); ax.set_ylabel('36 tissues'); ax.set_xlabel('Effect on coupling (adjusted)')
lab(ax, 'd', -0.08)
ax = fig.add_subplot(gs[1, 2]); order = ['0', '0-0.1', '0.1-0.5', '0.5-0.8', '>0.8']; g = Dz.groupby('p_coloc_bin').mean_r
med = [g.get_group(b).median() for b in order]; q1 = [g.get_group(b).quantile(0.25) for b in order]; q3 = [g.get_group(b).quantile(0.75) for b in order]
ax.errorbar(range(5), med, yerr=[np.subtract(med, q1), np.subtract(q3, med)], fmt='o-', color=OI['purple'], ms=4, capsize=2, lw=1.2)
ax.set_xticks(range(5)); ax.set_xticklabels(['none', '<0.1', '0.1–0.5', '0.5–0.8', '>0.8'], fontsize=5.4, rotation=30); ax.set_xlabel('Colocalisation score'); ax.set_ylabel('Adjacent-pair correlation')
lab(ax, 'e')
ax = fig.add_subplot(gs[1, 3]); PC['bin'] = pd.qcut(PC.pred_genetic_r, 8, duplicates='drop'); b = PC.groupby('bin', observed=True).agg(p=('pred_genetic_r', 'mean'), o=('obs_r', 'mean'), se=('obs_r', lambda x: x.std() / np.sqrt(len(x))))
ax.errorbar(b.p, b.o, yerr=1.96 * b.se, fmt='o', color=OI['blue'], ms=4, capsize=2); sl = np.polyfit(PC.pred_genetic_r, PC.obs_r, 1); xx = np.linspace(b.p.min(), b.p.max(), 10)
ax.set_xscale('symlog', linthresh=0.01); ax.set_xlabel('Predicted genetic covariance'); ax.set_ylabel('Observed coupling')
ES = pd.read_csv(OUTDIR + 'eqtl_law_v2_summary.csv'); es = ES[(ES.score_threshold == 0.1) & (ES.observed_scale == 'obs_r_int_pc15') & (ES.prediction == 'pred_r')].iloc[0]
ax.text(0.05, 0.9, f'r = {es.pearson_r:.2f}', transform=ax.transAxes, fontsize=5.6, va='top'); lab(ax, 'f')
# f: tissue x tissue heatmap of the effect of sharing a variant; the diagonal (same tissue) dominates each row
SP['diff'] = SP.b_own - SP.b_other
Mx = SP.pivot_table(index='coupling_tissue', columns='eqtl_tissue', values='diff', aggfunc='mean')
tis = Mx.mean(1).sort_values(ascending=False).index.tolist(); Mx = Mx.reindex(index=tis, columns=tis)
v = np.nanpercentile(np.abs(Mx.values), 98)
ax = fig.add_subplot(gs[2, 0:4]); im = ax.imshow(Mx.values, cmap='RdBu_r', vmin=-v, vmax=v, aspect='auto', interpolation='nearest')
ax.set_xticks(range(len(tis))); ax.set_xticklabels([tissue_label(t) for t in tis], rotation=90, fontsize=4.0)
ax.set_yticks(range(len(tis))); ax.set_yticklabels([tissue_label(t) for t in tis], fontsize=4.0)
ax.set_xlabel('Other tissue in which the variant is shared', fontsize=6); ax.set_ylabel('Tissue in which coupling is measured', fontsize=6)
ax.tick_params(length=1, pad=1)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02); cb.set_label('Own-tissue minus other-tissue effect', fontsize=5.2); cb.ax.tick_params(labelsize=5)
ax.text(1.0, 1.02, f'red: own tissue stronger ({100 * (SP["diff"] > 0).mean():.0f}% of 1,118 tissue pairs)', transform=ax.transAxes, ha='right', va='bottom', fontsize=5)
lab(ax, 'g')
save(fig, 'Fig3_genetic_component'); print('Fig3_genetic_component ok')
