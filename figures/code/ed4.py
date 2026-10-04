"""Extended Data Fig. 4: heatmaps summarising the atlas, the laws and the perturbation rules."""
import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import A, OI, OUTDIR, W, lab, np, os, pd, plt, save, sys, tissue_label
import schem
R = pd.read_csv(A + 'atlas_results.csv'); L = pd.read_csv(OUTDIR + 'isochore_v2.csv'); L = L.rename(columns={**{f'held_pred_L{k}': f'pred_L{k}' for k in (1, 2, 5, 10, 20, 30)}, **{f'held_obs_L{k}': f'obs_gc_component_L{k}' for k in (1, 2, 5, 10, 20, 30)}}); E = pd.read_csv(A + 'eqtl_cis_test.csv')
fig = plt.figure(figsize=(W, W * 0.95)); gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.55, height_ratios=[1, 1.5])

# a: per-tissue layer budget (tissues x layers)
R = R.sort_values('cis_L1_minusGC_tech', ascending=False).reset_index(drop=True)
cols = ['landscape_share', 'eta2_GCslope_batch', 'domain_L10_30_raw', 'domain_minusGC_tech', 'cis_L1_minusGC_tech', 'cis_decay_length_genes']
names = ['landscape\nshare', 'GC batch\nη²', 'domain,\nnaive', 'domain,\n−GC', 'cis,\nadjacent', 'decay\n(genes)']
M = R[cols].copy()
for c in cols: M[c] = (M[c] - M[c].min()) / (M[c].max() - M[c].min())
ax = fig.add_subplot(gs[0, :])
im = ax.imshow(M.T.values, aspect='auto', cmap='magma')
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=5.2)
ax.set_xticks(range(len(R))); ax.set_xticklabels([tissue_label(t) for t in R.tissue], fontsize=4.4, rotation=90)
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.ax.tick_params(labelsize=5); cb.set_label('scaled within row', fontsize=5.2)
for s in ax.spines.values(): s.set_visible(False)
lab(ax, 'a', -0.05)

# b: isochore law, predicted vs observed by distance (tissues x distances), ratio obs/pred
dists = [1, 2, 5, 10, 20, 30]
Lr = L.copy(); Lr['tissue'] = Lr.tissue.map(tissue_label)
Mb = np.column_stack([Lr[f'obs_gc_component_L{d}'] / Lr[f'pred_L{d}'] for d in dists])
order = np.argsort(np.nanmedian(Mb, 1))[::-1]
ax = fig.add_subplot(gs[1, 0])
im = ax.imshow(Mb[order], aspect='auto', cmap='RdBu_r', vmin=0.6, vmax=1.4)
ax.set_xticks(range(len(dists))); ax.set_xticklabels([f'{d}' for d in dists], fontsize=5.4); ax.set_xlabel('Distance (genes)', fontsize=6)
ax.set_yticks(range(len(Lr))); ax.set_yticklabels(Lr.tissue.values[order], fontsize=4.0)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02); cb.ax.tick_params(labelsize=5); cb.set_label('observed / predicted (held-out)', fontsize=5.2)
ax.set_title('Isochore law holds at every distance', fontsize=6, loc='left', pad=4)
for s in ax.spines.values(): s.set_visible(False)
lab(ax, 'b', -0.42)

# c: intervention rules as a rule table heatmap
rules = ['shared element\n(HBG1/2 edit)', 'shared enhancer\n(BET inhibitor)', 'shared variant\n(eQTL)', 'disease programme\n(fibrosis)', 'promoter silencing\n(CRISPRi)', 'protein, trans\n(BCL11A edit)']
obs = ['propagates with\ncoupling', 'distance only\n(<50 kb)', 'no cis effect']
G = np.array([[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
ax = fig.add_subplot(gs[1, 1])
im = ax.imshow(G, aspect='auto', cmap='Blues', vmin=0, vmax=1.6)
ax.set_xticks(range(3)); ax.set_xticklabels(obs, fontsize=5.0); ax.set_yticks(range(len(rules))); ax.set_yticklabels(rules, fontsize=5.0)
for i in range(G.shape[0]):
    j = int(np.argmax(G[i])); ax.plot(j, i, marker='o', ms=5, mfc=OI['blue'], mec='white', mew=0.8)
ax.set_xlabel('Observed response of neighbours', fontsize=6); ax.set_ylabel('Intervention acts on', fontsize=6, labelpad=26)
for s in ax.spines.values(): s.set_visible(False)
lab(ax, 'c', -0.42)
save(fig, 'ExtData_Fig4_heatmaps'); print('ed4 ok')
