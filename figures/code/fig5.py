import sys; sys.path.insert(0, '/home/claude/figs'); from style import *
import pyannotables as pa
from scipy import stats
G = pa.tables()['homo_sapiens-GRCh38-ensembl100']; G = G[~G.index.duplicated()]
BY = pd.read_csv('/mnt/user-data/outputs/NAR_piloto/bystander_liver_t_stage.csv'); ED = pd.read_csv('/home/claude/edit/edit_effects.csv', index_col=0)
WB = pd.read_csv(A + 'pairs_whole_blood.csv.gz'); DR = pd.read_csv('/home/claude/slam/drug_cis_results_thr3.csv'); CR = pd.read_csv('/home/claude/perturb/crispri_pairs.csv.gz')
import schem
fig = plt.figure(figsize=(W, W * 0.9)); gs = fig.add_gridspec(3, 6, hspace=0.6, wspace=1.4)
def quint(d, a, b):
    q = pd.qcut(d.r, 5, labels=False); return [d.r[q == i].mean() for i in range(5)], [stats.pearsonr(d[a][q == i], d[b][q == i])[0] for i in range(5)]
ax = fig.add_subplot(gs[0, 0:3]); schem.perturbation_rules(ax); lab(ax, 'a', -0.04)
ax = fig.add_subplot(gs[0, 3:6]); x, y = quint(BY, 't1', 't2'); ax.plot(x, y, 'o-', color=OI['purple'], ms=4, lw=1.2); ax.plot([-0.3, 0.6], [-0.3, 0.6], 'k:', lw=0.6)
ax.axhline(0, color='0.8', lw=0.5); ax.set_xlabel('Coupling in healthy liver (GTEx)'); ax.set_ylabel('Concordance of fibrosis effects\n(437 biopsies)'); lab(ax, 'b')
ax = fig.add_subplot(gs[1, 0:2]); pos = ((G.Start + G.End) / 2).reindex(ED.index); ch = G.Chromosome.astype(str).reindex(ED.index)
for (c, p0, col, tcol, nm) in [('11', 5.25e6, OI['red'], 't_HBG', 'HBG1/2 promoter edit'), ('2', 60.45e6, OI['blue'], 't_BCL11A', 'BCL11A enhancer edit (Casgevy)')]:
    w = (ch == c) & ((pos - p0).abs() < 1.5e6); xx = (pos[w] - p0) / 1e6; yy = ED.loc[w, tcol]
    ax.scatter(xx, yy, s=7, color=col, alpha=0.75, lw=0, label=nm)
    for gid in ED.index[w][np.abs(yy.values) > 4]: ax.annotate(ED.at[gid, 'sym'], (xx[gid], yy[gid]), fontsize=4.8, xytext=(2, 1), textcoords='offset points', color=col)
ax.axhline(0, color='0.7', lw=0.5); ax.set_xlabel('Distance from edited site (Mb)'); ax.set_ylabel('Response (t)'); ax.legend(loc='upper right', fontsize=5.2); lab(ax, 'c')
ax = fig.add_subplot(gs[1, 2:4]); P = WB[WB.dist > 0]
for tcol, col, nm in [('t_HBG', OI['red'], 'HBG1/2'), ('t_BCL11A', OI['blue'], 'BCL11A')]:
    d = P.assign(t1=ED[tcol].reindex(P.g1).values, t2=ED[tcol].reindex(P.g2).values).dropna(); x, y = quint(d, 't1', 't2'); ax.plot(x, y, 'o-', color=col, ms=4, lw=1.2, label=nm)
ax.axhline(0, color='0.8', lw=0.5); ax.set_xlabel('Coupling in healthy blood (GTEx)'); ax.set_ylabel('Concordance of edit responses'); ax.legend(loc='upper left'); lab(ax, 'd')
ax = fig.add_subplot(gs[1, 4:6]); CR = CR[np.isfinite(CR.resp) & np.isfinite(CR.r)]; CR['resp'] = CR.resp.clip(-20, 20); C2 = CR[(CR.type == 'cis') & (CR.kd > 2)].copy()
C2['db'] = pd.cut(C2.dist, [-1, 1e4, 5e4, 2e5, 5e5, 1e6], labels=['<10 kb', '10–50 kb', '50–200 kb', '0.2–0.5 Mb', '0.5–1 Mb']); qs = C2.r.quantile([1 / 3, 2 / 3]).values
C2['rt'] = pd.cut(C2.r, [-1, qs[0], qs[1], 1], labels=['low', 'mid', 'high']); tab = C2.pivot_table(index='db', columns='rt', values='resp', aggfunc='mean', observed=True)
xx = np.arange(len(tab))
for j, (c, col) in enumerate(zip(['low', 'mid', 'high'], ['0.75', OI['sky'], OI['blue']])): ax.bar(xx + (j - 1) * 0.26, tab[c], 0.26, color=col, label=f'{c} coupling')
ax.axhline(0, color='k', lw=0.6); ax.set_xticks(xx); ax.set_xticklabels(tab.index, fontsize=5.6, rotation=25); ax.set_ylabel('Neighbour response (robust z)'); ax.legend(loc='lower right', fontsize=5.6)
lab(ax, 'e', -0.12)
ax = fig.add_subplot(gs[2, 0:4]); DR['grp'] = np.where(DR.perturbation.str.contains('JQ1'), 'BET inhibitor (JQ1)', np.where(DR.perturbation.str.contains('BRD4'), 'BRD4 degradation', np.where(DR['class'] == 'signalling', 'signalling inhibitor', 'CDK9 inhibitor')))
cmap = {'BET inhibitor (JQ1)': OI['green'], 'BRD4 degradation': OI['sky'], 'CDK9 inhibitor': OI['orange'], 'signalling inhibitor': '0.55'}
d = DR.sort_values(['grp', 'slope_on_baseline_coupling']); y = np.arange(len(d))
ax.barh(y, d.slope_on_baseline_coupling, color=[cmap[g] for g in d.grp], height=0.7); ax.axvline(0, color='k', lw=0.6)
ax.set_yticks(y); ax.set_yticklabels([f'{a} · {b}' for a, b in zip(d.cell, d.perturbation)], fontsize=5); ax.set_xlabel('Propagation to coupled neighbours (slope)')
from matplotlib.patches import Patch; ax.legend(handles=[Patch(color=v, label=k) for k, v in cmap.items()], loc='upper right', fontsize=5.4); lab(ax, 'f', -0.55)
save(fig, 'Fig5_perturbations'); print('Fig5_perturbations ok')
