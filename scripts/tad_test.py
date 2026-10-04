"""Does the cis layer respect TAD boundaries? Adjacent gene pairs within one TAD vs across a boundary, at matched distance,
using the tissue's own TAD partition (McArthur & Capra 20-bin landscape, hg19; genes placed with GRCh37 coordinates)."""
import glob, os, numpy as np, pandas as pd, pyannotables as pa, statsmodels.formula.api as smf
TD = glob.glob('/home/claude/repo/data/external/TAD-full/*/data/20binsTADlandscape')[0]; CHR = [str(i) for i in range(1, 23)] + ['X']
def tads(name):
    d = f'{TD}/{name}/'
    b6 = pd.read_csv(d + f'bin_6_{name}.bed', sep='\t', header=None, names=['chr', 's', 'e']); b15 = pd.read_csv(d + f'bin_15_{name}.bed', sep='\t', header=None, names=['chr', 's', 'e'])
    b6['bin'] = b6.e - b6.s + 1; b6['exp_end'] = b6.s + 10 * b6.bin - 1; ends = {c: set(g.e) for c, g in b15.groupby('chr')}
    b6['end'] = [next((e for e in range(r.exp_end - 3, r.exp_end + 4) if e in ends.get(r.chr, set())), np.nan) for _, r in b6.iterrows()]
    T = b6.dropna(subset=['end']).rename(columns={'s': 'start'}); T['chr'] = T.chr.str.replace('chr', ''); T = T[T.chr.isin(CHR)].sort_values(['chr', 'start']).reset_index(drop=True); T['tad_id'] = np.arange(len(T)); return T
G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()][['Chromosome', 'Start', 'End']]; G37.columns = ['chr', 'start', 'end']; G37['chr'] = G37.chr.astype(str); G37 = G37[G37.chr.isin(CHR)]; G37['mid'] = (G37.start + G37.end) / 2
def assign(T):
    out = pd.Series(-1, index=G37.index)
    for c, g in T.groupby('chr'):
        idx = G37.index[G37.chr == c]; mm = G37.loc[idx, 'mid'].values; st_ = g.start.values; en = g.end.values.astype(int); ids = g.tad_id.values
        j = np.searchsorted(st_, mm, side='right') - 1; ok = (j >= 0) & (mm <= en[np.clip(j, 0, len(en) - 1)]); out.loc[idx[ok]] = ids[j[ok]]
    return out
MATCH = {'liver': 'Liver_leung2015', 'adrenal_gland': 'adrenal_schmitt2016', 'artery_aorta': 'aorta_leung2015', 'bladder': 'bladder_schmitt2016', 'brain_cortex': 'cortex_DLPFC_schmitt2016',
         'brain_frontal_cortex_ba9': 'cortex_DLPFC_schmitt2016', 'heart_left_ventricle': 'leftVentricle_leung2015', 'lung': 'lung_schmitt2016', 'pancreas': 'pancreas_schmitt2016',
         'muscle_skeletal': 'psoasMuscle_schmitt2016', 'small_intestine_terminal_ileum': 'smallBowel_schmitt2016', 'spleen': 'spleen_schmitt2016',
         'cells_ebv_transformed_lymphocytes': 'GM12878_lymphoblastoid_Lieberman', 'cells_cultured_fibroblasts': 'IMR90_fetalLungFibroblast_Lieberman'}
TADA = {n: assign(tads(n)) for n in sorted(set(MATCH.values()))}
norm = lambda s: s.lower().replace('-', '_'); BINS = [0, 1e3, 5e3, 2e4, 1e5, 5e5, np.inf]
P = {norm(os.path.basename(f)[6:-7]): pd.read_csv(f) for f in glob.glob('/home/claude/atlas/pairs_*.csv.gz')}
S = {norm(os.path.basename(f)[7:-7]): pd.read_csv(f) for f in glob.glob('/home/claude/atlas/shares_*.csv.gz')}
rows, spec = [], []
for t, tn in MATCH.items():
    d = P[t].merge(S[t][['g1', 'g2', 'egene1', 'egene2', 'share_any', 'share_same']], on=['g1', 'g2']); d = d[d.dist > 0].copy(); d['bin'] = pd.cut(d.dist, BINS).astype(str)
    def fit(tadname, sub=None):
        a = TADA[tadname]; x = d.copy(); x['t1'] = a.reindex(x.g1).values; x['t2'] = a.reindex(x.g2).values; x = x[(x.t1 >= 0) & (x.t2 >= 0)]
        x['same_tad'] = (x.t1 == x.t2).astype(int)
        if sub is not None: x = x[sub(x)]
        m = smf.ols('r ~ same_tad + C(bin) + C(orientation)', data=x).fit(); return m.params['same_tad'], m.pvalues['same_tad'], len(x), x.same_tad.mean()
    b_all, p_all, n_all, f_same = fit(tn)
    b_ng, p_ng, n_ng, _ = fit(tn, sub=lambda x: ~x.share_any)                      # no shared eQTL: non-genetic component
    others = [fit(o)[0] for o in TADA if o != tn]
    rows.append({'tissue': t, 'tad_map': tn, 'pairs': n_all, 'frac_same_tad': f_same, 'b_same_tad': b_all, 'p': p_all, 'b_same_tad_no_shared_eQTL': b_ng, 'p_no_shared_eQTL': p_ng,
                 'b_same_tad_other_maps_median': np.median(others), 'own_map_rank_among_maps': 1 + sum(o > b_all for o in others), 'n_maps': 1 + len(others)})
R = pd.DataFrame(rows); R.to_csv('/home/claude/atlas/tad_cis_test.csv', index=False); pd.set_option('display.width', 250); print(R.round(4).to_string(index=False))
print(f"\nsame-TAD effect > 0 in {(R.b_same_tad > 0).sum()}/{len(R)} (P<0.05 in {(R.p < 0.05).sum()}); median {R.b_same_tad.median():.3f}; without shared eQTL median {R.b_same_tad_no_shared_eQTL.median():.3f}; own map > median of other maps in {(R.b_same_tad > R.b_same_tad_other_maps_median).sum()}/{len(R)}")
