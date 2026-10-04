"""TAD partitions and gene-to-TAD assignment, shared by tad_test.py, fric1.py and tad_bootstrap.py.
TADs: McArthur & Capra 20-bin TAD landscape (hg19); genes placed by GRCh37 midpoint."""
import glob, numpy as np, pandas as pd, pyannotables as pa
from poslayers.config import DATA

CHR = [str(i) for i in range(1, 23)] + ['X']
_found = glob.glob(DATA + 'TAD-full/*/data/20binsTADlandscape')
TD = _found[0] if _found else DATA + 'TAD-full/MISSING/data/20binsTADlandscape'
MATCH = {'liver': 'Liver_leung2015', 'adrenal_gland': 'adrenal_schmitt2016', 'artery_aorta': 'aorta_leung2015', 'bladder': 'bladder_schmitt2016',
         'brain_cortex': 'cortex_DLPFC_schmitt2016', 'brain_frontal_cortex_ba9': 'cortex_DLPFC_schmitt2016', 'heart_left_ventricle': 'leftVentricle_leung2015',
         'lung': 'lung_schmitt2016', 'pancreas': 'pancreas_schmitt2016', 'muscle_skeletal': 'psoasMuscle_schmitt2016',
         'small_intestine_terminal_ileum': 'smallBowel_schmitt2016', 'spleen': 'spleen_schmitt2016',
         'cells_ebv_transformed_lymphocytes': 'GM12878_lymphoblastoid_Lieberman', 'cells_cultured_fibroblasts': 'IMR90_fetalLungFibroblast_Lieberman'}

G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()][['Chromosome', 'Start', 'End']]
G37.columns = ['chr', 'start', 'end']; G37['chr'] = G37.chr.astype(str); G37 = G37[G37.chr.isin(CHR)].copy(); G37['mid'] = (G37.start + G37.end) / 2


def tads(name):
    """TAD intervals for one cell type or tissue from the 20-bin landscape (bin 6 starts closed by a bin-15 end)."""
    d = f'{TD}/{name}/'
    b6 = pd.read_csv(d + f'bin_6_{name}.bed', sep='\t', header=None, names=['chr', 's', 'e'])
    b15 = pd.read_csv(d + f'bin_15_{name}.bed', sep='\t', header=None, names=['chr', 's', 'e'])
    b6['bin'] = b6.e - b6.s + 1; b6['exp_end'] = b6.s + 10 * b6.bin - 1; ends = {c: set(g.e) for c, g in b15.groupby('chr')}
    b6['end'] = [next((e for e in range(r.exp_end - 3, r.exp_end + 4) if e in ends.get(r.chr, set())), np.nan) for _, r in b6.iterrows()]
    T = b6.dropna(subset=['end']).rename(columns={'s': 'start'}); T['chr'] = T.chr.str.replace('chr', '')
    T = T[T.chr.isin(CHR)].sort_values(['chr', 'start']).reset_index(drop=True); T['tad_id'] = np.arange(len(T))
    return T


def assign(T):
    """Series gene_id -> TAD id (-1 when the gene midpoint is outside every TAD)."""
    out = pd.Series(-1, index=G37.index)
    for c, g in T.groupby('chr'):
        idx = G37.index[G37.chr == c]; mm = G37.loc[idx, 'mid'].values; st_ = g.start.values; en = g.end.values.astype(int); ids = g.tad_id.values
        j = np.searchsorted(st_, mm, side='right') - 1; ok = (j >= 0) & (mm <= en[np.clip(j, 0, len(en) - 1)]); out.loc[idx[ok]] = ids[j[ok]]
    return out
