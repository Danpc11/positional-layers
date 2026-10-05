import os, sys
from poslayers.config import OUTDIR
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.environ.get('LIVER_SPECTRA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Liver_Spectra')), 'scripts'))  # liver pipeline: github.com/Danpc11/Liver_Spectra
from common import *
X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl'))
s = M.index[M.estadio != 'Control']; grid = load_grid(); N = chrom_lengths(grid)
U = pd.read_csv(tab('Table_S2a_consensus_spectrum_universal_peaks.csv')); U['chr'] = U.chr.astype(str)
def per_chrom_signal(vals, g, n):
    x = np.zeros(n); x[g.grid_index.values - 1] = vals; return x
def pgram(x): n = len(x); return np.abs(np.fft.rfft(x)[1:n // 2 + 1]) ** 2
def whiten(P, n):
    k = np.arange(1, len(P) + 1); lx = np.log10(k / n); b = np.polyfit(lx, np.log10(P + 1e-12), 1); w = P / 10 ** (b[0] * lx + b[1]); return w / w.mean()
EUL = 0.5772156649
rows = []; stat_tot = dyn_tot = 0.0
for c in CHR:
    g = keep[keep.chr == c].sort_values('grid_index'); n = N[c]
    Y = A.loc[g.gene_id, s].values                      # genes x samples, batch-corrected log expression
    Yc = Y - Y.mean(0, keepdims=True)                   # centred across genes within each sample (as in the pipeline)
    mu = Yc.mean(1); D = Yc - mu[:, None]               # static landscape and per-sample deviations
    stat_tot += np.sum(mu ** 2) * Y.shape[1]; dyn_tot += np.sum(D ** 2)
    Wfull = np.array([whiten(pgram(per_chrom_signal(Yc[:, i], g, n)), n) for i in range(Y.shape[1])])
    Wmu = whiten(pgram(per_chrom_signal(mu, g, n)), n)
    Wdyn = np.array([whiten(pgram(per_chrom_signal(D[:, i], g, n)), n) for i in range(Y.shape[1])])
    cons = np.exp(np.log(Wfull).mean(0) + EUL); cons_dyn = np.exp(np.log(Wdyn).mean(0) + EUL)
    for k in range(len(cons)):
        rows.append({'chr': c, 'k': k + 1, 'period': n / (k + 1), 'cons_full': cons[k], 'w_static': Wmu[k], 'cons_dynamic': cons_dyn[k],
                     'frac_full_gt3': np.mean(Wfull[:, k] > 3), 'frac_dyn_gt3': np.mean(Wdyn[:, k] > 3)})
R = pd.DataFrame(rows).merge(U[['chr', 'k', 'universal']], on=['chr', 'k'], how='left')
R.to_csv(OUTDIR + 'p1_spectra.csv', index=False)
u = R.universal.fillna(False).astype(bool)
print('frequencies:', len(R), '| universal peaks (pipeline):', int(u.sum()))
print('Parseval share of total positional variance: static landscape %.1f%%, per-sample deviations %.1f%%' % (100 * stat_tot / (stat_tot + dyn_tot), 100 * dyn_tot / (stat_tot + dyn_tot)))
print('r(log consensus spectrum, log spectrum of the static landscape) = %.3f' % np.corrcoef(np.log(R.cons_full), np.log(R.w_static))[0, 1])
print('universal peaks with static-landscape power > 3x: %d / %d' % ((R.w_static[u] > 3).sum(), u.sum()))
print('frequencies universal in the DYNAMIC spectra (>3x in >=90%% of biopsies): %d' % (R.frac_dyn_gt3 >= 0.9).sum())
print('r(log consensus, log dynamic consensus) = %.3f' % np.corrcoef(np.log(R.cons_full), np.log(R.cons_dynamic))[0, 1])
