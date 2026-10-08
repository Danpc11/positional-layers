"""Extended Data Fig. 1: validation of the spectral decomposition and of the GC correction (panels moved from Fig. 1)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import OI, OUTDIR, W, lab, save
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from theory_sim import simulate
P1 = pd.read_csv(OUTDIR + 'p1_spectra.csv'); S = pd.read_csv(OUTDIR + 'sim_decay_recovery.csv')
fig, axs = plt.subplots(1, 3, figsize=(W, W * 0.34)); plt.subplots_adjust(wspace=0.55)
# a: the identity holds to numerical precision in a simulated genome
Xs, chrs, GC, mu_, z = simulate(seed=1, land_sd=4); c0 = chrs == 0; Xc = Xs[:, c0]; nn = Xc.shape[1]; mbar = Xc.mean(0); Dv = Xc - mbar
lhs = np.mean([np.abs(np.fft.fft(x)) ** 2 for x in Xc], 0); Cl = np.array([np.mean(np.sum(Dv * np.roll(Dv, -L, axis=1), 1)) for L in range(nn)])
land = np.abs(np.fft.fft(mbar)) ** 2; cov = np.real(np.fft.fft(Cl)); k = slice(1, nn // 2)
ax = axs[0]; ax.loglog(lhs[k], (land + cov)[k], '.', ms=1.4, color='0.35', rasterized=True); lim = [lhs[k].min(), lhs[k].max()]; ax.plot(lim, lim, color=OI['red'], lw=0.7)
ax.set_xlabel('Mean single-sample spectrum'); ax.set_ylabel('Tissue-average + covariance terms'); ax.text(0.95, 0.06, f'tissue average: {np.sum(land) / np.sum(lhs):.0%} of power', transform=ax.transAxes, fontsize=5.6, ha='right'); lab(ax, 'a')
# b: in 427 liver biopsies the consensus spectrum is the spectrum of the tissue-average profile
ax = axs[1]; ax.loglog(P1.w_static, P1.cons_full, '.', ms=1, color='0.45', rasterized=True); u = P1.universal.fillna(False).astype(bool)
ax.loglog(P1.w_static[u], P1.cons_full[u], '.', ms=3, color=OI['red'], label='186 universal peaks')
ax.set_xlabel('Spectrum of the tissue-average profile'); ax.set_ylabel('Consensus spectrum'); ax.legend(loc='lower right', markerscale=2, fontsize=5.5)
ax.text(0.05, 0.9, f"r = {np.corrcoef(np.log(P1.w_static), np.log(P1.cons_full))[0, 1]:.3f}", transform=ax.transAxes, fontsize=6); lab(ax, 'b')
# c: GC correction restores the true decay length in simulations
ax = axs[2]; g = S.groupby('true_lambda')[['naive', 'corrected']].median()
ax.plot(g.index, g.naive, 'o-', color=OI['red'], ms=3, lw=1, label='without GC correction'); ax.plot(g.index, g.corrected, 's-', color=OI['blue'], ms=3, lw=1, label='GC-corrected')
ax.plot([1, 13], [1, 13], 'k:', lw=0.7, label='identity'); ax.set_xlim(0, 13.5); ax.set_ylim(0, 17); ax.set_xlabel('True decay length (genes)'); ax.set_ylabel('Estimated decay length (genes)'); ax.legend(loc='upper left', fontsize=5.5)
lab(ax, 'c')
save(fig, 'ExtData_Fig1_validation'); print('ExtData_Fig1_validation ok')
