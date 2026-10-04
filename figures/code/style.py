import os, sys
from poslayers.config import DATA, OUTDIR, FIGDIR
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np, pandas as pd
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'], 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7.5,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.2, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42, 'savefig.dpi': 300})
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'yellow': '#F0E442', 'grey': '#7F7F7F', 'black': '#000000'}
W = 180 / 25.4
OUT = FIGDIR
def lab(ax, s, dx=-0.16, dy=1.24): ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', ha='left') if False else ax.text(dx, dy, s.lower(), transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', ha='left')
def save(fig, name): fig.savefig(OUT + name + '.pdf', bbox_inches='tight'); fig.savefig(OUT + name + '.png', dpi=220, bbox_inches='tight'); plt.close(fig)
A = OUTDIR
NAMES = {'cells_ebv-transformed_lymphocytes': 'EBV lymphocytes', 'cells_ebv_transformed_lymphocytes': 'EBV lymphocytes', 'cells_cultured_fibroblasts': 'Fibroblasts', 'esophagus_gastroesophageal_junction': 'Oesophagus, GEJ',
         'esophagus_mucosa': 'Oesophagus, mucosa', 'esophagus_muscularis': 'Oesophagus, muscularis', 'skin_sun_exposed_lower_leg': 'Skin', 'brain_frontal_cortex_ba9': 'Brain, frontal cortex',
         'brain_caudate_basal_ganglia': 'Brain, caudate', 'brain_cerebellar_hemisphere': 'Brain, cerebellum', 'small_intestine_terminal_ileum': 'Small intestine', 'adipose_visceral_omentum': 'Adipose, visceral',
         'adipose_subcutaneous': 'Adipose, subcutaneous', 'heart_left_ventricle': 'Heart, ventricle', 'heart_atrial_appendage': 'Heart, atrium', 'colon_transverse': 'Colon, transverse', 'colon_sigmoid': 'Colon, sigmoid'}
def tissue_label(t):
    t = t.replace('-', '_') if t not in NAMES else t
    return NAMES.get(t, NAMES.get(t.replace('_', '-'), t.replace('_', ' ').capitalize()))
