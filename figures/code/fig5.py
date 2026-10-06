"""Fig. 5: layers of positional expression variation and the shared-source model (schematics drawn in code)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import OI, W, lab, save
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyBboxPatch, Rectangle, Polygon

GENE, GREY, DARK, PURPLE = OI['blue'], '0.55', '0.15', '#7b5ea7'
rng = np.random.default_rng(3)
fig = plt.figure(figsize=(W, W * 0.6)); gs = fig.add_gridspec(2, 2, hspace=0.12, wspace=0.06, height_ratios=[1, 1.05])

def gene(ax, x0, y, w=10, h=4, color=GENE, tss=True):
    ax.add_patch(Rectangle((x0, y - h / 2), w, h, color=color, lw=0, zorder=3))
    if tss:                                                          # transcription start: bent arrow
        ax.plot([x0, x0, x0 + 4.5], [y + h / 2, y + h / 2 + 3, y + h / 2 + 3], color=color, lw=0.9, zorder=3, solid_capstyle='butt')
        ax.add_patch(Polygon([[x0 + 4.5, y + h / 2 + 4.3], [x0 + 6.5, y + h / 2 + 3], [x0 + 4.5, y + h / 2 + 1.7]], color=color, lw=0, zorder=3))
def chrom(ax, x0, x1, y): ax.plot([x0, x1], [y, y], color=GREY, lw=0.9, zorder=1)
def arrow(ax, xy0, xy1, color, lw=1.0, ms=6, ls='-'):
    ax.annotate('', xy=xy1, xytext=xy0, arrowprops=dict(arrowstyle='-|>', color=color, lw=lw, mutation_scale=ms, linestyle=ls, shrinkA=0, shrinkB=0))

# ---------------------------------------------------------------- a: the four layers, in the order they are removed
ax = fig.add_subplot(gs[0, :]); ax.set_xlim(0, 400); ax.set_ylim(0, 100); ax.axis('off')
ax.text(200, 99, 'Separating layers of expression variation along a chromosome', ha='center', va='top', fontsize=7.5, fontweight='bold')
t = np.linspace(0, 1, 140)
avg = 0.45 + 0.22 * np.exp(-((t - 0.3) / 0.12) ** 2) + 0.16 * np.exp(-((t - 0.72) / 0.1) ** 2)
gc = 0.25 + 0.5 * np.exp(-((t - 0.32) / 0.14) ** 2) + 0.35 * np.exp(-((t - 0.75) / 0.12) ** 2)
cn = np.where((t > 0.38) & (t < 0.72), 0.75, np.where(t >= 0.72, 0.45, 0.2))
shared = np.convolve(rng.normal(size=t.size + 8), np.ones(9) / 9, 'valid')[: t.size]
cols = [('1', 'Tissue average', 'Subtract the tissue average\n(mean landscape)'),
        ('2', 'GC-associated layer', 'Remove the per-sample GC trend\n(expression tracks gene GC)'),
        ('3', 'Copy number (tumours)', 'Remove each gene\'s dosage\n(expression tracks copy number)'),
        ('4', 'Cis layer', 'Measure coupling between\nneighbouring genes')]
for k, (num, title, action) in enumerate(cols):
    x0 = k * 100 + 9; wpl = 80; y0, hpl = 34, 40
    ax.plot([x0, x0, x0 + wpl], [y0 + hpl, y0, y0], color='0.6', lw=0.7)                         # mini axes
    X = x0 + 2 + t * (wpl - 4); Y = lambda v: y0 + 3 + v * (hpl - 6)
    if k == 0:
        ax.plot(X, Y(np.clip(avg + 0.09 * rng.normal(size=t.size), 0, 1)), color=GENE, lw=0.8)
        ax.plot(X, Y(avg), color=DARK, lw=1.2, ls=(0, (3, 1.5)))
        ax.text(x0 + wpl, y0 + hpl + 1, 'one sample', color=GENE, fontsize=5.5, ha='right', va='bottom')
        ax.text(x0 + wpl, y0 + hpl - 5, 'tissue average', color=DARK, fontsize=5.5, ha='right', va='bottom')
    elif k == 1:
        ax.plot(X, Y(gc), color=PURPLE, lw=1.3)
        ax.plot(X, Y(np.clip(0.7 * gc + 0.1 + 0.05 * rng.normal(size=t.size), 0, 1)), color='0.45', lw=0.7)
        ax.text(x0 + wpl, y0 + hpl + 1, 'gene GC content', color=PURPLE, fontsize=5.5, ha='right', va='bottom')
        ax.text(x0 + wpl, y0 + hpl - 5, 'expression', color='0.4', fontsize=5.5, ha='right', va='bottom')
    elif k == 2:
        ax.plot(X, Y(cn), color=OI['green'], lw=1.3, drawstyle='steps-mid')
        ax.plot(X, Y(np.clip(0.8 * cn + 0.08 + 0.05 * rng.normal(size=t.size), 0, 1)), color='0.45', lw=0.7)
        ax.text(x0 + wpl, y0 + hpl + 1, 'copy number', color=OI['green'], fontsize=5.5, ha='right', va='bottom')
        ax.text(x0 + wpl, y0 + hpl - 5, 'expression', color='0.4', fontsize=5.5, ha='right', va='bottom')
    else:
        g1 = 0.62 + 0.25 * shared / np.abs(shared).max() + 0.04 * rng.normal(size=t.size)
        g2 = 0.30 + 0.25 * shared / np.abs(shared).max() + 0.04 * rng.normal(size=t.size)
        ax.plot(X, Y(np.clip(g1, 0, 1)), color=GENE, lw=0.8); ax.plot(X, Y(np.clip(g2, 0, 1)), color=PURPLE, lw=0.8)
        ax.text(x0 + wpl, y0 + hpl + 1, 'neighbouring genes vary together', color=DARK, fontsize=5.5, ha='right', va='bottom')
    ax.text(x0, 86, f'{num}', fontsize=7, fontweight='bold', color='white', ha='center', va='center',
            bbox=dict(boxstyle='circle,pad=0.25', fc=DARK, ec='none'))
    ax.text(x0 + 5, 86, title, fontsize=7, fontweight='bold', ha='left', va='center')
    ax.text(x0 + wpl / 2, 26, action, fontsize=6, ha='center', va='top', linespacing=1.3)
    if k < 3: arrow(ax, (x0 + wpl + 3, y0 + hpl / 2), (x0 + wpl + 15, y0 + hpl / 2), DARK, lw=1.2, ms=8)
lab(ax, 'a', -0.0)

# ---------------------------------------------------------------- b: one variant, two genes
ax = fig.add_subplot(gs[1, 0]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')
ax.text(50, 99, 'A shared variant sets the sign of coupling', ha='center', va='top', fontsize=7.5, fontweight='bold')
for y, same in ((63, True), (25, False)):
    chrom(ax, 3, 60, y); gene(ax, 8, y); gene(ax, 42, y)
    ax.text(13, y - 5, 'Gene 1', fontsize=5.5, ha='center', va='top'); ax.text(47, y - 5, 'Gene 2', fontsize=5.5, ha='center', va='top')
    vx, vy = 30, y + 17
    ax.plot([vx, 14], [vy, y + 8.5], color=OI['orange'], lw=0.9, ls=(0, (2.5, 1.5))); ax.plot([vx, 46], [vy, y + 8.5], color=OI['orange'], lw=0.9, ls=(0, (2.5, 1.5)))
    ax.add_patch(Polygon([[vx, vy + 2.6], [vx + 2.2, vy], [vx, vy - 2.6], [vx - 2.2, vy]], color=OI['orange'], lw=0, zorder=4))
    ax.text(vx, vy + 4, 'variant', fontsize=5.5, ha='center', va='bottom', color=OI['orange'])
    arrow(ax, (21, y - 4), (21, y + 5), OI['green'], lw=1.3, ms=7)                                  # effect on gene 1: up, beside the gene
    if same: arrow(ax, (55, y - 4), (55, y + 5), OI['green'], lw=1.3, ms=7)                         # effect on gene 2
    else: arrow(ax, (55, y + 5), (55, y - 4), OI['red'], lw=1.3, ms=7)
    ax.text(64, y + 7, 'Same-direction effects' if same else 'Opposite-direction effects', fontsize=6.5, fontweight='bold', va='center')
    ax.text(64, y + 1, 'positive contribution to coupling' if same else 'negative contribution to coupling', fontsize=6,
            color=OI['green'] if same else OI['red'], va='center')
ax.text(50, 4, 'contribution of the variant  =  2p(1 − p) · β$_1$ · β$_2$', fontsize=6, ha='center', va='center', color='0.3')
lab(ax, 'b', -0.0)

# ---------------------------------------------------------------- c: architecture constrains which genes share sources
ax = fig.add_subplot(gs[1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')
ax.text(50, 99, 'Architecture constrains coupling', ha='center', va='top', fontsize=7.5, fontweight='bold')
def coupling_arc(x0, x1, y, h, strong):
    ax.add_patch(Arc(((x0 + x1) / 2, y), x1 - x0, h, theta1=0, theta2=180, color=DARK if strong else '0.55',
                     lw=2.2 if strong else 0.9, ls='-' if strong else (0, (2, 1.5)), zorder=2))
y = 60
ax.add_patch(FancyBboxPatch((6, y - 7), 47, 14, boxstyle='round,pad=0,rounding_size=2.5', fc='#eaf3fb', ec=OI['sky'], lw=0.9, zorder=0))
chrom(ax, 3, 97, y); gene(ax, 10, y); gene(ax, 40, y); gene(ax, 76, y)
ax.plot([63, 63], [y - 9, y + 20], color=PURPLE, lw=1.0, ls=(0, (2.5, 1.5)))
coupling_arc(15, 45, y + 6, 22, True); coupling_arc(45, 81, y + 6, 22, False)
ax.text(29.5, y - 10, 'within a domain', fontsize=6, ha='center', va='top'); ax.text(63, y - 10, 'boundary', fontsize=5.5, ha='center', va='top', color=PURPLE)
ax.text(81, y - 10, 'across a boundary', fontsize=6, ha='center', va='top')
y = 22
chrom(ax, 3, 97, y); gene(ax, 10, y); gene(ax, 44, y, color=OI['green']); gene(ax, 76, y)
for dx in (0, 3, 6):                                                         # transcripts leaving the active gene
    ax.plot([47 + dx, 49 + dx], [y + 10, y + 13.5], color=OI['green'], lw=0.8)
coupling_arc(15, 81, y + 6, 26, False)
ax.text(49, y - 5, 'active gene', fontsize=5.5, ha='center', va='top', color=OI['green'])
ax.text(50, y - 11, 'an active gene between two genes: weaker coupling', fontsize=6, ha='center', va='top')
ax.plot([66, 72], [94 - 6, 94 - 6], color=DARK, lw=2.2); ax.text(73, 88, 'strong', fontsize=5.5, va='center')
ax.plot([83, 89], [88, 88], color='0.55', lw=0.9, ls=(0, (2, 1.5))); ax.text(90, 88, 'weak', fontsize=5.5, va='center')
ax.text(65, 88, 'coupling:', fontsize=5.5, va='center', ha='right')
lab(ax, 'c', -0.0)
save(fig, 'Fig5_model'); print('Fig5_model ok')
