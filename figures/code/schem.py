"""Conceptual schematics drawn with matplotlib, Nature style (sans, 5.5-6.5 pt, Okabe-Ito)."""
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Arc, Circle, Rectangle
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'grey': '#7F7F7F'}
def _clean(ax): ax.set_xticks([]); ax.set_yticks([]); [ax.spines[s].set_visible(False) for s in ax.spines]
def decomposition(ax):
    """x(s,j) = landscape + technical GC term + cis term + noise, as stacked mini-profiles."""
    _clean(ax); rng = np.random.default_rng(2); n = 70; x = np.arange(n)
    mu = rng.gamma(1.3, 1, n); gc = np.convolve(rng.normal(0, 1, n + 20), np.ones(9) / 9, 'same')[10:10 + n] * 1.4; cis = np.convolve(rng.normal(0, 1, n + 10), [0.5, 1, 0.5], 'same')[5:5 + n] * 0.9
    rows = [('sample profile  x(s, j)', mu + gc + cis + rng.normal(0, 0.3, n), 'k', ''), ('landscape  μ_j', mu, '0.45', 'tissue mean · not positional'), ('technical  b_s · g_j', gc, OI['red'], 'per-sample GC bias × isochores'), ('cis  c(s, j)', cis, OI['blue'], 'shared regulation · carries gene order')]
    for i, (name, y, col, note) in enumerate(rows):
        off = -i * 4.2; ax.plot(x, 0.55 * (y - y.mean()) + off, color=col, lw=0.9); ax.text(-2, off, name, ha='right', va='center', fontsize=5.8, fontweight='bold' if i == 0 else 'normal')
        if note: ax.text(n + 2, off, note, ha='left', va='center', fontsize=5.3, color=col if col != '0.45' else '0.35')
        if i > 0: ax.text(-16, off + 2.1, '=' if i == 1 else '+', ha='center', va='center', fontsize=7, color='0.3')
    ax.set_xlim(-34, n + 70); ax.set_ylim(-15.5, 2.3); ax.text(n / 2, -15.2, 'gene order along a chromosome →', ha='center', fontsize=5.5, color='0.35')
def two_scales(ax):
    _clean(ax); y = 0; ax.plot([0, 100], [y, y], color='k', lw=1.2)
    for i, xg in enumerate(range(6, 100, 11)): ax.add_patch(FancyArrowPatch((xg - 3, y + 0.6), (xg + 3, y + 0.6), arrowstyle='-|>', mutation_scale=5, color=OI['blue'] if i < 2 else OI['purple'], lw=2))
    ax.annotate('', xy=(3, 2.4), xytext=(20, 2.4), arrowprops=dict(arrowstyle='<->', lw=0.7, color=OI['blue'])); ax.text(11.5, 3.0, 'short range\n~1 gene', ha='center', fontsize=5.5, color=OI['blue'])
    ax.annotate('', xy=(25, -2.2), xytext=(97, -2.2), arrowprops=dict(arrowstyle='<->', lw=0.7, color=OI['purple'])); ax.text(61, -3.1, 'domain scale · 7–24 genes (~1–4 Mb)', ha='center', va='top', fontsize=5.5, color=OI['purple'])
    ax.add_patch(Arc((61, y), 60, 7.5, theta1=0, theta2=180, color=OI['purple'], lw=0.9, ls='--')); ax.set_xlim(-3, 103); ax.set_ylim(-6, 7)
def shared_variant(ax):
    _clean(ax); 
    for k, (y0, s2, lab, col) in enumerate([(5, +1, 'same-direction effects → positive coupling', OI['blue']), (0, -1, 'opposite effects → negative coupling', OI['red'])]):
        ax.plot([0, 60], [y0, y0], color='k', lw=1); ax.add_patch(Rectangle((28, y0 - 0.35), 4, 0.7, color=OI['orange'])); ax.plot(30, y0 + 1.3, marker='*', ms=7, color=OI['orange'])
        for xg, s in [(12, 1), (48, s2)]:
            ax.add_patch(FancyArrowPatch((xg - 4, y0 + 0.45), (xg + 4, y0 + 0.45), arrowstyle='-|>', mutation_scale=5, color='0.3', lw=1.8))
            ax.annotate('', xy=(xg, y0 + 2.9 if s > 0 else y0 + 0.9), xytext=(xg, y0 + 0.9 if s > 0 else y0 + 2.9), arrowprops=dict(arrowstyle='-|>', lw=1, color=OI['green'] if s > 0 else OI['red']))
        ax.plot([30, 12], [y0 + 1.3, y0 + 0.9], color=OI['orange'], lw=0.6, ls=':'); ax.plot([30, 48], [y0 + 1.3, y0 + 0.9], color=OI['orange'], lw=0.6, ls=':')
        ax.text(64, y0 + 1.2, lab, fontsize=5.5, va='center', color=col)
    ax.text(30, 8.6, 'shared regulatory variant', ha='center', fontsize=5.5, color=OI['orange']); ax.set_xlim(-2, 125); ax.set_ylim(-1.5, 9.5)
def saturation(ax):
    c = np.logspace(-1.5, 1.5, 200); K = 1; th = c / (c + K); ax.plot(c, th, color=OI['purple'], lw=1.3); ax.set_xscale('log')
    ax.set_xlabel('contact  c', fontsize=5.8); ax.set_ylabel('occupancy  θ = c/(c+K)', fontsize=5.8); ax.set_ylim(0, 1.08)
    ax.annotate('low expression\nk → 1', xy=(0.08, 0.075), xytext=(0.04, 0.55), fontsize=5.5, color=OI['blue'], arrowprops=dict(arrowstyle='-|>', lw=0.6, color=OI['blue']))
    ax.annotate('high expression\nk → 0', xy=(12, 0.92), xytext=(1.6, 0.35), fontsize=5.5, color=OI['red'], arrowprops=dict(arrowstyle='-|>', lw=0.6, color=OI['red']))
    ax.text(0.04, 1.06, 'elasticity = 1 − θ\n(hypothesis)', fontsize=6, va='top', fontweight='bold'); ax.tick_params(labelsize=5.5)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
def perturbation_rules(ax):
    _clean(ax); rules = [('shared element hit\n(HBG1/2 edit · BET inhibitor · variant · disease)', 'neighbours move with coupling', OI['green'], True), ('silences one promoter\n(CRISPRi)', 'only genes < 50 kb, by distance', OI['orange'], False), ('acts through a protein\n(BCL11A edit · kinase inhibitor)', 'no cis propagation', '0.55', False)]
    for i, (what, out, col, shared) in enumerate(rules):
        y = 9 - i * 4.3; ax.plot([2, 40], [y, y], color='k', lw=0.9)
        for xg in (10, 32): ax.add_patch(FancyArrowPatch((xg - 3.5, y + 0.45), (xg + 3.5, y + 0.45), arrowstyle='-|>', mutation_scale=5, color='0.3', lw=1.6))
        if shared: ax.add_patch(Rectangle((19.5, y - 0.35), 3, 0.7, color=OI['red'])); ax.plot([21, 10], [y + 1.1, y + 0.9], color=OI['red'], lw=0.6, ls=':'); ax.plot([21, 32], [y + 1.1, y + 0.9], color=OI['red'], lw=0.6, ls=':'); ax.plot(21, y + 1.6, marker='v', color=OI['red'], ms=5)
        elif i == 1: ax.plot(7, y + 1.6, marker='v', color=OI['orange'], ms=5)
        else: ax.plot(10, y + 2.2, marker='o', mfc='white', color='0.4', ms=4); ax.annotate('', xy=(70, y + 2.2), xytext=(12, y + 2.2), arrowprops=dict(arrowstyle='-|>', lw=0.6, color='0.4', ls='--'))
        ax.text(44, y + 0.9, what, fontsize=5.3, va='center'); ax.text(44, y - 1.1, out, fontsize=5.6, va='center', fontweight='bold', color=col)
    ax.set_xlim(0, 100); ax.set_ylim(-2.5, 12)

# ---------------- Heatmaps ----------------
def heat(ax, M, xt, yt, xlab, ylab, cmap='viridis', vmin=None, vmax=None, fmt='{:.2f}', cbar_label=None, annot=True, fs=5.0):
    import matplotlib.pyplot as plt
    im = ax.imshow(M, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(xt))); ax.set_xticklabels(xt, fontsize=5.4, rotation=30, ha='right')
    ax.set_yticks(range(len(yt))); ax.set_yticklabels(yt, fontsize=5.4)
    ax.set_xlabel(xlab, fontsize=6); ax.set_ylabel(ylab, fontsize=6)
    if annot:
        import numpy as _np
        lo, hi = _np.nanmin(M), _np.nanmax(M)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if _np.isnan(M[i, j]): continue
                rel = (M[i, j] - lo) / (hi - lo + 1e-12)
                ax.text(j, i, fmt.format(M[i, j]), ha='center', va='center', fontsize=fs, color='white' if rel < 0.45 else 'black')
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=5); 
    if cbar_label: cb.set_label(cbar_label, fontsize=5.4)
    for s in ax.spines.values(): s.set_visible(False)
    return im

def layers_diagram(ax):
    """Fig 6a: the four layers, what each is, and the operation that removes it."""
    _clean(ax)
    rows = [('landscape  μ_j', 'mean level of each gene', 'subtract the gene mean', '0.45', 'no'),
            ('GC-associated  b_s·g_j', 'GC bias × isochores\n(may include biology)', 'regress out GC within sample', OI['red'], 'mostly no'),
            ('dosage  d(s,j)', 'copy-number segments', 'regress out continuous CN', OI['orange'], 'no'),
            ('cis  c(s,j)', 'local covariance; consistent with\nshared regulation', '— retained for analysis', OI['blue'], 'in part')]
    ax.text(0, 4.6, 'layer', fontsize=5.8, fontweight='bold'); ax.text(27, 4.6, 'what it is', fontsize=5.8, fontweight='bold')
    ax.text(62, 4.6, 'how to remove it', fontsize=5.8, fontweight='bold'); ax.text(100, 4.6, 'genome\norganisation?', fontsize=5.8, fontweight='bold', ha='center')
    for i, (lay, what, how, col, org) in enumerate(rows):
        y = 3.4 - i * 1.15
        ax.add_patch(FancyBboxPatch((-1.5, y - 0.33), 118, 0.66, boxstyle='round,pad=0.1', fc=('#eef5fb' if org == 'in part' else '0.97'), ec='none', zorder=0))
        ax.text(0, y, lay, fontsize=5.6, va='center', color=col, fontweight='bold' if org == 'in part' else 'normal')
        ax.text(27, y, what, fontsize=5.4, va='center'); ax.text(62, y, how, fontsize=5.4, va='center', color='0.3')
        ax.text(100, y, org, fontsize=5.6, va='center', ha='center', fontweight='bold' if org == 'in part' else 'normal', color=OI['green'] if org == 'in part' else '0.55')
    ax.set_xlim(-3, 120); ax.set_ylim(-1.2, 5.4)

def mechanism_diagram(ax):
    """Fig 6b: the cis layer as shared upstream regulation, two components."""
    _clean(ax)
    ax.plot([2, 96], [0, 0], color='k', lw=1.3)
    for xg, col in [(12, OI['blue']), (26, OI['blue']), (58, OI['purple']), (86, OI['purple'])]:
        ax.add_patch(FancyArrowPatch((xg - 4, 0.6), (xg + 4, 0.6), arrowstyle='-|>', mutation_scale=6, color=col, lw=2.4))
    ax.add_patch(Rectangle((17, -0.4), 4, 0.8, color=OI['orange'])); ax.text(19, -1.3, 'shared element', fontsize=5.3, ha='center', color=OI['orange'])
    ax.plot(19, 1.7, marker='*', ms=8, color=OI['orange']); ax.text(19, 2.4, 'variant', fontsize=5.1, ha='center', color=OI['orange'])
    for xg in (12, 26): ax.plot([19, xg], [1.5, 1.0], color=OI['orange'], lw=0.7, ls=':')
    ax.add_patch(Rectangle((70, -0.4), 4, 0.8, color=OI['green'])); ax.text(72, -1.3, 'shared enhancer', fontsize=5.3, ha='center', color=OI['green'])
    ax.add_patch(Arc((65, 0), 42, 7.0, theta1=0, theta2=180, color=OI['purple'], lw=1.1, ls='--'))
    ax.text(65, 4.0, 'coupling ~ contact^k, k ≈ 0.5', fontsize=5.4, ha='center', color=OI['purple'])
    for xb, lab in [(40, 'stable\nboundary'), (96, '')]:
        if lab: ax.plot([xb, xb], [-1.0, 1.6], color='0.25', lw=1.6); ax.text(xb, 2.3, lab, fontsize=5.3, ha='center', color='0.25')
    ax.annotate('', xy=(12, -2.3), xytext=(26, -2.3), arrowprops=dict(arrowstyle='<->', lw=0.7, color=OI['blue']))
    ax.text(19, -3.0, 'adjacent\ncontact-independent', fontsize=5.3, ha='center', va='top', color=OI['blue'])
    ax.annotate('', xy=(58, -2.3), xytext=(86, -2.3), arrowprops=dict(arrowstyle='<->', lw=0.7, color=OI['purple']))
    ax.text(72, -3.0, 'domain scale\ncontact-dependent', fontsize=5.3, ha='center', va='top', color=OI['purple'])
    ax.text(50, 6.2, 'lower with STAG2 mutation', fontsize=5.4, ha='center', color=OI['red'])
    ax.set_xlim(0, 100); ax.set_ylim(-6.2, 7.0)
