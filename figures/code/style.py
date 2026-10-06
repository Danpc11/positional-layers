from poslayers.config import OUTDIR, FIGDIR
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'], 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7.5,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.2, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42, 'savefig.dpi': 600})
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'yellow': '#F0E442', 'grey': '#7F7F7F', 'black': '#000000'}
W = 180 / 25.4
OUT = FIGDIR
_PANELS = {}


def lab(ax, s, *_, **__):
    """Register a panel letter. Its position is set in save(), after layout, at the OUTER left edge of the panel
    (tick labels and axis labels included) and on the panel's top edge, aligned across panels of the same row."""
    _PANELS.setdefault(id(ax.figure), []).append((ax, s.lower()))


def _data_boxes(ax, r):
    """Display-space points and boxes occupied by data in ax: line vertices (densified), scatter offsets, bars, error bars."""
    import numpy as _np
    from matplotlib.collections import LineCollection, PathCollection, PolyCollection
    pts, boxes = [], []
    for ln in ax.get_lines():
        if not ln.get_visible(): continue
        xy = _np.column_stack([ln.get_xdata(orig=False), ln.get_ydata(orig=False)]).astype(float)
        if len(xy) == 0: continue
        if len(xy) > 1 and ln.get_linestyle() not in ('None', '', ' '):
            t = _np.linspace(0, 1, 12)[:, None]
            xy = _np.vstack([a + t * (b - a) for a, b in zip(xy[:-1], xy[1:])] + [xy[-1:]])
        pts.append(ln.get_transform().transform(xy))
    for c in ax.collections:
        if not c.get_visible(): continue
        if isinstance(c, LineCollection):
            for seg in c.get_segments():
                if len(seg) > 1:
                    t = _np.linspace(0, 1, 12)[:, None]; seg = _np.vstack([a + t * (b - a) for a, b in zip(seg[:-1], seg[1:])])
                pts.append(c.get_transform().transform(seg))
        elif isinstance(c, PathCollection) and len(c.get_offsets()):
            pts.append(c.get_offset_transform().transform(_np.asarray(c.get_offsets(), float)))
        elif isinstance(c, PolyCollection):
            try: boxes.append(c.get_window_extent(r))
            except Exception: pass
    for p in ax.patches:
        if p.get_visible() and p.get_width() if hasattr(p, 'get_width') else False:
            boxes.append(p.get_window_extent(r))
    return (_np.vstack(pts) if pts else _np.empty((0, 2))), boxes


def _overlaps(bb, pts, boxes):
    import numpy as _np
    if len(pts) and _np.any((pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1) & (pts[:, 1] >= bb.y0) & (pts[:, 1] <= bb.y1)): return True
    return any(bb.overlaps(b) for b in boxes)


def _fix_legends(fig):
    """Move any legend that covers data to the first free corner of its axes, or above the axes if none is free.
    A legend with the attribute keep_position = True was placed by hand and is not moved."""
    r = fig.canvas.get_renderer()
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is None or not leg.get_visible() or getattr(leg, 'keep_position', False): continue   # placed by hand: leave it
        if not leg.get_window_extent(r).overlaps(ax.get_window_extent(r)): continue      # already outside the axes
        pts, boxes = _data_boxes(ax, r)
        if not _overlaps(leg.get_window_extent(r), pts, boxes): continue
        if getattr(ax, 'name', '') == 'polar': continue
        handles = list(getattr(leg, 'legend_handles', None) or getattr(leg, 'legendHandles', []))
        labels = [t.get_text() for t in leg.get_texts()]
        if not handles or len(handles) != len(labels): continue
        fs = leg.get_texts()[0].get_fontsize() if leg.get_texts() else 5.5
        for loc in ('upper left', 'upper right', 'lower left', 'lower right', 'center left', 'center right', 'upper center', 'lower center'):
            new = ax.legend(handles, labels, loc=loc, fontsize=fs, frameon=False, handlelength=1.4, borderaxespad=0.3)
            fig.canvas.draw()
            if not _overlaps(new.get_window_extent(r), pts, boxes): break
        else:
            ax.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.0, 1.01), fontsize=fs, frameon=False, handlelength=1.4,
                      ncol=1, borderaxespad=0.0)
            fig.canvas.draw()


def _place_letters(fig):
    items = _PANELS.pop(id(fig), [])
    if not items: return
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    boxes = []
    for ax, s in items:
        tb = ax.get_tightbbox(r); ext = [tb] + [a.get_window_extent(r) for a in (ax.yaxis.label, ax.xaxis.label, ax.title) if a.get_text()]
        X0 = min(e.x0 for e in ext); Y1 = max(e.y1 for e in ext)
        leg = ax.get_legend()
        if leg is not None and leg.get_visible(): Y1 = max(Y1, leg.get_window_extent(r).y1)
        (x0, _), (_, y1) = inv.transform([[X0, 0], [0, Y1]])
        boxes.append([ax, s, x0, y1, ax.get_position().y1])
    # panels whose AXES tops are within 3% of the figure height form a row; every letter in a row sits above the highest
    # element of that row (tick labels, titles and legends placed above the axes included)
    boxes.sort(key=lambda b: -b[4]); rows = []
    for b in boxes:
        if rows and abs(rows[-1][0][4] - b[4]) < 0.03: rows[-1].append(b)
        else: rows.append([b])
    for row in rows:
        top = max(b[3] for b in row)
        for ax, s, x0, _, _ in row:
            fig.text(x0, top + 0.004, s, fontsize=9, fontweight='bold', ha='left', va='bottom')


def save(fig, name):
    """Save at exactly 180 mm width (Nature double column) without rescaling text: the figure is scaled (width and height together,
    keeping its aspect ratio) until the tight bounding box is W wide, so font sizes in the code are the final printed sizes. Raster elements at 600 ppi."""
    fig.canvas.draw(); _fix_legends(fig); _place_letters(fig)
    for _ in range(6):
        fig.canvas.draw(); bb = fig.get_tightbbox(fig.canvas.get_renderer()).padded(0.02)
        if abs(bb.width - W) < 0.002: break
        w0 = fig.get_figwidth(); w1 = w0 + (W - bb.width)
        fig.set_size_inches(w1, fig.get_figheight() * w1 / w0, forward=True)          # keep the aspect ratio of the layout
        fig.canvas.draw(); _place_letters(fig)
    fig.savefig(OUT + name + '.pdf', bbox_inches='tight', pad_inches=0.02, dpi=600)
    fig.savefig(OUT + name + '.png', dpi=300, bbox_inches='tight', pad_inches=0.02); plt.close(fig)
A = OUTDIR
NAMES = {'cells_ebv-transformed_lymphocytes': 'EBV lymphocytes', 'cells_ebv_transformed_lymphocytes': 'EBV lymphocytes', 'cells_cultured_fibroblasts': 'Fibroblasts', 'esophagus_gastroesophageal_junction': 'Oesophagus, GEJ',
         'esophagus_mucosa': 'Oesophagus, mucosa', 'esophagus_muscularis': 'Oesophagus, muscularis', 'skin_sun_exposed_lower_leg': 'Skin', 'brain_frontal_cortex_ba9': 'Brain, frontal cortex',
         'brain_caudate_basal_ganglia': 'Brain, caudate', 'brain_cerebellar_hemisphere': 'Brain, cerebellum', 'small_intestine_terminal_ileum': 'Small intestine', 'adipose_visceral_omentum': 'Adipose, visceral',
         'adipose_subcutaneous': 'Adipose, subcutaneous', 'heart_left_ventricle': 'Heart, ventricle', 'heart_atrial_appendage': 'Heart, atrium', 'colon_transverse': 'Colon, transverse', 'colon_sigmoid': 'Colon, sigmoid'}
def tissue_label(t):
    t = t.replace('-', '_') if t not in NAMES else t
    return NAMES.get(t, NAMES.get(t.replace('_', '-'), t.replace('_', ' ').capitalize()))


def placeholder(ax, title, description, fontsize=5.6, key=None):
    """Schematic panel. If figures/schematics/<key>.png exists, the drawing is placed in the panel at its own aspect ratio;
    otherwise a reserved frame with a description of the intended content is drawn (dashed frame, light fill)."""
    import os
    png = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'schematics', f'{key}.png') if key else None
    if png and os.path.exists(png):
        img = plt.imread(png); ax.imshow(img, interpolation='lanczos'); ax.set_anchor('NW'); ax.axis('off'); return
    import textwrap
    from matplotlib.patches import FancyBboxPatch
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    ax.add_patch(FancyBboxPatch((0.01, 0.02), 0.98, 0.96, boxstyle='round,pad=0,rounding_size=0.03', transform=ax.transAxes,
                                fc='#f4f4f4', ec='0.55', lw=0.7, ls=(0, (3, 2)), zorder=0))
    width_in = ax.get_position().width * ax.figure.get_figwidth() * 0.9
    n = max(14, int(width_in * 72 / (fontsize * 0.55)))
    body = '\n'.join(textwrap.fill(p.strip(), n) for p in description.replace('\n', ' ').split('|')) if '|' in description else textwrap.fill(description.replace('\n', ' '), n)
    ax.text(0.5, 0.92, textwrap.fill(title, max(12, int(n * 0.85))), transform=ax.transAxes, ha='center', va='top', fontsize=fontsize + 0.8, fontweight='bold', color='0.25')
    ax.text(0.5, 0.70, body, transform=ax.transAxes, ha='center', va='top', fontsize=fontsize, color='0.35', linespacing=1.35)
