"""Where inputs, intermediate tables and figures live. Every script imports these three strings.

    POSLAYERS_DATA     raw public inputs (default: data/)        -- see data/MANIFEST.md
    POSLAYERS_RESULTS  every table the scripts produce (default: results/)
    POSLAYERS_FIGS     figure files (default: figures/output/)

Each value ends with '/', so scripts write DATA + 'gtex/...' or OUTDIR + f'pairs_{t}.csv.gz'.
"""
import os

def _dir(var, default):
    p = os.path.abspath(os.environ.get(var, default)).rstrip('/') + '/'
    return p

DATA = _dir('POSLAYERS_DATA', 'data')
OUTDIR = _dir('POSLAYERS_RESULTS', 'results')
FIGDIR = _dir('POSLAYERS_FIGS', 'figures/output')
for _d in (OUTDIR, FIGDIR):
    os.makedirs(_d, exist_ok=True)
