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


def upsert_csv(df, path, key):
    """Write rows to a results table so that re-running a script REPLACES its earlier rows instead of appending to them.

    Rows already in `path` whose `key` matches a row of `df` are dropped; the result must have one row per key."""
    import pandas as pd
    out = df
    if os.path.exists(path):
        prev = pd.read_csv(path)
        if len(prev):
            merged = prev.merge(df[key].drop_duplicates(), on=key, how='left', indicator=True)
            out = pd.concat([prev[(merged['_merge'] == 'left_only').values], df], ignore_index=True)
    if out.duplicated(key).any():
        raise ValueError(f'{path}: more than one row per {key}')
    out.to_csv(path, index=False)
    return out


# Re-running a script recomputes everything by default. Set POSLAYERS_RESUME=1 to reuse units (tissues, replicates) that a
# previous run already wrote; only do this when data, filters and code are unchanged since that run.
RESUME = os.environ.get('POSLAYERS_RESUME', '0') == '1'
