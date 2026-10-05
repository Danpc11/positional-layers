"""Extracts Hi-C contacts between gene promoters (TSS) up to 2 Mb apart, streaming only the needed regions
from the public Rao et al. 2014 in situ Hi-C maps (hg19). Nothing large is downloaded.

Usage:
    python scripts/make_tss_table.py     # writes DATA/hic/genes_hg19_tss.tsv
    python scripts/hic_extract.py
Output: DATA/hic/hic_contacts_<cell>.csv.gz and DATA/hic/hic_expected_<cell>.csv.gz (a few tens of MB each).
"""
import os
from poslayers.config import DATA
import numpy as np
import pandas as pd
import hicstraw
os.makedirs(DATA + 'hic', exist_ok=True)
MAPS = {'GM12878': 'https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic',
        'IMR90': 'https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic'}
RES, MAXD, WIN = 25000, 2_000_000, 6_000_000          # bin size, maximum distance, window streamed at a time
G = pd.read_csv(DATA + 'hic/genes_hg19_tss.tsv', sep='\t', dtype={'chr': str})
for name, url in MAPS.items():
    print(f'== {name}', flush=True); hic = hicstraw.HiCFile(url)
    lengths = {c.name: c.length for c in hic.getChromosomes()}
    rows, expected = [], []
    for c, gc in G.groupby('chr'):
        if c not in lengths: continue
        try: mzd = hic.getMatrixZoomData(c, c, 'observed', 'KR', 'BP', RES)
        except Exception as e: print('  skip', c, e); continue
        band = {}
        for start in range(0, lengths[c], WIN - MAXD):
            end = min(start + WIN, lengths[c])
            for r in mzd.getRecords(start, end, start, end):
                d = r.binY - r.binX
                if 0 <= d <= MAXD and np.isfinite(r.counts): band[(r.binX, r.binY)] = r.counts
            if end >= lengths[c]: break
        # expected contact per distance (zeros included): sum of counts / number of bin pairs at that distance
        nb = lengths[c] // RES + 1; s = {}
        for (x, y), v in band.items(): s[y - x] = s.get(y - x, 0.0) + v
        for d in range(0, MAXD // RES + 1): expected.append((c, d * RES, s.get(d * RES, 0.0) / max(nb - d, 1)))
        gc = gc.sort_values('tss'); b = (gc.tss.values // RES) * RES; ids = gc.gene_id.values
        for i in range(len(ids)):                                  # each gene with every gene up to 2 Mb downstream
            j = i + 1
            while j < len(ids) and gc.tss.values[j] - gc.tss.values[i] <= MAXD:
                x, y = min(b[i], b[j]), max(b[i], b[j])
                rows.append((ids[i], ids[j], c, int(gc.tss.values[j] - gc.tss.values[i]), band.get((x, y), 0.0), int(y - x)))
                j += 1
        print(f'  chr{c}: {len(band):,} contacts in band, {len(rows):,} gene pairs so far', flush=True)
    pd.DataFrame(rows, columns=['g1', 'g2', 'chr', 'tss_distance', 'contact_KR', 'bin_distance']).to_csv(DATA + f'hic/hic_contacts_{name}.csv.gz', index=False)
    pd.DataFrame(expected, columns=['chr', 'bin_distance', 'expected_KR']).to_csv(DATA + f'hic/hic_expected_{name}.csv.gz', index=False)
    print(f'  written hic_contacts_{name}.csv.gz', flush=True)
