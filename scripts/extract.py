"""Stream the GDC-PANCAN HTSeq matrix (log2(count+1)) from the zip and keep protein-coding genes for selected cohorts."""
import sys, subprocess, numpy as np, pandas as pd
cohorts = sys.argv[1:]
BP = pd.read_csv('GDC-PANCAN_basic_phenotype.tsv', sep='\t')
pcol = [c for c in BP.columns if 'project' in c.lower()][0]; BP['proj'] = BP[pcol].astype(str).str.replace('TCGA-', '')
BM = pd.read_csv('/mnt/user-data/uploads/mart_export__1_.txt', sep='\t', low_memory=False)
pc = set(BM.loc[BM['Gene type'] == 'protein_coding', 'Gene stable ID'])
p = subprocess.Popen(['unzip', '-p', 'GDC-PANCAN.htseq_counts.tsv.zip', 'GDC-PANCAN.htseq_counts.tsv'], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
hdr = p.stdout.readline().rstrip('\n').split('\t')
samp2proj = dict(zip(BP['sample'], BP.proj))
cols = {c: [i for i, s in enumerate(hdr) if i > 0 and samp2proj.get(s) == c and s[13:15] in ('01', '03', '11')] for c in cohorts}
for c in cohorts: print(c, len(cols[c]), 'samples', flush=True)
allidx = np.array(sorted(set(i for v in cols.values() for i in v))); pos = {i: k for k, i in enumerate(allidx)}
genes, rows = [], []
for n, line in enumerate(p.stdout):
    g = line[:line.index('\t')].split('.')[0]
    if g not in pc: continue
    v = np.array(line.rstrip('\n').split('\t'), dtype=object)[allidx].astype(np.float32); genes.append(g); rows.append(v)
X = np.vstack(rows); print('genes kept', len(genes), flush=True)
for c in cohorts:
    sel = [pos[i] for i in cols[c]]
    np.savez_compressed(f'expr_{c}.npz', X=X[:, sel], genes=np.array(genes), samples=np.array([hdr[i] for i in cols[c]]))
print('done')
