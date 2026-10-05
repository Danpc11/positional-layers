import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import W, lab, save, placeholder
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(W, W * 0.55)); gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.12, height_ratios=[1, 1.2])
ax = fig.add_subplot(gs[0, :]); placeholder(ax, 'Schematic: what the expression spectrum contains', 'Left to right: tissue mean landscape (gene-order independent) | GC-associated layer (per-sample GC slope × gene GC) |\ncopy-number dosage (tumours) | cis covariance (shared regulatory sources). Only the last depends on genome organization.'); lab(ax, 'a', -0.02)
ax = fig.add_subplot(gs[1, 0]); placeholder(ax, 'Schematic: shared regulatory sources', 'Sources (enhancers, genotype, copy-number segments) vary across individuals.\nNearby genes that share a source covary; a shared eQTL is a natural perturbation\nof one source and sets the sign of the coupling.'); lab(ax, 'b', -0.04)
ax = fig.add_subplot(gs[1, 1]); placeholder(ax, 'Schematic: architecture sets which genes share sources', 'Coupling is higher within a domain and decreases across each\nconvergent CTCF loop anchor and across an active intervening gene.\nSTAG2 mutations lower coupling between adjacent genes.'); lab(ax, 'c', -0.04)
save(fig, 'Fig5_model'); print('Fig5_model ok')
