"""Fig. 5: synthesis of how chromosome organisation constrains co-expression (schematic, figures/schematics/fig5.png)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import save, placeholder
import matplotlib.pyplot as plt
img = plt.imread(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'schematics', 'fig5.png'))
w_mm = 89; h_mm = w_mm * img.shape[0] / img.shape[1]
fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4)); ax = fig.add_axes([0, 0, 1, 1])
placeholder(ax, 'Schematic: chromosome organisation and co-expression', 'Within a domain, across a boundary, active intervening gene.', key='fig5')
save(fig, 'Fig5_model_summary', width=89 / 25.4); print('Fig5_model_summary ok')
