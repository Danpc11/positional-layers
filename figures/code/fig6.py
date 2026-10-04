import os
DATA = os.environ.get('POSLAYERS_DATA', 'data').rstrip('/') + '/'
import sys; sys.path.insert(0, os.path.dirname(__file__))
from style import *
import schem
fig = plt.figure(figsize=(W, W * 0.62)); gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.18, height_ratios=[1, 1.25])
ax = fig.add_subplot(gs[0, :]); schem.layers_diagram(ax); lab(ax, 'a', -0.02)
ax = fig.add_subplot(gs[1, 0]); schem.mechanism_diagram(ax); lab(ax, 'b', -0.04)
ax = fig.add_subplot(gs[1, 1]); schem.perturbation_rules(ax); lab(ax, 'c', -0.04)
save(fig, 'Fig6_model'); print('fig6 ok')
