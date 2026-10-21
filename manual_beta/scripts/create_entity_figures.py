"""Geometric sketches of the 3D polyhedral and simplex hierarchies."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch

OUT = Path(__file__).resolve().parents[1] / 'figures'
BLUE, ORANGE, DARK = '#176b91', '#d27a30', '#203547'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42})

def polygon(ax, points, fill='#edf4f8', z=2):
    ax.add_patch(Polygon(points, closed=True, facecolor=fill, edgecolor=BLUE, linewidth=1.5, zorder=z))

def line(ax, points, color=BLUE, dashed=False, z=4):
    ax.plot(*np.asarray(points).T, color=color, lw=2, ls='--' if dashed else '-', zorder=z)

def arrow(ax, start, end, color=BLUE):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=12, color=color, lw=1.3))

def make(simplex=False):
    fig, ax = plt.subplots(figsize=(10, 2.7))
    ax.set(xlim=(-1.1, 9.1), ylim=(-.95, 1.6), aspect='equal')
    ax.axis('off')
    centers = [0, 2.65, 5.3, 7.95]
    if simplex:
        p = np.array([[-.75,-.25],[.65,-.25],[.32,1.05],[-.05,.2]])
        polygon(ax, p[[0,1,2]])
        for i in [0,1,2]:
            line(ax, p[[3,i]], dashed=True, z=3)
        polygon(ax, p[[0,1,2]], '#fae3cc', z=3.5)
        # Dashed edges meet at the hidden fourth vertex.
        for i in [0,1,2]:
            line(ax, p[[3,i]], dashed=True, z=4)
        face = p[[0,1,2]] + [centers[1],0]
        labels = ['Tetrahedron', 'Triangle', 'Edge', 'Vertex']
        detail = ['4 vertices', '3 vertices', '2 vertices', '1 vertex']
        title = 'Simplex<3, order>: boundary entities'
    else:
        p = np.array([[-.75,.05],[-.4,-.25],[.4,-.25],[.75,.05],[.4,.48],[-.4,.48]])
        back = p + [0,.65]
        polygon(ax, back)
        for i in range(3):
            polygon(ax, np.array([p[i],p[i+1],back[i+1],back[i]]), '#d9e9f1')
        polygon(ax, p, '#fae3cc', z=4)
        face = p + [centers[1],.15]
        labels = ['Region', 'Face', 'Edge', 'Vertex']
        detail = ['region ID', 'planar polygon', 'shared by faces', 'position']
        title = 'Polyhedron<3, order>: boundary entities'
    polygon(ax, face)
    line(ax, face[:2], ORANGE)
    edge = np.array([[-.65,.2],[.65,.2]]) + [centers[2],0]
    line(ax, edge)
    ax.scatter(*edge.T, s=25, color=BLUE, zorder=5)
    ax.scatter(*edge[1], s=40, color=ORANGE, zorder=6)
    ax.scatter(centers[3],.2,s=50,color=ORANGE,zorder=6)
    for i, x in enumerate(centers):
        ax.text(x,-.53, f'{labels[i]} ({3-i}-cell)', ha='center', color=DARK, weight='bold')
        ax.text(x,-.8, detail[i], ha='center', color=DARK, fontsize=10)
    for x in centers[:-1]:
        arrow(ax, (x+.92,.45), (x+1.66,.45))
        if not simplex:
            arrow(ax, (x+1.66,.08), (x+.92,.08), ORANGE)
    ax.text(3.95,1.4,title,ha='center',color=DARK,fontsize=12)
    name = 'simplex_complex' if simplex else 'cell_complex'
    for ext in ['pdf','png']:
        fig.savefig(OUT / f'{name}.{ext}', bbox_inches='tight', dpi=180)
    plt.close(fig)

make()
make(simplex=True)
