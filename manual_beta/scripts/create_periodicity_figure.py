"""Translation lattice, affine field and periodically tiled polycrystal schematic."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
OUT=Path(__file__).resolve().parents[1]/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
blue='#176b91'; orange='#d27a30'; dark='#203547'
B=np.array([[1,.35],[0,1.]])

def arrow(ax,a,b,color=blue):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,color=color,lw=1.4,zorder=10))

# Construct whole Voronoi grains for the central seed representatives.
# Their union is a jagged fundamental unit; translated units tile the plane.
# No grain is clipped at the lattice parallelogram.
seeds=np.array([[-.32,-.3],[.2,-.32],[-.2,.17],[.31,.18],[.04,.42]])
images=[(p+[i,j],k) for i in range(-1,2) for j in range(-1,2) for k,p in enumerate(seeds)]
square=np.array([[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]])
def clip(poly,n,c):
    result=[]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        da,db=a@n-c,b@n-c
        if da<=1e-10:result.append(a)
        if (da<0)!=(db<0):result.append(a+(b-a)*da/(da-db))
    return np.array(result)
grains=[]
for k,p in enumerate(seeds):
    poly=4*square.copy()
    for q,_ in images:
        if np.array_equal(p,q):continue
        poly=clip(poly,q-p,(q@q-p@p)/2)
        if len(poly)<3:break
    if len(poly)>=3:grains.append((poly,k))
# Shared grain edges are internal; retain only the unit's external perimeter.
from collections import Counter
edge_counts=Counter()
edge_points={}
for poly,_ in grains:
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        key=tuple(sorted((tuple(np.round(a,9)),tuple(np.round(b,9)))))
        edge_counts[key]+=1
        edge_points[key]=np.array([a,b])
perimeter=[edge_points[e] for e,n in edge_counts.items() if n==1]
# Whole grains cover one fundamental-cell area without artificial cut faces.
area=sum(abs(np.dot(p[:,0],np.roll(p[:,1],-1))-np.dot(p[:,1],np.roll(p[:,0],-1)))/2 for p,_ in grains)
assert np.isclose(area,1), area
assert all(n in (1,2) for n in edge_counts.values())
colors=plt.get_cmap('Set2').colors

def cell(ax,i,j,alpha=1):
    for poly,k in grains:
        ax.add_patch(Polygon((poly+[i,j])@B.T,fc=colors[k],ec='white',lw=.45,alpha=alpha,zorder=3))
    for edge in perimeter:
        ax.plot(*((edge+[i,j])@B.T).T,color=blue,lw=1.3 if i==j==0 else .5,alpha=alpha,zorder=6 if i==j==0 else 4)

def lattice(ax,filled):
    for i in range(-3,4):
        for j in range(-3,4):
            if filled:cell(ax,i,j,float(np.exp(-.28*(i*i+j*j))))
            pt=B@np.array([i,j]);ax.plot(*pt,'.',color=blue,ms=2,alpha=.35,zorder=5)
    for i in range(-3,4):
        for ends in [np.array([[i,-3],[i,3]]),np.array([[-3,i],[3,i]])]:
            ax.plot(*(ends@B.T).T,color=blue,lw=.5,alpha=.18,zorder=1)
    if not filled:cell(ax,0,0)
    ax.plot(0,0,'o',color=dark,ms=3,zorder=12)
    ax.annotate('(0, 0)',(0,0),xytext=(5,-14),textcoords='offset points',fontsize=9,zorder=15,
                bbox=dict(fc='white',ec='none',alpha=.85,pad=1))
    for v,label,offset in [(B[:,0],r'$\mathbf{a}_1$',(0,.12)),(B[:,1],r'$\mathbf{a}_2$',(-.25,0))]:
        arrow(ax,np.zeros(2),v);ax.text(*(v+offset),label,color=blue,zorder=15)
    ax.set(xlim=(-2.65,2.65),ylim=(-2.1,2.1),aspect='equal');ax.axis('off')

fig,axs=plt.subplots(1,3,figsize=(12,3.7))
ax,cx,dx=axs.flat
fieldfig,bx=plt.subplots(figsize=(6,3.5))
for i in range(-1,2):
    for j in range(-1,2):
        q=(np.array([[0,0],[1,0],[1,1],[0,1]])+[i,j])@B.T
        ax.add_patch(Polygon(q,fc='#edf4f8' if i==j==0 else 'white',ec=blue,alpha=1 if i==j==0 else .3))
for v,t in [(B[:,0],r'$\mathbf{a}_1$'),(B[:,1],r'$\mathbf{a}_2$')]:
    arrow(ax,np.zeros(2),v);ax.text(*(v*.6+[0,.1]),t,color=blue)
p=np.array([.2,.55]);ax.scatter(*p,c=orange);ax.scatter(*(p+B[:,0]),c=orange);arrow(ax,p,p+B[:,0]);ax.text(.65,.78,'equivalent points',ha='center',fontsize=9)
ax.set(xlim=(-.5,2),ylim=(-.4,1.5),aspect='equal');ax.axis('off')
x=np.linspace(0,1,300);v=.12*np.sin(2*np.pi*x)
bx.plot(x,v,label='periodic fluctuation',color=blue);bx.plot(x,.5*x+v,label='affine + periodic field',color=orange)
bx.scatter([0,1],[0,.5],c=orange);bx.set(xlabel='x / L',ylabel='Displacement (illustrative units)');bx.legend(frameon=False,fontsize=9);bx.grid(alpha=.15)
lattice(cx,False);lattice(dx,True)
for a,t in zip(axs.flat,['(a) Translation lattice','(b) Polycrystal unit at the origin','(c) Space filled by translated units']):a.set_title(t,fontsize=11,pad=12)
fig.tight_layout(h_pad=2.5,w_pad=2)
for ext in ['pdf','png']:fig.savefig(OUT/f'periodicity.{ext}',bbox_inches='tight',dpi=180)

bx.set_title('Affine jump + periodic fluctuation',fontsize=11)
fieldfig.tight_layout()
for ext in ['pdf','png']:fieldfig.savefig(OUT/f'affine_periodic_field.{ext}',bbox_inches='tight',dpi=180)
