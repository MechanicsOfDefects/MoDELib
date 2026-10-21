"""Reproducible manual illustrations; analytical curves are not simulation results."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
from mpl_toolkits.mplot3d.art3d import Line3DCollection
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
blue='#176b91'; orange='#d27a30'; dark='#203547'; pale='#edf4f8'
def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.png'),dpi=150,bbox_inches='tight')
    plt.close(fig)
def box(ax,x,y,w,h,text,color=pale):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02',fc=color,ec=blue,lw=1.3))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=10,color=dark)
def arrow(ax,a,b):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,color=blue,lw=1.5))
# Mesh illustration reads the actual repository fixture, no synthetic tetrahedra.
lines=(ROOT.parent/'Library/Meshes/unitCube24.msh').read_text().splitlines()
k=lines.index('$Nodes'); nn=int(lines[k+1]); pts={int(s.split()[0]):np.array(list(map(float,s.split()[1:4]))) for s in lines[k+2:k+2+nn]}
k=lines.index('$Elements'); ne=int(lines[k+1]); tets=[]
for s in lines[k+2:k+2+ne]:
    a=list(map(int,s.split()))
    if a[1]==4:tets.append(a[3+a[2]:])
edges=set()
for tet in tets:
    for i in range(4):
        for j in range(i):edges.add(tuple(sorted((tet[i],tet[j]))))
fig=plt.figure(figsize=(10,4.3))
for idx,title in enumerate(['Simplicial volume: 24 tetrahedra','Polyhedral boundary: cube faces']):
    ax=fig.add_subplot(1,2,idx+1,projection='3d')
    if idx==0: segs=[[pts[a],pts[b]] for a,b in edges]
    else:
        corners=[v for v in pts.values() if np.all(np.abs(v)==.5)]
        segs=[[a,b] for i,a in enumerate(corners) for b in corners[:i] if np.sum(a!=b)==1]
    ax.add_collection3d(Line3DCollection(segs,colors=blue,linewidths=.8 if idx==0 else 2))
    ax.set(xlim=(-.55,.55),ylim=(-.55,.55),zlim=(-.55,.55),title=title)
    ax.set_box_aspect((1,1,1));ax.set_axis_off();ax.view_init(22,32)
fig.text(.5,.02,'Same domain; volume interpolation and boundary traversal use different representations.',ha='center',fontsize=10)
save(fig,'mesh_views')
fig,ax=plt.subplots(figsize=(10,4));ax.axis('off');ax.set(xlim=(0,10),ylim=(0,4))
box(ax,.1,2.7,2.6,.9,'Material file\nscales, crystal type, laws')
box(ax,3.65,2.7,2.7,.9,'Polycrystal\nshared material + grains')
box(ax,7.2,2.7,2.6,.9,'Mesh regions\ngeometric partition')
arrow(ax,(2.7,3.15),(3.65,3.15));arrow(ax,(7.2,3.15),(6.35,3.15))
for x,txt in [(1.4,'Grain 1: C2G1\nlattice + slip systems'),(6,'Grain 2: C2G2\nlattice + slip systems')]:box(ax,x,.7,2.7,1,txt)
arrow(ax,(4.4,2.7),(2.75,1.7));arrow(ax,(5.6,2.7),(7.35,1.7));arrow(ax,(4.1,1.2),(6,1.2));ax.text(5.05,.55,'grain boundary',ha='center',fontsize=9)
save(fig,'polycrystal')
fig,(ax,bx)=plt.subplots(1,2,figsize=(10,3.6));x=np.linspace(0,1,300)
for y,l,c in [(2*(x-.5)*(x-1),'N₀',blue),(4*x*(1-x),'N₁/₂',orange),(2*x*(x-.5),'N₁','#668f52')]:ax.plot(x,y,label=l,color=c)
ax.set(title='Quadratic Lagrange edge basis',xlabel='Edge coordinate',ylabel='Shape function');ax.legend(frameon=False);ax.grid(alpha=.15)
bx.set(xlim=(0,5),ylim=(0,4));bx.axis('off')
for y,t in [(3,'Mesh + LagrangeElement'),(1.8,'TrialFunction + weak forms'),(.6,'WeakProblem: constraints + solve')]:box(bx,.2,y,4.6,.7,t)
arrow(bx,(2.5,3),(2.5,2.5));arrow(bx,(2.5,1.8),(2.5,1.3));fig.tight_layout();save(fig,'fem')
fig,ax=plt.subplots(figsize=(10,4.2));ax.set(xlim=(0,10),ylim=(0,4.5));ax.axis('off')
box(ax,.3,3.4,9.4,.7,'DefectiveCrystalBase: mesh, polycrystal, periodic lattice, optional FEM')
for x,t in [(.3,'DDD'),(2.7,'Elastic deformation'),(5.1,'Inclusions'),(7.5,'Cluster dynamics')]:
    box(ax,x,1.8,2.2,.8,t);arrow(ax,(x+1.1,3.4),(x+1.1,2.6))
box(ax,.3,.25,9.4,.8,'MicrostructureContainer: solve / reSolve → minimum dt → output → update')
for x in [.3,2.7,5.1,7.5]:arrow(ax,(x+1.1,1.8),(x+1.1,1.05))
save(fig,'multiphysics')

fig=plt.figure(figsize=(6,4));ax=fig.add_subplot(111,projection='3d');t=np.linspace(0,2*np.pi,120)
ax.plot(np.cos(t),np.sin(t),0*t,color=blue);ax.plot(np.cos(t),0*t,np.sin(t),color=blue,alpha=.3)
v=np.array([.55,.45,.7]);v/=np.linalg.norm(v)
for p,l in [(np.array([0,0,1.3]),'e₃ = R / R'),(v,'κ(θ, φ)')]:ax.quiver(0,0,0,*p,color=orange);ax.text(*p,l,fontsize=10)
ax.set_axis_off();ax.set_box_aspect((1,1,1));ax.view_init(20,30);save(fig,'fourier_sphere')
# Standalone plots place each illustration next to the physics that uses it.
fig,ax=plt.subplots(figsize=(7,2.8));a=np.logspace(-4,4,400)
ax.semilogx(a,1/(1+a),label='Applied stress weight',color=blue)
ax.semilogx(a,a/(1+a),label='Imposed strain weight',color=orange)
ax.set(xlabel='Scalar stiffness ratio α',ylabel='Weight');ax.legend(frameon=False,fontsize=10);ax.grid(alpha=.15);fig.tight_layout();save(fig,'controller_weights')
fig,ax=plt.subplots(figsize=(7,2.8));x=np.linspace(-3,3,500)
ax.plot(x,np.sqrt(x*x+.3**2),color=orange,label='a = 0.3');ax.plot(x,np.abs(x),color=blue,ls='--',label='a = 0')
ax.set(xlabel='Signed separation along a line',ylabel='Regularized distance Rₐ');ax.legend(frameon=False);ax.grid(alpha=.15);fig.tight_layout();save(fig,'regularized_distance')

# Regenerate the polycrystalline comparison used in Figure 1.1.
import runpy
runpy.run_path(str(ROOT / "scripts/create_polycrystal_figure.py"))

# Geometric entity hierarchies.
runpy.run_path(str(ROOT / "scripts/create_entity_figures.py"))

runpy.run_path(str(ROOT / "scripts/create_periodicity_figure.py"))
