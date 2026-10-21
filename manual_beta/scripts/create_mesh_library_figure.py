"""Render actual library mesh exterior triangles, colored by region ID."""
from pathlib import Path
from collections import defaultdict
from itertools import combinations
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
root=Path(__file__).resolve().parents[1]
items=[('unitCube24.msh','(a) Cube'),('unitCylinderAR2_2828.msh','(b) Cylinder'),('bicrystal_48.msh','(c) Bicrystal'),('poly10_periodicXYZ.msh','(d) Periodic polycrystal')]
fig=plt.figure(figsize=(10,7.4));plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
for panel,(name,title) in enumerate(items,1):
    lines=(root.parent/'Library/Meshes'/name).read_text().splitlines()
    k=lines.index('$Nodes');n=int(lines[k+1]);pts={int(t[0]):np.array(t[1:4],float) for t in (l.split() for l in lines[k+2:k+2+n])}
    k=lines.index('$Elements');tets=[]
    for l in lines[k+2:k+2+int(lines[k+1])]:
        a=list(map(int,l.split()))
        if a[1] in (4,11):tets.append((a[3+a[2]:3+a[2]+4],a[4]))
    faces=defaultdict(list)
    for tet,gid in tets:
        for f in combinations(tet,3):faces[tuple(sorted(f))].append(gid)
    ids=sorted(set(g for _,g in tets));palette=plt.get_cmap('tab10')
    xyz=np.array(list(pts.values()));center=(xyz.min(0)+xyz.max(0))/2;scale=np.max(np.ptp(xyz,axis=0))
    triangles=[];colors=[]
    for f,regions in faces.items():
        if len(regions)==1:
            triangles.append([(pts[v]-center)/scale for v in f]);colors.append(palette(ids.index(regions[0])%10))
    ax=fig.add_subplot(2,2,panel,projection='3d')
    ax.add_collection3d(Poly3DCollection(triangles,facecolors=colors,edgecolors='#253e50',linewidths=.22))
    ax.set(xlim=(-.55,.55),ylim=(-.55,.55),zlim=(-.55,.55));ax.set_box_aspect((1,1,1),zoom=1.2)
    ax.view_init(22,-55);ax.set_proj_type('ortho');ax.set_axis_off();ax.set_title(title,fontsize=12,pad=0)
    ax.text2D(.5,-.01,name,transform=ax.transAxes,ha='center',fontsize=9)
    ax.text2D(.5,-.08,f'{len(ids)} region'+('s' if len(ids)>1 else '')+f' | {len(tets):,} tetrahedra',transform=ax.transAxes,ha='center',fontsize=9,color='#203547')
    print(name,len(ids),len(tets))
fig.subplots_adjust(left=0,right=1,top=.95,bottom=.08,wspace=.02,hspace=.25)
for ext in ['pdf','png']:fig.savefig(root/'figures'/f'mesh_library_examples.{ext}',bbox_inches='tight',dpi=170)
