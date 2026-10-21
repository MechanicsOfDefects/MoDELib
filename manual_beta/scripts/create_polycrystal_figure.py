"""Figure 1.1: the same ten grains with and without coplanar triangle edges."""
from pathlib import Path
from collections import defaultdict
from itertools import combinations
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
root=Path(__file__).resolve().parents[1]
source=root.parent/'Library/Meshes/polycrystalCube_10_periodicX.msh'
s=source.read_text().splitlines();k=s.index('$Nodes');n=int(s[k+1])
points={int(t[0]):np.array(t[1:4],float) for t in (x.split() for x in s[k+2:k+2+n])}
k=s.index('$Elements');grains=defaultdict(list)
for line in s[k+2:k+2+int(s[k+1])]:
    a=list(map(int,line.split()))
    if a[1] in (4,11):grains[a[3]].append(a[3+a[2]:3+a[2]+4])
assert len(grains)==10
allp=np.array(list(points.values()));center=(allp.max(0)+allp.min(0))/2;scale=np.max(np.ptp(allp,axis=0))
points={k:(p-center)/scale for k,p in points.items()}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
fig=plt.figure(figsize=(10,4.8));colors=plt.get_cmap('tab10')
for panel,title in enumerate(['Simplicial mesh: volume triangulation','Polyhedral mesh: grains and grain boundaries']):
    ax=fig.add_subplot(1,2,panel+1,projection='3d')
    polygons=[]; facecolors=[]
    for j,(gid,tets) in enumerate(sorted(grains.items())):
        faces=defaultdict(int)
        for tet in tets:
            for f in combinations(tet,3):faces[tuple(sorted(f))]+=1
        boundary=[f for f,count in faces.items() if count==1]
        ids=set(v for f in boundary for v in f)
        shift=.22*np.mean([points[v] for v in ids],axis=0)
        tris=np.array([[points[v]+shift for v in f] for f in boundary])
        color=colors(j)
        if panel==0:
            polygons.extend(tris);facecolors.extend([color]*len(tris))
        else:
            planes=defaultdict(list)
            for f in boundary:
                p=np.array([points[v] for v in f]);normal=np.cross(p[1]-p[0],p[2]-p[0]);normal/=np.linalg.norm(normal)
                if normal[np.argmax(np.abs(normal))]<0:normal=-normal
                key=tuple(np.round(np.r_[normal,np.dot(normal,p[0])],6));planes[key].append(f)
            for fs in planes.values():
                counts=defaultdict(int)
                for f in fs:
                    for e in combinations(f,2):counts[tuple(sorted(e))]+=1
                adjacency=defaultdict(list)
                for (a,b),count in counts.items():
                    if count==1:adjacency[a].append(b);adjacency[b].append(a)
                assert all(len(v)==2 for v in adjacency.values())
                first=next(iter(adjacency));loop=[first];prev=None;cur=first
                while True:
                    nxt=next(v for v in adjacency[cur] if v!=prev)
                    if nxt==first:break
                    loop.append(nxt);prev,cur=cur,nxt
                assert len(loop)==len(adjacency)
                polygons.append([points[v]+shift for v in loop]);facecolors.append(color)
    ax.add_collection3d(Poly3DCollection(polygons,facecolors=facecolors,edgecolors='#344450',linewidths=.24 if panel==0 else .65,alpha=1,zsort='average'))
    ax.set(xlim=(-.67,.67),ylim=(-.67,.67),zlim=(-.67,.67),title=title)
    ax.set_box_aspect((1,1,1),zoom=1.5);ax.view_init(22,-55);ax.set_proj_type('ortho');ax.set_axis_off()
fig.subplots_adjust(left=0,right=1,top=.9,bottom=.12,wspace=.02)
fig.text(.5,.055,'10 grains • Matching colors identify the same grain • Grains separated for visibility',ha='center',fontsize=10,color='#203547')
for ext in ['pdf','png']:fig.savefig(root/'figures'/('polycrystal_mesh_views.'+ext),dpi=180,bbox_inches='tight')
print('Rendered',len(grains),'grains and',sum(map(len,grains.values())),'tetrahedra from',source.name)
