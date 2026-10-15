import numpy as np
import sys, os
from matplotlib import pyplot as plt
sys.path.append('../../python/')

from readF import *
from readFile import *
from readEVL import *
import pathlib

# -- Experimental Data Carpenter Rogerson -- #
EXPIa=np.array([[0.1245,1.3249],[0.3034,1.1033],[0.6301,1.163],[0.8557,1.0831],
    [1.3613,1.4861],[1.5946,1.3451],[1.7191,1.4055],[2.5514,1.3048],[2.7381,1.2645],
    [3.7105,1.6474],[4.4417,2.0101],[4.6283,2.5945],[5.7174,2.9572],
    [6.0596,2.9773],[6.1141,3.1184],[6.6197,3.4207],[7.2887,4.3073],[8.1677,5.1537]])  #IODIDE PURITY a # x dpa y strain units e-4
EXPIc=np.array([[ 0.1400,-0.1461],[0.3034,-0.4081],[0.6379, -0.4081],[1.2368,-0.7103],
    [1.5013,-0.8111],[1.6180,-0.6902],[2.1936,-0.6902],[2.6059,-1.3149],
    [3.1037,-1.1537],[3.5860,-1.2544],[4.2083,-1.0126]]) #IODIDE PURITY c # x dpa y strain units e-4
EXPZa=np.array([[0.1245,0.5189],[0.2723,0.7406],[0.4745,0.8816],[0.6845,1.1234],[1.0346,1.4660],[1.2757,1.5063],
    [1.5013,1.2645],[1.6024,1.0428],[1.9758,1.4660],[2.3492,1.2645],[3.1271,1.7884],
    [3.3293,1.8690],[3.6482,2.1310],[3.9594,2.1713],[4.4261,2.3325],[4.6283,2.1713],
    [5.6940,2.4332],[5.9507,2.9370], [6.4719,2.8766],[6.6430,3.1385],[6.7753,3.1788]])  #Zone Refined a # x dpa y strain units e-4
EXPZc=np.array([[0.1478,-0.4282],[0.2800,-0.2267],[0.6145,-0.7909],[0.8634,-0.7506],[1.3613,-0.7506],
    [1.4935,-0.6297],[1.8124,-0.7708],[2.1625,-0.8917],[2.3647,-0.7103],[3.4693,-1.0126]]) #Zone Refined c # x dpa y strain units e-4

# ------------------------------------------------------------------------ #
#                          Main Inputs
# ------------------------------------------------------------------------ #
folder = str(pathlib.Path().resolve()) # simulation folder
textfont = 18
axisfont = 16
interiorDistance_nm = 300.0 # nodes farther than this from the grain boundary are "interior"
profileStep = -1            # step of the profiles and histograms (-1: the last one)

# ------------------------------------------------------------------------ #
#                          Material and units
# ------------------------------------------------------------------------ #
print(f"Processing {folder}")
matfile = os.path.join(folder, 'inputFiles', 'Zr_CD4opt.txt')
DoseRate = get_scalar(matfile, 'doseRate_dpaPerSec')
mu0_SI = get_scalar(matfile, 'mu0_SI')
rho_SI = get_scalar(matfile, 'rho_SI')
b_SI = get_scalar(matfile, 'b_SI')
t_dd2SI = b_SI / np.sqrt(mu0_SI / rho_SI) # time unit b/cs in seconds
mSize = len(get_vector(matfile, 'mobileSpeciesVector'))
polarity = get_vector(matfile, 'immobileSpeciesVector')
nF = len(polarity)
burgers = get_matrix(matfile, 'immobileSpeciesBurgers') # columns, in lattice components
mMin = get_vector(matfile, 'shapeTransitionMinSize')
mMax = get_vector(matfile, 'shapeTransitionMaxSize')
w0 = get_scalar(matfile, 'shapeTransitionWidth')

# HEX lattice of MoDELib, in units of b: a1, a2, c (ideal c/a)
lattice = np.array([[1.0, 0.5, 0.0], [0.0, np.sqrt(3.0)/2.0, 0.0], [0.0, 0.0, np.sqrt(8.0/3.0)]])
omega = np.linalg.det(lattice) / 2.0 # atomic volume, in b^3 (two atoms per cell)
bMag = np.linalg.norm(lattice @ burgers, axis=0) # Burgers vector magnitudes, in b
names = ['c' if abs(burgers[2,k]) > 0 else 'a%d' % k for k in range(nF)]
mobileNames = ['v', 'i', '2i', '3i'][:mSize]

def diameter_nm(n, c):
    # mean diameter of the clusters: twice the effective radius R=(1-S)R_pyr+S*R_loop (report Eq. 30)
    # n: number density in b^-3 (as in the evl files), c: defects per atom
    m = np.where(n > 0, c / np.maximum(n, 1e-300) / omega, np.nan) # defects per cluster (nan where there are no loops)
    S = 1.0 / (1.0 + np.exp(-(m - 0.5*(mMin + mMax)) / (w0 * (mMax - mMin))))
    R = (1.0 - S) * (m * omega / np.sqrt(8.0))**(1.0/3.0) + S * np.sqrt(m * omega / (np.pi * bMag))
    return 2.0 * R * b_SI * 1e9

# ------------------------------------------------------------------------ #
#                          Read the output
# ------------------------------------------------------------------------ #
F, Flabels = readFfile(folder)
runIDs = getFarray(F, Flabels, 'runID').astype(int)
dpa = getFarray(F, Flabels, 'time [b/cs]') * t_dd2SI * DoseRate

X = readCDNodes(folder) # positions of the nodes, in b
lo, hi = X.min(axis=0), X.max(axis=0)
distance_nm = np.minimum(X - lo, hi - X).min(axis=1) * b_SI * 1e9 # distance from the nearest face (grain boundary)
interior = distance_nm > interiorDistance_nm

mobile, density, content = [], [], [] # interior means at each step
for runID in runIDs:
    cd = readEVLtxt(os.path.join(folder, 'evl', 'evl_%d' % runID)).clusterDynamics
    mobile.append(cd[interior, :mSize].mean(axis=0))
    density.append(cd[interior, mSize:mSize+nF].mean(axis=0) / b_SI**3) # [m^-3]
    content.append(cd[interior, mSize+nF:].mean(axis=0))
mobile, density, content = np.array(mobile), np.array(density), np.array(content)
diameter = diameter_nm(density * b_SI**3, content) # from the mean density and content

# ------------------------------------------------------------------------ #
#                          Figure 1: growth strain
# ------------------------------------------------------------------------ #
b11 = getFarray(F, Flabels, 'betaP_11')
b22 = getFarray(F, Flabels, 'betaP_22')
b33 = getFarray(F, Flabels, 'betaP_33')
swelling = (b11 + b22 + b33) / 3

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(dpa, b11*100, label=r'$\beta_{11}$', color='red')
ax.plot(dpa, b22*100, label=r'$\beta_{22}$', color='green', linestyle='--')
ax.plot(dpa, b33*100, label=r'$\beta_{33}$', color='blue')
ax.plot(dpa, swelling*100, label=r'Swelling ($\mathrm{tr}(\beta^P)/3$)', color='gray')
ax.scatter(EXPIa[:,0],EXPIa[:,1]*10**-2,color="firebrick",marker="s",s=80,label=r'Experiment: Iodide Purity')
ax.scatter(EXPZa[:,0],EXPZa[:,1]*10**-2,color="orange",marker="s",s=80,label=r'Experiment: Zone Refined Purity')
ax.scatter(EXPIc[:,0],EXPIc[:,1]*10**-2,color="firebrick",marker="s",s=80)
ax.scatter(EXPZc[:,0],EXPZc[:,1]*10**-2,color="orange",marker="s",s=80)
ax.set_xlabel('Irradiation Dose [dpa]', fontsize=textfont)
ax.set_ylabel('Strain [%]', fontsize=textfont)
ax.legend(fontsize=axisfont)
ax.grid(True, linestyle=':')
ax.set_xlim(left=0)
ax.axhline(0, color='k', linewidth=0.5)
ax.tick_params(axis='both', which='major', direction='in', length=6, labelsize=axisfont)
fig.tight_layout()
fig.savefig(os.path.join(folder, 'StrainComponentsWithSwelling.pdf'))

# ------------------------------------------------------------------------ #
#                          Figure 2: interior means against dose
# ------------------------------------------------------------------------ #
fig2, ax2 = plt.subplots(2, 2, figsize=(14, 10))
for k in range(mSize):
    ax2[0,0].semilogy(dpa, mobile[:,k], label=mobileNames[k])
ax2[0,0].set_ylabel('Mobile content [defects/atom]', fontsize=textfont)
for k in range(nF):
    ax2[0,1].semilogy(dpa, density[:,k], label=names[k])
    ax2[1,0].plot(dpa, diameter[:,k], label=names[k])
    ax2[1,1].plot(dpa, content[:,k], label=names[k])
ax2[0,1].set_ylabel(r'Loop density [m$^{-3}$]', fontsize=textfont)
ax2[1,0].set_ylabel('Mean loop diameter [nm]', fontsize=textfont)
ax2[1,1].set_ylabel('Defects in loops [defects/atom]', fontsize=textfont)
for a in ax2.flat:
    a.set_xlabel('Irradiation Dose [dpa]', fontsize=textfont)
    a.legend(fontsize=axisfont)
    a.grid(True, linestyle=':')
    a.tick_params(axis='both', which='major', direction='in', length=6, labelsize=axisfont)
fig2.suptitle('Interior of the grain (farther than %g nm from the boundary)' % interiorDistance_nm, fontsize=textfont)
fig2.tight_layout()
fig2.savefig(os.path.join(folder, 'InteriorMeans_vs_DPA.pdf'))

# ------------------------------------------------------------------------ #
#                          Figure 3: profiles from the grain boundary
# ------------------------------------------------------------------------ #
runID = runIDs[profileStep]
cd = readEVLtxt(os.path.join(folder, 'evl', 'evl_%d' % runID)).clusterDynamics
edges = np.linspace(0.0, distance_nm.max(), 31)
centers = 0.5 * (edges[1:] + edges[:-1])
bins = np.clip(np.digitize(distance_nm, edges) - 1, 0, len(centers) - 1)
def binMean(v): # mean of each column over the nodes of each distance bin
    return np.array([v[bins == j].mean(axis=0) if np.sum(bins == j) >= 20 else np.full(v.shape[1], np.nan) for j in range(len(centers))])
pMobile = binMean(cd[:, :mSize])
pDensity = binMean(cd[:, mSize:mSize+nF])
pContent = binMean(cd[:, mSize+nF:])
pDiameter = diameter_nm(pDensity, pContent) # from the binned density and content, not the mean of nodal diameters

fig3, ax3 = plt.subplots(2, 2, figsize=(14, 10))
for k in range(mSize):
    ax3[0,0].semilogy(centers, pMobile[:,k], label=mobileNames[k])
ax3[0,0].set_ylabel('Mobile content [defects/atom]', fontsize=textfont)
ax3[0,0].set_ylim(bottom=1e-16) # the boundary value (thermal equilibrium) is far lower
for k in range(nF):
    ax3[0,1].semilogy(centers, pDensity[:,k] / b_SI**3, label=names[k])
    ax3[1,0].plot(centers, pDiameter[:,k], label=names[k])
    ax3[1,1].plot(centers, pContent[:,k], label=names[k])
ax3[0,1].set_ylabel(r'Loop density [m$^{-3}$]', fontsize=textfont)
ax3[1,0].set_ylabel('Mean loop diameter [nm]', fontsize=textfont)
ax3[1,1].set_ylabel('Defects in loops [defects/atom]', fontsize=textfont)
for a in ax3.flat:
    a.set_xlabel('Distance from the grain boundary [nm]', fontsize=textfont)
    a.legend(fontsize=axisfont)
    a.grid(True, linestyle=':')
    a.tick_params(axis='both', which='major', direction='in', length=6, labelsize=axisfont)
fig3.suptitle('Profiles at %.2f dpa' % dpa[profileStep], fontsize=textfont)
fig3.tight_layout()
fig3.savefig(os.path.join(folder, 'Profiles_vs_Distance.pdf'))

# ------------------------------------------------------------------------ #
#                          Figure 4: distribution of the loop diameters
# ------------------------------------------------------------------------ #
# Each node carries one mean size per family: the spread is spatial (boundary vs interior), not a size distribution.
# The mesh is finer near the boundary, so the near-boundary sizes count more nodes.
nodalDiameter = diameter_nm(cd[:, mSize:mSize+nF], cd[:, mSize+nF:])
fig4, ax4 = plt.subplots(1, 2, figsize=(14, 6))
for j, k in enumerate([0, 1]):
    d = nodalDiameter[:,k]
    ax4[j].hist(d[np.isfinite(d)], bins=50, color='steelblue')
    ax4[j].set_xlabel('Mean loop diameter at the node [nm]', fontsize=textfont)
    ax4[j].set_ylabel('Number of nodes', fontsize=textfont)
    ax4[j].set_title('%s loops, %.2f dpa' % (names[k], dpa[profileStep]), fontsize=textfont)
    ax4[j].tick_params(axis='both', which='major', direction='in', length=6, labelsize=axisfont)
fig4.tight_layout()
fig4.savefig(os.path.join(folder, 'DiameterDistribution.pdf'))

print('Saved StrainComponentsWithSwelling.pdf, InteriorMeans_vs_DPA.pdf, Profiles_vs_Distance.pdf, DiameterDistribution.pdf')
plt.show()
