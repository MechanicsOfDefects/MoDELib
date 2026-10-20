import csv
import os
import sys
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'python'))
from readF import readFfile, getFarray
from readFile import get_scalar

# Preserve the original Carpenter-Rogerson data, exported with every test run.
folder = str(Path(sys.argv[1]).resolve())
series = {name: [] for name in ('EXPIa', 'EXPIc', 'EXPZa', 'EXPZc')}
with open(Path(folder) / 'experimental-data.csv') as stream:
    for row in csv.DictReader(stream):
        series[row['series']].append([float(row['dose_dpa']), float(row['strain_in_1e_minus_4'])])
EXPIa, EXPIc, EXPZa, EXPZc = (np.asarray(series[name]) for name in series)

# ------------------------------------------------------------------------ #
#                          Main Inputs
# ------------------------------------------------------------------------ #

textfont = 18
axisfont = 16

# ------------------------------------------------------------------------ #
#                          Main Processing Loop
# ------------------------------------------------------------------------ #
print(f"Processing {folder}")
matfile = os.path.join(folder, 'inputFiles', 'Zr4_Fitted.txt')
F, Flabels = readFfile(folder)
t = getFarray(F, Flabels, 'time [b/cs]')
DoseRate = get_scalar(matfile, 'doseRate_dpaPerSec')
mu0_SI = get_scalar(matfile, 'mu0_SI')
rho_SI = get_scalar(matfile, 'rho_SI')
b_SI = get_scalar(matfile, 'b_SI')
v_dd2SI = np.sqrt(mu0_SI / rho_SI)
t_dd2SI = b_SI / v_dd2SI
dpa = t * t_dd2SI * DoseRate  # <-- missing before

# Plotting strain components and swelling
b11 = getFarray(F, Flabels, 'betaP_11')
b22 = getFarray(F, Flabels, 'betaP_22')
b33 = getFarray(F, Flabels, 'betaP_33')
b12 = getFarray(F, Flabels, 'betaP_12')
b13 = getFarray(F, Flabels, 'betaP_13')
b23 = getFarray(F, Flabels, 'betaP_23')

swelling_trace_over_3 = []
for i in range(len(t)):
    betaP = np.array([
        [b11[i], b12[i], b13[i]],
        [b12[i], b22[i], b23[i]],
        [b13[i], b23[i], b33[i]],
    ])
    trace = np.trace(betaP)
    swelling_trace_over_3.append(trace / 3)
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(dpa+0.5, b11*100, label=r'$\beta_{11}$', color='red')
ax.plot(dpa+0.5, b22*100, label=r'$\beta_{22}$', color='green')
ax.plot(dpa+0.5, b33*100, label=r'$\beta_{33}$', color='blue')
ax.plot(dpa+0.5, swelling_trace_over_3, label=r'Swelling ($\mathrm{tr}(\beta^P)/3$)', color='gray', linestyle='-')


ax.set_xlim(0,9.2)
ax.set_ylim(-0.08,0.06)

# Experimental Data
ax.scatter(EXPIa[:,0],EXPIa[:,1]*10**-2,color="firebrick",marker="s",s=80,label=r'Experiment: Iodide Purity');
ax.scatter(EXPZa[:,0],EXPZa[:,1]*10**-2,color="orange",marker="s",s=80,label=r'Experiment: Zone Refined Purity');
ax.scatter(EXPIc[:,0],EXPIc[:,1]*10**-2,color="firebrick",marker="s",s=80);
ax.scatter(EXPZc[:,0],EXPZc[:,1]*10**-2,color="orange",marker="s",s=80);


ax.set_xlabel('Irradiation Dose [dpa]', fontsize=textfont)
ax.set_ylabel('Strain [%]', fontsize=textfont)
ax.legend(fontsize=axisfont)
ax.grid(True, linestyle=':')
ax.set_xlim(left=0)
ax.axhline(0, color='k', linewidth=0.5)
ax.tick_params(axis='both', which='major', direction='in', length=6, labelsize=axisfont)
fig.tight_layout()



# # Plot dislocation density
# glissile_density = getFarray(F, Flabels, 'glissile density [m^-2]')
# sessile_density  = getFarray(F, Flabels, 'sessile density [m^-2]')
# fig2, ax2 = plt.subplots(figsize=(10, 6))
# ax2.plot(dpa, glissile_density+sessile_density, label='Dislocation Density', color='k', linewidth=2)
# ax2.set_xlabel('Irradiation Dose [dpa]', fontsize=textfont)
# ax2.set_ylabel('Dislocation Density [m$^{-2}$]', fontsize=textfont)
# ax2.legend(fontsize=axisfont)
# ax2.grid(True, linestyle=':')
# ax2.set_xlim(left=0)
# ax2.tick_params(axis='both', which='major', direction='in', length=6, labelsize=axisfont)
# fig2.tight_layout()

# Headless output for CTest.
# Save figures
save_path = os.path.join(folder, 'StrainComponentsWithSwelling.pdf')
fig.savefig(save_path)
fig.savefig(os.path.join(folder, "StrainComponentsWithSwelling.png"), dpi=150)
# save_path2 = os.path.join(folder, 'DislocationDensities_vs_DPA.pdf')
# fig2.savefig(save_path2)   
plt.close(fig)
print(f"Saved {save_path}")
# plt.close(fig2)
# print(f"Saved {save_path2}")

