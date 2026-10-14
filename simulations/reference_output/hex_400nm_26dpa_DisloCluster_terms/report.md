# 20261002_232731_4b9e666-dirty_hex400_DisloCluster_terms

- material: `Zr_CD4opt.txt` at 573 K, 1e-07 dpa/s
- geometry: `hex_400nm` (hexagonal), 400 x 346 x 653 nm, every face a grain boundary
- integration of the immobile species: **cvode**
- snapshots: 7, 0 .. 26 dpa
- finite-element nodes: 32317 (129268 mobile and 258536 immobile unknowns)
- status: completed; wall time 2347 s

## Verification

| Check | Expected | Result | |
|---|---|---|---|
| no negative field | minimum >= 0 | minimum 1.30e-70 | PASS |
| the three <a> families are equal without stress | relative difference < 1e-8 | 9.1e-13 | PASS |
| <c> density in the interior against production times lifetime | <= 1.493e+21 m^-3, equal while the loops do not coalesce | 1.481e+21 m^-3 at 26 dpa | PASS |
| <c> density at the first dose beyond 20 lifetimes | 1.493e+21 m^-3 within 5 % | 1.481e+21 m^-3 at 1 dpa | PASS |
| vacancies at the grain boundary at thermal equilibrium | exp(-Ef/kT) = 8.459e-15 | largest relative difference 1.1e-05 | PASS |

## Intervals

| dose from | dose to | steps | wall |
|---:|---:|---:|---:|
| 0 | 1 | 1 | 174 s |
| 1 | 6 | 5 | 467 s |
| 6 | 11 | 5 | 428 s |
| 11 | 16 | 5 | 421 s |
| 16 | 21 | 5 | 426 s |
| 21 | 26 | 5 | 426 s |

## Values in the interior of the grain

Means over the nodes farther than 80 % of the half size from the grain boundary.

| dose [dpa] | Cv | Ci | C2i | C3i | N_c | N_a1 | N_a2 | N_a3 | C_c | C_a1 | C_a2 | C_a3 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 5.252e-06 | 6.359e-10 | 1.640e-11 | 1.059e-11 | 2.959e-05 | 2.959e-05 | 2.959e-05 | 2.959e-05 | 1.000e-33 | 1.000e-33 | 1.000e-33 | 1.000e-33 |
| 1 | 4.877e-06 | 5.883e-10 | 1.026e-11 | 6.644e-12 | 1.481e+21 | 1.701e+21 | 1.701e+21 | 1.701e+21 | 1.923e-03 | 3.887e-05 | 3.887e-05 | 3.887e-05 |
| 6 | 4.862e-06 | 5.871e-10 | 1.007e-11 | 6.522e-12 | 1.481e+21 | 1.798e+21 | 1.798e+21 | 1.798e+21 | 2.170e-03 | 3.730e-05 | 3.730e-05 | 3.730e-05 |
| 11 | 4.862e-06 | 5.869e-10 | 1.002e-11 | 6.494e-12 | 1.481e+21 | 1.801e+21 | 1.801e+21 | 1.801e+21 | 2.174e-03 | 3.729e-05 | 3.729e-05 | 3.729e-05 |
| 16 | 4.861e-06 | 5.868e-10 | 9.995e-12 | 6.475e-12 | 1.481e+21 | 1.801e+21 | 1.801e+21 | 1.801e+21 | 2.176e-03 | 3.727e-05 | 3.727e-05 | 3.727e-05 |
| 21 | 4.861e-06 | 5.868e-10 | 9.974e-12 | 6.462e-12 | 1.481e+21 | 1.801e+21 | 1.801e+21 | 1.801e+21 | 2.176e-03 | 3.726e-05 | 3.726e-05 | 3.726e-05 |
| 26 | 4.860e-06 | 5.868e-10 | 9.959e-12 | 6.453e-12 | 1.481e+21 | 1.800e+21 | 1.800e+21 | 1.800e+21 | 2.176e-03 | 3.725e-05 | 3.725e-05 | 3.725e-05 |

Mobile species in clusters per atom, number densities N in m^-3, contents C in defects per atom.

## Comparison with the reference run

Interior values of this run and of the reference run of the same geometry, temperature and dose rate: DisloCluster ZrMicro/output/20260806_153853_15c465f_hex400 (commit 15c465f), table provenance_comparison.md, column 'hexagonal 400 nm'.
Interior values (the tenth of the nodes with the highest vacancy concentration) of the DisloCluster run of the hexagonal crystal: MoDELib3 of 2026-08-06 with the material file of that date (four loop families, atomic volume 1.2e-29 m^3), steps of 1 dpa with 20 sub-steps for the loops. This run uses the two terms of that code that differ from the report (loopNucleationPerReaction=1: one loop per clustering reaction; coalescenceClimbBurgers=0: climb speed of coalescence without the Burgers-vector magnitude), so the two runs solve the same equations and the ratios measure the numerical differences.

| quantity | dose [dpa] | MoDELib | DisloCluster | ratio |
|---|---:|---:|---:|---:|
| <c> loop density [m^-3] | 1 | 1.482e+21 | 1.481e+21 | 1.00 |
| <c> loop density [m^-3] | 6 | 1.481e+21 | 1.482e+21 | 1.00 |
| <c> loop density [m^-3] | 11 | 1.481e+21 | 1.482e+21 | 1.00 |
| <c> loop density [m^-3] | 16 | 1.481e+21 | 1.482e+21 | 1.00 |
| <c> loop density [m^-3] | 21 | 1.481e+21 | 1.482e+21 | 1.00 |
| <c> loop density [m^-3] | 26 | 1.481e+21 | 1.482e+21 | 1.00 |
| <a>1 loop density [m^-3] | 1 | 1.634e+21 | 1.336e+21 | 1.22 |
| <a>1 loop density [m^-3] | 6 | 1.757e+21 | 1.670e+21 | 1.05 |
| <a>1 loop density [m^-3] | 11 | 1.765e+21 | 1.675e+21 | 1.05 |
| <a>1 loop density [m^-3] | 16 | 1.768e+21 | 1.677e+21 | 1.05 |
| <a>1 loop density [m^-3] | 21 | 1.769e+21 | 1.677e+21 | 1.05 |
| <a>1 loop density [m^-3] | 26 | 1.768e+21 | 1.677e+21 | 1.05 |
| <c> loop diameter [nm] | 1 | 5.604e+01 | 5.588e+01 | 1.00 |
| <c> loop diameter [nm] | 6 | 5.919e+01 | 5.779e+01 | 1.02 |
| <c> loop diameter [nm] | 11 | 5.931e+01 | 5.790e+01 | 1.02 |
| <c> loop diameter [nm] | 16 | 5.935e+01 | 5.794e+01 | 1.02 |
| <c> loop diameter [nm] | 21 | 5.937e+01 | 5.796e+01 | 1.02 |
| <c> loop diameter [nm] | 26 | 5.936e+01 | 5.796e+01 | 1.02 |
| <a>1 loop diameter [nm] | 1 | 9.499e+00 | 9.280e+00 | 1.02 |
| <a>1 loop diameter [nm] | 6 | 9.100e+00 | 9.310e+00 | 0.98 |
| <a>1 loop diameter [nm] | 11 | 9.082e+00 | 9.300e+00 | 0.98 |
| <a>1 loop diameter [nm] | 16 | 9.075e+00 | 9.290e+00 | 0.98 |
| <a>1 loop diameter [nm] | 21 | 9.073e+00 | 9.290e+00 | 0.98 |
| <a>1 loop diameter [nm] | 26 | 9.074e+00 | 9.290e+00 | 0.98 |
| vacancies [m^-3] | 1 | 3.952e+23 | 4.300e+23 | 0.92 |
| vacancies [m^-3] | 6 | 3.933e+23 | 3.962e+23 | 0.99 |
| vacancies [m^-3] | 11 | 3.928e+23 | 3.959e+23 | 0.99 |
| vacancies [m^-3] | 16 | 3.923e+23 | 3.956e+23 | 0.99 |
| vacancies [m^-3] | 21 | 3.919e+23 | 3.954e+23 | 0.99 |
| vacancies [m^-3] | 26 | 3.916e+23 | 3.952e+23 | 0.99 |
| interstitials [m^-3] | 1 | 4.832e+19 | 5.268e+19 | 0.92 |
| interstitials [m^-3] | 6 | 4.804e+19 | 4.893e+19 | 0.98 |
| interstitials [m^-3] | 11 | 4.794e+19 | 4.886e+19 | 0.98 |
| interstitials [m^-3] | 16 | 4.788e+19 | 4.882e+19 | 0.98 |
| interstitials [m^-3] | 21 | 4.784e+19 | 4.879e+19 | 0.98 |
| interstitials [m^-3] | 26 | 4.781e+19 | 4.878e+19 | 0.98 |
| di-interstitials [m^-3] | 1 | 7.262e+17 | 1.214e+18 | 0.60 |
| di-interstitials [m^-3] | 6 | 7.072e+17 | 7.811e+17 | 0.91 |
| di-interstitials [m^-3] | 11 | 7.007e+17 | 7.745e+17 | 0.90 |
| di-interstitials [m^-3] | 16 | 6.963e+17 | 7.702e+17 | 0.90 |
| di-interstitials [m^-3] | 21 | 6.931e+17 | 7.671e+17 | 0.90 |
| di-interstitials [m^-3] | 26 | 6.907e+17 | 7.648e+17 | 0.90 |
| tri-interstitials [m^-3] | 1 | 4.711e+17 | 7.587e+17 | 0.62 |
| tri-interstitials [m^-3] | 6 | 4.589e+17 | 4.914e+17 | 0.93 |
| tri-interstitials [m^-3] | 11 | 4.547e+17 | 4.874e+17 | 0.93 |
| tri-interstitials [m^-3] | 16 | 4.519e+17 | 4.847e+17 | 0.93 |
| tri-interstitials [m^-3] | 21 | 4.499e+17 | 4.828e+17 | 0.93 |
| tri-interstitials [m^-3] | 26 | 4.483e+17 | 4.814e+17 | 0.93 |

Figure: `figures/comparison_with_reference.png`.

## Files

| Path | Contents |
|---|---|
| `provenance.json`, `provenance.md` | settings, input files and executables with their SHA-256, commit, machine, times |
| `config.json`, `inputFiles/` | the configuration and the input files as run |
| `evl/` | the fields at the snapshot doses (`evl_<n>.txt`) and the positions of the finite-element nodes (`cdNodes.txt`) |
| `F/` | the average plastic distortion (growth strain) at every step |
| `logs/` | the output of the executables, one file per dose interval |
| `volume_means.csv`, `summary.json` | volume means at the snapshot doses; the summary and the checks |
| `figures/` | volume means against dose, mobile profile, comparison with DisloCluster |
| `gb/` | each field against the distance from the grain boundary, one curve per dose |
| `3d/` | each field on two mid-planes of the grain, at each snapshot dose; `loops_<family>_<dose>dpa.png`: the discrete loops of a family |
| `volume_average/` | volume means against dose: point defects, loop densities (also by region) and diameters (also against the measured ranges of the D1/M1 report), growth strain and its rate, hardening estimate, retained defects and their shares, fluxes to the sinks, balance of v and i; `additional_analysis.csv`, `growth_strain.csv` |
| `boundary_flux.md` | balance of each mobile species at each dose, and the share that leaves through the grain boundary |
| `gb/denuded_zone.png` | number density and mean diameter of the loops against the distance to the nearest face, at the last dose |
| `size_spectrum/` | distribution of the loop diameters over the grain, its interior and its boundary shell, at the last dose |
| `discrete_loops/` | discrete loops drawn from the fields at up to four doses: `loops_<dose>dpa.csv`, `manifest.json`, three-dimensional views of all the loops and of each family |
| `tem_slices/` | the discrete loops of a foil projected along [0001], [01-10] and [2-1-10]; one image at the last dose and a montage over the doses |
| `movies/` | the mid-plane figures of four fields as animated images, one frame per snapshot dose |

Not produced, because the model of this run does not carry the quantities: the second moments of the size distributions (`q_*` figures of DisloCluster) and the accumulated conservation error. Each family has one loop size at each point; the size spectra and the discrete loops show how that size varies over the grain.
