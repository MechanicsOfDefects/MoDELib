# 20261002_205527_4b9e666-dirty_hex_400nm_hex400

- material: `Zr_CD4opt.txt` at 573 K, 1e-07 dpa/s
- geometry: `hex_400nm` (hexagonal), 400 x 346 x 653 nm, every face a grain boundary
- integration of the immobile species: **cvode**
- snapshots: 7, 0 .. 26 dpa
- finite-element nodes: 32317 (129268 mobile and 258536 immobile unknowns)
- status: completed; wall time 2468 s

## Verification

| Check | Expected | Result | |
|---|---|---|---|
| no negative field | minimum >= 0 | minimum 1.30e-70 | PASS |
| the three <a> families are equal without stress | relative difference < 1e-8 | 6.6e-13 | PASS |
| <c> density in the interior against production times lifetime | <= 1.493e+21 m^-3, equal while the loops do not coalesce | 1.483e+21 m^-3 at 26 dpa | PASS |
| <c> density at the first dose beyond 20 lifetimes | 1.493e+21 m^-3 within 5 % | 1.483e+21 m^-3 at 1 dpa | PASS |
| vacancies at the grain boundary at thermal equilibrium | exp(-Ef/kT) = 8.459e-15 | largest relative difference 1.1e-05 | PASS |

## Intervals

| dose from | dose to | steps | wall |
|---:|---:|---:|---:|
| 0 | 1 | 1 | 181 s |
| 1 | 6 | 5 | 477 s |
| 6 | 11 | 5 | 460 s |
| 11 | 16 | 5 | 457 s |
| 16 | 21 | 5 | 447 s |
| 21 | 26 | 5 | 442 s |

## Values in the interior of the grain

Means over the nodes farther than 80 % of the half size from the grain boundary.

| dose [dpa] | Cv | Ci | C2i | C3i | N_c | N_a1 | N_a2 | N_a3 | C_c | C_a1 | C_a2 | C_a3 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 5.252e-06 | 6.359e-10 | 1.640e-11 | 1.059e-11 | 2.959e-05 | 2.959e-05 | 2.959e-05 | 2.959e-05 | 1.000e-33 | 1.000e-33 | 1.000e-33 | 1.000e-33 |
| 1 | 4.882e-06 | 5.986e-10 | 1.079e-11 | 6.984e-12 | 1.483e+21 | 1.120e+21 | 1.120e+21 | 1.120e+21 | 2.708e-03 | 2.650e-05 | 2.650e-05 | 2.650e-05 |
| 6 | 4.858e-06 | 5.946e-10 | 1.040e-11 | 6.737e-12 | 1.483e+21 | 1.301e+21 | 1.301e+21 | 1.301e+21 | 2.946e-03 | 2.910e-05 | 2.910e-05 | 2.910e-05 |
| 11 | 4.857e-06 | 5.944e-10 | 1.035e-11 | 6.704e-12 | 1.483e+21 | 1.305e+21 | 1.305e+21 | 1.305e+21 | 2.953e-03 | 2.913e-05 | 2.913e-05 | 2.913e-05 |
| 16 | 4.856e-06 | 5.942e-10 | 1.032e-11 | 6.683e-12 | 1.483e+21 | 1.307e+21 | 1.307e+21 | 1.307e+21 | 2.956e-03 | 2.914e-05 | 2.914e-05 | 2.914e-05 |
| 21 | 4.856e-06 | 5.942e-10 | 1.030e-11 | 6.668e-12 | 1.483e+21 | 1.308e+21 | 1.308e+21 | 1.308e+21 | 2.958e-03 | 2.915e-05 | 2.915e-05 | 2.915e-05 |
| 26 | 4.855e-06 | 5.942e-10 | 1.028e-11 | 6.657e-12 | 1.483e+21 | 1.308e+21 | 1.308e+21 | 1.308e+21 | 2.958e-03 | 2.915e-05 | 2.915e-05 | 2.915e-05 |

Mobile species in clusters per atom, number densities N in m^-3, contents C in defects per atom.

## Comparison with the reference run

Interior values of this run and of the reference run of the same geometry, temperature and dose rate: DisloCluster ZrMicro/output/20260806_153853_15c465f_hex400 (commit 15c465f), table provenance_comparison.md, column 'hexagonal 400 nm'.
Interior values (the tenth of the nodes with the highest vacancy concentration) of the DisloCluster run of the hexagonal crystal: MoDELib3 of 2026-08-06 with the material file of that date (four loop families, atomic volume 1.2e-29 m^3), steps of 1 dpa with 20 sub-steps for the loops. Two terms of that run differ from this one. (1) Nucleation by clustering: that run created one loop for each reaction between mobile species whose product is not mobile, holding the 4 or 5 interstitials of the reaction. This run nucleates the loops from the cascades only, as the report states for its figures; Eq. 37 of the report, which divides the clustered content by the nucleation size, would add 0.3 %. The one-loop-per-reaction source raises the <a> loop density by a factor 1.26. (2) Climb speed of coalescence: that code omits the magnitude of the Burgers vector (1.63 b for the <c> loops), which makes the coalescence of the <c> loops with the network faster and their diameter smaller by a factor 0.83. With these two terms written as in that code and the mobile concentrations of that run, the rate equations of MoDELib give its loop densities and diameters to 0.1 %: the expected ratios of this run to that one are 0.79 to 0.81 for the <a> density and 1.19 to 1.20 for the <c> diameter.

| quantity | dose [dpa] | MoDELib | DisloCluster | ratio |
|---|---:|---:|---:|---:|
| <c> loop density [m^-3] | 1 | 1.483e+21 | 1.481e+21 | 1.00 |
| <c> loop density [m^-3] | 6 | 1.483e+21 | 1.482e+21 | 1.00 |
| <c> loop density [m^-3] | 11 | 1.483e+21 | 1.482e+21 | 1.00 |
| <c> loop density [m^-3] | 16 | 1.483e+21 | 1.482e+21 | 1.00 |
| <c> loop density [m^-3] | 21 | 1.483e+21 | 1.482e+21 | 1.00 |
| <c> loop density [m^-3] | 26 | 1.483e+21 | 1.482e+21 | 1.00 |
| <a>1 loop density [m^-3] | 1 | 1.159e+21 | 1.336e+21 | 0.87 |
| <a>1 loop density [m^-3] | 6 | 1.344e+21 | 1.670e+21 | 0.81 |
| <a>1 loop density [m^-3] | 11 | 1.354e+21 | 1.675e+21 | 0.81 |
| <a>1 loop density [m^-3] | 16 | 1.359e+21 | 1.677e+21 | 0.81 |
| <a>1 loop density [m^-3] | 21 | 1.361e+21 | 1.677e+21 | 0.81 |
| <a>1 loop density [m^-3] | 26 | 1.362e+21 | 1.677e+21 | 0.81 |
| <c> loop diameter [nm] | 1 | 6.641e+01 | 5.588e+01 | 1.19 |
| <c> loop diameter [nm] | 6 | 6.929e+01 | 5.779e+01 | 1.20 |
| <c> loop diameter [nm] | 11 | 6.945e+01 | 5.790e+01 | 1.20 |
| <c> loop diameter [nm] | 16 | 6.952e+01 | 5.794e+01 | 1.20 |
| <c> loop diameter [nm] | 21 | 6.955e+01 | 5.796e+01 | 1.20 |
| <c> loop diameter [nm] | 26 | 6.955e+01 | 5.796e+01 | 1.20 |
| <a>1 loop diameter [nm] | 1 | 9.639e+00 | 9.280e+00 | 1.04 |
| <a>1 loop diameter [nm] | 6 | 9.375e+00 | 9.310e+00 | 1.01 |
| <a>1 loop diameter [nm] | 11 | 9.355e+00 | 9.300e+00 | 1.01 |
| <a>1 loop diameter [nm] | 16 | 9.347e+00 | 9.290e+00 | 1.01 |
| <a>1 loop diameter [nm] | 21 | 9.344e+00 | 9.290e+00 | 1.01 |
| <a>1 loop diameter [nm] | 26 | 9.343e+00 | 9.290e+00 | 1.01 |
| vacancies [m^-3] | 1 | 3.952e+23 | 4.300e+23 | 0.92 |
| vacancies [m^-3] | 6 | 3.926e+23 | 3.962e+23 | 0.99 |
| vacancies [m^-3] | 11 | 3.921e+23 | 3.959e+23 | 0.99 |
| vacancies [m^-3] | 16 | 3.916e+23 | 3.956e+23 | 0.99 |
| vacancies [m^-3] | 21 | 3.912e+23 | 3.954e+23 | 0.99 |
| vacancies [m^-3] | 26 | 3.909e+23 | 3.952e+23 | 0.99 |
| interstitials [m^-3] | 1 | 4.893e+19 | 5.268e+19 | 0.93 |
| interstitials [m^-3] | 6 | 4.848e+19 | 4.893e+19 | 0.99 |
| interstitials [m^-3] | 11 | 4.838e+19 | 4.886e+19 | 0.99 |
| interstitials [m^-3] | 16 | 4.831e+19 | 4.882e+19 | 0.99 |
| interstitials [m^-3] | 21 | 4.826e+19 | 4.879e+19 | 0.99 |
| interstitials [m^-3] | 26 | 4.823e+19 | 4.878e+19 | 0.99 |
| di-interstitials [m^-3] | 1 | 7.534e+17 | 1.214e+18 | 0.62 |
| di-interstitials [m^-3] | 6 | 7.237e+17 | 7.811e+17 | 0.93 |
| di-interstitials [m^-3] | 11 | 7.168e+17 | 7.745e+17 | 0.93 |
| di-interstitials [m^-3] | 16 | 7.121e+17 | 7.702e+17 | 0.92 |
| di-interstitials [m^-3] | 21 | 7.087e+17 | 7.671e+17 | 0.92 |
| di-interstitials [m^-3] | 26 | 7.061e+17 | 7.648e+17 | 0.92 |
| tri-interstitials [m^-3] | 1 | 4.884e+17 | 7.587e+17 | 0.64 |
| tri-interstitials [m^-3] | 6 | 4.694e+17 | 4.914e+17 | 0.96 |
| tri-interstitials [m^-3] | 11 | 4.650e+17 | 4.874e+17 | 0.95 |
| tri-interstitials [m^-3] | 16 | 4.619e+17 | 4.847e+17 | 0.95 |
| tri-interstitials [m^-3] | 21 | 4.598e+17 | 4.828e+17 | 0.95 |
| tri-interstitials [m^-3] | 26 | 4.581e+17 | 4.814e+17 | 0.95 |

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
