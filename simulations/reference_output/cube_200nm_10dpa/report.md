# 20261002_215830_4b9e666-dirty_cube_200nm_reference

- material: `Zr_CD4opt.txt` at 573 K, 1e-07 dpa/s
- geometry: `cube_200nm` (cubic), 200 x 200 x 200 nm, every face a grain boundary
- integration of the immobile species: **cvode**
- snapshots: 9, 0 .. 10 dpa
- finite-element nodes: 47828 (191312 mobile and 382624 immobile unknowns)
- status: completed; wall time 2977 s

## Verification

| Check | Expected | Result | |
|---|---|---|---|
| no negative field | minimum >= 0 | minimum 1.30e-70 | PASS |
| the three <a> families are equal without stress | relative difference < 1e-8 | 4.9e-13 | PASS |
| <c> density in the interior against production times lifetime | <= 1.493e+21 m^-3, equal while the loops do not coalesce | 1.484e+21 m^-3 at 10 dpa | PASS |
| <c> density at the first dose beyond 20 lifetimes | 1.493e+21 m^-3 within 5 % | 1.492e+21 m^-3 at 0.1 dpa | PASS |
| vacancies at the grain boundary at thermal equilibrium | exp(-Ef/kT) = 8.459e-15 | largest relative difference 1.1e-05 | PASS |

## Intervals

| dose from | dose to | steps | wall |
|---:|---:|---:|---:|
| 0 | 0.0001 | 3 | 413 s |
| 0.0001 | 0.001 | 3 | 379 s |
| 0.001 | 0.01 | 3 | 374 s |
| 0.01 | 0.1 | 3 | 350 s |
| 0.1 | 1 | 3 | 346 s |
| 1 | 2 | 3 | 365 s |
| 2 | 5 | 3 | 352 s |
| 5 | 10 | 3 | 365 s |

## Values in the interior of the grain

Means over the nodes farther than 80 % of the half size from the grain boundary.

| dose [dpa] | Cv | Ci | C2i | C3i | N_c | N_a1 | N_a2 | N_a3 | C_c | C_a1 | C_a2 | C_a3 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 4.860e-06 | 6.158e-10 | 7.771e-12 | 5.046e-12 | 2.959e-05 | 2.959e-05 | 2.959e-05 | 2.959e-05 | 1.000e-33 | 1.000e-33 | 1.000e-33 | 1.000e-33 |
| 0.0001 | 4.861e-06 | 6.154e-10 | 7.762e-12 | 5.041e-12 | 1.006e+20 | 1.960e+19 | 1.960e+19 | 1.960e+19 | 4.846e-07 | 7.174e-08 | 7.174e-08 | 7.174e-08 |
| 0.001 | 4.865e-06 | 6.124e-10 | 7.691e-12 | 4.996e-12 | 7.499e+20 | 1.936e+20 | 1.936e+20 | 1.936e+20 | 3.737e-06 | 8.179e-07 | 8.179e-07 | 8.179e-07 |
| 0.01 | 4.779e-06 | 5.909e-10 | 7.009e-12 | 4.559e-12 | 1.492e+21 | 1.493e+21 | 1.493e+21 | 1.493e+21 | 1.280e-05 | 1.687e-05 | 1.687e-05 | 1.687e-05 |
| 0.1 | 4.693e-06 | 5.793e-10 | 6.581e-12 | 4.285e-12 | 1.492e+21 | 1.439e+21 | 1.439e+21 | 1.439e+21 | 2.109e-04 | 3.052e-05 | 3.052e-05 | 3.052e-05 |
| 1 | 4.529e-06 | 5.694e-10 | 6.156e-12 | 4.010e-12 | 1.484e+21 | 1.390e+21 | 1.390e+21 | 1.390e+21 | 2.897e-03 | 3.115e-05 | 3.115e-05 | 3.115e-05 |
| 2 | 4.526e-06 | 5.692e-10 | 6.143e-12 | 4.002e-12 | 1.484e+21 | 1.386e+21 | 1.386e+21 | 1.386e+21 | 2.927e-03 | 3.119e-05 | 3.119e-05 | 3.119e-05 |
| 5 | 4.525e-06 | 5.689e-10 | 6.131e-12 | 3.994e-12 | 1.484e+21 | 1.388e+21 | 1.388e+21 | 1.388e+21 | 2.929e-03 | 3.121e-05 | 3.121e-05 | 3.121e-05 |
| 10 | 4.523e-06 | 5.685e-10 | 6.114e-12 | 3.983e-12 | 1.484e+21 | 1.391e+21 | 1.391e+21 | 1.391e+21 | 2.934e-03 | 3.123e-05 | 3.123e-05 | 3.123e-05 |

Mobile species in clusters per atom, number densities N in m^-3, contents C in defects per atom.

## Comparison with the reference run

Volume means of this run and of the reference run of the same geometry, temperature and dose rate: 20261001_200157_a70f2f3_cube_200nm_reference.
The two runs do not solve the same model. MoDELib integrates the equations of the D1/M1 report with four loop families and the parameters of Zr_CD4opt.txt. ENNDS integrates the later model of DisloCluster, with nine family slots, second moments of the size distributions and refitted parameters, and its vacancy loops are split into faulted, perfect and pyramid populations, summed here. The ratios show how far the two parameter sets are from each other; they are not a test of the code.

| quantity | dose [dpa] | MoDELib | ENNDS | ratio |
|---|---:|---:|---:|---:|
| vacancies [m^-3] | 0.0001 | 2.088e+23 | 1.282e+23 | 1.63 |
| vacancies [m^-3] | 0.001 | 2.089e+23 | 1.279e+23 | 1.63 |
| vacancies [m^-3] | 0.01 | 2.055e+23 | 1.275e+23 | 1.61 |
| vacancies [m^-3] | 0.1 | 2.003e+23 | 1.275e+23 | 1.57 |
| vacancies [m^-3] | 1 | 1.951e+23 | 1.274e+23 | 1.53 |
| vacancies [m^-3] | 2 | 1.945e+23 | 1.266e+23 | 1.54 |
| vacancies [m^-3] | 5 | 1.937e+23 | 1.257e+23 | 1.54 |
| vacancies [m^-3] | 10 | 1.927e+23 | 1.236e+23 | 1.56 |
| interstitials [m^-3] | 0.0001 | 2.673e+19 | 2.127e+20 | 0.13 |
| interstitials [m^-3] | 0.001 | 2.664e+19 | 2.127e+20 | 0.13 |
| interstitials [m^-3] | 0.01 | 2.588e+19 | 2.124e+20 | 0.12 |
| interstitials [m^-3] | 0.1 | 2.512e+19 | 2.124e+20 | 0.12 |
| interstitials [m^-3] | 1 | 2.476e+19 | 2.122e+20 | 0.12 |
| interstitials [m^-3] | 2 | 2.468e+19 | 2.112e+20 | 0.12 |
| interstitials [m^-3] | 5 | 2.457e+19 | 2.101e+20 | 0.12 |
| interstitials [m^-3] | 10 | 2.441e+19 | 2.077e+20 | 0.12 |
| <c> vacancy loops [m^-3] | 0.0001 | 1.006e+20 | 7.816e+19 | 1.29 |
| <c> vacancy loops [m^-3] | 0.001 | 7.499e+20 | 3.137e+20 | 2.39 |
| <c> vacancy loops [m^-3] | 0.01 | 1.492e+21 | 6.858e+20 | 2.17 |
| <c> vacancy loops [m^-3] | 0.1 | 1.493e+21 | 3.095e+21 | 0.48 |
| <c> vacancy loops [m^-3] | 1 | 1.489e+21 | 1.780e+22 | 0.08 |
| <c> vacancy loops [m^-3] | 2 | 1.488e+21 | 3.337e+22 | 0.04 |
| <c> vacancy loops [m^-3] | 5 | 1.488e+21 | 8.010e+22 | 0.02 |
| <c> vacancy loops [m^-3] | 10 | 1.488e+21 | 1.580e+23 | 0.01 |
| <a>1 interstitial loops [m^-3] | 0.0001 | 1.957e+19 | 1.726e+19 | 1.13 |
| <a>1 interstitial loops [m^-3] | 0.001 | 1.945e+20 | 6.142e+19 | 3.17 |
| <a>1 interstitial loops [m^-3] | 0.01 | 1.750e+21 | 1.101e+20 | 15.89 |
| <a>1 interstitial loops [m^-3] | 0.1 | 5.054e+21 | 5.052e+20 | 10.00 |
| <a>1 interstitial loops [m^-3] | 1 | 1.676e+22 | 4.453e+21 | 3.76 |
| <a>1 interstitial loops [m^-3] | 2 | 3.026e+22 | 8.841e+21 | 3.42 |
| <a>1 interstitial loops [m^-3] | 5 | 7.076e+22 | 2.200e+22 | 3.22 |
| <a>1 interstitial loops [m^-3] | 10 | 1.383e+23 | 4.394e+22 | 3.15 |
| content of the <c> loops [per atom] | 0.0001 | 4.838e-07 | 2.595e-07 | 1.86 |
| content of the <c> loops [per atom] | 0.001 | 3.670e-06 | 6.574e-07 | 5.58 |
| content of the <c> loops [per atom] | 0.01 | 9.794e-06 | 2.374e-06 | 4.12 |
| content of the <c> loops [per atom] | 0.1 | 8.568e-05 | 1.775e-05 | 4.83 |
| content of the <c> loops [per atom] | 1 | 1.835e-03 | 1.495e-04 | 12.27 |
| content of the <c> loops [per atom] | 2 | 2.343e-03 | 2.936e-04 | 7.98 |
| content of the <c> loops [per atom] | 5 | 2.656e-03 | 7.229e-04 | 3.67 |
| content of the <c> loops [per atom] | 10 | 2.737e-03 | 1.428e-03 | 1.92 |
| content of the <a>1 loops [per atom] | 0.0001 | 7.105e-08 | 1.493e-07 | 0.48 |
| content of the <a>1 loops [per atom] | 0.001 | 7.617e-07 | 2.865e-06 | 0.27 |
| content of the <a>1 loops [per atom] | 0.01 | 1.257e-05 | 4.606e-06 | 2.73 |
| content of the <a>1 loops [per atom] | 0.1 | 7.338e-05 | 7.554e-06 | 9.71 |
| content of the <a>1 loops [per atom] | 1 | 1.235e-04 | 3.512e-05 | 3.52 |
| content of the <a>1 loops [per atom] | 2 | 1.725e-04 | 6.579e-05 | 2.62 |
| content of the <a>1 loops [per atom] | 5 | 3.188e-04 | 1.577e-04 | 2.02 |
| content of the <a>1 loops [per atom] | 10 | 5.627e-04 | 3.109e-04 | 1.81 |

Figure: `figures/comparison_with_ENNDS.png`.

## Files

| Path | Contents |
|---|---|
| `provenance.json`, `provenance.md` | settings, input files and executables with their SHA-256, commit, machine, times |
| `config.json`, `inputFiles/` | the configuration and the input files as run |
| `evl/` | the fields at the snapshot doses (`evl_<n>.txt`) and the positions of the finite-element nodes (`cdNodes.txt`) |
| `F/` | the average plastic distortion (growth strain) at every step |
| `logs/` | the output of the executables, one file per dose interval |
| `volume_means.csv`, `summary.json` | volume means at the snapshot doses; the summary and the checks |
| `figures/` | volume means against dose, mobile profile, comparison with ENNDS |
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
