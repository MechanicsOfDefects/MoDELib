# Simulations: zirconium spatially resolved cluster dynamics

[`run_simulation.ipynb`](run_simulation.ipynb) sets up and runs a simulation of irradiated zirconium in one grain whose surface is a grain boundary, and gathers its output. The mobile species are solved by MoDELib at steady state on the finite-element mesh; the loop populations are integrated at every node by CVODE. The equations are those of [`docs/reports/GW_Phase4_D1M1.pdf`](../docs/reports/GW_Phase4_D1M1.pdf), with the parameters of `Library/Materials/Zr_CD4opt.txt`.

The sections of the notebook are:

1. **Environment.** Finds the MoDELib checkout from the notebook's own folder.
2. **Build.** Configures and builds `DDomp` and `microstructureGenerator` (CMake, SUNDIALS).
3. **The application.** The geometries, presets and materials available.
4. **Control parameters.** Edit these.
5. **Material and its overrides.** Edit these.
6. **Validate.** Unknown keys and inconsistent settings raise here, before anything runs.
7. **Run.**
8. **Gather the output.**

The logic behind the notebook is in [`simulation_driver.py`](simulation_driver.py) (build, stage, run) and [`postprocessing.py`](postprocessing.py) (volume means, figures, checks, report).

**Requirements.** Python with NumPy, SciPy, Matplotlib and Jupyter. A C++20 compiler, CMake, Eigen, FFTW and SUNDIALS 7 where the executables are built: on Windows that is WSL, which the notebook calls; on Linux and macOS it is the machine itself. SUNDIALS is looked for in `$SUNDIALS_ROOT`, `~/sundials-install`, a `Libraries/sundials-7.1.1` folder beside the repository, and the system locations.

**Presets.** `PRESET` in §4 names a configuration:

| Preset | Doses [dpa] | Steps | Purpose |
|---|---|---|---|
| `reference` | 1e-4, 1e-3, 1e-2, 0.1, 1, 2, 5, 10 | 3 per interval | the 200 nm cube to 10 dpa, at the doses of the ENNDS reference run |
| `hex400` | 1, 6, 11, 16, 21, 26 | 1 dpa each | the hexagonal crystal of the D1/M1 report (its Figure 27), with the dose steps of the DisloCluster run that produced it |
| `hex400dc` | as `hex400` | as `hex400` | the same crystal with the two terms of the DisloCluster code that differ from the report (`loopNucleationPerReaction=1`, `coalescenceClimbBurgers=0`, clustering nucleation on): a check of the numerical accuracy against the DisloCluster run |
| `smoke` | 1e-4 | 3 | a check of the build and of the workflow |

## The verification cases

There are two: a cube of 200 nm compared with the ENNDS reference run, and the hexagonal crystal of the D1/M1 report compared with the DisloCluster run behind its figures.

### Cube of 200 nm (`reference`)

| | |
|---|---|
| Geometry | [`geometry/cube_200nm`](geometry/cube_200nm): a cube of 200 nm, 6762 mesh nodes with a refined layer of 40 nm at the surface; 47 828 finite-element nodes (quadratic tetrahedra) |
| Boundary | every face is a grain boundary: the mobile species are at thermal equilibrium there |
| Conditions | 573 K, 1e-7 dpa/s, no applied stress |
| Material | `Zr_CD4opt.txt`: four mobile species (v, i, 2i, 3i), four loop families (c, a1, a2, a3) |
| Initial state | no loops |

The report of a run (`report.md`) holds three kinds of results.

- **Checks** against values known independently of the run: no negative field; the three prismatic families equal without stress; the vacancy concentration at the grain boundary equal to its thermal equilibrium; the density of the vacancy loops against production times lifetime.
- **Values in the interior of the grain** at each snapshot dose.
- **A comparison with the ENNDS reference run** of the same geometry ([`reference/ENNDS_cube_200nm_10dpa`](reference/ENNDS_cube_200nm_10dpa)). The two codes do not solve the same model: MoDELib integrates the equations of the D1/M1 report with four loop families, ENNDS the later model of DisloCluster with nine family slots, second moments and refitted parameters. The comparison shows how far apart the two descriptions are at each dose. It is not a test that one code reproduces the other.

### Hexagonal crystal of 400 nm (`hex400`)

| | |
|---|---|
| Geometry | [`geometry/hex_400nm`](geometry/hex_400nm): a hexagonal prism along the c axis, 400 nm corner to corner, 346.4 nm flat to flat, 653.2 nm tall. The mesh is that of the report (its Figure 24, right): 4672 mesh nodes, 20 726 tetrahedra, 32 317 finite-element nodes |
| Boundary, conditions, material, initial state | as for the cube |
| Steps | 26 steps of 1 dpa, as in the DisloCluster run; snapshots at 1, 6, 11, 16, 21 and 26 dpa |

This case reproduces Figure 27 of the report: `figures/loop_populations_6dpa.png` shows the three prismatic families and the basal family at 6 dpa, each with its content on the basal mid-plane and some of its loops drawn in their habit plane with magnified radii (7 times for the prismatic loops, 1.5 times for the basal loops).

The report of the run compares the interior of the grain, at each snapshot dose, with the values of the DisloCluster run ([`reference/DisloCluster_hex_400nm_26dpa`](reference/DisloCluster_hex_400nm_26dpa)): loop densities and diameters of the four families and the concentrations of the four mobile species. That run used the same atomic volume and fitted parameters, so the comparison is quantitative. The mobile concentrations, the basal density and the prismatic diameter agree within a few percent after the first dose. Two quantities differ, and both differences are traced to a term of the DisloCluster code. Its prismatic density is 1.24 times higher, because it created one loop for each reaction between mobile species whose product is not mobile; this run nucleates the loops from the cascades only, as the report states for its figures (and Eq. 37 of the report, which divides the clustered content by the nucleation size, would add 0.3 %). Its basal diameter is 0.83 times smaller, because its climb speed of coalescence omitted the magnitude of the Burgers vector. With these two terms written as in that code, the rate equations of MoDELib reproduce its interior loop densities and diameters to 0.1 %.

**Numerical accuracy against DisloCluster (`hex400dc`).** With the two terms of the DisloCluster code switched on, the two codes solve the same equations. From 6 dpa onward the interior values of the two runs agree:

| Quantity | MoDELib / DisloCluster, 6 to 26 dpa |
|---|---|
| `<c>` loop density | 1.00 |
| `<a>` loop density | 1.05 |
| `<c>` loop diameter | 1.02 |
| `<a>` loop diameter | 0.98 |
| vacancies, interstitials | 0.98 to 0.99 |
| di- and tri-interstitials | 0.90 to 0.93 |

At 1 dpa, the end of the first step, the differences are larger (`<a>` density 1.22, di- and tri-interstitials 0.6): the two codes do not integrate the first step of 1 dpa in the same way (MoDELib with CVODE, DisloCluster with 20 explicit sub-steps). The deficit of di- and tri-interstitials is the same with the report's terms (`hex400`), so it comes from the solution of the mobile species, not from the loops; its origin has not been traced. The independent implementation of the rate equations (`testsPy/clusterDynamics/referenceRates.py`), with the same two terms and the mobile concentrations of the DisloCluster run, gives its interior loop densities and diameters to 0.1 %.

The rate equations themselves are tested against an independent implementation in [`testsPy/clusterDynamics`](../testsPy/clusterDynamics).

## Output

Each run writes `output/<UTC YYYYMMDD_HHMMSS>_<git commit>[-dirty]_<tag>/`:

| Path | Contents |
|---|---|
| `provenance.json`, `provenance.md` | every setting and override; input files and executables with their SHA-256; MoDELib commit and whether the working tree was modified; machine; start, finish, wall time, status |
| `config.json`, `inputFiles/` | the configuration, and the input files as run |
| `run_record.json`, `logs/` | the dose intervals with their wall times; the output of the executables, one file per interval |
| `evl/` | the fields at the snapshot doses (`evl_<n>.txt`, full precision) and the positions of the finite-element nodes (`cdNodes.txt`) |
| `F/` | the average plastic distortion (growth strain) at every step |
| `volume_means.csv`, `summary.json` | volume means at the snapshot doses; the summary and the checks |
| `report.md` | checks, intervals, interior values, comparison with the reference |
| `figures/` | volume means against dose; mobile species against the distance from the grain boundary; comparison with the reference run; `loop_populations_<dose>dpa.png`: the loop families in the grain, as in Figure 27 of the D1/M1 report |
| `gb/` | each field against the distance from the grain boundary, one curve per dose; `denuded_zone.png`: density and diameter of the loops against the distance to the nearest face |
| `3d/` | each field on two mid-planes of the grain, at each snapshot dose, on a logarithmic color scale; the discrete loops of the `c` and `a1` families |
| `volume_average/` | volume means against dose: point defects, loop densities (also by region), loop diameters (also against the measured ranges quoted in the D1/M1 report), growth strain and its rate, a hardening estimate, the retained defects and their shares, the fluxes to the sinks, the balance of the vacancies and of the interstitials; `additional_analysis.csv`, `growth_strain.csv` |
| `boundary_flux.md` | the balance of each mobile species at each dose, and the share that leaves through the grain boundary |
| `size_spectrum/` | distribution of the loop diameters over the grain, its interior and its boundary shell |
| `discrete_loops/` | discrete loops drawn from the fields at up to four doses: tables, `manifest.json`, three-dimensional views |
| `tem_slices/` | the discrete loops of a foil of 100 nm projected along [0001], [01-10] and [2-1-10] |
| `movies/` | the mid-plane figures of four fields as animated images |

The folders follow those of the DisloCluster and ENNDS runs. Two of their figure families are not produced, because the model of MoDELib does not carry the quantities: the second moments of the size distributions (`q_*`) and the accumulated conservation error. Each loop family has one size at each point, so the size spectra and the discrete loops show how that size varies over the grain, not a distribution at a point. The figures are written by [`postprocessing.py`](postprocessing.py) and [`postprocessing_extended.py`](postprocessing_extended.py).

Run folders are not tracked by git. Only `output/README.md` is. The three verification runs are kept in [`reference_output/`](reference_output), without their fields and logs.

A run can be continued, or used as the initial state of another study: `evl/evl_<n>.txt` holds the fields of the mobile and immobile species, and a run that finds it as `evl/evl_0.txt` in its own folder starts from them.

## Running headless

```bash
MODELIB_SIMULATION_PRESET=smoke jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=-1 --output-dir <scratch> simulations/run_simulation.ipynb
```

Set `PYTHONIOENCODING=utf-8` first on Windows.

## Source of the geometry and of the reference

The mesh and `case.json` of `geometry/cube_200nm`, and the volume means of `reference/ENNDS_cube_200nm_10dpa`, are from the ENNDS repository (<https://github.com/Ghoniem/ENNDS>, version 3.1.0, `application_fields/radiation_damage/zirconium_spatial_cluster_dynamics`), included here with the permission of its owner.

The mesh of `geometry/hex_400nm` and the interior values of `reference/DisloCluster_hex_400nm_26dpa` are from the DisloCluster repository of the same owner: the mesh `Gmsh/meshes/hexagonal_400nm_15000el_o1_bl_65b0ac58.msh`, and the table `provenance_comparison.md` of the run `ZrMicro/output/20260806_153853_15c465f_hex400`.
