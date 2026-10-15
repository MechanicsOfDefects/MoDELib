# Mechanics of Defects Evolution Library (MoDELib)

**A C++ library for the discrete mechanics of crystal defects.**

MoDELib evolves defects that are represented individually inside a crystal, and computes their interactions from continuum elasticity. One object, the `DefectiveCrystal`, holds every *physics* that a simulation enables:

- **dislocation dynamics:** a network of discrete dislocation loops that glide, climb, react and nucleate;
- **cluster dynamics:** mobile point-defect species that diffuse and react on the finite-element mesh, and immobile loop populations that grow by absorbing them;
- **elastic deformation:** a uniform deformation applied by a load controller, or a finite-element solution on the mesh;
- **inclusions:** Eshelby inclusions, spherical or polyhedral.

The library separates three things that a simulation needs:

- **the crystal:** its structure, slip systems, mobility laws and domain mesh;
- **the microstructure:** which defects are present, and where;
- **the algorithm:** how forces become velocities, and how the configuration is updated.

The governing rule follows from that separation:

> Each physics is a *microstructure* with the same duties: solve, propose a time step, update its configuration, and report its fields. The crystal adds the fields and takes the smallest time step.

**Languages:** C++20 for the library and its tools. Python for input generation, post-processing and the optional `pyMoDELib` bindings (pybind11).\
**Version:** 2.0.0, with the cluster-dynamics changes of [section 1.5](#15-cluster-dynamics-since-200). Input files written for 1.0.0 must be converted: see [section 9](#9-migrating-input-files-from-100).\
**License:** GNU GPL v2, as stated in the header of every source file. The repository has no separate `LICENSE` file.\
**Primary reference:** G. Po, M. S. Mohamed, T. Crosby, C. Erel, A. El-Azab, N. Ghoniem, *Recent progress in discrete dislocation dynamics and its applications to micro plasticity*, JOM 66 (2014) 2108–2120.

## Contents

1. [What is new in 2.0.0](#1-what-is-new-in-200)
2. [Capabilities](#2-capabilities)
3. [Methodology](#3-methodology)
4. [Repository structure](#4-repository-structure)
5. [Installation](#5-installation)
6. [Quick start](#6-quick-start)
7. [Tutorials](#7-tutorials)
8. [Outputs](#8-outputs)
9. [Migrating input files from 1.0.0](#9-migrating-input-files-from-100)
10. [Tests and validation](#10-tests-and-validation)
11. [Status and known issues](#11-status-and-known-issues)
12. [Documentation](#12-documentation)
13. [How to cite](#13-how-to-cite)
14. [Repository conventions](#14-repository-conventions)
15. [License and contact](#15-license-and-contact)

---

## 1. What is new in 2.0.0

Version 2.0.0 adds a spatially resolved model of irradiation-induced loop populations, replaces the mobility laws by composable ones, and reorganizes the input files. The sections named in each item give the details.

### 1.1 Cluster dynamics with immobile species

In 1.0.0 cluster dynamics solved for the mobile point-defect species only. In 2.0.0 the mobile species feed **immobile loop populations** that are fields on the same finite-element mesh ([section 3.6](#36-cluster-dynamics)).

- **Eight immobile fields.** For each of four loop families, `c`, `a1`, `a2` and `a3`, the code carries a number density $N_j$ and a defect content $c_j$ (the point defects stored in that family, per lattice site). The mean cluster size is $n_j = c_j/(N_j\,\Omega)$.
- **Sink strengths that follow the cluster shape.** A cluster is treated as a bi-pyramid below `shapeTransitionMinSize` defects and as a dislocation loop above `shapeTransitionMaxSize`, with a sigmoidal transition between the two. The sink strength, the capture bias and the eigenstrain all follow that transition.
- **Anisotropic capture.** The bias of a loop family combines a drift factor (`immobileBias`, one value per mobile species) and the diffusion-anisotropy factor of the mobile species (the ratio of its diffusion coefficients along and normal to the `c` axis).
- **Network dislocations as sinks.** `dislocationSinks_SI` gives a dislocation density per family, in addition to the loops.
- **Growth strain.** Cluster dynamics now contributes to the average plastic distortion of the crystal: the eigenstrain of the loop populations, and the relaxation volume of the dissolved defects. In 1.0.0 that contribution was zero.
- **Dissociation of mobile clusters**, with a binding energy per species (`mobileSpeciesBindingEnergy_eV`).
- **New files:** `FirstOrderReaction.h`, `ImmobileSinkRate.h`, `SpatialODESolver.h`, the material file `Zr4_Fitted.txt`, and the tutorial `spatialCDtest`.
- **Species counts are build options.** The numbers of mobile and immobile species are compile-time constants. They are now set by the CMake options `MODELIB_CD_MSIZE` and `MODELIB_CD_ISIZE` (defaults 4 and 8) instead of being edited in `ClusterDynamicsParameters.h` ([section 5](#5-installation)).

### 1.2 Physics selection

The four switches `useDislocations`, `useClusterDynamics`, `useInclusions` and `useElasticDeformation` of `DD.txt` are replaced by one line in a new input file, `DefectiveCrystal.txt`:

```
physics=DislocationDynamics ClusterDynamics ElasticDeformation;
```

The run controls that are not specific to dislocations (`Nsteps`, `dtMax`, `outputFrequency`, `outputBinary`, `useFEM`, the restart step and the periodic image settings) moved to the same file. A physics can also be re-solved several times within a step (`maxResolveSteps`), which the coupled climb and diffusion problem uses. See [section 3.1](#31-the-defective-crystal).

### 1.3 Composable mobility laws

The per-crystal mobility classes (FCC, BCC, HEX basal, HEX prismatic, HEX pyramidal) are removed. A mobility law is now declared **per slip system** in the material file and assembled from building blocks ([section 3.3](#33-mobility-laws)):

- `viscousDrag` and `kinkPair`, the two elementary laws;
- `edgeScrew`, which blends an edge law and a screw law;
- `interpolated`, which interpolates any number of laws over the character angle (linear, logarithmic or Lagrange);
- `python`, a law written in a Python module.

Slip systems are named by their crystallography, for example `a/2<111>{110}` instead of `full<111>{110}`.

### 1.4 Other additions

- **Cubic fluorite crystal** with the slip systems `a/2<110>{111}`, `a/2<110>{110}` and `a/2<110>{100}`, and a `UO2.txt` material file.
- **New material files:** `FeCrAl_Fe.txt`, `UO2.txt`, `Zr4_Fitted.txt`, `Zr_BMD19.txt`, `Zr_CD3_BMD19.txt`.
- **Stacking-fault noise from MD** (`MDStackingFaultNoise`): a correlated stacking-fault energy field on the glide plane, which acts on partial dislocations. Correlation data for Al–Mg (5, 10 and 15 at.%) and Fe–Cr are in `Library/GlidePlaneNoise/`.
- **Gamma surfaces of HEX prismatic planes**, given as wave vectors and energy points in the material file.
- **Dislocation nucleation** in the bulk (`bulkNucleationModel`), as a first model ([section 11](#11-status-and-known-issues)).
- **Nodal velocity constraints** (`nodalVelocityConstraints`): directions removed from the velocity of every node.
- **Dislocation density per slip system** in the output (`outputDislocationDensityPerSlipSystem`).
- **Microstructure generation** moved from `DislocationMicrostructure/` to `MicrostructureGeneration/`. Loop densities and radii are given in SI units (`targetDensity_SI`, `radiusDistributionMean_SI`); Frank loops can be restricted to grains and planes.
- **One-dimensional interpolants** (`LinearInterpolant`, `LagrangeInterpolant`) with periodic extrapolation.
- **New meshes:** `unitCube_15K`, `unitCube_70K`, `disk2d`, and periodic polycrystal cubes of 10 and 50 grains.
- **`pyMoDELib`:** bindings for the mobility laws, the noise generators, gamma surfaces, slip systems, lattices, the stress of a straight segment and the interpolants.
- **`DDqt`:** a slip-system tab with a color per slip system, segments and slipped areas coloured by slip system, stacking faults coloured by misfit energy, animated playback, and face labels and periodic faces on the mesh.
- **Tutorials:** `annealing`, `dipoles`, `spatialCDtest` and `uniformLoadController` are new; `dipoleNoise` is rewritten; `crossSlip` and `polycrystal` are removed.
- **Scripts:** `lib/` (now in `python/`) holds readers of the output files; `testsPy/` holds verification scripts.
- **Build:** portable `CMakeLists.txt` with options for the optional parts ([section 5](#5-installation)).

### 1.5 Cluster dynamics since 2.0.0

The immobile loop populations now obey the full rate equations of the D1/M1 report [`docs/reports/GW_Phase4_D1M1.pdf`](docs/reports/GW_Phase4_D1M1.pdf), integrated at every finite-element node by CVODE ([section 3.6](#36-cluster-dynamics)). The explicit integrator of 2.0.0 (`FirstOrderReaction.h`, `ImmobileSinkRate.h`, `SpatialODESolver.h`), the tutorial `spatialCDtest` and the input `immobileIntegrator` are removed; `Zr4_Fitted.txt` is kept for reference only and does not run with the current code.

- **Rate equations of the report:** nucleation in cascades and by clustering of mobile species, absorption with a size gate, dissolution of the vacancy embryos and thermal emission, like-loop and loop-network coalescence, and stress-oriented nucleation. New classes `ImmobileRateEquations`, `ImmobileSpeciesSolver`, `LoopSinkReaction` and `ImmobileClusterSource`; the mobile and immobile solves exchange their fields in memory.
- **SUNDIALS 7 is required** for cluster dynamics with immobile species.
- **Cluster geometry used consistently:** the capture length blends the sink strength of a sphere and of a loop with the shape sigmoid; coalescence uses the effective radius and a climb speed from the net growth of the family; the eigenstrain is one atomic volume per defect in both shapes, isotropic for bi-pyramids and along the Burgers vector for loops.
- **Inputs:** `immobileSpeciesBias` (capture bias of each family for each mobile species; `discreteDislocationBias` is used by the discrete segments only), `clusterDissolutionRatePrefactor` and `clusterDissolutionEnergy`, `clusterBindingEnergy_eV`, `cLL`, `cLN`, `kappaLL`, `kappaLN`. The size gate needs no input. `atomicVolume_SI`, `r_min`, `loopSinkScale` and `pyramidRelaxationVolume_SI` are no longer read.
- **Atomic volume** of the lattice (`Polycrystal::Omega`, 2.39e-29 m³ for the HEX lattice of `Zr_CD4opt.txt`).
- **Material file** `Zr_CD4opt.txt` and tutorial `spatialCDcvode`, with a plotting script.
- **Conversion of the immobile fields to discrete loops** at `clusterDiscretizationTime` ([section 3.6](#36-cluster-dynamics)).
- **Any number of loop families, in any crystal:** `MODELIB_CD_ISIZE` is any even number; the diffusion-anisotropy bias distinguishes the families by their Burgers vector.
- **Starting a study from another one.** The cluster-dynamics block of `evl_<n>.txt` is written at full precision, and the positions of the finite-element nodes in `evl/cdNodes.txt`. A run that finds a configuration file as `evl/evl_0.txt` starts from its fields.
- **Documentation on the Wiki**: at each push to `master` a workflow converts the Doxygen output to Wiki pages ([section 12](#12-documentation)).

## 2. Capabilities

| Capability | Where | Notes |
|---|---|---|
| Dislocation network topology | `include/LoopNetwork/`, `include/DislocationDynamics/` | loops, loop nodes and loop links; network nodes and segments shared between loops |
| Elastic fields of dislocations | `DislocationDynamicsBase/StraightDislocationSegment.h`, `StressStraight.h` | piecewise-straight segments; `_MODEL_NON_SINGULAR_DD_` selects classical (0), Cai (1, default) or Lazar (2) fields; `coreSize` in `DD.txt` |
| Glide solver | `GalerkinGlideSolver.h` | nodal velocities by Galerkin projection of the line velocity; periodic constraints by a null-space projection |
| Glide solver in Python | `PyGlideSolver.h`, `python/MLglideSolver.py` | `glideSolverType=pybind11`; the module path is `pyModuleName` in `DD.txt` |
| Climb solver | `GalerkinClimbSolver.h` | `climbSolverType=Galerkin`; needs `ClusterDynamics` in `physics`, since climb is driven by the point-defect concentrations |
| Mobility laws | `include/DislocationMobilities/` | per slip system: `viscousDrag`, `kinkPair`, `edgeScrew`, `interpolated`, `python` |
| Crystal structures and slip systems | `include/PolycrystallineMaterials/` | FCC, BCC, HEX and cubic fluorite; names in [section 3.3](#33-mobility-laws) |
| Gamma surfaces and second phases | `GammaSurface.h`, `SecondPhase.h` | FCC {111} and HEX basal from three energies; HEX prismatic from points; second phases `L12` (FCC), `sigma` and `chi` (BCC) |
| Polycrystals | `Grain.h`, `GrainBoundary.h`, `Polycrystal.h` | grains are mesh regions, each with its rotation `C2G<n>` |
| Cross slip | `CrossSlipModels.h` | `crossSlipModel`: 0 none, 1 deterministic (largest glide Peach–Koehler force) |
| Junctions and remeshing | `DislocationJunctionFormation.h`, `DislocationNetworkRemesh.h`, `DislocationNodeContraction.h` | segment lengths kept between `Lmin` and `Lmax` |
| Nucleation | `DislocationNucleation.h` | `bulkNucleationModel=1` inserts shear loops in mesh elements; a first model |
| Periodic domains | `GlidePlanes/PeriodicGlidePlane.h` | periodic faces chosen per mesh face (`periodicFaceIDs`); image sums set by `periodicImageSize`; `EwaldLengthFactor` |
| Finite domains | `include/FEM/`, `ElasticDeformationFEM.h` | `useFEM=1`; direct solvers (CHOLMOD, UMFPACK) or iterative solvers |
| Uniform load control | `DislocationDynamicsBase/UniformController.h` | stress control, strain control, or mixed, per Voigt component ([section 3.5](#35-load-control)) |
| Glide-plane noise | `GlidePlanes/AnalyticalSolidSolutionNoise.h`, `MDSolidSolutionNoise.h`, `MDStackingFaultNoise.h` | stress noise of a solid solution and stacking-fault energy noise, sampled on a grid (FFTW) from analytical or MD correlations |
| Stochastic force | `DD.txt` | `useStochasticForce`, `stochasticForceSeed` |
| Inclusions | `SphericalInclusion.h`, `PolyhedronInclusion.h` | eigendistortion, velocity reduction factor and second phase per inclusion |
| Cluster dynamics | `include/ClusterDynamics/` | mobile species: anisotropic diffusion, production from a dose rate, dissociation, second-order reactions, sinks. Immobile species: loop families with number density and defect content, nucleation, conversion to discrete loops. Growth strain ([section 3.6](#36-cluster-dynamics)) |
| Microstructure generation | `include/MicrostructureGeneration/` | shear, prismatic and Frank loops; stacking-fault tetrahedra; periodic dipoles; a planar loop from its nodes; spherical and polyhedral inclusions. By target density or individually |
| Time stepping | `DDtimeStepper.h` | `fixed` or `adaptive`; velocity filter; subcycling when `subcyclingBins` has more than one value |
| Lattice arithmetic | `include/Lattices/` | lattice and reciprocal vectors, rational directions, LLL reduction, CSL and DSCL of bicrystals |
| Visualization | `tools/DDqt/` | Qt 6 and VTK viewer of configurations, slip systems, glide planes, inclusions, quadrature data, fields and the mesh |

## 3. Methodology

### 3.1 The defective crystal

`DefectiveCrystal` is a container of physics. The `physics` line of `DefectiveCrystal.txt` lists the ones that exist in a run:

| Name in `physics` | Physics added |
|---|---|
| `DislocationDynamics` | dislocation network |
| `ClusterDynamics` | mobile and immobile point-defect species |
| `ElasticDeformation` | uniform elastic deformation with load control, or the finite-element elastic solution |
| `InclusionMicrostructure` | inclusions |

The names are separated by spaces and are not case sensitive. They are created, solved and written in the order given. At least one is required: the value `none` of the template must be replaced. `useFEM` enables the finite-element mesh for the physics that can use it.

One step of `DefectiveCrystal::runSingleStep` does the following, in this order:

1. **Solve.** Each physics solves for its own rates.
2. **Re-solve.** The solve is repeated `maxResolveSteps` times, so that coupled physics (climb and diffusion) can exchange their fields within the step.
3. **Time step.** Each physics proposes a time step; the smallest one is used, capped by `dtMax`.
4. **Output.** Every `outputFrequency` steps, the configuration and the averaged quantities are written.
5. **Update.** Each physics advances its configuration by the time step.

For the dislocation network, cross slip is carried out at the end of the solve, before the nodes move. The update then moves the nodes and executes the discrete events: remeshing, junction formation, remeshing again, and nucleation.

### 3.2 Dislocation glide

Each segment carries quadrature points. At a point with unit tangent $\boldsymbol{\xi}$ and Burgers vector $\mathbf{b}$, the Peach–Koehler force per unit length is

$$\mathbf{f}_{PK} = (\boldsymbol{\sigma}\,\mathbf{b}) \times \boldsymbol{\xi},$$

where $\boldsymbol{\sigma}$ is the sum of the stresses of all enabled physics. The total force adds the stacking-fault force and a line-tension force scaled by `alphaLineTension`:

$$\mathbf{f} = \mathbf{f}_{PK} + \mathbf{f}_{SF} + \mathbf{f}_{LT}.$$

The mobility law of the slip system turns the force and the local stress into a glide velocity $\mathbf{v}$ at that point. Nodal velocities $\mathbf{V}$ then follow from a Galerkin projection with the segment shape functions $\mathbf{N}$:

$$\left[\int \mathbf{N}^\top \mathbf{N}\, dL\right] \mathbf{V} = \int \mathbf{N}^\top \mathbf{v}\, dL .$$

With adaptive time stepping, $dt = \min(\texttt{dxMax}/v_{max},\ \texttt{dtMax})$.

### 3.3 Mobility laws

**Slip systems.** `enabledSlipSystems` in the material file lists the slip systems, by these names:

| `crystalStructure` | Slip-system names |
|---|---|
| `FCC` | `a/2<110>{111}` (full), `a/6<112>{111}` (Shockley partials), `a/3<112>{111}` (Kear) |
| `BCC` | `a/2<111>{110}`, `a<100>{110}` |
| `HEX` | `<a>{basal}`, `<Shockley>{basal}`, `<a>{prismatic}`, `<a>{pyramidal}` |
| `CubicFluorite` | `a/2<110>{111}`, `a/2<110>{110}`, `a/2<110>{100}` |

A name that is not in this table creates no slip system, without an error.

**One law per slip system.** Each enabled slip system needs a line `mobility_<slip system>=<type> <arguments>;`. The type is one of:

| Type | Arguments | Velocity |
|---|---|---|
| `edgeScrew` | two labels: edge, screw | $v = v_{screw}\cos^2\theta + v_{edge}\sin^2\theta$, with $\theta$ the angle between $\mathbf{b}$ and $\boldsymbol{\xi}$ |
| `interpolated` | `linear`, `logarithmic` or `lagrange`, then labels ending in `_<angle in degrees>` | interpolation over the character angle, periodic with period $\pi$; `logarithmic` interpolates $\log v$ |
| `python` | path of a Python module | the module's class `MobilitySolver` returns the velocity from `velocityPy(stress, b, xi, n, T, dL, dt)`; needs the pybind11 build |

Each label names a further line, `mobility_<label>=...`, with an elementary law:

| Elementary law | Arguments, in this order | Meaning |
|---|---|---|
| `viscousDrag` | `B0_SI B1_SI` | $v = \tau b / B$, with $B = B_0 + B_1 T$ in Pa·s |
| `kinkPair` | `h w B0_SI B1_SI Bk_SI dH0_eV p q Tf tauC_SI a0 a1 a2 a3 angle` | the kink-pair law of Po et al. (2016): kink height `h` and width `w` in units of $b$; drag $B_0 + B_1 T$ and kink drag `Bk`; kink-pair enthalpy `dH0` with exponents `p`, `q`; athermal temperature `Tf` as a fraction of $T_m$; Peierls stress `tauC`; non-Schmid coefficients `a0`–`a3`; the angle, in degrees, between the glide plane and the conjugate plane of the non-Schmid terms |

The law of tungsten shows the pattern:

```
enabledSlipSystems=a/2<111>{110};
mobility_a/2<111>{110}=edgeScrew a/2<111>{110}e a/2<111>{110}s;
mobility_a/2<111>{110}e= viscousDrag 4.26e-04 0.87e-06;
mobility_a/2<111>{110}s= kinkPair 0.94280904158 25 9.8e-4 0.0 8.3e-05 1.63 0.86 1.69 0.8 2.03e9 1.50 1.15 2.32 4.29 60.0;
```

`UO2.txt` shows the `interpolated` type, with four laws at 0°, 30°, 90° and 150°.

### 3.4 Elastic fields

The stress of the network is the sum over straight segments, in the non-singular form selected at compile time ([section 2](#2-capabilities)) with core width `coreSize`. In a finite domain, the infinite-medium fields are corrected by a finite-element solution on the mesh (the superposition principle). In a periodic domain, the fields of `periodicImageSize` images are summed along each periodic shift vector.

The derivations are in `manual_beta/MoDELib_manual.pdf`: eigendistortion theory, discrete loops in anisotropic and isotropic media, the solid angle as a line integral, and the fields of straight segments.

### 3.5 Load control

`ElasticDeformation.txt` gives an applied stress $\boldsymbol{\sigma}_0$, an applied strain $\boldsymbol{\varepsilon}_0$, their rates, and a machine stiffness ratio $\alpha_i$ per Voigt component (order 11, 22, 33, 12, 23, 13). With $D_i = \alpha_i C_{ii}$, material stiffness $\mathbf{C}$ and plastic strain $\boldsymbol{\varepsilon}_p$, the stress in the volume is

$$\boldsymbol{\sigma} = \mathbf{C}\,(\mathbf{C}+\mathbf{D})^{-1}\left[\boldsymbol{\sigma}_0 + \mathbf{D}\,(\boldsymbol{\varepsilon}_0 - \boldsymbol{\varepsilon}_p)\right].$$

| `stiffnessRatio` component | Result |
|---|---|
| $\alpha_i = 0$ | stress control: $\sigma_i = \sigma_{0i}$ |
| $\alpha_i \to \infty$ (for example `1e20`) | strain control: $\varepsilon_i = \varepsilon_{0i}$ |

The plastic strain is the sum over all physics, so the growth strain of cluster dynamics enters the load control as dislocation slip does. With `useElasticDeformationFEM=1` the uniform controller is replaced by the finite-element elastic solution.

### 3.6 Cluster dynamics

Cluster dynamics is solved on the finite-element mesh (quadratic tetrahedra) when `useFEM=1` and `useClusterDynamicsFEM=1`, in a domain that is not periodic. The equations are those of the D1/M1 report [`docs/reports/GW_Phase4_D1M1.pdf`](docs/reports/GW_Phase4_D1M1.pdf), with the changes listed below.

**Fields.**

- *Mobile species*, $m = 1 \dots$ `MODELIB_CD_MSIZE`. `mobileSpeciesVector` gives the signed point-defect content of each cluster: $-1$ for a vacancy, $+1$, $+2$, $+3$ for interstitial clusters. The stored field $c_m$ is the defect content per atomic site.
- *Immobile species*, $k \in \{c, a_1, a_2, a_3\}$: the clusters per atom $n_k$ and the defect content per atom $c_k$ of each loop family, with the Burgers vector `immobileSpeciesBurgers` and the sign $p_k$ of `immobileSpeciesVector` ($-1$ vacancy type, $+1$ interstitial type). The mean size is $m_k = c_k/n_k$. The `evl` files hold the number density per unit volume, $n_k/\Omega$, in units of $b^{-3}$.
- $\Omega$ is the atomic volume of the lattice (`Polycrystal::Omega`, $\sqrt2/2\,b^3$ for ideal HEX).

**Mobile species.** Each step solves the quasi-steady balance (report, Eq. 1) for the mobile contents, with:

- **diffusion** $\mathbf{D}_m = \mathbf{D}_{0m}\exp(-\mathbf{E}^m_m/k_BT)$, component by component in the crystal frame, so that diffusion can be anisotropic; $\bar D_m = (\det \mathbf{D}_m)^{1/3}$;
- **production** $G f_m$: dose rate times `mobileSpeciesCascadeFractions`$_m$; the rest of the cascade defects is born in the immobile families;
- **dissociation** of a cluster into a smaller cluster and a single defect, with the binding energy `mobileSpeciesBindingEnergy_eV`;
- **second-order reactions** between the pairs listed in `reactionPrefactorMap`, with rate constant $p\,4\pi(r_a+r_b)(\bar D_a+\bar D_b)/\Omega$;
- **network sinks** `mobileDislocationSinkEfficiency`$_m\times$`dislocationSinkDensity_SI`$\times\bar D_m$;
- **loop sinks** $\sum_k \Gamma_{km}$ (below), and the source of single vacancies and interstitials released by the immobile families, $A^c + A^{emit}$ (report, Eq. 55).

The reactions make the problem nonlinear; it is solved by Newton iterations to a relative change of $10^{-5}$. On the boundary, the concentration is the equilibrium value under the local stress. What the immobile families gain or lose through absorption, dissolution and emission is what the mobile species lose or gain.

**Cluster shape.** A cluster of $m$ defects is a bi-pyramid when small and a loop when large, with the loop fraction

$$S_k(m) = \frac{1}{1+\exp\left[-(m-m_{0k})/w_k\right]},\qquad m_{0k} = \tfrac{1}{2}(\texttt{shapeTransitionMinSize}_k+\texttt{shapeTransitionMaxSize}_k),\qquad w_k = \texttt{shapeTransitionWidth}\,(\texttt{shapeTransitionMaxSize}_k-\texttt{shapeTransitionMinSize}_k),$$

the bi-pyramid and loop radii $R^{pyr} = (m\Omega/\sqrt{8})^{1/3}$ and $R^{loop} = \sqrt{m\Omega/(\pi|\mathbf{b}_k|)}$, the effective radius $R_k = (1-S_k)R^{pyr}_k + S_kR^{loop}_k$, and the capture length $L_k = (1-S_k)\,4\pi R^{pyr}_k + S_k\,2\pi R^{loop}_k$, which blends the sink strength of a sphere and of a loop (as in MoDELib3).

**Immobile species.** At each node (report, Eqs. 36 and 40),

$$\frac{dn_k}{dt} = \frac{s^{nuc}_k + G^L_k}{n^{nuc}_k} - a^n_k - \Delta N_k ,\qquad \frac{dc_k}{dt} = \sum_m p_k p_m\,\Gamma_{km} + G^L_k + s^{nuc}_k - A^c_k - A^{emit}_k - \Delta C_k .$$

| Term | Expression | Input |
|---|---|---|
| Absorption | $\Gamma_{km} = g_{km}\,\bar D_m\,Z_{km} L_k n_k c_m/\Omega$; like defects add to the content, opposite defects remove from it | |
| Bias | $Z_{km} = Z^0_{km} Z^{DAD}_k(p_m)$, $p_m = (D_{m,33}/D_{m,11})^{1/6}$: $p_m$ for basal loops, $(p_m+p_m^{-2})/2$ for prismatic loops | `immobileSpeciesBias` (one row per immobile family) |
| Size gate | $g_{km} = 1$ for like defects; for opposite defects a smooth step in $m_k$, closed at the smallest immobile cluster of the family's sign (the largest mobile cluster of that sign plus one: 2 for vacancy families, 4 for interstitial families here) and fully open at twice that size. The same factor multiplies the loss of the mobile species | none (from `mobileSpeciesVector`) |
| Cascade nucleation | $G^L_k = G(\varepsilon_{iL}\chi^i_k + \varepsilon_{vL}\chi^v_k)$, with $\varepsilon_{iL}$, $\varepsilon_{vL}$ the fractions of the interstitials and vacancies that `mobileSpeciesCascadeFractions` leaves out | |
| Clustering nucleation | $s^{nuc}_k = R^{2,i}\chi^i_k + R^{2,v}\chi^v_k$: the content of the reactions between mobile species whose product is larger than the largest, or smaller than the smallest, mobile cluster | `reactionPrefactorMap` (zero entries switch it off) |
| Number of new clusters | content divided by $n^{nuc}_k$ | `nucleatedClusterSize` |
| Orientation of the nucleation | $\chi_k \propto \exp(\boldsymbol\sigma:\boldsymbol\Omega^k/k_BT)$, normalized over the families of the same sign, with the formation-volume tensor of a cluster of $n^{nuc}$ defects; equal fractions without stress | |
| Dissolution | $a^n = \nu_a n$, $A^c = \nu_a \min(n^{nuc}n, c)$, $\nu_a = \nu_a^0\exp(-E_a/k_BT)$. The clusters that dissolve are counted as embryos of the nucleated size: grown loops keep their excess content | `clusterDissolutionRatePrefactor`, `clusterDissolutionEnergy` |
| Thermal emission of vacancies | detailed balance with the absorption of vacancies, $A^{emit}_k = \bar D_v Z_{kv} L_k n_k/\Omega\;e^{-E_b/k_BT}$, with a constant $E_b$ per vacancy family | `clusterBindingEnergy_eV` |
| Coalescence | $\Delta N = (\nu_{LL}\phi_{LL} + \nu_{LN}\phi_{LN})\,n$, $\Delta C = \nu_{LN}\phi_{LN}\,c$, with $\phi_{LL} = 1-e^{-\kappa_{LL}\frac43\pi R_{LL}^3N}$, $\phi_{LN} = 1-e^{-\kappa_{LN}\pi R_{LN}^2\rho_N}$, $\nu_{LL} = c_{LL}v_{abs}N^{1/3}$, $\nu_{LN} = c_{LN}v_{abs}\sqrt{\rho_N}$, $N = n/\Omega$, and $R_{LL}$, $R_{LN}$ the effective radius capped at $N^{-1/3}$ and $\rho_N^{-1/2}$ | `cLL`, `cLN`, `kappaLL`, `kappaLN`, `dislocationSinkDensity_SI` |
| Climb speed of coalescence | $v_{abs} = (dR_k/dm)\,\dot G^{net}_k/n_k$, the rate of change of the effective radius, with $dR_k/dm$ including the change of shape and $\dot G^{net}_k = \max(0, \sum_m p_kp_m\Gamma_{km} - A^{emit}_k)$ the net growth of the family. A shrinking or empty family does not coalesce | |

The network density $\rho_N$ is a constant. `captureZoneFactor` is not read: the two overlap coefficients and the two rate coefficients are given directly.

**Numerical method.** With the mobile field frozen over the step, the equations of a node involve the state of that node only, so the mesh gives as many independent systems of `MODELIB_CD_ISIZE` equations as it has nodes. They are stiff: at 573 K the lifetime of the vacancy clusters is four hours and a step is days. Each is integrated by CVODE (SUNDIALS, variable-order BDF) with a dense Newton iteration, relative tolerance $10^{-8}$; the nodes are distributed over the OpenMP threads. The time step `dtMax` only sets how often the mobile species are solved again. A build without SUNDIALS cannot run cluster dynamics with immobile species.

**Growth strain.** The average plastic distortion of cluster dynamics is the volume average of

$$\boldsymbol{\beta}^P = \sum_j s_j c_j \left[(1-S_j)\,\frac{\mathbf{I}}{3} + S_j\,\hat{\mathbf{b}}_j\otimes\hat{\mathbf{b}}_j\right],$$

with $c_j$ the defect content (defects per atom): a cluster changes the volume by one atomic volume per defect whatever its shape, isotropically as a bi-pyramid and along its Burgers vector as a loop (Christiaen et al., Acta Mater. 179 (2019) 93), so the sigmoid changes only the anisotropy. To this is added a relaxation part, the same expression times `immobileSpeciesRelRelaxVol` for each family (0.1 for `c`, 0.01 for `a` in `Zr_CD4opt.txt`), and the relaxation volumes of the mobile defects (`mobileSpeciesRelRelaxVol`). The input `pyramidRelaxationVolume_SI` is no longer read. Vacancy loops on the basal plane and interstitial loops on the prismatic planes therefore produce the anisotropic shape change known as irradiation growth.

**Conversion to discrete loops.** A field of clusters cannot glide, intersect a dislocation or be seen as an obstacle; a discrete loop can. At the first update after `clusterDiscretizationTime`, with `DislocationDynamics` in `physics` and a network that holds no loop yet, each family is converted as follows.

1. *Number.* The number of clusters in the crystal is not a parameter: it is the integral of the number density, $\mathcal{N}_j = \int N_j\,dV$. The number of discrete loops is $\mathcal{N}_j/\Lambda$, rounded, where $\Lambda$ is `clusterDiscretizationFactor`, the number of clusters that one discrete loop stands for. $\Lambda = 1$ gives the population of the field; a larger value makes the population affordable for dislocation dynamics.
2. *Position.* Each loop is placed in a mesh element drawn with probability $\mathcal{N}_e/\mathcal{N}_j$, the share of the clusters that the element holds, and at a uniformly random point of that element. The loops therefore follow the spatial variation of the field.
3. *Size.* A loop holds the defects of the clusters it stands for, $m = \Lambda\,\mathcal{C}_e/\mathcal{N}_e$, with $\mathcal{C}_e = \int_e c_j/\Omega\,dV$. A prismatic loop of radius $r$ and Burgers vector $\mathbf{b}_j$ stores $\pi r^2 |\mathbf{b}_j|/\Omega$ defects, so

$$r = \sqrt{\frac{m\,\Omega}{\pi\,|\mathbf{b}_j|}} .$$

   Rounding the number of loops changes the total by a fraction of a loop. All radii are then multiplied by one factor, so that the loops hold exactly the defects of the field: $\sum_i \pi r_i^2|\mathbf{b}_j|/\Omega = \int c_j/\Omega\,dV$.
4. *Coalescence.* Two loops of the same family overlap when the distance of their centers is less than the sum of their radii. They are merged only if the centers also lie within one spacing of the habit planes along the normal: loops on different parallel planes cannot become one loop without climbing. The merged loop has the area of the two, $r = \sqrt{r_1^2+r_2^2}$, which conserves the defects, and its center is their area-weighted centroid. Merging at constant total area cannot remove overlap: with $N$ loops in a volume $V$, $r \propto N^{-1/2}$ and the spacing $d = (V/N)^{1/3} \propto N^{-1/3}$, so the ratio $2r/d \propto N^{-1/6}$ grows with every merge. A pass that would make the mean diameter larger than the mean spacing is therefore rejected, and the family is reported as saturated.
5. *Shape and crystallography.* A loop is a regular polygon of 12 sides in the lattice plane nearest to its center, with normal along $\mathbf{b}_j$ (a pure prismatic loop). A polygon of circumradius $R$ has the area $\tfrac{n}{2}R^2\sin(2\pi/n)$, 4.5 % less than the disc for $n=12$. The stored defects are proportional to the area, so the circumradius is $R = r\sqrt{2\pi/(n\sin(2\pi/n))} = 1.023\,r$. The Burgers vector is the lattice vector `immobileSpeciesBurgers`, along the area normal for a vacancy loop and opposite to it for an interstitial loop. The loops are sessile: they move by climb.
6. *Loops that are not created.* A loop smaller than `minimumLoopSize` is not created. A loop with a node outside its grain is moved toward the center of the grain in steps of 10 % of the distance, and is not created if it still does not fit after 30 steps. If the created loops hold the fraction $\varphi$ of the defects of the family, the fields $N_j$ and $c_j$ are both multiplied by $1-\varphi$. The same factor on both leaves the cluster size of what remains unchanged.

For each family the code prints three totals before and after the conversion: the number of loops, the stored defects, and the sink strength $\sum 2\pi r$. The defects are conserved by construction. The number is conserved for $\Lambda=1$, up to rounding and coalescence. The sink strength is not conserved when clusters are lumped: at constant area it scales as $\Lambda^{-1/2}$.

Two consequences should guide the choice of the conversion time:

- The eigenstrain of a discrete loop is that of a loop, $\mathbf{b}\otimes\hat{\mathbf{b}}\,\pi r^2/V$. A family that is still in the bi-pyramid range of the sigmoid has an isotropic eigenstrain in the fields, and its conversion changes the growth strain. Set `minimumLoopSize` to the loop radius at `shapeTransitionMaxSize` to leave such a family in the fields.
- After the conversion the loops evolve by dislocation dynamics, and the time step is the smaller of `dtMax` and the step that the network proposes. The `dtMax` of a run of cluster dynamics alone (days) is far larger than what the climb of discrete loops needs.

**Coupling with dislocation dynamics.** Climb needs cluster dynamics: the climb velocity of a segment follows from the difference between the concentration in equilibrium with its climb force and the concentration that cluster dynamics gives at that point. When the network holds discrete loops, the loop sink term is left out of the mobile equation and the immobile fields are frozen; this applies to the share of a family that a conversion left in the fields.

**Time step.** Cluster dynamics proposes `dtMax`. A run of cluster dynamics alone uses `dtMax` as its step.

### 3.7 Units

Input files give material constants in SI units, in variables whose names end in `_SI`. Internally, and in the output files, quantities are dimensionless:

| Quantity | Unit |
|---|---|
| length | Burgers vector magnitude $b$ |
| stress | shear modulus $\mu = \mu_0 + \mu_1 T$ |
| speed | shear-wave speed $c_s = \sqrt{\mu/\rho}$ |
| time | $b/c_s$ |
| energy | $\mu b^3$ |

The values of these units are printed when a simulation starts. For zirconium ($b = 0.3233$ nm, $\mu = 33$ GPa, $\rho = 6520$ kg/m³), one time unit is $1.437\times10^{-13}$ s.

### 3.8 The input files

A simulation folder holds `inputFiles/`, `evl/` and `F/`. These files in `inputFiles/` define the run:

| File | Declares | Read when |
|---|---|---|
| `DefectiveCrystal.txt` | the `physics`, `useFEM`, `Nsteps`, `dtMax`, `maxResolveSteps`, the restart step, the periodic images, the output frequency and format | always |
| `polycrystal.txt` | the material file, temperature, mesh file, crystal rotations `C2G<n>`, mesh mapping `F` and `X0`, and periodic faces | always |
| the material file, for example `W.txt` | crystal structure, slip systems and their mobility laws, elastic constants, gamma surfaces, noise, second phases, cluster-dynamics parameters | always |
| `initialMicrostructure.txt` | a list of microstructure files, each declaring one family of initial defects | by `microstructureGenerator` |
| `DD.txt` | the numerical controls of dislocation dynamics | with `DislocationDynamics` |
| `ElasticDeformation.txt` | the applied stress and strain, their rates, the stiffness ratio, `useElasticDeformationFEM` | with `ElasticDeformation` |
| `ClusterDynamics.txt` | `useClusterDynamicsFEM` | with `ClusterDynamics` |

Mesh nodes $\mathbf{X}$ are mapped to $\mathbf{x} = \mathbf{F}(\mathbf{X} - \mathbf{X}_0)$, so one unit mesh serves domains of any size and shape.

## 4. Repository structure

```
MoDELib/
├── include/                      # headers (397 files)
│   ├── DislocationDynamics/      # DefectiveCrystal, network, solvers, cross slip, junctions, remeshing, nucleation
│   ├── DislocationDynamicsBase/  # microstructure interfaces, segment fields, inclusions, load controller,
│   │                             # cluster-dynamics parameters
│   ├── DislocationDynamicsIO/    # configuration (evl) and auxiliary (ddAux) files, parameters
│   ├── MicrostructureGeneration/ # generators and specifications of initial microstructures
│   ├── DislocationMobilities/    # mobility laws and their selector
│   ├── PolycrystallineMaterials/ # crystals, slip systems, grains, grain boundaries, gamma surfaces
│   ├── GlidePlanes/              # glide planes, periodic glide planes, noise
│   ├── ClusterDynamics/          # mobile and immobile species on the FEM mesh, reactions, sink rates
│   ├── ElasticDeformation/  InclusionMicrostructure/  DiscreteCrackMechanics/
│   ├── LoopNetwork/              # the generic loop-network topology
│   ├── FEM/  Mesh/  Quadrature/  # finite elements, simplicial meshes (Gmsh), Gauss–Legendre rules
│   ├── Lattices/  Lattices_NEW/  # lattice arithmetic
│   ├── Geometry/  Math/  IO/  Utilities/
│   ├── ParticleInteraction/  MPI/  DDvtk/
├── src/                          # implementations (168 files), one folder per header folder
├── tools/
│   ├── DDomp/                    # the simulation driver
│   ├── MicrostructureGenerator/  # writes the initial configuration
│   ├── DDqt/                     # Qt 6 / VTK viewer
│   ├── pyMoDELib/                # pybind11 module
│   ├── DDconverter/              # binary configuration files to text
│   ├── DDsegments/               # reads segments from an OVITO DXA .vtk file
│   └── QuadratureGeneration/     # generates the Gauss–Legendre tables
├── python/                       # modlibUtils.py (input and output helpers), readers of the evl, ddAux and F files
│                                 # (readEVL.py, readAUX.py, readF.py, readFile.py, getLoops.py), examples of the
│                                 # bindings, MLglideSolver.py, DDMobilityPy.py
├── Library/                      # templates copied into each simulation
│   ├── DefectiveCrystal/         # DefectiveCrystal.txt
│   ├── DislocationDynamics/      # DD.txt
│   ├── ElasticDeformation/       # ElasticDeformation.txt
│   ├── ClusterDynamics/          # ClusterDynamics.txt
│   ├── Materials/                # Al, AlMg5, AlMg10, AlMg15, Co, Cu, Ni (FCC); W, FeCrAl_Fe (BCC);
│   │                             # Zr, Zr_CD2, Zr_CD4, Zr_BMD19, Zr_CD3_BMD19, Zr_CD4opt, Zr4_Fitted (HEX); UO2 (fluorite)
│   ├── Meshes/                   # Gmsh meshes: unit cubes, cylinders, a disk, bicrystal, polycrystals
│   ├── Microstructures/          # one template per defect family
│   └── GlidePlaneNoise/          # noise templates and MD correlation data
├── tutorials/                    # annealing, dipoleNoise, dipoles, spatialCDcvode, uniformLoadController
├── testsPy/                      # verification scripts (section 10)
├── docs/                         # reports and manuals (the D1/M1 report of the cluster-dynamics model)
├── manual_beta/                  # the theory manual (tex, pdf)
├── doxygen/                      # Doxyfile and pages of the API documentation
├── scripts/                      # shell scripts for cleaning and for image and video conversion
├── .github/workflows/            # workflow.yml: Doxygen pages; container.yml: container image
├── CMakeLists.txt  Dockerfile  Contributors.txt  README.md  git_basic_commands.txt
```

## 5. Installation

### Dependencies

| Dependency | Needed for | Required |
|---|---|---|
| C++20 compiler with OpenMP | everything | yes |
| CMake 3.20 or later | the build | yes |
| Eigen 3 | linear algebra | yes |
| FFTW 3 | glide-plane noise | yes: configuration fails without it |
| Boost | headers | no (a warning if absent) |
| SuiteSparse (CHOLMOD, UMFPACK) | direct FEM solvers | no (a warning if absent) |
| Python 3 with NumPy | the tutorials' `generateInputFiles.py` | for the tutorials |
| Matplotlib | the tutorials' plotting scripts | for the plots |
| pybind11 | `pyMoDELib`, the Python glide solver and mobility | when `USE_PYBIND11=ON` (the default); skipped with a warning if absent |
| SUNDIALS 7 (CVODE) | integration of the immobile cluster-dynamics species | for cluster dynamics with immobile species (the build succeeds without it, with a warning) |
| Qt 6 and VTK 9.4.2 or later | `DDqt` | when `BUILD_DDQT=ON` (the default) |

### Build

```bash
git clone https://github.com/MechanicsOfDefects/MoDELib.git
cd MoDELib
cmake -S . -B build
cmake --build build -j
```

The build type is `Release`, with `-O3 -march=native -fopenmp`. Eigen is found through its CMake package, or in the usual MacPorts, Homebrew and Linux locations.

| Option | Default | Effect |
|---|---|---|
| `USE_PYBIND11` | `ON` | builds `pyMoDELib` and the Python solver hooks; needs Python 3 and pybind11 |
| `BUILD_TOOLS` | `ON` | builds `microstructureGenerator`, `DDomp`, `DDqt` and `pyMoDELib`. `DDconverter`, `DDsegments` and `QuadratureGeneration` are built separately |
| `BUILD_DDQT` | `ON` | builds the viewer; needs Qt 6 and VTK. Configuration fails if they are missing, so set it to `OFF` on a machine without them |
| `MODELIB_MARCH` | `native` | the value of `-march`. Use for example `x86-64-v2` for binaries that must run on other machines |
| `MODELIB_CD_MSIZE` | `4` | number of mobile cluster-dynamics species |
| `USE_SUNDIALS` | `ON` | uses SUNDIALS 7 (CVODE) for the immobile cluster-dynamics species when it is found; a warning otherwise, and cluster dynamics with immobile species cannot run. Pass `-DCMAKE_PREFIX_PATH=<SUNDIALS install>` if it is not in a system location |
| `MODELIB_CD_ISIZE` | `8` | number of immobile cluster-dynamics fields: two per loop family, so an even number. `0` for none |
| `EIGEN3_INCLUDE_DIRS` | found automatically | path of the Eigen headers |
| `CMAKE_PREFIX_PATH` | | where Qt 6 is installed, for example `~/Qt/6.7.3/macos/lib/cmake` |
| `Python_EXECUTABLE`, `pybind11_DIR` | | the interpreter and the pybind11 installation to use |

Common configurations:

```bash
# command-line tools only (no Python bindings, no viewer); works on a plain Debian or Ubuntu machine
cmake -S . -B build -DUSE_PYBIND11=OFF -DBUILD_DDQT=OFF

# library only
cmake -S . -B build -DBUILD_TOOLS=OFF -DUSE_PYBIND11=OFF

# macOS with MacPorts, clang 18 and Python 3.11
cmake -S . -B build -DCMAKE_CXX_COMPILER=clang++-mp-18 -DCMAKE_PREFIX_PATH=$HOME/Qt/6.7.3/macos/lib/cmake \
      -DPython_EXECUTABLE=/opt/local/bin/python3.11 -Dpybind11_DIR=$(/opt/local/bin/python3.11 -m pybind11 --cmakedir)
```

On Debian or Ubuntu the packages for the library and the command-line tools are `build-essential cmake libeigen3-dev libfftw3-dev libboost-dev libsuitesparse-dev`. That configuration builds with g++ 13.3 and CMake 3.28 on Ubuntu 24.04.

**The species counts must match the material file.** The numbers of cluster-dynamics species are fixed when the library is compiled, and the material file must give arrays of the same sizes. The defaults (4 mobile, 8 immobile) are those of `Zr_CD4opt.txt` and of the `spatialCDcvode` tutorial. The `annealing` tutorial uses `Zr_BMD19.txt`, with one mobile species and no immobile species, and needs its own build:

```bash
cmake -S . -B build_1_0 -DMODELIB_CD_MSIZE=1 -DMODELIB_CD_ISIZE=0 -DUSE_PYBIND11=OFF -DBUILD_DDQT=OFF
cmake --build build_1_0 -j
```

`MODELIB_CD_ISIZE` is twice the number of loop families of the material file: any even number, or 0. The diffusion-anisotropy bias of a family is that of a loop normal to the third lattice vector when its Burgers vector is along that vector (the basal loops of a HEX crystal), and that of a loop containing it otherwise. With isotropic diffusion the bias is 1 for every family, in any crystal.

### SUNDIALS

SUNDIALS provides CVODE, the integrator of the immobile cluster-dynamics species ([section 3.6](#36-cluster-dynamics)). It is required for cluster dynamics with immobile species: MoDELib builds without it, but such a run cannot proceed. Version 7 is required.

From a package manager, where the packaged version is 7 or later:

```bash
sudo apt-get install libsundials-dev     # Ubuntu 25.04 or later, Debian 13 or later. Ubuntu 24.04 has version 6, which is too old
brew install sundials                    # macOS with Homebrew
```

From source, in a folder of your own. No administrator rights are needed, and only CVODE and its serial vector are built:

```bash
git clone --depth 1 --branch v7.1.1 https://github.com/LLNL/sundials.git sundials-src
cmake -S sundials-src -B sundials-build -DCMAKE_INSTALL_PREFIX=$HOME/sundials-install -DBUILD_ARKODE=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_CVODES=OFF -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_INSTALL=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build sundials-build -j
cmake --install sundials-build
```

Then configure MoDELib with the install prefix, unless SUNDIALS is in a system location:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=$HOME/sundials-install
cmake --build build -j
```

The configuration prints `SUNDIALS_VERSION:7.1.1` when it is found, and a warning when it is not. To check a build, run two steps of the `spatialCDcvode` tutorial: the log of each step holds the line `immobile ODE solver (CVODE)`. `-DUSE_SUNDIALS=OFF` builds without it even when it is installed. On Windows, install it inside WSL, where MoDELib is built.

### Container image

Each release is also published as a container image with the command-line tools, `Library/`, `python/` and the tutorials:

```bash
docker pull ghcr.io/mechanicsofdefects/modelib:v2.2.0
docker run -it --rm ghcr.io/mechanicsofdefects/modelib:v2.2.0
```

The container starts in `/opt/MoDELib/tutorials`, with `DDomp` and `microstructureGenerator` on the `PATH`. They are built with the default species counts (4 and 8) and `-march=x86-64-v2`, so that the binaries do not depend on the build machine. `DDqt` and `pyMoDELib` are not in the image, and it is built without SUNDIALS: cluster dynamics with immobile species does not run there.

The image is built from the `Dockerfile` at the root of the repository. To build it locally:

```bash
docker build -t modelib .
```

The build ends with a smoke test: two steps of the `dipoleNoise` tutorial. For an image with other species counts, for example for the `annealing` tutorial, give the counts and the tests that apply to them:

```bash
docker build -t modelib:annealing --build-arg CD_MSIZE=1 --build-arg CD_ISIZE=0 --build-arg SMOKE_TESTS=dipoleNoise .
```

## 6. Quick start

A run has three steps: generate the input files, generate the initial microstructure, and evolve it.

```bash
cd tutorials/dipoleNoise
python3 generateInputFiles.py                                         # writes inputFiles/, creates evl/ and F/
../../build/tools/MicrostructureGenerator/microstructureGenerator .   # writes evl/evl_0.txt
../../build/tools/DDomp/DDomp .                                       # runs Nsteps steps
```

Both executables take the simulation folder as their only argument. The default is the current folder.

`generateInputFiles.py` shows the pattern every simulation follows. It copies a template from `Library/`, then changes the copy:

```python
import sys
sys.path.append("../../python/")
from modlibUtils import *          # also brings in os, shutil and numpy as np

for folder in ['evl', 'F', 'inputFiles']:
    os.makedirs(folder, exist_ok=True)

# the physics and the run controls
shutil.copy2('../../Library/DefectiveCrystal/DefectiveCrystal.txt', 'inputFiles/DefectiveCrystal.txt')
setInputVariable('inputFiles/DefectiveCrystal.txt', 'physics', 'DislocationDynamics ElasticDeformation')
setInputVariable('inputFiles/DefectiveCrystal.txt', 'useFEM', '0')
setInputVariable('inputFiles/DefectiveCrystal.txt', 'Nsteps', '10000')
setInputVariable('inputFiles/DefectiveCrystal.txt', 'outputFrequency', '100')
setInputVector('inputFiles/DefectiveCrystal.txt', 'periodicImageSize', np.array([2, 2, 2]), 'images in each direction')

# numerical controls of dislocation dynamics
shutil.copy2('../../Library/DislocationDynamics/DD.txt', 'inputFiles/DD.txt')
setInputVariable('inputFiles/DD.txt', 'glideSolverType', 'Galerkin')
setInputVariable('inputFiles/DD.txt', 'Lmin', '5')
setInputVariable('inputFiles/DD.txt', 'Lmax', '20')

# glide-plane noise
shutil.copy2('../../Library/GlidePlaneNoise/AnalyticalSolidSolutionNoise.txt', 'inputFiles/AnalyticalSolidSolutionNoise.txt')
setInputVariable('inputFiles/AnalyticalSolidSolutionNoise.txt', 'MSSS_SI', '0.45e18')

# material
shutil.copy2('../../Library/Materials/AlMg15.txt', 'inputFiles/AlMg15.txt')
setInputVariable('inputFiles/AlMg15.txt', 'glidePlaneNoise', 'AnalyticalSolidSolutionNoise.txt')
b_SI = getValueInFile('inputFiles/AlMg15.txt', 'b_SI')

# domain: a box mapped from the unit-cube mesh
shutil.copy2('../../Library/Meshes/unitCube24.msh', 'inputFiles/unitCube24.msh')
pf = PolyCrystalFile('AlMg15.txt')
pf.absoluteTemperature = 300
pf.meshFile = 'unitCube24.msh'
pf.alignToSlipSystem0 = 1                                      # orient the crystal and the box along slip system 0
pf.boxScaling = np.array([102.4e-9, 102.4e-9, 102.4e-9]) / b_SI   # edge lengths in units of b
pf.periodicFaceIDs = np.array([-1])                            # -1: every face periodic
pf.write('inputFiles')
```

The script goes on to set the applied strain rate in `ElasticDeformation.txt`, and to declare a periodic dipole in `periodicDipoleIndividual.txt`, which `initialMicrostructure.txt` then lists.

Three points about `setInputVariable`:

- it changes nothing, and gives no warning, when the variable is not in the file, so a misspelt or renamed variable is silently ignored;
- it replaces the old value everywhere in the line, the comment included;
- for a variable whose value is empty in the template (such as `nodalVelocityConstraints=;`), use `setInputVector`.

The same objects are available from Python through `pyMoDELib`:

```python
import sys, os
import numpy as np
sys.path.append("../build/tools/pyMoDELib")
import pyMoDELib

ddBase = pyMoDELib.DislocationDynamicsBase(os.path.abspath("../tutorials/uniformLoadController"))   # after generateInputFiles.py

generator = pyMoDELib.MicrostructureGenerator(ddBase)
spec = pyMoDELib.ShearLoopIndividualSpecification()
spec.slipSystemIDs = [0, -1]
spec.loopRadii = [27.0e-8, 27.0e-8]
spec.loopCenters = np.array([[200.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
spec.loopSides = [10, 10]
generator.addShearLoopIndividual(spec)

crystal = pyMoDELib.DefectiveCrystal(ddBase)
crystal.initializeConfiguration(generator.configIO)
crystal.runSteps()            # or runSingleStep(); displacement(points) and stress(points) return the fields
```

`python/defectFields.py` gives the complete example. The module also exposes the mobility laws (`DislocationMobilitySelector`), the noise generators, `GammaSurface`, `SlipSystem`, `SingleCrystalBase`, `StressStraight` and the interpolants; `testsPy/` shows their use.

## 7. Tutorials

Each tutorial is a folder with one `generateInputFiles.py`.

| Tutorial | Physics | Material | Domain | Initial microstructure | Exercises |
|---|---|---|---|---|---|
| `dipoleNoise` | dislocations, elastic deformation | AlMg15, FCC, `a/2<110>{111}` | about 100 nm box, periodic, 300 K | one periodic dipole | glide through analytical solid-solution noise under a strain rate (strain control); a nodal velocity constraint |
| `dipoles` | dislocations, elastic deformation | W, BCC, `a/2<111>{110}` | 500 × 100 × 500 nm box, periodic, 300 K | one periodic dipole | the `edgeScrew` law with a kink-pair screw under a stress ramp; elastic energy output |
| `uniformLoadController` | dislocations, elastic deformation | W, BCC, `a/2<111>{110}` | 1 µm box, periodic, 300 K | prismatic loops and interstitial Frank loops, $10^{13}$ m$^{-2}$ each | the uniform load controller under a constant stress of $0.01\,\mu$ |
| `spatialCDcvode` | cluster dynamics | Zr (`Zr_CD4opt.txt`), HEX | 1 µm box, not periodic, 573 K | none: the grain starts without loops | nucleation, growth, dissolution and coalescence of the loops, integrated by CVODE; needs a build with SUNDIALS |
| `annealing` | dislocations, cluster dynamics, elastic deformation | Zr (`Zr_BMD19.txt`), HEX | 1 µm box, not periodic, 700 K | vacancy and interstitial Frank loops, $10^{13}$ m$^{-2}$, radius 50 nm | loop annealing by climb coupled to vacancy diffusion, with the finite-element elastic and diffusion solutions |

Notes:

- **`spatialCDcvode`** takes 30 steps of 0.1 dpa ($10^{-7}$ dpa/s), 3 dpa in all. `plotClusterDynamics.py` plots the growth strain against dose with the measurements of Carpenter and Rogerson, the mobile concentrations and the loop populations against dose, their profiles against the distance from the grain boundary, and the distribution of the loop diameters; it needs Matplotlib.
- **`annealing`** needs a build with one mobile species and no immobile species ([section 5](#5-installation)). `plotLoops.py` plots the loop radii and the dislocation density against time; it needs `pyMoDELib` and Matplotlib.
- The simulation time of the dislocation tutorials is set by `Nsteps` in the script.

## 8. Outputs

A run writes into the simulation folder:

| File | Contents |
|---|---|
| `evl/evl_<runID>.txt` or `.bin` | the configuration: network nodes, loops, loop links, loop nodes, spherical and polyhedral inclusions, the finite-element displacement, and the cluster-dynamics fields |
| `evl/ddAux_<runID>.txt` or `.bin` | auxiliary data of the dislocation network: quadrature points with their stress, forces, velocity and point-defect concentrations (with `outputQuadraturePoints=1`) |
| `F/F_0.txt` | one row per output step: `runID`, time, `dt`, the average plastic distortion and its rate, then the columns added by each physics |
| `F/F_labels.txt` | the label of each column of `F_0.txt`, one per line |

Further points:

- **Text layout of `evl`.** Ten lines give the number of rows of each block; the blocks follow in the order of the table.
- **Finite-element nodes.** `evl/cdNodes.txt` is written once per run: the identifier and the position of each finite-element node. The rows of the cluster-dynamics block follow these nodes, which are not those of the mesh file since the elements are quadratic.
- **Starting from another study.** The cluster-dynamics block is written with all its digits. Copy `evl_<n>.txt` of a study to `evl/evl_0.txt` of a new folder, after `microstructureGenerator`, and the new run starts from its fields.
- **Cluster-dynamics fields.** The last block of `evl_<runID>.txt` has one row per finite-element node. Its columns are the mobile concentrations, then the number densities $N_c, N_{a1}, N_{a2}, N_{a3}$ in units of $b^{-3}$, then the defect contents $c_c, c_{a1}, c_{a2}, c_{a3}$. Cluster dynamics is written in text format only (`outputBinary=0`).
- **Columns of `F_0.txt`.** The crystal writes `runID`, `time`, `dt`, the nine components of `betaP`, its trace and norm, and the same for its rate. `DislocationDynamics` adds the glissile, sessile, boundary and grain-boundary densities in m$^{-2}$; the densities and the plastic distortion per slip system, and the elastic and core energies, when the corresponding switches of `DD.txt` are set. `ElasticDeformation` adds the strain and stress components (with the uniform controller). `ClusterDynamics` adds no column: its strain is in `betaP`.
- **Format and frequency.** `outputBinary` and `outputFrequency` are in `DefectiveCrystal.txt`. `tools/DDconverter` converts the binary files in `./evl` to text; it has its own `Makefile` and is not part of the CMake build.
- **Restart.** `startAtTimeStep=-1` restarts from the last step recorded in `F/F_0.txt`.
- **Reading in Python.** `python/readEVL.py` reads a configuration file `evl_<n>.txt` (network nodes, loops, links, loop nodes, inclusions, ellipsoidal inclusions, displacement and cluster-dynamics blocks) with a header of ten or eleven counts, and `evl/cdNodes.txt`; `python/readAUX.py` reads `ddAux_<n>.txt`; `python/readF.py` reads the `F` files (`readFfile`, `getFarray(F, labels, label)`); `python/readFile.py` reads values from the input files; `python/getLoops.py` gives the center and mean radius of each loop. `python/modlibUtils.py` has simpler versions of `readEVLtxt`, `readAUXtxt` and `readFfile`.
- **Viewing.** `DDqt` opens a simulation folder and displays its configurations.

## 9. Migrating input files from 1.0.0

Input files of 1.0.0 do not run unchanged. A missing variable stops the run with the message `File ... does not cointain line with format <name>=...;`, which names the variable to add. The simplest route is to start from the templates of `Library/` and carry the values over.

**Variables that left `DD.txt`:**

| In 1.0.0 (`DD.txt`) | In 2.0.0 |
|---|---|
| `useDislocations`, `useClusterDynamics`, `useInclusions`, `useElasticDeformation` | the `physics` list of `DefectiveCrystal.txt` |
| `useCracks` | removed |
| `useFEM`, `Nsteps`, `dtMax`, `startAtTimeStep`, `outputFrequency`, `outputBinary` | `DefectiveCrystal.txt` |
| `periodicImageSize`, `EwaldLengthFactor` | `DefectiveCrystal.txt` |
| `useElasticDeformationFEM`, `inertiaReliefPenaltyFactor` | `ElasticDeformation.txt` |
| `useClusterDynamicsFEM` | `ClusterDynamics.txt` |
| `use_stochasticForce` | renamed `useStochasticForce` |
| `useSubCycling` | removed: subcycling is on when `subcyclingBins` has more than one value |
| `stepsBetweenBVPupdates`, `use_directSolver_FEM`, `solverTolerance` | removed |

**New variables:** `physics` and `maxResolveSteps` in `DefectiveCrystal.txt`; `bulkNucleationModel`, `surfaceNucleationModel`, `outputDislocationDensityPerSlipSystem` and `nodalVelocityConstraints` in `DD.txt`.

**Changed defaults in the templates:** `Nsteps` 250 → 10; `dtMax` 1e25 → 1e35; `periodicImageSize` `1 1 1` → `0 0 0`; `useElasticDeformationFEM` 1 → 0; `useClusterDynamicsFEM` 1 → 0; `subcyclingBins` `1 2 5 10 50 100` → `1`.

**Material file:**

| In 1.0.0 | In 2.0.0 |
|---|---|
| `enabledSlipSystems=full<111>{110}` and `full<100>{110}` (BCC) | `a/2<111>{110}`, `a<100>{110}` |
| `full`, `Shockley`, `Kear` (FCC) | `a/2<110>{111}`, `a/6<112>{111}`, `a/3<112>{111}` |
| `fullBasal`, `ShockleyBasal`, `fullPrismatic`, `fullPyramidal` (HEX) | `<a>{basal}`, `<Shockley>{basal}`, `<a>{prismatic}`, `<a>{pyramidal}` |
| `dislocationMobilityType...` and the variables `B0e_SI`, `B1e_SI`, `B0s_SI`, `B1s_SI`, `Bk_SI`, `dH0_eV`, `p`, `q`, `Tf`, `tauC_SI`, `a0`–`a3` | one `mobility_<slip system>` line per slip system ([section 3.3](#33-mobility-laws)); the same numbers become the arguments of `viscousDrag` (edge) and `kinkPair` (screw), with the kink height, the kink width and the conjugate angle added |
| `ISF_SI`, `USF_SI`, `MSF_SI` (HEX) | `basalISF_SI`, `basalUSF_SI`, `basalMSF_SI`; a prismatic system also needs `prismaticWaveVectors` and `prismaticGammaSurfacePoints` |
| `dOmegav`, `Ufv_eV`, `Umv_eV`, `D0v_SI` | the cluster-dynamics block (`mobileSpeciesVector`, ...) |
| `correlationFile_L`, `correlationFile_T` (MD solid-solution noise file) | `correlationFile_xz`, `correlationFile_yz` |

An old slip-system name gives no error: the crystal is simply created without slip systems.

**Microstructure files:** in the density specifications of shear and Frank loops, `targetDensity`, `radiusDistributionMean` and `radiusDistributionStd` are now `targetDensity_SI`, `radiusDistributionMean_SI` and `radiusDistributionStd_SI`. Frank loops gained `allowedGrainIDs` and `allowedPlaneIDs`.

**Cluster dynamics.** With immobile species, the material file needs the cluster-dynamics variables of `Zr_CD4opt.txt`, which is the reference for the 4 + 8 species model ([section 1.5](#15-cluster-dynamics-since-200) lists the keys that changed).

**Source code.** Code that includes headers of `DislocationMicrostructure/` must use `MicrostructureGeneration/`. `DislocationNetwork<dim,corder>` is now `DislocationNetwork<dim>`.

## 10. Tests and validation

The repository has no automated test suite. What it has:

- **`testsPy/`**, scripts that are run by hand from their folder and produce plots:

  | Folder | Checks | Needs |
  |---|---|---|
  | `dislocationMobility` | velocity against character angle and stress, for every material of `Library/Materials` | `pyMoDELib` |
  | `gammaSurface` | gamma surfaces and their cuts, for the planes and second phases of W | `pyMoDELib` |
  | `interpolation1D` | linear and Lagrange interpolants with periodic extrapolation | `pyMoDELib` |
  | `stressStraight` | timing of the stress of a straight segment | `pyMoDELib` |
  | `analyticalSolidSolutionCorrelations`, `mdSolidSolutionCorrelations`, `mdSolidSolutionNoiseSample`, `mdStackingFaultCorrelations`, `mdStackingFaultNoiseSample` | sampled noise against its input correlation | `pyMoDELib` with FFTW, SciPy |
  | `periodicFields`, `periodicEnergy` | fields and energy of periodic dipoles against the number of images | `pyMoDELib`; the executables |
  | `clusterDynamics` | the rate equations of the immobile clusters against an independent implementation (`referenceRates.py`), a restart, and stress-oriented nucleation; prints PASS or FAIL for each check | the executables built with SUNDIALS, NumPy, SciPy |

  The plots use Matplotlib, several of them with LaTeX labels.
- **The tutorials:** five complete simulations ([section 7](#7-tutorials)).
- **Continuous integration:** `.github/workflows/workflow.yml` builds the Doxygen pages on each push to `master` and publishes them; `.github/workflows/doxygen-wiki.yml` converts them to the pages of the Wiki. `.github/workflows/container.yml` builds the container image when a release is published; that build compiles the library and the two command-line tools, and runs two steps of the `dipoleNoise` tutorial.

### Validation of the immobile-species integration

The three cases of `testsPy/clusterDynamics` use `tutorials/spatialCDcvode` (1 µm box, 573 K, `Zr_CD4opt.txt`, steps of 0.1 dpa). `referenceRates.py` is an implementation of the equations of the report written separately from the C++ code and integrated with SciPy (BDF, relative tolerance $10^{-10}$).

| Check | Case | Expected | Result |
|---|---|---|---|
| Fields after a step, against the reference equations | `cvode`, steps 0 and 2, 61 nodes spread over the mesh, the 8 immobile fields | equal | largest relative difference $3.7\times10^{-5}$ and $1.2\times10^{-5}$ |
| Saturation of the vacancy-loop density | `cvode`, interior nodes after 0.3 dpa | at most $G\varepsilon_{vL}\tau/n^{nuc} = 1.7916\times10^{-8}$ per atom, the value without coalescence | $1.7856\times10^{-8}$ |
| The three prismatic families without stress | `cvode` | equal | relative difference $4\times10^{-14}$ |
| No negative value | `cvode` | | smallest value of the fields $10^{-70}$ |
| A new study started from a configuration file | `cvodeRestart`: two steps from `evl_2.txt` of `cvode` | the fields of the uninterrupted run | largest difference $1.2\times10^{-10}$ of the largest value of each field |
| Nucleation oriented by stress | `cvodeStress`: 2 MPa along $a_1$ | fractions 0.4915, 0.2543, 0.2543 of the interstitial families; fields equal to the reference | largest relative difference $3.7\times10^{-5}$ |

These checks establish that the code integrates its rate equations. They do not compare the model with measurements. The numbers above were obtained before the model changes of [section 1.5](#15-cluster-dynamics-since-200) (capture length, coalescence, gate, atomic volume); `referenceRates.py` follows the current equations, and `runCases.sh` should be run again.

**What these tests do not cover.** The conversion to discrete loops has no test case since the removal of `spatialCDtest`. All cases are single crystals of HEX zirconium.

## 11. Status and known issues

| Component | Status |
|---|---|
| Dislocation glide, junctions, remeshing | in use; exercised by the tutorials |
| Periodic and finite domains, inclusions | in use |
| Glide-plane noise | in use; `dipoleNoise` tutorial. The templates `MDSolidSolution.txt` and `MDStackingFault.txt` are placeholders: set `type`, the correlation files, and three values for `gridSize` and `gridSpacing_SI`, as the `testsPy` scripts do |
| Cluster dynamics, mobile and immobile species | in use; `spatialCDcvode` (CVODE, kinetics of the D1/M1 report). `Zr_CD4opt.txt` holds the parameters of the report, valid at 573 K and without anisotropic diffusion. The network density is constant. The emission model of the vacancy loops has a constant binding energy. Grain boundaries are not sinks for the loops. The diffusion-anisotropy bias is that of a HEX crystal |
| Climb | Galerkin climb solver coupled to the cluster-dynamics concentrations; `annealing` tutorial |
| Species counts | fixed at compile time (`MODELIB_CD_MSIZE`, `MODELIB_CD_ISIZE`). Two builds are needed to run all the tutorials. Making the counts run-time parameters is the open item |
| Material files | all are in the 2.0.0 format. `Zr_CD4opt.txt` has the variables of the immobile species; `Zr4_Fitted.txt` is the parameter set of the removed explicit integrator, kept for reference, and does not run. `Zr_CD2.txt`, `Zr_CD3_BMD19.txt` and `Zr_CD4.txt` have two, three and four mobile species and need a build with `MODELIB_CD_ISIZE=0` and the matching `MODELIB_CD_MSIZE`. The mobility numbers of `UO2.txt` are those of W and are placeholders |
| Conversion of loop fields to discrete loops | in use ([section 3.6](#36-cluster-dynamics)). One conversion per run, of all families at once. Every cluster of a converted family becomes a loop, also in the bi-pyramid range unless `minimumLoopSize` excludes it. The relaxation-volume strain of a converted family (`immobileSpeciesRelRelaxVol`) is not carried by the discrete loops. No tutorial or test case exercises it |
| Nucleation | `bulkNucleationModel=1` is a first model: it inserts a shear loop of radius $50\,b$ in every mesh element where a resolved shear stress is below $0.05\,\mu$. `surfaceNucleationModel` is read and has no effect |
| Cross slip | model 1 (deterministic) is in use. Model 2 (thermally activated, HEX) has not been updated: it relies on the names of the 1.0.0 mobility classes and on material variables that no material file has |
| Mobility type `none` | not usable: an enabled slip system needs a law |
| Kear partials (FCC) | `a/3<112>{111}` reads the mobility line of `a/6<112>{111}` |
| Polycrystals | `PolyCrystalFile.write` in `modlibUtils.py` writes the rotation of grain 1 only; add `C2G2`, ... by hand for a multi-grain mesh. Grain-boundary transmission is not connected to the solver |
| Cracks | the `CircularCrack` and `CrackMesh` microstructure types are accepted and generate nothing; there is no crack physics |
| Python examples | `python/modelibPy11.py` still refers to the removed `crossSlip` tutorial and to the old density names. `python/DDMobilityPy.py` and `python/MLglideSolver.py` do not define the classes (`MobilitySolver`, `PIGNN`) that the C++ code calls |
| `testsPy` | the scripts that build a simulation folder (`periodicFields`, `periodicEnergy` and the five noise tests) do not yet copy `DefectiveCrystal.txt`, which 2.0.0 requires |
| `DDqt` | the immobile cluster-dynamics fields are not plotted |
| Sources outside the build | `PeriodicGlidePlane_NEW.h/.cpp`, `GalerkinGlideSolver_lumped.cpp`, `Lattices_NEW/`, `Math/Rational_NEW.h` and `MPI/` are in the repository and are not compiled |
| Doxygen pages | the main page of the GitHub Pages site still describes the earlier header-only MODEL library. The Wiki has a home page for 2.2.0 (`doxygen/dox/Home.dox`), which replaces the Wiki's Home at each push to `master` |
| Container image | built without SUNDIALS: cluster dynamics with immobile species does not run in it |
| Manual | `manual_beta/`: the theory chapters are written; the section on running MoDELib is referenced but not yet written. It does not describe the 2.0.0 input files |

## 12. Documentation

- **API reference (Doxygen):** <https://mechanicsofdefects.github.io/MoDELib/>, built from `doxygen/` on each push to `master`.
- **Theory manual:** `manual_beta/MoDELib_manual.pdf`. Chapters: the elastic theory of discrete dislocations; discrete dislocation dynamics in MoDELib (network topology, the uniform load controller); a review of tensor calculus; the elastic fields of piecewise-straight loops.
- **Input reference:** the templates in `Library/` carry a comment on each variable. `Library/DefectiveCrystal/DefectiveCrystal.txt` and `Library/DislocationDynamics/DD.txt` list the controls; `Library/Materials/Zr_CD4opt.txt` lists the cluster-dynamics variables.
- **Wiki:** <https://github.com/MechanicsOfDefects/MoDELib/wiki>. `.github/workflows/doxygen-wiki.yml` converts the Doxygen output to Markdown (doxybook2) and pushes it to the Wiki of the repository at each push to `master`: the page `Home` (`doxygen/dox/Home.dox`), the other pages of `doxygen/dox`, one page per class and per file with its diagram, and a sidebar. It needs the Wiki enabled with a first page, and read and write permissions for workflows (Settings → Actions → General). To render the pages locally, run `doxygen Doxyfile` in `doxygen/` and open `doxygen/HTML/index.html`.
- **Report of the cluster-dynamics model:** `docs/reports/GW_Phase4_D1M1.pdf`.
- **Git notes:** `git_basic_commands.txt`.

## 13. How to cite

If you use MoDELib, please cite the method paper:

> G. Po, M. S. Mohamed, T. Crosby, C. Erel, A. El-Azab, N. Ghoniem, *Recent progress in discrete dislocation dynamics and its applications to micro plasticity*, JOM 66 (2014) 2108–2120.

```bibtex
@article{Po2014JOM,
  author  = {Po, Giacomo and Mohamed, Mamdouh S. and Crosby, Tamer and Erel, Can and El-Azab, Anter and Ghoniem, Nasr},
  title   = {Recent progress in discrete dislocation dynamics and its applications to micro plasticity},
  journal = {JOM},
  volume  = {66},
  pages   = {2108--2120},
  year    = {2014}
}
```

Cite also the papers behind the features you use:

- **Kink-pair mobility law:** G. Po, Y. Cui, D. Rivera, D. Cereceda, T. D. Swinburne, J. Marian, N. Ghoniem, *A phenomenological dislocation mobility law for bcc metals*, Acta Mater. 119 (2016) 123–135, [doi:10.1016/j.actamat.2016.08.016](https://doi.org/10.1016/j.actamat.2016.08.016).
- **Non-singular fields, option 1:** W. Cai, A. Arsenlis, C. R. Weinberger, V. V. Bulatov, *A non-singular continuum theory of dislocations*, J. Mech. Phys. Solids 54 (2006) 561–587.
- **Non-singular fields, option 2:** G. Po, M. Lazar, D. Seif, N. Ghoniem, *Singularity-free dislocation dynamics with strain gradient elasticity*, J. Mech. Phys. Solids 68 (2014) 161–178, [doi:10.1016/j.jmps.2014.03.005](https://doi.org/10.1016/j.jmps.2014.03.005).

## 14. Repository conventions

- **Templates are copied, not edited.** A simulation copies its input files from `Library/` into its own `inputFiles/` and changes the copies. The templates stay as the reference.
- **One folder per simulation.** `inputFiles/`, `evl/` and `F/` sit together, and every tool takes that folder as its argument.
- **Input syntax.** `name=value;` with `#` starting a comment. Vectors and matrices are space-separated values, with one matrix row per line, ended by one `;`.
- **SI at the boundary only.** Input variables in SI units end in `_SI` (or name their unit, as in `dH0_eV`). Everything inside the code is in units of $b$, $\mu$ and $c_s$.
- **Headers and sources mirror each other.** A class declared in `include/<Module>/` is implemented in `src/<Module>/`.
- **Every source file states the license** in its first lines.
- **Options are read where they are used.** A class reads its own parameters from the input files with `TextFileParser`. Every variable that is read is required.

## 15. License and contact

MoDELib is distributed without any warranty under the GNU General Public License v2, as stated in each source file.

MoDELib was created by Giacomo Po, whose copyright notices in the tools date from 2012. The contributors are listed in `Contributors.txt` (38 names).

Questions and bug reports: please open an issue at <https://github.com/MechanicsOfDefects/MoDELib/issues>.
