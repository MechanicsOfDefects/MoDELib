# Mechanics of Defects Evolution Library (MoDELib)

**A C++ library for the discrete mechanics of crystal defects.**

MoDELib evolves defects that are represented individually inside a crystal, and computes their interactions from continuum elasticity. One object, the `DefectiveCrystal`, holds every kind of defect that a simulation enables:

- a network of discrete dislocation loops;
- Eshelby inclusions, spherical or polyhedral;
- a uniform elastic deformation applied by a load controller, with an optional finite-element correction;
- mobile point-defect species diffusing and reacting on the finite-element mesh (cluster dynamics);
- cracks (interface only; see [section 9](#9-status)).

The library separates three things that a simulation needs:

- **the crystal:** its structure, slip systems, mobility law and domain mesh;
- **the microstructure:** which defects are present, and where;
- **the algorithm:** how forces become velocities, and how the configuration is updated.

The governing rule follows from that separation:

> Each defect type is a *microstructure* with the same four duties: solve, propose a time step, update its configuration, and report its fields. The crystal adds the fields and takes the smallest time step.

**Languages:** C++20 for the library and its tools. Python for input generation, post-processing and the optional `pyMoDELib` bindings (pybind11).\
**Version:** 1.0.0.\
**License:** GNU GPL v2, as stated in the header of every source file. The repository has no separate `LICENSE` file.\
**Primary reference:** G. Po, M. S. Mohamed, T. Crosby, C. Erel, A. El-Azab, N. Ghoniem, *Recent progress in discrete dislocation dynamics and its applications to micro plasticity*, JOM 66 (2014) 2108–2120.

## Contents

1. [Capabilities](#1-capabilities)
2. [Methodology](#2-methodology)
3. [Repository structure](#3-repository-structure)
4. [Installation](#4-installation)
5. [Quick start](#5-quick-start)
6. [Tutorials](#6-tutorials)
7. [Outputs](#7-outputs)
8. [Tests and validation](#8-tests-and-validation)
9. [Status](#9-status)
10. [Documentation](#10-documentation)
11. [How to cite](#11-how-to-cite)
12. [Repository conventions](#12-repository-conventions)
13. [License and contact](#13-license-and-contact)

---

## 1. Capabilities

| Capability | Where | Notes |
|---|---|---|
| Dislocation network topology | `include/LoopNetwork/`, `include/DislocationDynamics/` | loops, loop nodes and loop links; network nodes and segments shared between loops |
| Elastic fields of dislocations | `DislocationDynamicsBase/StraightDislocationSegment.h`, `StressStraight.h` | piecewise-straight segments; `_MODEL_NON_SINGULAR_DD_` selects classical (0), Cai (1, default) or Lazar (2) fields; `coreSize` in `DD.txt` |
| Glide solver | `GalerkinGlideSolver.h` | nodal velocities by Galerkin projection of the line velocity; lumped matrix; periodic constraints by a null-space projection |
| Glide solver in Python | `PyGlideSolver.h`, `python/MLglideSolver.py` | `glideSolverType=pybind11`; the module path is `pyModuleName` in `DD.txt` |
| Climb solver | `GalerkinClimbSolver.h` | `climbSolverType=Galerkin`; climb steps are taken when the plastic distortion rate falls below `glideEquilibriumRate` |
| Mobility laws | `include/DislocationMobilities/` | FCC, BCC, HEX basal, HEX prismatic, HEX pyramidal; a Python law (`DislocationMobilityPy.h`); grain-boundary mobility |
| Crystal structures and slip systems | `include/PolycrystallineMaterials/` | FCC (`full`, `Shockley`, `Kear`), BCC (`full<111>{110}`, `full<100>{110}`), HEX (`fullBasal`, `ShockleyBasal`, `fullPrismatic`, `fullPyramidal`) |
| Polycrystals | `Grain.h`, `GrainBoundary.h`, `Polycrystal.h` | grains are mesh regions; `grainBoundaryTransmissionModel` in `DD.txt` |
| Cross slip | `CrossSlipModels.h` | `crossSlipModel`: 0 none, 1 deterministic (largest glide Peach–Koehler force), 2 thermally activated (Escaig; HEX only) |
| Junctions and remeshing | `DislocationJunctionFormation.h`, `DislocationNetworkRemesh.h`, `DislocationNodeContraction.h` | segment lengths kept between `Lmin` and `Lmax` |
| Periodic domains | `GlidePlanes/PeriodicGlidePlane.h` | periodic faces chosen per mesh face (`periodicFaceIDs`); image sums set by `periodicImageSize`; `EwaldLengthFactor` |
| Finite domains | `include/FEM/`, `DislocationDynamicsBase/ElasticDeformationFEM.h` | `useFEM=1`; direct solvers (CHOLMOD, UMFPACK) or an iterative solver |
| Uniform load control | `DislocationDynamicsBase/UniformController.h` | stress control, strain control, or mixed, per Voigt component (section 2.4) |
| Glide-plane noise | `GlidePlanes/AnalyticalSolidSolutionNoise.h`, `MDSolidSolutionNoise.h` | solid-solution stress noise sampled on a grid (FFTW); analytical correlations, or correlations from MD (Mo–Nb–Ti files in `Library/GlidePlaneNoise/`) |
| Stochastic force | `DD.txt` | `use_stochasticForce`, `stochasticForceSeed` |
| Inclusions | `SphericalInclusion.h`, `PolyhedronInclusion.h` | eigendistortion, velocity reduction factor and second phase per inclusion |
| Cluster dynamics | `include/ClusterDynamics/` | mobile species on the finite-element mesh: anisotropic diffusion, second-order reactions, production from a dose rate, sinks at discrete dislocations with an elastic bias |
| Microstructure generation | `include/DislocationMicrostructure/` | shear, prismatic and Frank loops; stacking-fault tetrahedra; periodic dipoles; a planar loop from its nodes; inclusions; dislocation lines from an OVITO DXA `.vtk` file. Each by target density or individually |
| Time stepping | `DDtimeStepper.h` | `fixed` or `adaptive`; velocity filter; optional subcycling |
| Lattice arithmetic | `include/Lattices/` | lattice and reciprocal vectors, rational directions, LLL reduction, CSL and DSCL of bicrystals |
| Visualization | `tools/DDqt/` | Qt6 and VTK viewer of configurations, glide planes, inclusions, quadrature data and the mesh |

## 2. Methodology

### 2.1 The defective crystal

`DefectiveCrystal` is a container of microstructures. `DD.txt` decides which ones exist in a run:

| Switch in `DD.txt` | Microstructure added |
|---|---|
| `useInclusions` | inclusion microstructure |
| `useElasticDeformation` | uniform elastic deformation and load control |
| `useClusterDynamics` | mobile point-defect species |
| `useDislocations` | dislocation network |

`useFEM` enables the finite-element mesh for the microstructures that can use it.

One step of `DefectiveCrystal::runSingleStep` does the following, in this order:

1. **Solve.** Each microstructure solves for its own rates.
2. **Time step.** Each microstructure proposes a time step; the smallest one is used, capped by `dtMax`.
3. **Output.** Every `outputFrequency` steps, the configuration and the averaged quantities are written.
4. **Update.** Each microstructure advances its configuration by the time step.

For the dislocation network, cross slip is carried out at the end of the solve, before the nodes move. The update then moves the nodes and executes the remaining discrete events: remeshing, junction formation, and remeshing again.

### 2.2 Dislocation glide

Each segment carries quadrature points. At a point with unit tangent $\boldsymbol{\xi}$ and Burgers vector $\mathbf{b}$, the Peach–Koehler force per unit length is

$$\mathbf{f}_{PK} = (\boldsymbol{\sigma}\,\mathbf{b}) \times \boldsymbol{\xi},$$

where $\boldsymbol{\sigma}$ is the sum of the stresses of all enabled microstructures. The total force adds the stacking-fault force and a line-tension force scaled by `alphaLineTension`:

$$\mathbf{f} = \mathbf{f}_{PK} + \mathbf{f}_{SF} + \mathbf{f}_{LT}.$$

The mobility law of the slip system turns the force and the local stress into a glide velocity $\mathbf{v}$ at that point. Nodal velocities $\mathbf{V}$ then follow from a Galerkin projection with the segment shape functions $\mathbf{N}$:

$$\left[\int \mathbf{N}^\top \mathbf{N}\, dL\right] \mathbf{V} = \int \mathbf{N}^\top \mathbf{v}\, dL .$$

With adaptive time stepping, $dt = \min(\texttt{dxMax}/v_{max},\ \texttt{dtMax})$.

### 2.3 Elastic fields

The stress of the network is the sum over straight segments, in the non-singular form selected at compile time (section 1) with core width `coreSize`. In a finite domain, the infinite-medium fields are corrected by a finite-element solution on the mesh (the superposition principle). In a periodic domain, the fields of `periodicImageSize` images are summed along each periodic shift vector.

The derivations are in `manual_beta/MoDELib_manual.pdf`: eigendistortion theory, discrete loops in anisotropic and isotropic media, the solid angle as a line integral, and the fields of straight segments.

### 2.4 Load control

`ElasticDeformation.txt` gives an applied stress $\boldsymbol{\sigma}_0$, an applied strain $\boldsymbol{\varepsilon}_0$, their rates, and a machine stiffness ratio $\alpha_i$ per Voigt component (order 11, 22, 33, 12, 23, 13). With $D_i = \alpha_i C_{ii}$, material stiffness $\mathbf{C}$ and plastic strain $\boldsymbol{\varepsilon}_p$, the stress in the volume is

$$\boldsymbol{\sigma} = \mathbf{C}\,(\mathbf{C}+\mathbf{D})^{-1}\left[\boldsymbol{\sigma}_0 + \mathbf{D}\,(\boldsymbol{\varepsilon}_0 - \boldsymbol{\varepsilon}_p)\right].$$

| `stiffnessRatio` component | Result |
|---|---|
| $\alpha_i = 0$ | stress control: $\sigma_i = \sigma_{0i}$ |
| $\alpha_i \to \infty$ (for example `1e20`) | strain control: $\varepsilon_i = \varepsilon_{0i}$ |

### 2.5 Units

Input files give material constants in SI units, in variables whose names end in `_SI`. Internally, and in the output files, quantities are dimensionless:

| Quantity | Unit |
|---|---|
| length | Burgers vector magnitude $b$ |
| stress | shear modulus $\mu = \mu_0 + \mu_1 T$ |
| speed | shear-wave speed $c_s = \sqrt{\mu/\rho}$ |
| time | $b/c_s$ |
| energy | $\mu b^3$ |

The values of these units are printed when a simulation starts.

### 2.6 The four input files

A simulation folder holds `inputFiles/`, `evl/` and `F/`. Four files in `inputFiles/` define the run:

| File | Declares |
|---|---|
| `DD.txt` | which microstructures are enabled, and every numerical control |
| `polycrystal.txt` | the material file, temperature, mesh file, crystal rotation `C2G1`, mesh mapping `F` and `X0`, and periodic faces |
| the material file, for example `W.txt` | crystal structure, enabled slip systems, elastic constants, mobility parameters, stacking-fault energies, noise, second phases, cluster-dynamics parameters |
| `initialMicrostructure.txt` | a list of microstructure files, each declaring one family of initial defects |

Mesh nodes $\mathbf{X}$ are mapped to $\mathbf{x} = \mathbf{F}(\mathbf{X} - \mathbf{X}_0)$, so one unit mesh serves domains of any size and shape.

## 3. Repository structure

```
MoDELib/
├── include/                      # headers (377 files)
│   ├── DislocationDynamics/      # DefectiveCrystal, network, solvers, cross slip, junctions, remeshing
│   ├── DislocationDynamicsBase/  # microstructure interfaces, segment fields, inclusions, load controller
│   ├── DislocationDynamicsIO/    # configuration (evl) and auxiliary (ddAux) files, parameters
│   ├── DislocationMicrostructure/# generators and specifications of initial microstructures
│   ├── DislocationMobilities/    # mobility laws
│   ├── PolycrystallineMaterials/ # lattices, slip systems, grains, grain boundaries, gamma surfaces
│   ├── GlidePlanes/              # glide planes, periodic glide planes, noise
│   ├── ClusterDynamics/          # mobile species on the FEM mesh, reactions
│   ├── ElasticDeformation/  InclusionMicrostructure/  DiscreteCrackMechanics/
│   ├── LoopNetwork/              # the generic loop-network topology
│   ├── FEM/  Mesh/  Quadrature/  # finite elements, simplicial meshes (Gmsh), Gauss–Legendre rules
│   ├── Lattices/  Lattices_NEW/  # lattice arithmetic
│   ├── Geometry/  Math/  IO/  Utilities/
│   ├── ParticleInteraction/  MPI/  DDvtk/
├── src/                          # implementations (160 files), one folder per header folder
├── tools/
│   ├── DDomp/                    # the simulation driver
│   ├── MicrostructureGenerator/  # writes the initial configuration
│   ├── DDqt/                     # Qt6/VTK viewer
│   ├── pyMoDELib/                # pybind11 module
│   ├── DDconverter/              # binary configuration files to text
│   ├── DDsegments/               # reads segments from an OVITO DXA .vtk file
│   └── QuadratureGeneration/     # generates the Gauss–Legendre tables
├── python/                       # modlibUtils.py (input and output helpers), examples of the bindings,
│                                 # MLglideSolver.py, DDMobilityPy.py
├── Library/                      # templates copied into each simulation
│   ├── DislocationDynamics/      # DD.txt
│   ├── ElasticDeformation/       # ElasticDeformation.txt
│   ├── Materials/                # Al, AlMg5, AlMg10, AlMg15, Co, Cu, Ni, W, Zr, Zr_CD2, Zr_CD4
│   ├── Meshes/                   # Gmsh meshes: unit cubes, cylinders, bicrystal, polycrystals
│   ├── Microstructures/          # one template per defect family
│   └── GlidePlaneNoise/          # noise templates and MD correlation data
├── tutorials/                    # crossSlip, dipoleNoise, polycrystal
├── tests/periodicFields/         # fields of periodic dipoles against the number of images
├── manual_beta/                  # the theory manual (tex, pdf)
├── doxygen/                      # Doxyfile and pages of the API documentation
├── scripts/                      # shell scripts for cleaning and for image and video conversion
├── .github/workflows/            # workflow.yml: Doxygen pages; container.yml: container image
├── CMakeLists.txt  Dockerfile  Contributors.txt  README.md  git_basic_commands.txt
```

## 4. Installation

### Dependencies

| Dependency | Needed for | Required |
|---|---|---|
| C++20 compiler with OpenMP | everything | yes |
| CMake 3.20 or later | the build (`cmake_path` is used) | yes |
| Eigen 3 | linear algebra | yes |
| FFTW 3 | glide-plane noise | yes: configuration fails without it |
| Boost | headers | no (a warning if absent) |
| SuiteSparse (CHOLMOD, UMFPACK) | direct FEM solvers | no (a warning if absent) |
| Python 3 with NumPy | the tutorials' `generateInputFiles.py` | for the tutorials |
| pybind11 | `pyMoDELib`, the Python glide solver and mobility | when `USE_PYBIND11=ON` (the default) |
| Qt 6 and VTK | `DDqt` | when `BUILD_TOOLS=ON` (the default) |

### Build

```bash
git clone https://github.com/MechanicsOfDefects/MoDELib.git
cd MoDELib
```

Before configuring, set the Eigen path in `CMakeLists.txt` to the one on your machine. The path is assigned inside the file, so `-DEIGEN3_INCLUDE_DIRS=...` on the command line has no effect:

```cmake
set(EIGEN3_INCLUDE_DIRS /opt/local/include/eigen3)   # MacPorts; on Debian/Ubuntu: /usr/include/eigen3
```

Then:

```bash
cmake -S . -B build
cmake --build build -j
```

The build type is `Release`, with `-Ofast -march=native -fopenmp`. With the Eigen path set and the `DDqt` line removed (see below), the library, `microstructureGenerator` and `DDomp` build with g++ 13.3 and CMake 3.28 on Ubuntu 24.04.

| Option | Default | Effect |
|---|---|---|
| `USE_PYBIND11` | `ON` | builds `pyMoDELib` and the Python solver hooks; needs Python and pybind11 |
| `BUILD_TOOLS` | `ON` | builds `microstructureGenerator`, `DDomp`, `DDqt` and (with pybind11) `pyMoDELib`. `DDconverter`, `DDsegments` and `QuadratureGeneration` are built separately |

Two points to know:

- **Tools without Qt and VTK.** `tools/CMakeLists.txt` adds `DDqt` whenever tools are built, and `DDqt` requires Qt 6 and VTK. To build the command-line tools on a machine without them, comment out the `DDqt` line in that file.
- **Library only.** `cmake -S . -B build -DBUILD_TOOLS=OFF -DUSE_PYBIND11=OFF` builds `libMoDELib` alone.

On Debian or Ubuntu the packages for the library and the command-line tools are `libeigen3-dev libfftw3-dev libboost-dev libsuitesparse-dev`. pybind11, Qt 6 and VTK (built with Qt support) are needed only for the optional parts.

### Container image

Each release is also published as a container image with the command-line tools, `Library/`, `python/` and the tutorials:

```bash
docker pull ghcr.io/mechanicsofdefects/modelib:v1.0.0
docker run -it --rm ghcr.io/mechanicsofdefects/modelib:v1.0.0
```

The container starts in `/opt/MoDELib/tutorials`, with `DDomp` and `microstructureGenerator` on the `PATH`. `DDqt` and `pyMoDELib` are not in the image. The image is built from the `Dockerfile` at the root of the repository, which applies the two changes above and replaces `-march=native` by `-march=x86-64-v2`, so that the binaries do not depend on the build machine. To build it locally: `docker build -t modelib .`

## 5. Quick start

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

# numerical controls
shutil.copy2('../../Library/DislocationDynamics/DD.txt', 'inputFiles/DD.txt')
setInputVariable('inputFiles/DD.txt', 'useFEM', '0')
setInputVariable('inputFiles/DD.txt', 'glideSolverType', 'Galerkin')
setInputVariable('inputFiles/DD.txt', 'Nsteps', '1000')

# glide-plane noise
shutil.copy2('../../Library/GlidePlaneNoise/AnalyticalSolidSolutionNoise.txt', 'inputFiles/AnalyticalSolidSolutionNoise.txt')
setInputVariable('inputFiles/AnalyticalSolidSolutionNoise.txt', 'MSSS_SI', '0.45e18')

# material
shutil.copy2('../../Library/Materials/AlMg15.txt', 'inputFiles/AlMg15.txt')
setInputVariable('inputFiles/AlMg15.txt', 'enabledSlipSystems', 'Shockley')
setInputVariable('inputFiles/AlMg15.txt', 'glidePlaneNoise', 'AnalyticalSolidSolutionNoise.txt')
b_SI = getValueInFile('inputFiles/AlMg15.txt', 'b_SI')

# domain: a 50 nm box mapped from the unit-cube mesh
shutil.copy2('../../Library/Meshes/unitCube24.msh', 'inputFiles/unitCube24.msh')
pf = PolyCrystalFile('AlMg15.txt')
pf.absoluteTemperature = 300
pf.meshFile = 'unitCube24.msh'
pf.grain1globalX1 = np.array([0, 1, 1])                        # crystal direction along global x1
pf.grain1globalX3 = np.array([-1, 1, -1])                      # crystal direction along global x3
pf.boxEdges = np.array([[0, 1, 1], [2, 1, -1], [-1, 1, -1]])   # row i is the direction of box edge i
pf.boxScaling = np.array([50e-9, 50e-9, 50e-9]) / b_SI         # edge lengths in units of b
pf.periodicFaceIDs = np.array([-1])                            # -1: every face periodic
pf.write('inputFiles')
```

The script goes on to set the applied stress in `ElasticDeformation.txt`, and to declare two periodic dipoles in `periodicDipoleIndividual.txt`, which `initialMicrostructure.txt` then lists.

The same objects are available from Python through `pyMoDELib`:

```python
import sys, os
import numpy as np
sys.path.append("../build/tools/pyMoDELib")
import pyMoDELib

ddBase = pyMoDELib.DislocationDynamicsBase(os.path.abspath("../tutorials/crossSlip"))

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

`python/modelibPy11.py` and `python/defectFields.py` give the complete examples.

## 6. Tutorials

Each tutorial is a folder with one `generateInputFiles.py`.

| Tutorial | Material | Domain | Initial microstructure | Exercises |
|---|---|---|---|---|
| `crossSlip` | W, BCC, `full<111>{110}` | 1 µm box, periodic | prismatic loops, density $10^{13}$ m$^{-2}$, radius 100 nm | deterministic cross slip (`crossSlipModel=1`) |
| `dipoleNoise` | AlMg15, FCC, Shockley partials | 50 nm box, periodic | two periodic dipoles | analytical solid-solution noise; elastic energy per length |
| `polycrystal` | W, BCC, `full<111>{110}` | 1 µm cube of 10 grains, not periodic | shear loops in one grain and one slip system | a polycrystalline mesh (see section 9) |

All three apply a shear stress of $0.01\,\mu$ on the last Voigt component, and use the Galerkin glide solver with climb disabled.

## 7. Outputs

A run writes into the simulation folder:

| File | Contents |
|---|---|
| `evl/evl_<runID>.txt` or `.bin` | the configuration: network nodes, loops, loop nodes, loop links, spherical and polyhedral inclusions |
| `evl/ddAux_<runID>.txt` or `.bin` | auxiliary data: quadrature points, mesh nodes, periodic glide-plane patches (as enabled in `DD.txt`) |
| `F/F_0.txt` | one row per output step: `runID`, time, `dt`, the average plastic distortion and its rate, then the columns added by each microstructure |
| `F/F_labels.txt` | the label of each column of `F_0.txt`, one per line |

Further points:

- **Format.** `outputBinary` in `DD.txt` selects binary or text. `tools/DDconverter` converts the binary files in `./evl` to text; it has its own `Makefile` and is not part of the CMake build.
- **Frequency.** `outputFrequency` sets the number of steps between outputs.
- **Restart.** `startAtTimeStep=-1` restarts from the last step recorded in `F/F_0.txt`.
- **Reading in Python.** `modlibUtils.py` provides `readEVLtxt`, `readAUXtxt`, `readFfile` and `getFarray(F, labels, label)`.
- **Viewing.** `DDqt` opens a simulation folder and displays its configurations.

## 8. Tests and validation

The repository has no automated test suite. What it has:

- **`tests/periodicFields/test.py`:** builds periodic dipoles in copper through `pyMoDELib` and compares their fields for increasing numbers of periodic images. It needs the pybind11 build.
- **The tutorials:** three complete simulations (section 6).
- **Continuous integration:** `.github/workflows/workflow.yml` builds the Doxygen pages on each push to `master` and publishes them. `.github/workflows/container.yml` builds the container image when a release is published; that build compiles the library and the two command-line tools, and runs two steps of the `dipoleNoise` tutorial.

## 9. Status

| Component | Status |
|---|---|
| Dislocation glide, junctions, remeshing, cross slip | in use; exercised by the tutorials |
| Periodic and finite domains, polycrystals, inclusions | in use |
| Glide-plane noise | in use; `dipoleNoise` tutorial |
| Climb | Galerkin climb solver present; coupled to the cluster-dynamics concentrations; disabled in the tutorials |
| Cluster dynamics | present; the number of species is fixed at compile time in `ClusterDynamicsParameters.h`; parameters in `Zr_CD2.txt` and `Zr_CD4.txt` |
| Tutorials | `dipoleNoise` runs as committed. `polycrystal` stops in `microstructureGenerator`, which asks for `C2G2`: `PolyCrystalFile.write` in `modlibUtils.py` writes the rotation of grain 1 only |
| Cracks | `CrackSystem` is an interface that returns zero fields; `useCracks=0` in the template |
| `Lattices_NEW/` | a second lattice implementation beside `Lattices/`; not part of the build |
| `MPI/` | headers present; not on the build's include path |
| Doxygen pages | the main page still describes the earlier header-only MODEL library |
| Manual | `manual_beta/`: the theory chapters are written; the section on running MoDELib is referenced but not yet written |

## 10. Documentation

- **API reference (Doxygen):** <https://mechanicsofdefects.github.io/MoDELib/>, built from `doxygen/` on each push to `master`.
- **Theory manual:** `manual_beta/MoDELib_manual.pdf`. Chapters: the elastic theory of discrete dislocations; discrete dislocation dynamics in MoDELib (network topology, the uniform load controller); a review of tensor calculus; the elastic fields of piecewise-straight loops.
- **Input reference:** the templates in `Library/` carry a comment on each variable. `Library/DislocationDynamics/DD.txt` lists every numerical control.
- **Git notes:** `git_basic_commands.txt`.

## 11. How to cite

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

- **BCC mobility law:** G. Po, Y. Cui, D. Rivera, D. Cereceda, T. D. Swinburne, J. Marian, N. Ghoniem, *A phenomenological dislocation mobility law for bcc metals*, Acta Mater. 119 (2016) 123–135, [doi:10.1016/j.actamat.2016.08.016](https://doi.org/10.1016/j.actamat.2016.08.016).
- **Non-singular fields, option 1:** W. Cai, A. Arsenlis, C. R. Weinberger, V. V. Bulatov, *A non-singular continuum theory of dislocations*, J. Mech. Phys. Solids 54 (2006) 561–587.
- **Non-singular fields, option 2:** G. Po, M. Lazar, D. Seif, N. Ghoniem, *Singularity-free dislocation dynamics with strain gradient elasticity*, J. Mech. Phys. Solids 68 (2014) 161–178, [doi:10.1016/j.jmps.2014.03.005](https://doi.org/10.1016/j.jmps.2014.03.005).

## 12. Repository conventions

- **Templates are copied, not edited.** A simulation copies its input files from `Library/` into its own `inputFiles/` and changes the copies. The templates stay as the reference.
- **One folder per simulation.** `inputFiles/`, `evl/` and `F/` sit together, and every tool takes that folder as its argument.
- **Input syntax.** `name=value;` with `#` starting a comment. Vectors and matrices are space-separated values, with one matrix row per line, ended by one `;`.
- **SI at the boundary only.** Input variables in SI units end in `_SI` (or name their unit, as in `dH0_eV`). Everything inside the code is in units of $b$, $\mu$ and $c_s$.
- **Headers and sources mirror each other.** A class declared in `include/<Module>/` is implemented in `src/<Module>/`.
- **Every source file states the license** in its first lines.
- **Options are read where they are used.** A class reads its own parameters from the input files with `TextFileParser`.

## 13. License and contact

MoDELib is distributed without any warranty under the GNU General Public License v2, as stated in each source file.

MoDELib was created by Giacomo Po, whose copyright notices in the tools date from 2012. The contributors are listed in `Contributors.txt` (37 names).

Questions and bug reports: please open an issue at <https://github.com/MechanicsOfDefects/MoDELib/issues>.
