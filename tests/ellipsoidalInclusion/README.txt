Ellipsoidal inclusion tests

These programs use the implementation, base class, and spherical inclusion
from this repository. No mock headers or copied implementation are used.
The production files are:
  include/InclusionMicrostructure/EllipsoidalInclusion.h
  src/InclusionMicrostructure/EllipsoidalInclusion.cpp

Focused build from the MoDELib2 root:
  cmake -S tests/ellipsoidalInclusion -B build/ellipsoidalInclusion \
    -DCMAKE_BUILD_TYPE=Release
  cmake --build build/ellipsoidalInclusion -j2
  ctest --test-dir build/ellipsoidalInclusion --output-on-failure

If dependencies are not found, supply:
  -DEIGEN3_INCLUDE_DIR=/directory/containing/Eigen
  -DELLIPSOIDAL_BOOST_INCLUDE_DIR=/directory/containing/boost
  -DMODELIB_EXTRA_INCLUDE_DIRS=/additional/dependency/include
Eigen3::Eigen is preferred when its package configuration is available.
The focused build avoids unrelated application, GUI, and Python targets.

Alternatively, in the normal MoDELib build, set BUILD_TESTS=ON
and build testEllipsoidalInclusion. The tests then link the MoDELib library.
This uses the repository's existing BUILD_TESTS option (default ON).

Programs
- testEllipsoidalInclusionIO: text record field order, round trips, stress
  reconstruction, defaults, and malformed/truncated input; registered with CTest.
- testEllipsoidalInclusion: 833 correctness checks, registered with CTest.
  Covers general strain sphere limits, rotations, traction, equilibrium,
  scaling, matrix-axis validation, cache refresh and concurrent cache reads,
  and the default 1e-3 tolerance plus explicit high-accuracy overrides.
- benchmarkEllipsoidalInclusion: stress timings against SphericalInclusion;
  explicitly uses 1e-10 for comparison with the original timing study.
- studyEllipsoidalTolerance [tolerance]: speed and accuracy study, with 1e-10
  as the timing baseline and 1e-12 as the accuracy reference. The optional
  trial tolerance defaults to 1e-4. To test the production default:
    ./build/ellipsoidalInclusion/studyEllipsoidalTolerance 1e-3

Benchmarks are not CTest pass/fail tests. Use Release builds for timings.
Tolerance-study errors are full-tensor relative Frobenius errors and are
sampled results, not rigorous bounds on individual stress components.

Production constructor:
  EllipsoidalInclusion(C, axes, eT, nu, mu, mobilityReduction, phaseID, sph,
                      quadratureTolerance=1e-3, quadratureMaxDepth=30)
Columns of axes are mutually orthogonal semiaxis vectors, not unit vectors.
Scope: 3D, equal isotropic inclusion/matrix moduli, infinite matrix, uniform
symmetric eigenstrain. DefectiveCrystalConfigIO supports text and binary records; inclusion
individual generation and simulation initialization are connected; DDqt visualization is connected. Stress uses cached interior values and adaptive
Gauss-Kronrod exterior integrals with safeguarded Newton for the lower bound.

See IMPLEMENTATION.txt for usage examples, mathematical formulation, and
limitations. reports/ preserves the original validation and performance
reports; its timings and package-relative paths refer to the original work.

EllipsoidalInclusionIO<3> follows the existing inclusion IO record interface:
construct from an inclusion, explicit fields, or std::stringstream; write
with operator<<. A record has 24 whitespace-separated fields:
  inclusionID C[3] axes[9, row-major] eT[9, row-major] mobilityReduction phaseID
The axes matrix still stores semiaxis vectors in its columns. No newline is
appended by operator<<. Check the input stream's fail() after parsing.
Geometry and material validation occur when constructing EllipsoidalInclusion.
As with spherical records, supply elastic constants and the second-phase
pointer when reconstructing the physical inclusion. Quadrature settings are
not serialized.

DefectiveCrystalConfigIO provides const and mutable ellipsoidalInclusions() accessors.
Both text and binary formats store the ellipsoidal count immediately after
spherical inclusions and before polyhedral inclusions; record blocks use the
same order. clear(), print(), and bin2txt() include the new records.
The header now has 11 counts: nodes, loops, loopLinks, loopNodes, spherical
inclusions, ellipsoidal inclusions, polyhedral inclusions, polyhedral nodes,
polyhedral edges, displacement rows, and CD rows.
This is a format change: old ten-count files require conversion. For old text
files, insert a line containing 0 after the spherical-inclusion count (line 5).
Binary files require conversion with a reader for their original format;
the existing native size_t/raw-record binary format remains ABI-dependent.
The binary reader also now consumes the displacement count that the writer
already emitted, keeping the record offset aligned.

- testDDconfigEllipsoidalIO checks mixed inclusion block/count ordering,
  text and binary round trips, binary-to-text conversion, trailing text
  matrices, const access, clear(), and empty configurations.

Individual specification
EllipsoidalInclusionIndividualSpecification describes exactly one inclusion.
It has default and filename constructors and is exposed in pyMoDELib with
the same writable fields as the C++ struct:
  center: 3-vector, simulation length units [b]
  crystalAxes_SI: 3x3 matrix, columns are orthogonal crystal-frame semiaxis vectors [m]
  crystalEigenDistortion: crystal-frame 1x9 row, dimensionless, ordered 11,12,13,21,22,23,31,32,33
  velocityReductionFactor: scalar (default 1)
  phaseID: integer (default 0)
  allowOverlap, allowOutside: booleans (default false)
Default geometry and eigen-distortion are zero; set geometry before use.
File constructors require all seven fields and check matrix dimensions via
TextFileParser. Geometry/material validation remains with EllipsoidalInclusion.
Use symmetric eigen-distortion for the current physical inclusion class.
The specification stores SI axes unchanged; generation computes C2G * crystalAxes_SI /
b_SI using the grain containing the center, while leaving the center in simulation units.

Example input: Library/Microstructures/ellipsoidalInclusionIndividual.txt.
Python example (after rebuilding pyMoDELib):
  spec = pyMoDELib.EllipsoidalInclusionIndividualSpecification()
  spec.center = np.array([200., 0., 0.])
  spec.crystalAxes_SI = np.diag([3e-8, 2e-8, 1e-8])
  spec.crystalEigenDistortion = np.array([[1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3]])
Generate from Python using:
  microstructureGenerator.addEllipsoidalInclusionIndividual(spec)
The C++ method has the same name. File dispatch accepts type=EllipsoidalInclusions
and style=Individual or individual. Density generation is also supported; see below.
The generator rotates crystalAxes_SI to the global frame and converts to simulation units, validates finite orthogonal
nonzero semiaxes and symmetric eigen-distortion, and inserts an IO record.
All insertion overloads allocate IDs above the maximum existing inclusion ID.
The center must always be inside the mesh to determine the grain orientation.
allowOutside=false additionally requires no contact with any
external boundary triangle, including cavity boundaries. Contact is rejected.
allowOverlap=false rejects contact/overlap with existing spheres and ellipsoids.
For polyhedra, overlapping axis-aligned bounding boxes are conservatively
rejected; this can reject otherwise valid placements. These checks are local
and do not consider periodic images. Rejected placements print a reason and
leave the configuration unchanged; malformed input throws. This follows the
existing generator's void add-method convention.
The existing spherical/polyhedral generators retain their own placement checks.
InclusionMicrostructure loads ellipsoids into the shared Eshelby container, so
stress is included in the existing stress/periodic-image sum. Simulation output
preserves ellipsoid records for restart. The grain at the center supplies the
second phase; phaseID must exist there. As for spheres, centers outside the mesh
are skipped by the loader even if generation used allowOutside=true. DDqt renders ellipsoids with their semiaxes and orientation. Set physics=InclusionMicrostructure (along with any other
required physics) in DefectiveCrystal.txt to enable the simulation component.
testEllipsoidalSpecification checks defaults, file parsing, matrix ordering,
unit conversion into a physical inclusion, and malformed matrix dimensions.

The placement test covers sphere limits, ellipsoid/sphere contact, rotation,
scaling, and boundary-face intersections (including degenerate triangles).

With the full-library build, testEllipsoidalGenerator additionally constructs a
real mesh/material fixture to test insertion, unit conversion, overlap/outside
flags, invalid input, mixed/sparse inclusion IDs, file dispatch, and output IO.

The full integration test also verifies mixed-shape loading, second-phase
metadata, stress superposition, binary restart, invalid phases, duplicate IDs,
and empty reinitialization.

DDqt displays ellipsoids using x=C+axes*u for unit-sphere surface coordinates u.
The existing inclusion visibility and opacity controls apply to ellipsoids.
The rendering uses simulation units, preserves all three semiaxis vectors,
and transforms surface normals for correct lighting. Ellipsoidal geometry is
replaced on each frame update, including frames with no ellipsoids.
When DDqt and BUILD_TESTS are enabled, ellipsoidalInclusionPolyData verifies
surface positions, semiaxis endpoints, rotations, signed/permuted axes,
sphere limits, and transformed normals without requiring a display.

All inclusion generators accept crystal-frame distortion components (row-major).
Ellipsoids, polyhedra, and spherical density use crystalEigenDistortion;
individual spheres use crystalEigenDistortions, one row per sphere.
Each tensor is converted to C2G * crystalEigenDistortion * C2G.transpose()
before insertion. Spheres use their center grain; polyhedra require all boundary
vertices in one grain, matching the simulation loader. DefectiveCrystalConfigIO stores global
tensors and restart loading does not rotate them again. Old input/property names
must be updated; no ambiguous global-frame aliases are provided.

Ellipsoidal density generation
Use EllipsoidalInclusionDensitySpecification and addEllipsoidalInclusionDensity
from C++ or Python, or type=EllipsoidalInclusions; style=density; in a file.
Library/Microstructures/ellipsoidalInclusionsDensity.txt contains all fields.
crystalAxes contains three orthogonal unit-vector columns (dimensionless).
axesLognormalDistributionMean_SI and axesLognormalDistributionStd_SI each contain
three values: arithmetic means >0 and standard deviations >=0 in meters.
The generator independently samples each semiaxis length. Internally,
sigma^2=log(1+(std/mean)^2), mu=log(mean)-sigma^2/2, and length=exp(mu+sigma*Z).
Zero standard deviation gives exactly the specified mean length. No reference
size or dimensionless multiplier is needed. Directions are preserved.
This replaces the old density crystalAxes_SI and lognormal M/S inputs.
The individual specification still accepts sized crystalAxes_SI vectors.
The accepted population is conditioned on containment, overlap, and the supported
aspect-ratio range, so its size distribution can differ from the proposal distribution.
targetDensity is number/m^3; ceil(targetDensity * meshVolume_SI) NEW inclusions are
requested. Zero density is a no-op. maxAttempts limits consecutive failed placements;
exhaustion throws with progress and retains already accepted records.
patternVector_SI is a crystal-frame plane normal times spacing in meters. Centers
are projected to planes n_global dot x = k * spacing through the global origin,
with n_global=C2G*normalized(patternVector_SI). Zero disables the pattern.
A projected center must remain in its original grain. Axes and eigen-distortion
are rotated using that grain. Individual placement flags/checks also apply.
The size sampler uses a default-seeded mt19937, like spherical generation;
center positions use the existing random lattice-point sampler.
