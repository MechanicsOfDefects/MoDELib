spatialCDtest: original surface-bounded ClusterDynamics FEM regression

Labels: simulation;long;ClusterDynamics. Excluded from the quick preset.
Run all simulation regressions with: ctest --preset simulation

From the repository root:
  cmake --build build --target testSpatialCD
  ctest --test-dir build -R '^spatialCDtest$' --output-on-failure

The original settings are preserved: unitCube_15K, no periodic faces, 553 K,
120 steps, dt=3.4e18 b/cs, original material and initial cluster populations,
no discrete microstructure. inputFiles freezes the generated settings.
All 120 strain/swelling samples are compared with the successful 2026-10-19
reference (relative 1e-4 plus absolute 1e-10 strain tolerance). The test also
checks finite values, run IDs, time/dose sampling, in-plane symmetry and shear.

Python 3, NumPy and Matplotlib are required. CMake Python3_EXECUTABLE selects
the interpreter; MODELIB_CD_MSIZE=4 and MODELIB_CD_ISIZE=8 are required.
The test is available with BUILD_TOOLS=OFF. Labels: integration, slow,
ClusterDynamics. Runtime is approximately eight minutes; timeout 30 minutes.

Each run saves a timestamped reports/ directory containing exact inputs and
hashes, F and evl outputs, log, all 60 experimental points with original units,
PDF/PNG plots, numerical CSV, reference and original tutorial PDF, result.json.
Reports are ignored by Git. Historical reports have been preserved.
Experimental data are for visual comparison, not an exact pass/fail target.
Plotting retains the +0.5 dpa shift and original swelling scaling convention.

The separate XYZ-periodic example is in tutorials/clusterDynamicsFEM.
