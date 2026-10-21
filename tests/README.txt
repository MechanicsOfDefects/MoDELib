CTest categories and presets
============================

From the repository root, using the standard build/ directory:
  cmake --preset default
  cmake --build --preset default --parallel 4
  ctest --preset quick

On a fresh build, supply the usual compiler and dependency options to cmake.
The default configure preset uses Unix Makefiles and does not override compiler,
physics species counts, assertions, or dependency paths in an existing cache.
Presets require CMake/CTest 3.20 or newer. Tests must be built before running them.

Presets:
  quick        Excludes simulation, long, and benchmark labels.
  numerical    Mathematical results and numerical algorithms.
  integration  IO, API ownership, topology, and viewer behavior.
  simulation   Complete simulation regressions, including spatialCDtest.
  benchmark    Explicit performance measurements; prints timing output.
  all          Every registered test, including long runs and benchmarks.

Examples:
  ctest --preset quick --parallel 4
  ctest --preset simulation
  ctest --preset benchmark
  ctest --preset all
  ctest --preset simulation -N

Labels may overlap. They control selection and timing summaries, not physics,
correctness checks, or default execution. Bare ctest still runs everything.
For a different build directory, use the labels directly:
  ctest --test-dir /path/to/build --output-on-failure -LE '^(simulation|long|benchmark)$'
  ctest --test-dir /path/to/build --output-on-failure -L '^numerical$'

spatialCDtest has labels simulation;long;ClusterDynamics. Its original settings,
experimental data, numerical reference checks, report location and timeout are
unchanged. It is registered only for the supported 4-mobile/8-immobile build.
A preset selecting no tests reports an error instead of a misleading success.
See spatialCDtest/README.txt for retained reports and comparison tolerances.

ellipsoidalStressBenchmark measures stress-evaluation timings and ratios. It runs
serially to avoid contention with other CTest tests. It has no machine-dependent
speed threshold; correctness is checked separately by ellipsoidalInclusion.
CTest retains its output in build/Testing/Temporary/LastTest.log; the benchmark
preset also prints it. Benchmark timings are not deterministic regression data.

Assign labels explicitly next to each new add_test call. Use numerical for
algorithm/reference checks, integration for interacting components, simulation
for complete runs, and benchmark for explicit timing experiments. Add long to
expensive tests regardless of their primary category. Optional GUI/Python tests
appear only when their targets are enabled.
