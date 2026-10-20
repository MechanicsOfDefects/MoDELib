ClusterDynamics FEM in a fully XYZ-periodic polycrystal

Uses poly10_periodicXYZ.msh, periodicFaceIDs=-1, and explicit identity C2G
for all ten grains. Internal and cross-grain periodic faces remain GB sinks;
same-grain periodic seams do not act as surface sinks.

Preserves the former spatialCDtest material, 553 K temperature, irradiation
rate, initial cluster populations, 120 steps, dt=3.4e18 b/cs and mesh scale
(boxScaling=3093). Only ClusterDynamics is enabled; no discrete defects.
The non-convex periodic cell has a different shape from the original cube.

From the repository root:
  python3 tutorials/clusterDynamicsFEM/generateInputFiles.py
  build/tools/MicrostructureGenerator/microstructureGenerator tutorials/clusterDynamicsFEM
  build/tools/DDomp/DDomp tutorials/clusterDynamicsFEM

Plot from this tutorial directory:
  python3 plotGrowthSwellingDensity.py
Use MPLBACKEND=Agg for a noninteractive run. The script saves
StrainComponentsWithSwelling.pdf and retains the original experimental data,
+0.5 dpa shift and swelling scaling convention. Density plotting is commented.
Requires the 4-mobile/8-immobile-species build and NumPy/Matplotlib for Python.

surfaceReference.pdf and results-2026-10-19-WeakProblem preserve the old
nonperiodic results for comparison. They are not results for this new mesh.
The original numerical regression remains in tests/spatialCDtest.
