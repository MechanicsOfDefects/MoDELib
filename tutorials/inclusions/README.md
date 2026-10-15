# Inclusions

This tutorial uses a nonperiodic tungsten cube of approximately 1 micrometer at
300 K, prismatic loops at `1e12 m^-2`, and uniform applied stress sigma_33 = 0.01
in simulation units. `InclusionMicrostructure` enables precipitate stress.
`C2G` is identity.

There are three separate `EllipsoidalInclusionDensitySpecification` inputs:

| File | Platelet base | Normal |
| --- | --- | --- |
| `ellipsoidalInclusionsDensity_100.txt` | (100) | [100] |
| `ellipsoidalInclusionsDensity_010.txt` | (010) | [010] |
| `ellipsoidalInclusionsDensity_001.txt` | (001) | [001] |

Each precipitate is an oblate ellipsoid with a circular 40 nm diameter and a
3 nm full thickness: its semiaxes are **20, 20, and 1.5 nm**. Size spread is zero,
so every platelet has these dimensions and the mean diameter is exactly 40 nm.
Independent nonzero in-plane size spreads would produce elliptical bases.
The third semiaxis is the platelet normal.

`patternVector_SI = [0, 0, 0]` disables positional planes: centers are randomly
placed throughout the mesh. The total target number density `20e19 m^-3` is split
equally among the families; each count is rounded up separately using the actual
mesh volume (67 per family in the default box). Overlap and boundary crossing
are disabled. The individual-inclusion example remains commented out.

Principal eigen-distortions remain `(0.01, -0.005, 0.02)` in each family's local
semiaxis frame, with the third component along its normal. The script expresses
them in crystal coordinates; the generator applies the grain's C2G tensor
rotation. Phase ID is 0 and the velocity reduction factor is 1.0.

## Run

With NumPy installed and MoDELib built in `build`, run from this directory:

```sh
python3 generateInputFiles.py
../../build/tools/MicrostructureGenerator/microstructureGenerator .
python3 verifyConfiguration.py
../../build/tools/DDqt/DDqt
```

In DDqt, open this simulation folder, load frame 0, and enable inclusions.
Rotate the camera to inspect the three platelet orientations.
To evolve the dislocations:

```sh
../../build/tools/DDomp/DDomp .
```

Adjust executable paths for your build. Both Python scripts accept
`--output-dir /path/to/simulation`. Verification checks the counts, circular
diameters, thickness, {100} normals, eigen-distortions, and disabled pattern.

## Modify the example

Edit `axesLognormalDistributionMean_SI`, `axesLognormalDistributionStd_SI`, `totalTargetDensity`, and `principalEigenDistortion` in
`generateInputFiles.py`; adjust verification expectations if changing dimensions.
The three density files are registered in `initialMicrostructure.txt` alongside
the prismatic-loop file. Old individual/density input files are unused.
Regeneration overwrites inputs but does not remove old simulation results; use
a fresh output folder for each run. Generated inputs and results are ignored
by Git.

The density input `crystalAxes` contains orthogonal unit-vector columns only;
it has no length scale. `axesLognormalDistributionMean_SI` and
`axesLognormalDistributionStd_SI` directly specify the arithmetic mean and
standard deviation of each semiaxis length in meters. The generator samples
those lengths independently and applies them to the unit vectors internally.
The mean values are `[20e-9, 20e-9, 1.5e-9]`; all standard deviations are zero.
