# mlMobility

Dislocation dynamics with a machine-learned mobility law. Three 11 nm prismatic loops on
slip systems 5, 6 and 11 in BCC Fe at 300 K, in the box of an MD simulation
(33.8 x 35.6 x 176.5 nm), under strain control along z at 2.5e7 1/s.

## The mobility

`Library/Materials/Fe_320.txt` selects a Python mobility for the a/2<111>{110} slip systems:

```
mobility_a/2<111>{110}=python ../../../python/MobilityLaws/Fe320_GPR/DDMobilityPy_Fe320.py;
```

The module evaluates a Gaussian-process regression model trained on MD velocities of Fe at
320 K. Its 8 features are the sorted eigenvalues of the stress, the absolute projections of
their eigenvectors on the line direction, the temperature and the character angle. The model
predicts sqrt(velocity); the module returns the velocity in m/s. A relative module path is
relative to the folder of the material file (`inputFiles/`).

Other Python mobilities are in `python/MobilityLaws`: two W-xRe screw-mobility components
(`WRe_GPR`, `WRe_Curtin`), selected by `Library/Materials/WRe_GPR.txt` and `WRe_Curtin.txt`.

## Settings that differ from the defaults

- `crossSlipAngle_deg=30`, so that near-screw dipole arms (20-30 degrees) can cross slip and annihilate.
- `alphaLineTension=4`, `Lmin=5`, `Lmax=25`, `quadPerLength=0.1`, `remeshFrequency=20`.
- `mu0_SI=49e9` to match the MD shear modulus.

## Requirements

MoDELib built with `USE_PYBIND11=ON`, and numpy, pandas, joblib and scikit-learn installed in
the Python that the build links against.

## Running

```
python3 generateInputFiles.py
../../build/tools/MicrostructureGenerator/microstructureGenerator .
../../build/tools/DDomp/DDomp .
```

Figures of the segments colored by character (screw/edge/mixed), densities, stress-strain
and movies: see `tests/LML`, for example `python3 simulationAnalysis.py ../../tutorials/mlMobility`.
