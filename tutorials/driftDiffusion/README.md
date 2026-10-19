# Drift-diffusion (diffusion bias)

Spatial cluster dynamics with stress-induced drift of the mobile species, coupled to the climb of discrete dislocation loops.
This branch (`feature/diffusion-bias`) ports the DiffusionBias code (theory and implementation notes in `report/main.tex`) onto MoDELib v2.

## Model

The FEM unknown of the mobile species is the activity W = c exp(E^e/k_B T), where E^e = -tr(sigma) msVector msRelRelaxVol omega / 3
is the hydrostatic interaction energy of each species with the local stress. The physical concentration is c = W exp(-E^e/k_B T).

* Diffusion tensors depend on the macroscopic external stress through the migration activation volumes
  (`mobileSpeciesMigrationRelVolume`, crystal frame, components [11 12 13 22 23 33]).
  The FEM operators, the reaction coefficients of each grain, and the analytical dislocation kernel (inverse and determinant of D)
  are rebuilt when the external stress changes.
* Reactions, sink absorption and the eigenstrain of the mobile species use the physical concentration c.
* Boundary values, dislocation-core values and restart input are converted from c to W with the nodal (patch-average) stress.
  Output files (evl cd matrix, quadrature points) and DDqt show the physical concentration c.
* `lumpingFactor` scales the lumped climb stiffness of the discrete dislocations.

With `mobileSpeciesRelRelaxVol=0` and `mobileSpeciesMigrationRelVolume=0` the model reduces to the original one.
Both new inputs are optional (defaults: zero migration volumes, `lumpingFactor=1`), so existing material files still work.

## Running the tutorial

Build with a single mobile species and no immobile species:

    cmake .. -DMODELIB_CD_MSIZE=1 -DMODELIB_CD_ISIZE=0
    make

then

    python3 generateInputFiles.py
    ../../build/tools/DDomp/DDomp .

The tutorial places a single prismatic <a> interstitial loop (radius 100 A) at the center of a 1 micron Zr box at 625 K,
with no irradiation and no applied stress, using `Library/Materials/Zr_driftDiffusion.txt`.

## Notes on the port

* The DiffusionBias `drift` tutorial is this one. Its `drift_E0` companion (same inputs, zero interaction energy) is `tutorials/driftDiffusion_E0`.
  They generated the loop with the `aLoop` generator, which is not in MoDELib v2; here the same loop is created with `PlanarLoop`.
* The `elasticBias` input of DiffusionBias is `discreteDislocationBias` in MoDELib v2.
* The climb stiffness of the discrete dislocations is only reassembled when the climb solver updates it, so it can lag a change of external stress.
* `FirstOrderReaction.h` (sink absorption by immobile clusters) still uses the zero-stress crystal-frame diffusion tensors.
