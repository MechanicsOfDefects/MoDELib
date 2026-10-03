# Provenance: 20261002_215830_4b9e666-dirty_cube_200nm_reference

- **Status:** completed
- **Started / finished (UTC):** 2026-10-02T21:58:30+00:00 / 2026-10-02T22:48:07+00:00 (2977.3 s)
- **MoDELib:** commit `4b9e666f31c9b51d5fe783a5d50d7beaa3f2dabe` (working tree modified)
- **Machine:** Windows-11-10.0.26200-SP0, AMD64, 24 CPUs; Python 3.14.3; executables run through WSL
- **SUNDIALS in the build:** yes

## Configuration

```json
{
 "CASE": {
  "geometry": "cube_200nm",
  "tag": null
 },
 "MATERIAL": {
  "material_file": "Zr_CD4opt.txt",
  "temperature_K": 573.0,
  "dose_rate_dpa_s": 1e-07,
  "overrides": {}
 },
 "MARCH": {
  "doses": [
   0.0001,
   0.001,
   0.01,
   0.1,
   1.0,
   2.0,
   5.0,
   10.0
  ],
  "steps_per_interval": 3,
  "step_dpa": null,
  "immobile_integrator": "cvode"
 },
 "RUN": {
  "threads": null,
  "keep_all_configurations": false
 },
 "OUTPUT": {
  "mesh_cut_figures": true,
  "figure_doses": null,
  "voronoi_samples": 2000000,
  "population_dose": null
 },
 "PRESET": "reference"
}
```

## Input files (SHA-256)

- `ClusterDynamics.txt` 3ef98cfb543e376da6de1870e2930ac12c6315d64422a3de2c4398d98afbd952
- `cube_200nm.msh` edf1e36f362dd0acfe5fa3deb8745387fa3a421ee03b803a1cf95b9fabc8b430
- `DefectiveCrystal.txt` 1adbbb98aec7f8c81fed57dc1c1910f65f28b3ac24cf1c40375d26c4f03e6314
- `initialMicrostructure.txt` e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
- `polycrystal.txt` f15744b14f67ae3d091ea6058e386724cd6594b105d6f4c23a912746cab3905e
- `Zr_CD4opt.txt` a38e0fa27ddd6f01020ea7d1c5c1ca236d06a7ce68a2f72a44a2a4ea0c676a60

## Executables (SHA-256)

- `DDomp` 9aefdaf608c02fbe42cc1b73dff7d9b4213724f4e9c688597b85b0ab89e561f5
- `microstructureGenerator` eb618cb9028ba1469d1ae157fda44437aaa46d7d05d7f817658dcdada8934d09
