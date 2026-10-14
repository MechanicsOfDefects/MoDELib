# Provenance: 20261002_232731_4b9e666-dirty_hex400_DisloCluster_terms

- **Status:** completed
- **Started / finished (UTC):** 2026-10-02T23:27:32+00:00 / 2026-10-03T00:06:38+00:00 (2346.6 s)
- **MoDELib:** commit `4b9e666f31c9b51d5fe783a5d50d7beaa3f2dabe` (working tree modified)
- **Machine:** Windows-11-10.0.26200-SP0, AMD64, 24 CPUs; Python 3.14.3; executables run through WSL
- **SUNDIALS in the build:** yes

## Configuration

```json
{
 "CASE": {
  "geometry": "hex_400nm",
  "tag": "hex400_DisloCluster_terms"
 },
 "MATERIAL": {
  "material_file": "Zr_CD4opt.txt",
  "temperature_K": 573.0,
  "dose_rate_dpa_s": 1e-07,
  "overrides": {
   "loopClusteringNucleation": "1",
   "loopNucleationPerReaction": "1",
   "coalescenceClimbBurgers": "0"
  }
 },
 "MARCH": {
  "doses": [
   1.0,
   6.0,
   11.0,
   16.0,
   21.0,
   26.0
  ],
  "steps_per_interval": 3,
  "step_dpa": 1.0,
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
  "population_dose": 6.0
 },
 "PRESET": "hex400dc"
}
```

## Input files (SHA-256)

- `ClusterDynamics.txt` 3ef98cfb543e376da6de1870e2930ac12c6315d64422a3de2c4398d98afbd952
- `DefectiveCrystal.txt` 2c9b9d929b2e071f372da9d1a264e9609c8a9cf14a041b09a09c2bce4d884f1b
- `hex_400nm.msh` 268415d6f1089b28ccfd226f9deb9d3605585afb96cf68e7edef8f5ae7aa1054
- `initialMicrostructure.txt` e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
- `polycrystal.txt` db1b7862b02f418477bc7490b9b6956737adac4f614fcb1f9359d982b4f8a2eb
- `Zr_CD4opt.txt` e05d6e569e3361e4a8bd7cb178859481200372830bb210178421157c29ba2ad7

## Executables (SHA-256)

- `DDomp` b0ff3cee071a4af41b3e28c4eac2dd8657d63d652bffbd38fa1068519180f569
- `microstructureGenerator` eb618cb9028ba1469d1ae157fda44437aaa46d7d05d7f817658dcdada8934d09
