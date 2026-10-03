# Reference output

The output of the three verification runs of MoDELib 2.2.0, made with `run_simulation.ipynb` and the presets of `simulation_driver.py`. Each folder has the layout described in [`../README.md`](../README.md#output), without `evl/` (the fields, 20 to 60 MB per run) and `logs/`. To obtain them, run the preset again.

| Folder | Preset | Case | Compared with |
|---|---|---|---|
| [`cube_200nm_10dpa`](cube_200nm_10dpa) | `reference` | cube of 200 nm, 573 K, 1e-7 dpa/s, to 10 dpa | the ENNDS run of the same grain (a later model, so the comparison is qualitative) |
| [`hex_400nm_26dpa`](hex_400nm_26dpa) | `hex400` | hexagonal crystal of 400 nm, to 26 dpa in steps of 1 dpa; the equations of the D1/M1 report | the DisloCluster run behind Figure 27 of the report |
| [`hex_400nm_26dpa_DisloCluster_terms`](hex_400nm_26dpa_DisloCluster_terms) | `hex400dc` | the same, with the two terms of the DisloCluster code that differ from the report | the same DisloCluster run: a check of numerical accuracy |

Start with `report.md` in each folder: the checks, the interior values and the comparison table. `figures/loop_populations_6dpa.png` in the two hexagonal folders is the view of Figure 27 of the report.

All three were run from commit 4b9e666 with the changes of 2.2.0 not yet committed (`-dirty` in the run names); `provenance.json` records the SHA-256 of every input file and executable.

