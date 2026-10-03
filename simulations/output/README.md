# Simulation output

One folder per run, named `<UTC YYYYMMDD_HHMMSS>_<git commit>[-dirty]_<tag>/`, written by [`../run_simulation.ipynb`](../run_simulation.ipynb). Each holds `provenance.json` and `provenance.md` beside the results; [`../README.md`](../README.md) lists the contents.

Run folders are not tracked by git; only this README is. To find runs from Python:

```python
import simulation_driver as driver
driver.find_runs()     # oldest first
driver.latest_run()
```
