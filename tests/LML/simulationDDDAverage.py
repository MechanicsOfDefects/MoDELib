# Averages the stress-strain and density curves of several runs (the simulation_summary.txt files
# written by simulationAnalysis.py) and saves the figures and CSVs to ./DDD_average_plots.
# Usage (from tests/LML):  python3 simulationDDDAverage.py [pattern]
# pattern is a glob for the output folders of simulationAnalysis.py, for example "./mlMobility*".
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    "figure.dpi": 120,
    "font.size": 11,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
})

plt.rcParams["text.usetex"] = False


# ============================================================
# User inputs
# ============================================================

# Option 1: automatically find run folders.
# This will find files like ./run1/simulation_summary.txt, ./run2/simulation_summary.txt, etc.
here = os.path.dirname(os.path.abspath(__file__))
runPattern = "./mlMobility*"  # output folders of simulationAnalysis.py, relative to tests/LML
if len(sys.argv) > 1:
    runPattern = sys.argv[1]
summary_files = sorted(glob.glob(os.path.join(here, runPattern, "simulation_summary.txt")))

# Option 2: manually list folders instead. Uncomment and edit if needed.
# ddd_roots = ["./run1", "./run2", "./run3"]
# summary_files = [os.path.join(r, "simulation_summary.txt") for r in ddd_roots]

# Plot/averaging settings
N_GRID = 600                 # number of strain points used for the average
MIN_RUNS_FOR_AVERAGE = 1     # use 2 or 3 if you do not want the avg where only one run exists
X_LIMIT = None               # example: (0.01, 10). Use None to auto-scale to data.
OUT_DIR = os.path.join(here, "DDD_average_plots")
os.makedirs(OUT_DIR, exist_ok=True)

print(f"Found {len(summary_files)} DDD summary files:")
for f in summary_files:
    print("  ", f)

if len(summary_files) == 0:
    raise FileNotFoundError("No simulation_summary.txt files found. Edit summary_files or ddd_roots above.")


# ============================================================
# Load and average functions
# ============================================================

COLS = {
    "strain": 2,       # percent strain, based on your current notebook
    "stress": 3,       # GPa
    "rho_total": 4,    # m^-2
    "rho_gliss": 5,
    "rho_sess": 6,
    "rho_screw": 9,
    "rho_edge": 10,
    "rho_mixed": 11,
}

def _unique_sorted_xy(x, y):
    """Sort by x and average y values if duplicate x values exist."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]

    order = np.argsort(x)
    x, y = x[order], y[order]

    xu, inv = np.unique(x, return_inverse=True)
    yu = np.zeros_like(xu, dtype=float)
    counts = np.zeros_like(xu, dtype=int)
    np.add.at(yu, inv, y)
    np.add.at(counts, inv, 1)
    yu = yu / counts
    return xu, yu

def load_dd_summary(path):
    """Load one simulation_summary.txt into a dictionary of arrays."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]

    run = {"path": path, "label": os.path.basename(os.path.dirname(path))}
    for key, col in COLS.items():
        if col < data.shape[1]:
            run[key] = data[:, col].astype(float)
        else:
            run[key] = np.full(data.shape[0], np.nan)

    # Make every y-array consistent with sorted/unique strain.
    x = run["strain"]
    for key in list(COLS.keys()):
        if key == "strain":
            continue
        xu, yu = _unique_sorted_xy(x, run[key])
        run[key + "_x"] = xu
        run[key + "_y"] = yu

    # Store a clean strain vector too.
    run["strain_clean"] = _unique_sorted_xy(x, x)[0]
    return run

def make_common_grid(runs, metric, n_grid=600):
    """Common strain grid spanning the union of all run strain ranges."""
    mins, maxs = [], []
    for r in runs:
        x = r[metric + "_x"]
        if len(x) > 1:
            mins.append(np.nanmin(x))
            maxs.append(np.nanmax(x))
    if len(mins) == 0:
        raise ValueError(f"No valid data found for {metric}")
    return np.linspace(min(mins), max(maxs), n_grid)

def interpolate_runs_to_grid(runs, metric, grid):
    """Interpolate each run. Outside each run's strain range stays NaN, not extrapolated."""
    Y = np.full((len(runs), len(grid)), np.nan)
    for i, r in enumerate(runs):
        x = r[metric + "_x"]
        y = r[metric + "_y"]
        if len(x) < 2:
            continue
        inside = (grid >= np.nanmin(x)) & (grid <= np.nanmax(x))
        Y[i, inside] = np.interp(grid[inside], x, y)
    return Y

def average_metric(runs, metric, n_grid=600, min_runs=1):
    grid = make_common_grid(runs, metric, n_grid=n_grid)
    Y = interpolate_runs_to_grid(runs, metric, grid)
    n_valid = np.sum(np.isfinite(Y), axis=0)

    avg = np.full(len(grid), np.nan)
    ok = n_valid >= min_runs
    avg[ok] = np.nanmean(Y[:, ok], axis=0)

    std = np.full(len(grid), np.nan)
    std[ok] = np.nanstd(Y[:, ok], axis=0)
    return grid, avg, std, n_valid, Y

runs = [load_dd_summary(f) for f in summary_files]
print(f"Loaded {len(runs)} runs.")


# ============================================================
# Plot helpers
# ============================================================

def plot_metric_average(ax, runs, metric, ylabel, *, n_grid=N_GRID, min_runs=MIN_RUNS_FOR_AVERAGE,
                        show_std=False, individual_color="0.65"):
    # Individual DDD runs
    for r in runs:
        ax.plot(r[metric + "_x"], r[metric + "_y"], color=individual_color, lw=1.2, alpha=0.45)

    # Average curve
    grid, avg, std, n_valid, Y = average_metric(runs, metric, n_grid=n_grid, min_runs=min_runs)
    ax.plot(grid, avg, color="black", lw=2.8, label=f"DDD average, n ≥ {min_runs}")

    # Optional +/- one standard deviation band
    if show_std:
        ok = np.isfinite(avg) & np.isfinite(std)
        ax.fill_between(grid[ok], avg[ok] - std[ok], avg[ok] + std[ok], color="black", alpha=0.15, lw=0)

    ax.set_xlabel(r"Strain (%)")
    ax.set_ylabel(ylabel)
    ax.grid(True, ls="--", alpha=0.25)
    if X_LIMIT is not None:
        ax.set_xlim(*X_LIMIT)
    return grid, avg, std, n_valid

def add_run_count_axis(ax, grid, n_valid):
    """Optional small right-axis showing how many runs contribute at each strain."""
    ax2 = ax.twinx()
    ax2.plot(grid, n_valid, color="0.25", lw=1.0, ls=":", alpha=0.7)
    ax2.set_ylabel("Runs in average")
    ax2.set_ylim(0, max(n_valid) + 0.5)
    return ax2


# ============================================================
# Stress-strain plot only
# ============================================================

fig, ax = plt.subplots(figsize=(7, 5))
grid, avg, std, n_valid = plot_metric_average(
    ax, runs, "stress", r"Stress (GPa)", show_std=False
)
add_run_count_axis(ax, grid, n_valid)
ax.legend(frameon=False, loc="best")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "DDD_stress_strain_average.pdf"), dpi=600)
fig.savefig(os.path.join(OUT_DIR, "DDD_stress_strain_average.png"), dpi=300)
plt.close(fig)


# ============================================================
# Stress + density components, same style as your original figure
# ============================================================

fig = plt.figure(figsize=(14, 7))
gs = GridSpec(2, 3, figure=fig, wspace=0.35, hspace=0.35)

axes = {
    "stress":    fig.add_subplot(gs[0, 0]),
    "rho_total": fig.add_subplot(gs[0, 1]),
    "rho_sess":  fig.add_subplot(gs[0, 2]),
    "rho_edge":  fig.add_subplot(gs[1, 0]),
    "rho_screw": fig.add_subplot(gs[1, 1]),
    "rho_mixed": fig.add_subplot(gs[1, 2]),
}

labels = {
    "stress": r"Stress (GPa)",
    "rho_total": r"$\rho_{\mathrm{total}}$ (m$^{-2}$)",
    "rho_sess": r"$\rho_{\mathrm{sessile}}$ (m$^{-2}$)",
    "rho_edge": r"$\rho_{\mathrm{edge}}$ (m$^{-2}$)",
    "rho_screw": r"$\rho_{\mathrm{screw}}$ (m$^{-2}$)",
    "rho_mixed": r"$\rho_{\mathrm{mixed}}$ (m$^{-2}$)",
}

for metric, ax in axes.items():
    plot_metric_average(ax, runs, metric, labels[metric], show_std=False)

axes["stress"].legend(frameon=False, loc="best")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "DDD_density_stress_average.pdf"), dpi=600)
fig.savefig(os.path.join(OUT_DIR, "DDD_density_stress_average.png"), dpi=300)
plt.close(fig)


# ============================================================
# Optional: save average curves to CSV
# ============================================================

metrics = ["stress", "rho_total", "rho_sess", "rho_edge", "rho_screw", "rho_mixed"]
for metric in metrics:
    grid, avg, std, n_valid, Y = average_metric(runs, metric, n_grid=N_GRID, min_runs=MIN_RUNS_FOR_AVERAGE)
    out = np.column_stack([grid, avg, std, n_valid])
    header = "strain_percent,average,std,n_runs"
    np.savetxt(os.path.join(OUT_DIR, f"DDD_average_{metric}.csv"), out, delimiter=",", header=header, comments="")

print(f"Saved figures and CSVs to: {OUT_DIR}")

