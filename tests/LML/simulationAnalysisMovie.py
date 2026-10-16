# Makes movies from the output of simulationAnalysis.py: microstructure + stress-strain + densities,
# stress-strain only, and microstructure only (mp4, or gif if ffmpeg is not available).
# Usage (from tests/LML):  python3 simulationAnalysisMovie.py [simPath]
# simPath is the simulation folder. If it is not given, the simPath below is used.
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation

# -----------------------------
# Inputs
# -----------------------------
here = os.path.dirname(os.path.abspath(__file__))
simPath = "../../tutorials/mlMobility"  # simulation folder, relative to tests/LML
if len(sys.argv) > 1:
    simPath = sys.argv[1]
root = os.path.join(here, os.path.basename(os.path.normpath(simPath)))  # output folder of simulationAnalysis.py
summary_txt = os.path.join(root, "simulation_summary.txt")
fig_dir = os.path.join(root, "figures")
out_mp4 = os.path.join(root, "movie.mp4")     # output
out_gif = os.path.join(root, "movie.gif")     # fallback

# -----------------------------
# Load summary (no header)
# cols:
# 0 evl
# 1 time_s
# 2 strain_pct
# 3 stress_GPa
# 4 rho_total
# 5 rho_gliss
# 6 rho_sess
# 7 rho_boundary
# 8 rho_grainB
# 9 rho_screw
# 10 rho_edge
# 11 rho_mixed
# -----------------------------
data = np.loadtxt(summary_txt)
evl = data[:, 0].astype(int)
time_s = data[:, 1]
strain = data[:, 2]         # already %
stress = data[:, 3]         # already GPa
rho_total = data[:, 4]      # m^-2
rho_gliss = data[:, 5]
rho_sess  = data[:, 6]
rho_screw = data[:, 9]
rho_edge  = data[:, 10]
rho_mixed = data[:, 11]

# Sort by strain (or by time if you prefer)
order = np.argsort(strain)
evl = evl[order]
time_s = time_s[order]
strain = strain[order]
stress = stress[order]
rho_total = rho_total[order]
rho_gliss = rho_gliss[order]
rho_sess  = rho_sess[order]
rho_screw = rho_screw[order]
rho_edge  = rho_edge[order]
rho_mixed = rho_mixed[order]

# Pre-build image paths (must match your naming)
img_paths = [os.path.join(fig_dir, f"evl_{k:08d}.png") for k in evl]

# Filter out frames missing images (just in case)
keep = [i for i,p in enumerate(img_paths) if os.path.exists(p)]
if len(keep) != len(img_paths):
    print(f"[warn] Missing {len(img_paths)-len(keep)} images; filtering frames.")
    evl = evl[keep]
    time_s = time_s[keep]
    strain = strain[keep]
    stress = stress[keep]
    rho_total = rho_total[keep]
    rho_gliss = rho_gliss[keep]
    rho_sess  = rho_sess[keep]
    rho_screw = rho_screw[keep]
    rho_edge  = rho_edge[keep]
    rho_mixed = rho_mixed[keep]
    img_paths = [img_paths[i] for i in keep]

nframes = len(evl)
if nframes < 2:
    raise RuntimeError("Not enough frames to animate.")

# -----------------------------
# Figure layout
# left column = microstructure spanning both rows
# right side = 2 rows x 3 cols plots
# -----------------------------
plt.rcParams["figure.dpi"] = 120
fig = plt.figure(figsize=(14, 7))
gs = GridSpec(2, 4, figure=fig, width_ratios=[1.35, 1.0, 1.0, 1.0], wspace=0.35, hspace=0.35)

ax_img  = fig.add_subplot(gs[0, 0]) 
ax_ss   = fig.add_subplot(gs[0, 1]) 
ax_tot  = fig.add_subplot(gs[0, 2]) 
ax_ses  = fig.add_subplot(gs[0, 3]) 
ax_glis = fig.add_subplot(gs[1, 0]) 
ax_mix  = fig.add_subplot(gs[1, 1]) 
ax_edg  = fig.add_subplot(gs[1, 2]) 
ax_scr  = fig.add_subplot(gs[1, 3]) 


# -----------------------------
# Styling helpers
# -----------------------------
def setup_curve_axes(ax, y, ylabel):
    ax.plot(strain, y, color="k", alpha=0.25, lw=2.0)  # full curve (black, semi-opaque)
    ax.set_xlabel("Strain (%)")
    ax.set_ylabel(ylabel)
    ax.grid(True, ls="--", alpha=0.25)

# Use a colormap for the evolving curve
cmap = plt.cm.jet
norm = mpl.colors.Normalize(vmin=0, vmax=nframes - 1)

# Stress-strain (full curve)
setup_curve_axes(ax_ss, stress, "Stress (GPa)")

# Density curves (full curve)
setup_curve_axes(ax_scr, rho_screw, r"$\rho_{\mathrm{screw}}$ (m$^{-2}$)")
setup_curve_axes(ax_mix, rho_mixed, r"$\rho_{\mathrm{mixed}}$ (m$^{-2}$)")
setup_curve_axes(ax_tot, rho_total, r"$\rho_{\mathrm{total}}$ (m$^{-2}$)")
setup_curve_axes(ax_edg, rho_edge,  r"$\rho_{\mathrm{edge}}$ (m$^{-2}$)")
setup_curve_axes(ax_ses, rho_sess,  r"$\rho_{\mathrm{sessile}}$ (m$^{-2}$)")
setup_curve_axes(ax_glis, rho_gliss,  r"$\rho_{\mathrm{glissile}}$ (m$^{-2}$)")

# (optional) log-scale for densities if they span decades:
# for ax in [ax_scr, ax_mix, ax_tot, ax_edg, ax_ses]:
#     ax.set_yscale("log")

# Image axis setup
ax_img.set_axis_off()
img0 = plt.imread(img_paths[0])
im_artist = ax_img.imshow(img0)
title_artist = ax_img.text(
    0.02, 0.98, "", transform=ax_img.transAxes,
    va="top", ha="left", color="w",
    bbox=dict(facecolor="k", alpha=0.45, edgecolor="none", pad=4)
)

# Animated (colored) evolving lines (start empty)
line_ss,  = ax_ss.plot([], [], lw=2.5)
line_scr, = ax_scr.plot([], [], lw=2.5)
line_mix, = ax_mix.plot([], [], lw=2.5)
line_tot, = ax_tot.plot([], [], lw=2.5)
line_edg, = ax_edg.plot([], [], lw=2.5)
line_ses, = ax_ses.plot([], [], lw=2.5)
line_glis, = ax_glis.plot([], [], lw=2.5)

# A moving marker at current point (optional but helpful)
mk_ss,  = ax_ss.plot([], [], marker="o", ms=6)
mk_scr, = ax_scr.plot([], [], marker="o", ms=6)
mk_mix, = ax_mix.plot([], [], marker="o", ms=6)
mk_tot, = ax_tot.plot([], [], marker="o", ms=6)
mk_edg, = ax_edg.plot([], [], marker="o", ms=6)
mk_ses, = ax_ses.plot([], [], marker="o", ms=6)
mk_glis, = ax_glis.plot([], [], marker="o", ms=6)

def update(i):
    color = cmap(norm(i))

    # Update microstructure image
    img = plt.imread(img_paths[i])
    im_artist.set_data(img)
    title_artist.set_text(f"strain = {strain[i]:.3f} %")

    # Update evolving traces up to i
    xs = strain[:i+1]

    line_ss.set_data(xs, stress[:i+1]);      line_ss.set_color(color)
    line_scr.set_data(xs, rho_screw[:i+1]);  line_scr.set_color(color)
    line_mix.set_data(xs, rho_mixed[:i+1]);  line_mix.set_color(color)
    line_tot.set_data(xs, rho_total[:i+1]);  line_tot.set_color(color)
    line_edg.set_data(xs, rho_edge[:i+1]);   line_edg.set_color(color)
    line_ses.set_data(xs, rho_sess[:i+1]);   line_ses.set_color(color)
    line_glis.set_data(xs, rho_gliss[:i+1]); line_glis.set_color(color)

    # Markers at current point
    mk_ss.set_data([strain[i]], [stress[i]]);       mk_ss.set_color(color)
    mk_scr.set_data([strain[i]], [rho_screw[i]]);   mk_scr.set_color(color)
    mk_mix.set_data([strain[i]], [rho_mixed[i]]);   mk_mix.set_color(color)
    mk_tot.set_data([strain[i]], [rho_total[i]]);   mk_tot.set_color(color)
    mk_edg.set_data([strain[i]], [rho_edge[i]]);    mk_edg.set_color(color)
    mk_ses.set_data([strain[i]], [rho_sess[i]]);    mk_ses.set_color(color)
    mk_glis.set_data([strain[i]], [rho_gliss[i]]);  mk_glis.set_color(color)

    return (im_artist, title_artist,
            line_ss, line_scr, line_mix, line_tot, line_edg, line_ses,
            mk_ss, mk_scr, mk_mix, mk_tot, mk_edg, mk_ses, mk_glis)

anim = FuncAnimation(fig, update, frames=nframes, interval=60, blit=False)

# -----------------------------
# Save movie
# -----------------------------
# Try MP4 (ffmpeg). If ffmpeg isn't available, save GIF.
try:
    writer = mpl.animation.FFMpegWriter(fps=20, bitrate=3000)
    anim.save(out_mp4, writer=writer)
    print(f"Saved: {out_mp4}")
except Exception as e:
    print(f"[warn] MP4 save failed ({e}). Saving GIF instead...")
    anim.save(out_gif, writer="pillow", fps=20)
    print(f"Saved: {out_gif}")

plt.close(fig)


# -----------------------------
# Save stress-strain panel only
# -----------------------------
out_ss_mp4 = os.path.join(root, "movie_stress_strain.mp4")
out_ss_gif = os.path.join(root, "movie_stress_strain.gif")

fig_ss, ax_ss_only = plt.subplots(figsize=(6.5, 5), dpi=120)

# Full stress-strain curve in background
ax_ss_only.plot(strain, stress, color="k", alpha=0.25, lw=2.0)
ax_ss_only.set_xlabel("Strain (%)")
ax_ss_only.set_ylabel("Stress (GPa)")
ax_ss_only.grid(True, ls="--", alpha=0.25)

# Optional: lock limits so the movie does not autoscale
ax_ss_only.set_xlim(np.min(strain), np.max(strain))

ymin = np.min(stress)
ymax = np.max(stress)
ypad = 0.05 * (ymax - ymin) if ymax > ymin else 0.1
ax_ss_only.set_ylim(ymin - ypad, ymax + ypad)

# Animated curve and marker
line_ss_only, = ax_ss_only.plot([], [], lw=2.5)
mk_ss_only, = ax_ss_only.plot([], [], marker="o", ms=6)

title_ss_only = ax_ss_only.set_title("Stress-Strain")

def update_ss_only(i):
    color = cmap(norm(i))

    xs = strain[:i+1]
    ys = stress[:i+1]

    line_ss_only.set_data(xs, ys)
    line_ss_only.set_color(color)

    mk_ss_only.set_data([strain[i]], [stress[i]])
    mk_ss_only.set_color(color)

    title_ss_only.set_text(
        f"Stress-Strain | strain = {strain[i]:.3f} %, stress = {stress[i]:.3f} GPa"
    )

    return line_ss_only, mk_ss_only, title_ss_only

anim_ss = FuncAnimation(
    fig_ss,
    update_ss_only,
    frames=nframes,
    interval=60,
    blit=False
)

try:
    writer = mpl.animation.FFMpegWriter(fps=20, bitrate=3000)
    anim_ss.save(out_ss_mp4, writer=writer)
    print(f"Saved: {out_ss_mp4}")
except Exception as e:
    print(f"[warn] stress-strain MP4 save failed ({e}). Saving GIF instead...")
    anim_ss.save(out_ss_gif, writer="pillow", fps=20)
    print(f"Saved: {out_ss_gif}")

plt.close(fig_ss)



# -----------------------------
# Save microstructure panel only
# -----------------------------
out_micro_mp4 = os.path.join(root, "movie_microstructure.mp4")
out_micro_gif = os.path.join(root, "movie_microstructure.gif")
fig_micro, ax_micro = plt.subplots(figsize=(7, 7), dpi=120)

ax_micro.set_axis_off()

img0_micro = plt.imread(img_paths[0])
im_micro_artist = ax_micro.imshow(img0_micro)

title_micro_artist = ax_micro.text(
    0.02, 0.98, "",
    transform=ax_micro.transAxes,
    va="top", ha="left",
    color="w",
    bbox=dict(facecolor="k", alpha=0.45, edgecolor="none", pad=4)
)

def update_micro(i):
    img = plt.imread(img_paths[i])
    im_micro_artist.set_data(img)

    title_micro_artist.set_text(
        f"strain = {strain[i]:.3f} %"
    )

    return im_micro_artist, title_micro_artist

anim_micro = FuncAnimation(
    fig_micro,
    update_micro,
    frames=nframes,
    interval=60,
    blit=False
)

try:
    writer = mpl.animation.FFMpegWriter(fps=5, bitrate=3000)
    anim_micro.save(out_micro_mp4, writer=writer)
    print(f"Saved: {out_micro_mp4}")
except Exception as e:
    print(f"[warn] microstructure MP4 save failed ({e}). Saving GIF instead...")
    anim_micro.save(out_micro_gif, writer="pillow", fps=20)
    print(f"Saved: {out_micro_gif}")

plt.close(fig_micro)