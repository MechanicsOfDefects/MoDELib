import sys, string, os, re, glob
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


# -----------------------------
# Return ordered evl numbers from evl directory
# -----------------------------
def list_evl_steps(evl_dir):
    # robust: accept evl files like "evl_40000.txt" or "40000.txt" etc.
    files = glob.glob(os.path.join(evl_dir, "*.txt"))
    steps = []
    for f in files:
        base = os.path.basename(f)
        m = re.search(r"(\d+)\.txt$", base)
        if m:
            steps.append(int(m.group(1)))
    return sorted(set(steps))

# -----------------------------
# Draw mesh bounding box
# -----------------------------
def draw_box(ax, xmin, xmax, lw=2.0, alpha=0.4):
    X0, X1 = np.array(xmin, float), np.array(xmax, float)
    xs = [X0[0], X1[0]]
    ys = [X0[1], X1[1]]
    zs = [X0[2], X1[2]]
    edges = [
        ([xs[0], xs[1]], [ys[0], ys[0]], [zs[0], zs[0]]),
        ([xs[0], xs[1]], [ys[1], ys[1]], [zs[0], zs[0]]),
        ([xs[0], xs[1]], [ys[0], ys[0]], [zs[1], zs[1]]),
        ([xs[0], xs[1]], [ys[1], ys[1]], [zs[1], zs[1]]),

        ([xs[0], xs[0]], [ys[0], ys[1]], [zs[0], zs[0]]),
        ([xs[1], xs[1]], [ys[0], ys[1]], [zs[0], zs[0]]),
        ([xs[0], xs[0]], [ys[0], ys[1]], [zs[1], zs[1]]),
        ([xs[1], xs[1]], [ys[0], ys[1]], [zs[1], zs[1]]),

        ([xs[0], xs[0]], [ys[0], ys[0]], [zs[0], zs[1]]),
        ([xs[1], xs[1]], [ys[0], ys[0]], [zs[0], zs[1]]),
        ([xs[0], xs[0]], [ys[1], ys[1]], [zs[0], zs[1]]),
        ([xs[1], xs[1]], [ys[1], ys[1]], [zs[0], zs[1]]),
    ]
    for ex, ey, ez in edges:
        ax.plot(ex, ey, ez, linewidth=lw, alpha=alpha, color="k")


# -----------------------------
#   Returns character angle in degrees in [0, 180), or None if angle is undefined (no slip system).
# -----------------------------
def angle_atan2_from_link(link, eps=1e-12):

    # Burgers vector
    b = np.array(link.burgers(), dtype=float)
    nb = np.linalg.norm(b)
    if nb < eps:
        return None
    b_unit = b / nb

    # Line direction
    t = np.array(link.chord(), dtype=float)
    nt = np.linalg.norm(t)
    if nt < eps:
        return None
    t_unit = t / nt

    # Slip plane normal (required for signed angle)
    ss = link.slipSystem()
    if ss is None:
        return None   # sessile / undefined plane
    n = np.array(ss.unitNormal, dtype=float)
    nn = np.linalg.norm(n)
    if nn < eps:
        return None
    n_unit = n / nn

    # atan2 formulation
    cos_theta = np.dot(b_unit, t_unit)
    sin_theta = np.dot(np.cross(b_unit, t_unit), n_unit)
    theta = np.degrees(np.arctan2(sin_theta, cos_theta))
    theta = (theta % 180 + 180) % 180
    return theta


# -----------------------------
#   Returns Network Lengths.
# -----------------------------
def network_length_from_networklinks(DN,screw_dev=5.0,edge_lo=80.0,edge_hi=100.0):
    bulkGliss = 0.0
    bulkGlissEdge = 0.0
    bulkGlissScrew = 0.0
    bulkGlissMixed = 0.0
    bulkSess  = 0.0
    boundary  = 0.0
    grainB    = 0.0
    skipped_zeroB = 0
    skipped_null  = 0
    skipped_theta = 0

    for k in DN.networkLinks().keys():
        link = DN.networkLinks().getRef(k)
        if link is None:
            skipped_null += 1
            continue
        if link.hasZeroBurgers():
            skipped_zeroB += 1
            continue

        L = np.linalg.norm(np.array(link.chord(), dtype=float))  # chord() is a VectorDim

        if link.isBoundarySegment():
            boundary += L
        elif link.isGrainBoundarySegment():
            grainB += L
        else:
            # bulk
            if link.isSessile():
                bulkSess += L
            else:
                bulkGliss += L
                theta = angle_atan2_from_link(link)
                if theta is None:
                    skipped_theta += 1
                    print("theta is None for a bulk glissile link (unexpected).")
                    continue
                    # raise RuntimeError("theta is None for a bulk glissile link (unexpected).")
                if (theta <= screw_dev) or (theta >= 180.0 - screw_dev):
                    bulkGlissScrew += L
                elif (edge_lo <= theta <= edge_hi):
                    bulkGlissEdge += L
                else:
                    bulkGlissMixed += L

    return (bulkGliss, bulkSess, boundary, grainB, bulkGlissScrew, bulkGlissEdge, bulkGlissMixed)


def dislocationSegments_from_Network(DN):
    node_pos = {}
    for k in DN.networkNodes().keys():
        n = DN.networkNodes().getRef(k)
        node_pos[n.networkID()] = np.array(n.position(), float)

    # Collect network link segments
    segments = []
    skipped = 0
    for lk in DN.networkLinks().keys():
        link = DN.networkLinks().getRef(lk)
        theta = angle_atan2_from_link(link)
        seg_type = "glissile" if link.isGlissile() else "sessile"

        if link.hasZeroBurgers():
            continue
        if link.isBoundarySegment():
            continue

        sid = link.source.networkID()
        tid = link.sink.networkID()
        a = node_pos.get(sid, None)
        b = node_pos.get(tid, None)
        if a is None or b is None:
            skipped += 1
            continue

        segments.append((lk, a, b, seg_type, theta))
    return segments



# -----------------------------
#   Returns the material file of a simulation, read from inputFiles/polycrystal.txt
# -----------------------------
def material_file(simPath):
    with open(os.path.join(simPath, 'inputFiles', 'polycrystal.txt')) as f:
        for line in f:
            m = re.match(r'\s*materialFile\s*=\s*([^;]+);', line)
            if m:
                return os.path.join(simPath, 'inputFiles', m.group(1).strip())
    raise RuntimeError("materialFile not found in "+os.path.join(simPath, 'inputFiles', 'polycrystal.txt'))
