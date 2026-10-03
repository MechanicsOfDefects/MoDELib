"""Gathers the output of a run of simulation_driver: volume means, figures and the report.

The fields are read from the cluster-dynamics block of the configuration files evl_<n>.txt, whose
rows follow the finite-element nodes of evl/cdNodes.txt. Units in the figures and tables:
  mobile species      clusters per atom (the code stores defects per atom: C = c/|m|)
  number densities    loops per m^3
  defect contents     defects per atom
"""

import csv
import json
from pathlib import Path

import numpy as np

kB_eV = 8.617333262e-5


def format_dose(dose):
    return ('%g' % dose).replace('.', 'p')


class Run:
    """The fields of a run at its snapshot doses."""

    def __init__(self, run_dir):
        self.dir = Path(run_dir)
        self.provenance = json.loads((self.dir / 'provenance.json').read_text(encoding='utf-8'))
        self.record = json.loads((self.dir / 'run_record.json').read_text(encoding='utf-8'))
        derived = self.provenance['derived']
        self.b_SI = derived['b_SI']
        self.omega_SI = derived['omega_SI']
        self.omega = self.omega_SI / self.b_SI ** 3
        self.size_b = derived['size_b']
        self.ms = np.array(derived['mobile_species'], float)
        self.nF = derived['families']
        self.mSize = len(self.ms)
        nodes = np.loadtxt(self.dir / 'evl' / 'cdNodes.txt')
        self.nodes = nodes[np.argsort(nodes[:, 0]), 1:4]  # in units of b
        self.snapshots = sorted((int(k), v) for k, v in self.record['snapshots'].items())
        self.doses = [dose for _, dose in self.snapshots]
        self.fields = [self.read_fields(step) for step, _ in self.snapshots]
        self.mobile_names = ['Cv' if m == -1 else ('Ci' if m == 1 else 'C%di' % m) for m in self.ms.astype(int)]
        self.family_names = ['c', 'a1', 'a2', 'a3'][:self.nF] if self.nF == 4 else ['f%d' % k for k in range(self.nF)]
        self.lo = self.nodes.min(0)
        self.hi = self.nodes.max(0)
        self.geometry = derived.get('geometry_type', 'cubic')
        self.distance = self.depth(self.nodes) * self.b_SI * 1e9  # [nm] from the grain boundary
        self._weights = None

    def depth(self, points):
        """Distance of points (in b) to the surface of the grain; negative outside.
        cubic: the bounding box. hexagonal: a prism along z with a vertex on +x, inscribed in the bounding box."""
        points = np.atleast_2d(points)
        d = np.minimum(points - self.lo, self.hi - points)
        if self.geometry != 'hexagonal':
            return d.min(1)
        center = 0.5 * (self.lo + self.hi)
        apothem = 0.5 * (self.hi[1] - self.lo[1])
        angles = np.radians(30.0 + 60.0 * np.arange(6))
        normals = np.stack([np.cos(angles), np.sin(angles)], 1)
        side = apothem - ((points[:, :2] - center[:2]) @ normals.T).max(1)
        return np.minimum(side, d[:, 2])

    def volume_b3(self):
        e = self.hi - self.lo
        if self.geometry == 'hexagonal':
            return 1.5 * np.sqrt(3.0) * (0.5 * e[0]) ** 2 * e[2]
        return float(np.prod(e))

    def outline(self):
        """The edges of the grain, as pairs of points in nm."""
        nm = self.b_SI * 1e9
        lo, hi = self.lo * nm, self.hi * nm
        if self.geometry == 'hexagonal':
            center = 0.5 * (lo + hi)
            R = 0.5 * (hi[0] - lo[0])
            ring = [np.array([center[0] + R * np.cos(a), center[1] + R * np.sin(a)]) for a in np.radians(60.0 * np.arange(6))]
            edges = []
            for i in range(6):
                a, b = ring[i], ring[(i + 1) % 6]
                edges += [(np.array([*a, lo[2]]), np.array([*b, lo[2]])), (np.array([*a, hi[2]]), np.array([*b, hi[2]])), (np.array([*a, lo[2]]), np.array([*a, hi[2]]))]
            return edges
        edges = []
        for axis in range(3):
            others = [i for i in range(3) if i != axis]
            for u in (lo[others[0]], hi[others[0]]):
                for v in (lo[others[1]], hi[others[1]]):
                    p0, p1 = np.array(lo, float), np.array(hi, float)
                    p0[others[0]] = p1[others[0]] = u
                    p0[others[1]] = p1[others[1]] = v
                    edges.append((p0, p1))
        return edges

    def read_fields(self, step):
        lines = (self.dir / 'evl' / ('evl_%d.txt' % step)).read_text().split('\n')
        n = int(lines[9])
        rows = [l for l in lines[10:] if l.strip()][-n:]
        return np.array([[float(x) for x in r.split()] for r in rows])

    def weights(self, samples=2000000, seed=0):
        """The share of the volume nearest to each node, by Monte Carlo sampling of the box."""
        if self._weights is None:
            from scipy.spatial import cKDTree
            rng = np.random.default_rng(seed)
            tree = cKDTree(self.nodes)
            counts = np.zeros(len(self.nodes))
            for _ in range(max(samples // 500000, 1)):
                points = self.lo + (self.hi - self.lo) * rng.random((500000, 3))
                points = points[self.depth(points) >= 0.0]  # inside the grain
                counts += np.bincount(tree.query(points)[1], minlength=len(self.nodes))
            self._weights = counts / counts.sum()
        return self._weights

    def quantity(self, name, index):
        """A field in the units of the figures. name: a mobile species, N_<family> or C_<family>."""
        F = self.fields[index]
        if name in self.mobile_names:
            m = self.mobile_names.index(name)
            return F[:, m] / abs(self.ms[m])
        kind, family = name.split('_')
        k = self.family_names.index(family)
        if kind == 'N':
            return F[:, self.mSize + k] / self.b_SI ** 3
        if kind == 'C':
            return F[:, self.mSize + self.nF + k]
        if kind == 'm':  # defects per loop
            N = F[:, self.mSize + k] * self.omega
            return np.where(N > 1e-30, F[:, self.mSize + self.nF + k] / np.maximum(N, 1e-300), 0.0)
        raise KeyError(name)

    def quantity_names(self):
        return self.mobile_names + ['N_' + f for f in self.family_names] + ['C_' + f for f in self.family_names]

    def interior(self, fraction=0.5):
        """Nodes farther from the grain boundary than the given fraction of the half size."""
        return self.distance > fraction * self.distance.max()


LABELS = {'Cv': r'$C_v$', 'Ci': r'$C_i$', 'C2i': r'$C_{2i}$', 'C3i': r'$C_{3i}$'}


def label(name):
    if name in LABELS:
        return LABELS[name]
    kind, family = name.split('_')
    family = r'\langle c\rangle' if family == 'c' else r'\langle a\rangle_%s' % family[1:]
    return {'N': r'$N_{%s}$ [m$^{-3}$]', 'C': r'$c_{%s}$', 'm': r'$m_{%s}$'}[kind] % family


def volume_means(run, samples=2000000):
    w = run.weights(samples)
    names = run.quantity_names()
    table = np.array([[np.sum(w * run.quantity(name, i)) for name in names] for i in range(len(run.doses))])
    with open(run.dir / 'volume_means.csv', 'w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['dose_dpa'] + names)
        for dose, row in zip(run.doses, table):
            writer.writerow([repr(dose)] + [repr(float(x)) for x in row])
    return names, table


def figure_volume_means(run, names, table, plt):
    out = run.dir / 'figures'
    out.mkdir(exist_ok=True)
    doses = np.array(run.doses)
    x = np.where(doses > 0, doses, doses[doses > 0].min() / 3.0)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    groups = [(run.mobile_names, 'volume mean (clusters per atom)'), (['N_' + f for f in run.family_names], r'volume mean number density [m$^{-3}$]'),
              (['C_' + f for f in run.family_names], 'volume mean content (defects per atom)')]
    for ax, (group, ylabel) in zip(axes, groups):
        for name in group:
            y = table[:, names.index(name)]
            ax.loglog(x[y > 0], y[y > 0], 'o-', label=label(name).split(' [')[0])
        ax.set_xlabel('dose (dpa)')
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / 'volume_means.png', dpi=150)
    plt.close(fig)


def profile(run, values, bins=40):
    """Mean of a field against the distance from the grain boundary."""
    edges = np.linspace(0.0, run.distance.max() * (1 + 1e-9), bins + 1)
    which = np.digitize(run.distance, edges) - 1
    x, y = [], []
    for b in range(bins):
        selection = which == b
        if selection.any():
            x.append(run.distance[selection].mean())
            y.append(values[selection].mean())
    return np.array(x), np.array(y)


def figures_grain_boundary(run, plt):
    out = run.dir / 'gb'
    out.mkdir(exist_ok=True)
    cmap = plt.get_cmap('coolwarm')
    indices = [i for i, d in enumerate(run.doses) if d > 0]
    names = run.mobile_names + ['N_' + run.family_names[0], 'C_' + run.family_names[0]]
    if run.nF > 1:
        names += ['N_' + run.family_names[1], 'C_' + run.family_names[1]]
    for name in names:
        fig, ax = plt.subplots(figsize=(6.5, 4.4))
        for j, i in enumerate(indices):
            x, y = profile(run, run.quantity(name, i))
            ok = y > 0
            ax.semilogy(x[ok], y[ok], lw=2, color=cmap(j / max(len(indices) - 1, 1)), label='%g dpa' % run.doses[i])
        ax.set_xlabel('distance from grain boundary, $x$ [nm]')
        ax.set_ylabel(label(name))
        ax.set_title(label(name).split(' [')[0] + ' vs distance from the GB')
        ax.grid(True, which='both', alpha=0.3)
        ax.legend(frameon=False, fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(out / ('gb_%s.png' % name), dpi=150)
        plt.close(fig)
    # the mobile species at the last dose, in one figure
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    for name in run.mobile_names:
        x, y = profile(run, run.quantity(name, len(run.doses) - 1))
        ax.semilogy(x[y > 0], y[y > 0], lw=2, label=label(name))
    ax.set_xlabel('distance from grain boundary, $x$ [nm]')
    ax.set_ylabel('clusters per atom')
    ax.set_title('mobile species at %g dpa' % run.doses[-1])
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    (run.dir / 'figures').mkdir(exist_ok=True)
    fig.savefig(run.dir / 'figures' / 'mobile_profile.png', dpi=150)
    plt.close(fig)


def figures_cuts(run, plt, doses=None, resolution=90):
    """Each field on two orthogonal mid-planes of the grain, on a logarithmic color scale."""
    from matplotlib import cm, colors
    from scipy.interpolate import LinearNDInterpolator
    out = run.dir / '3d'
    out.mkdir(exist_ok=True)
    nm = run.b_SI * 1e9
    lo, hi = run.lo * nm, run.hi * nm
    mid = 0.5 * (lo + hi)
    points = run.nodes * nm
    grid = np.linspace(0.0, 1.0, resolution)
    A, B = np.meshgrid(grid, grid)
    planes = []
    for axis in (2, 0):  # horizontal mid-plane, then a vertical one
        other = [a for a in range(3) if a != axis]
        P = np.zeros((resolution, resolution, 3))
        P[..., axis] = mid[axis]
        P[..., other[0]] = lo[other[0]] + A * (hi[other[0]] - lo[other[0]])
        P[..., other[1]] = lo[other[1]] + B * (hi[other[1]] - lo[other[1]])
        planes.append(P)
    indices = [i for i, d in enumerate(run.doses) if d > 0 and (doses is None or any(abs(d - x) <= 1e-9 * max(d, x) for x in doses))]
    for name in run.quantity_names():
        if name.split('_')[-1] in run.family_names[2:]:
            continue  # a2 and a3 are equal to a1 without applied stress
        for i in indices:
            values = run.quantity(name, i)
            positive = values[values > 0]
            if positive.size == 0:
                continue
            vmax = positive.max()
            vmin = max(positive.min(), vmax * 1e-4)
            norm = colors.LogNorm(vmin=vmin, vmax=vmax) if vmax > vmin * 1.0001 else colors.Normalize(vmin=0.9 * vmin, vmax=1.1 * vmax)
            interpolator = LinearNDInterpolator(points, np.log(np.maximum(values, vmin)))
            fig = plt.figure(figsize=(6.4, 5.0))
            ax = fig.add_subplot(111, projection='3d')
            for P in planes:
                V = np.exp(interpolator(P.reshape(-1, 3)).reshape(resolution, resolution))
                face = cm.jet(norm(np.nan_to_num(V, nan=vmin)))
                outside = np.isnan(V) | (run.depth(P.reshape(-1, 3) / nm).reshape(resolution, resolution) < 0.0)
                face[outside, 3] = 0.0  # nothing is drawn outside the grain
                ax.plot_surface(P[..., 0], P[..., 1], P[..., 2], facecolors=face, rstride=1, cstride=1, shade=False, antialiased=False, linewidth=0)
            for p0, p1 in run.outline():
                ax.plot(*zip(p0, p1), color='0.35', lw=0.8)
            ax.set_axis_off()
            ax.set_box_aspect(tuple(hi - lo))
            ax.view_init(elev=22, azim=-58)
            ax.set_title('%s   %g dpa' % (label(name).split(' [')[0], run.doses[i]), fontsize=14)
            ax.text2D(0.02, 0.05, 'axes: [2-1-10], [01-10], [0001]', transform=ax.transAxes, fontsize=7)
            bar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.jet), ax=ax, shrink=0.55, pad=0.08)
            bar.ax.set_title(label(name) + ('\nlog' if isinstance(norm, colors.LogNorm) else ''), fontsize=8)
            fig.savefig(out / ('%s_%sdpa.png' % (name, format_dose(run.doses[i]))), dpi=130)
            plt.close(fig)


def checks(run, names, table):
    """Verification of the run against values known independently. Returns rows (check, expected, result, verdict)."""
    rows = []
    material = next(p for p in (run.dir / 'inputFiles').iterdir() if p.name == run.provenance['configuration']['MATERIAL']['material_file']).read_text()

    def value(key):
        import re
        return [float(x) for x in re.search(r'^' + re.escape(key) + r' *=([^;]*);', material, re.M).group(1).split()]

    last = len(run.doses) - 1
    F = run.fields[last]
    rows.append(('no negative field', 'minimum >= 0', 'minimum %.2e' % F.min(), F.min() >= 0.0))
    if run.nF == 4:
        symmetry = max(np.abs(F[:, run.mSize + k] / np.maximum(F[:, run.mSize + 1], 1e-300) - 1.0).max() for k in (2, 3))
        rows.append(('the three <a> families are equal without stress', 'relative difference < 1e-8', '%.1e' % symmetry, symmetry < 1e-8))
    try:
        T = run.provenance['configuration']['MATERIAL']['temperature_K']
        time_unit = run.provenance['derived']['time_unit_s']
        tau = value('tau0_vLoop_SI')[0] * np.exp(value('Ea_vLoop_eV')[0] / (kB_eV * T))  # [s]
        fractions = np.array(value('mobileSpeciesCascadeFractions'))
        epsV = 1.0 - fractions[run.ms < 0].sum()
        G = run.provenance['configuration']['MATERIAL']['dose_rate_dpa_s'] * value('mobileSpeciesSurvivingEfficiency')[0]
        expected = G * epsV / value('nNuc')[0] * tau / run.omega_SI  # [m^-3]
        dose_needed = 20.0 * tau * G
        if run.doses[last] > dose_needed and run.provenance['configuration']['MARCH']['immobile_integrator'] == 'cvode':
            interior = run.interior(0.8)
            measured = run.quantity('N_c', last)[interior]
            # without coalescence the density saturates at production times lifetime; coalescence lowers it
            rows.append(('<c> density in the interior against production times lifetime', '<= %.3e m^-3, equal while the loops do not coalesce' % expected,
                         '%.3e m^-3 at %g dpa' % (measured.mean(), run.doses[last]), measured.mean() <= expected * 1.001))
            first = next(i for i, d in enumerate(run.doses) if d > dose_needed)
            early = run.quantity('N_c', first)[interior].mean()
            rows.append(('<c> density at the first dose beyond 20 lifetimes', '%.3e m^-3 within 5 %%' % expected, '%.3e m^-3 at %g dpa' % (early, run.doses[first]), abs(early / expected - 1.0) < 0.05))
    except (AttributeError, StopIteration, ValueError):
        pass
    cv = run.quantity('Cv', last)
    boundary = run.distance < 1e-6
    T = run.provenance['configuration']['MATERIAL']['temperature_K']
    equilibrium = np.exp(-value('mobileSpeciesEnergyFormation_eV')[0] / (kB_eV * T))
    error = np.abs(cv[boundary] / equilibrium - 1.0).max()
    rows.append(('vacancies at the grain boundary at thermal equilibrium', 'exp(-Ef/kT) = %.3e' % equilibrium, 'largest relative difference %.1e' % error, error < 1e-4))
    return rows


def reference_comparison(run, names, table, plt):
    """Volume means against those of the ENNDS run of the same geometry, where one is stored. Returns table rows."""
    reference_dir = Path(__file__).resolve().parent / 'reference' / ('ENNDS_%s_10dpa' % run.provenance['configuration']['CASE']['geometry'])
    if not (reference_dir / 'volume_means.csv').is_file():
        return interior_comparison(run, plt)
    meta = json.loads((reference_dir / 'reference.json').read_text())
    with open(reference_dir / 'volume_means.csv') as stream:
        data = list(csv.DictReader(stream))
    ref_dose = np.array([float(r['dose_dpa']) for r in data])
    omega_ref = meta['atomic_volume_m3']

    def ref(column, per_volume=True):
        return np.array([float(r[column]) for r in data]) / (omega_ref if per_volume else 1.0)

    pairs = [  # (name here, label, reference values in the same units)
        ('Cv', 'vacancies [m^-3]', ref('Cv')), ('Ci', 'interstitials [m^-3]', ref('Ci')),
        ('N_c', '<c> vacancy loops [m^-3]', ref('n_c_f') + ref('n_c_p') + ref('n_c_0')),
        ('N_a1', '<a>1 interstitial loops [m^-3]', ref('n_a1')),
        ('C_c', 'content of the <c> loops [per atom]', ref('c_c_f', False) + ref('c_c_p', False) + ref('c_c_0', False)),
        ('C_a1', 'content of the <a>1 loops [per atom]', ref('c_a1', False)),
    ]
    rows = []
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5))
    for ax, (name, text, reference) in zip(axes.ravel(), pairs):
        here = table[:, names.index(name)].copy()
        if name in run.mobile_names:
            here = here / run.omega_SI  # clusters per atom -> per m^3
        d = np.array(run.doses)
        ax.loglog(d[d > 0], here[d > 0], 'o-', label='MoDELib (this run)')
        ax.loglog(ref_dose[ref_dose > 0], reference[ref_dose > 0], 's--', label='ENNDS reference')
        ax.set_xlabel('dose (dpa)')
        ax.set_ylabel(text)
        ax.grid(True, which='both', alpha=0.3)
        ax.legend(frameon=False, fontsize=8)
        for dose, v in zip(d, here):
            match = np.where(np.abs(ref_dose - dose) <= 1e-9 * max(dose, 1e-30))[0]
            if dose > 0 and match.size:
                rows.append((text, dose, v, reference[match[0]], v / reference[match[0]] if reference[match[0]] > 0 else float('nan')))
    fig.suptitle('Volume means: MoDELib (D1/M1 model, 4 loop families) and the ENNDS reference (later model, 9 families)', fontsize=10)
    fig.tight_layout()
    (run.dir / 'figures').mkdir(exist_ok=True)
    fig.savefig(run.dir / 'figures' / 'comparison_with_ENNDS.png', dpi=150)
    plt.close(fig)
    return rows, meta


def interior_comparison(run, plt):
    """Interior values against those of a reference run that stores them (reference/<name>/reference.json with 'interior')."""
    geometry = run.provenance['configuration']['CASE']['geometry']
    folder = next((p for p in (Path(__file__).resolve().parent / 'reference').iterdir()
                   if p.is_dir() and geometry.replace('_', '') in p.name.replace('_', '') and (p / 'reference.json').is_file()), None)
    if folder is None:
        return None
    meta = json.loads((folder / 'reference.json').read_text())
    if 'interior' not in meta:
        return None
    import postprocessing_extended
    kin = postprocessing_extended.Kinetics(run)
    quantities = [('N_c', '<c> loop density [m^-3]'), ('N_a1', '<a>1 loop density [m^-3]'), ('d_c', '<c> loop diameter [nm]'), ('d_a1', '<a>1 loop diameter [nm]'),
                  ('Cv', 'vacancies [m^-3]'), ('Ci', 'interstitials [m^-3]'), ('C2i', 'di-interstitials [m^-3]'), ('C3i', 'tri-interstitials [m^-3]')]
    here = {q: [] for q, _ in quantities}
    doses = []
    for i, dose in enumerate(run.doses):
        key = '%g' % dose
        if key not in meta['interior']:
            continue
        doses.append(dose)
        # the interior of the reference run: the tenth of the nodes with the highest vacancy concentration
        vacancies = run.quantity(run.mobile_names[0], i)
        interior = vacancies >= np.percentile(vacancies, 90.0)
        n, c, size, R, d = kin.sizes(i)
        for q, _ in quantities:
            if q.startswith('N_'):
                here[q].append(float(run.quantity(q, i)[interior].mean()))
            elif q.startswith('d_'):
                k = run.family_names.index(q[2:])
                N = run.quantity('N_' + q[2:], i)[interior]
                here[q].append(float(np.sum(N * d[interior, k]) / max(N.sum(), 1e-300)))
            else:
                here[q].append(float(run.quantity(q, i)[interior].mean() / run.omega_SI))
    if not doses:
        return None
    rows = []
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5))
    for ax, (q, text) in zip(axes.ravel(), quantities):
        ref = [meta['interior']['%g' % dose][q] for dose in doses]
        ax.plot(doses, here[q], 'o-', label='MoDELib (this run)')
        ax.plot(doses, ref, 's--', label='DisloCluster reference')
        ax.set_xlabel('dose (dpa)'); ax.set_ylabel(text); ax.grid(True, alpha=0.3); ax.legend(frameon=False, fontsize=8)
        if not q.startswith('d_'):
            ax.set_yscale('log')
        rows += [(text, dose, a, b, a / b if b else float('nan')) for dose, a, b in zip(doses, here[q], ref)]
    fig.suptitle('Interior of the grain: this run and the reference run of the D1/M1 report', fontsize=10)
    fig.tight_layout()
    (run.dir / 'figures').mkdir(exist_ok=True)
    fig.savefig(run.dir / 'figures' / 'comparison_with_reference.png', dpi=150)
    plt.close(fig)
    return rows, meta


def gather(run_dir, cuts=True, figure_doses=None, samples=2000000, progress=print):
    """Computes the volume means, draws the figures and writes summary.json and report.md."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    run = Run(run_dir)
    names, table = volume_means(run, samples)
    progress('volume means written')
    figure_volume_means(run, names, table, plt)
    figures_grain_boundary(run, plt)
    progress('profiles written')
    if cuts:
        figures_cuts(run, plt, figure_doses)
        progress('mid-plane cuts written')
    extended, manifest = [], None
    try:
        import postprocessing_extended
        extended, manifest = postprocessing_extended.gather_extended(run, plt, samples, progress, population_dose=run.provenance['configuration']['OUTPUT'].get('population_dose'))
    except Exception as error:  # the report is written in any case
        import traceback
        (run.dir / 'logs').mkdir(exist_ok=True)
        (run.dir / 'logs' / 'postprocessing_extended_error.txt').write_text(traceback.format_exc(), encoding='utf-8')
        progress('extended figures not completed: %s (see logs/postprocessing_extended_error.txt)' % error)
    verification = checks(run, names, table)
    comparison = reference_comparison(run, names, table, plt)

    interior = run.interior(0.8)
    last = len(run.doses) - 1
    summary = {
        'run': run.dir.name, 'doses': run.doses, 'nodes': int(len(run.nodes)), 'wall_s': run.record.get('wall_s'),
        'volume_means_at_last_dose': {name: float(table[last, names.index(name)]) for name in names},
        'interior_means_at_last_dose': {name: float(run.quantity(name, last)[interior].mean()) for name in names},
        'checks': [{'check': c, 'expected': e, 'result': r, 'passed': bool(ok)} for c, e, r, ok in verification],
    }
    (run.dir / 'summary.json').write_text(json.dumps(summary, indent=1), encoding='utf-8')

    config = run.provenance['configuration']
    lines = ['# %s' % run.dir.name, '',
             '- material: `%s` at %g K, %g dpa/s' % (config['MATERIAL']['material_file'], config['MATERIAL']['temperature_K'], config['MATERIAL']['dose_rate_dpa_s']),
             '- geometry: `%s` (%s), %s nm, every face a grain boundary' % (config['CASE']['geometry'], run.geometry, ' x '.join('%.0f' % (e * run.b_SI * 1e9) for e in (run.hi - run.lo))),
             '- integration of the immobile species: **%s**' % config['MARCH']['immobile_integrator'],
             '- snapshots: %d, %g .. %g dpa' % (len(run.doses), run.doses[0], run.doses[-1]),
             '- finite-element nodes: %d (%d mobile and %d immobile unknowns)' % (len(run.nodes), len(run.nodes) * run.mSize, len(run.nodes) * 2 * run.nF),
             '- status: %s; wall time %.0f s' % (run.record.get('status'), run.record.get('wall_s') or 0.0), '',
             '## Verification', '', '| Check | Expected | Result | |', '|---|---|---|---|']
    lines += ['| %s | %s | %s | %s |' % (c, e, r, 'PASS' if ok else '**FAIL**') for c, e, r, ok in verification]
    lines += ['', '## Intervals', '', '| dose from | dose to | steps | wall |', '|---:|---:|---:|---:|']
    lines += ['| %g | %g | %d | %.0f s |' % (i['dose_from'], i['dose_to'], i['steps'], i['wall_s']) for i in run.record.get('intervals', [])]
    lines += ['', '## Values in the interior of the grain', '',
              'Means over the nodes farther than 80 % of the half size from the grain boundary.', '',
              '| dose [dpa] | ' + ' | '.join(names) + ' |', '|---:|' + '---:|' * len(names)]
    for i, dose in enumerate(run.doses):
        lines.append('| %g | ' % dose + ' | '.join('%.3e' % run.quantity(name, i)[interior].mean() for name in names) + ' |')
    lines += ['', 'Mobile species in clusters per atom, number densities N in m^-3, contents C in defects per atom.', '']
    if comparison:
        rows, meta = comparison
        lines += ['## Comparison with the reference run', '',
                  '%s of this run and of the reference run of the same geometry, temperature and dose rate: %s.' % ('Interior values' if 'interior' in meta else 'Volume means', meta['run']),
                  meta['note'], '',
                  '| quantity | dose [dpa] | MoDELib | %s | ratio |' % ('DisloCluster' if 'interior' in meta else 'ENNDS'), '|---|---:|---:|---:|---:|']
        lines += ['| %s | %g | %.3e | %.3e | %.2f |' % row for row in rows]
        lines += ['', 'Figure: `figures/%s`.' % ('comparison_with_reference.png' if 'interior' in meta else 'comparison_with_ENNDS.png'), '']
    lines += ['## Files', '',
              '| Path | Contents |', '|---|---|',
              '| `provenance.json`, `provenance.md` | settings, input files and executables with their SHA-256, commit, machine, times |',
              '| `config.json`, `inputFiles/` | the configuration and the input files as run |',
              '| `evl/` | the fields at the snapshot doses (`evl_<n>.txt`) and the positions of the finite-element nodes (`cdNodes.txt`) |',
              '| `F/` | the average plastic distortion (growth strain) at every step |',
              '| `logs/` | the output of the executables, one file per dose interval |',
              '| `volume_means.csv`, `summary.json` | volume means at the snapshot doses; the summary and the checks |',
              '| `figures/` | volume means against dose, mobile profile%s |' % ((', comparison with DisloCluster' if 'interior' in comparison[1] else ', comparison with ENNDS') if comparison else ''),
              '| `gb/` | each field against the distance from the grain boundary, one curve per dose |',
              '| `3d/` | each field on two mid-planes of the grain, at each snapshot dose; `loops_<family>_<dose>dpa.png`: the discrete loops of a family |']
    extra = {'volume_average/': 'volume means against dose: point defects, loop densities (also by region) and diameters (also against the measured ranges of the D1/M1 report), growth strain and its rate, hardening estimate, retained defects and their shares, fluxes to the sinks, balance of v and i; `additional_analysis.csv`, `growth_strain.csv`',
             'boundary_flux.md': 'balance of each mobile species at each dose, and the share that leaves through the grain boundary',
             'gb/denuded_zone.png': 'number density and mean diameter of the loops against the distance to the nearest face, at the last dose',
             'size_spectrum/': 'distribution of the loop diameters over the grain, its interior and its boundary shell, at the last dose',
             'discrete_loops/': 'discrete loops drawn from the fields at up to four doses: `loops_<dose>dpa.csv`, `manifest.json`, three-dimensional views of all the loops and of each family',
             'tem_slices/': 'the discrete loops of a foil projected along [0001], [01-10] and [2-1-10]; one image at the last dose and a montage over the doses',
             'movies/': 'the mid-plane figures of four fields as animated images, one frame per snapshot dose'}
    lines += ['| `%s` | %s |' % (key, extra[key]) for key in extended if key in extra]
    lines += ['', 'Not produced, because the model of this run does not carry the quantities: the second moments of the size distributions (`q_*` figures of DisloCluster) and the accumulated conservation error. Each family has one loop size at each point; the size spectra and the discrete loops show how that size varies over the grain.']
    (run.dir / 'report.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    progress('report written: %s' % (run.dir / 'report.md'))
    return summary
