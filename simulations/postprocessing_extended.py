"""The figures of a run beyond those of postprocessing.py, after the output folders of DisloCluster:

  gb/denuded_zone.png          number density and mean diameter of the loops against the distance from the grain boundary
  volume_average/              volume means against dose: point defects, loop densities and sizes, growth strain and its rate,
                               hardening, stored defects, fluxes to the sinks, balance of each point defect
  size_spectrum/               distribution of the loop diameters over the grain, its interior and its boundary shell
  discrete_loops/              discrete loops drawn from the fields: tables, a manifest and three-dimensional views
  3d/loops_<family>_<dose>.png the discrete loops of one family in the grain
  tem_slices/                  the discrete loops of a foil projected along three zone axes
  movies/                      each field on the mid-planes of the grain, one frame per snapshot dose
  boundary_flux.md             the balance of the mobile species and what leaves through the grain boundary

Everything is computed from the fields of the run and from its material file. The model carries the number
density and the content of each family, hence ONE loop size per family at each point: the size distributions
and the discrete loops reflect the variation of that size over the grain, not a distribution at a point.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'testsPy' / 'clusterDynamics'))

FAMILY_COLORS = {'c': '#1f4fbf', 'a1': '#c62828', 'a2': '#2e7d32', 'a3': '#e6b800'}
FAMILY_LABELS = {'c': r'$\langle c\rangle$ basal, vacancy', 'a1': r'$\langle a\rangle_1$ prismatic, interstitial',
                 'a2': r'$\langle a\rangle_2$ prismatic, interstitial', 'a3': r'$\langle a\rangle_3$ prismatic, interstitial'}
# Dispersed-barrier estimate of the hardening: nominal Taylor factor and obstacle strengths
TAYLOR_FACTOR = 3.06
OBSTACLE_STRENGTH = {'c': 0.25, 'a': 0.25}
# Ranges quoted in the D1/M1 report (section 4.1) for irradiated zirconium near 573 K
MEASURED = {'N_a': (1e22, 3e22), 'd_a': (9.0, 10.0), 'd_c': (95.0, 150.0), 'N_c': None}


def format_dose(dose):
    return ('%g' % dose).replace('.', 'p')


class Kinetics:
    """The rates of the D1/M1 model evaluated on the fields of a run (vectorized over the nodes)."""

    def __init__(self, run):
        import referenceRates
        self.run = run
        config = run.provenance['configuration']
        material = run.dir / 'inputFiles' / config['MATERIAL']['material_file']
        self.model = referenceRates.Model(str(material), config['MATERIAL']['temperature_K'])
        self.time_unit = run.provenance['derived']['time_unit_s']
        m = self.model
        text = material.read_text()
        self.fractions = np.array(referenceRates.value(str(material), 'mobileSpeciesCascadeFractions'))
        self.network = np.array(referenceRates.value(str(material), 'otherSinks_SI')) * m.b_SI ** 2
        nM = len(m.ms)
        r = np.where(np.abs(m.ms) == 1, (3.0 * m.omega / 4.0 / np.pi) ** (1.0 / 3.0), np.sqrt(np.abs(m.ms) * m.omega / np.pi))
        self.reactions = []  # (a, b, K) of every active reaction between mobile species
        for a, b, prefactor in referenceRates.matrix(str(material), 'reactionPrefactorMap', nM * (nM + 1) // 2):
            if prefactor > 0:
                ia, ib = list(m.ms).index(a), list(m.ms).index(b)
                self.reactions.append((ia, ib, prefactor * 4.0 * np.pi * (r[ia] + r[ib]) * (m.Dbar[ia] + m.Dbar[ib]) / m.omega))

    def sizes(self, index):
        """Defects per loop, effective radius and loop diameter [nm] of each family at each node."""
        run, m = self.run, self.model
        F = run.fields[index]
        n = F[:, run.mSize:run.mSize + run.nF] * m.omega
        c = F[:, run.mSize + run.nF:]
        present = (n > 1e-30) & (c > 1e-30)
        size = np.where(present, c / np.where(present, n, 1.0), 0.0)
        R = np.zeros_like(size)
        d = np.zeros_like(size)
        for k in range(run.nF):
            S = m.sigmoid(k, size[:, k])
            loop = np.sqrt(size[:, k] * m.omega / (np.pi * m.bMag[k]))
            R[:, k] = (1.0 - S) * (size[:, k] * m.omega / np.sqrt(8.0)) ** (1.0 / 3.0) + S * loop
            d[:, k] = 2.0 * R[:, k] * m.b_SI * 1e9
        return n, c, size, R, d

    def absorption(self, index):
        """Gamma[node, family, species]: defects of each mobile species absorbed by each family, per atom and per second."""
        run, m = self.run, self.model
        n, c, size, R, d = self.sizes(index)
        cM = np.maximum(run.fields[index][:, :run.mSize], 0.0)
        gamma = np.zeros((len(n), run.nF, run.mSize))
        for k in range(run.nF):
            x = (R[:, k] / m.rmin[k] - 1.0) / 0.3 if m.rmin[k] > 0 else np.ones(len(n))
            gate = np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, x * x * (3.0 - 2.0 * x)))
            for s in range(run.mSize):
                g = 1.0 if m.pF[k] * m.pM[s] > 0 else gate
                gamma[:, k, s] = g * m.Dbar[s] * 2.0 * np.pi * m.Z[k, s] * R[:, k] * n[:, k] / m.omega * cM[:, s]
        return gamma / self.time_unit

    def balance(self, index, w):
        """Volume means of the channels of each mobile species, in defects per atom per second."""
        run, m = self.run, self.model
        cM = np.maximum(run.fields[index][:, :run.mSize], 0.0)
        C = cM / np.abs(m.ms)
        gamma = self.absorption(index)
        out = {}
        n, c, size, R, d = self.sizes(index)
        release = np.zeros(len(cM))
        for k in range(run.nF):
            if m.pF[k] < 0:
                release += c[:, k] / m.tau + m.Dbar[m.iv] * 2.0 * np.pi * m.Z[k, m.iv] * R[:, k] * n[:, k] / m.omega * np.exp(-m.Eb[m.iv] / m.kT)
        for s in range(run.mSize):
            reaction = np.zeros(len(cM))
            for ia, ib, K in self.reactions:
                if s in (ia, ib):
                    reaction -= abs(m.ms[s]) * K * C[:, ia] * C[:, ib]          # consumed
                elif abs(m.ms[ia] + m.ms[ib] - m.ms[s]) < 1e-9:
                    reaction += abs(m.ms[s]) * (0.5 if ia == ib else 1.0) * K * C[:, ia] * C[:, ib]  # produced
            channels = {
                'production': m.G / self.time_unit * self.fractions[s],
                'reactions': float(np.sum(w * reaction)) / self.time_unit,
                'loops': -float(np.sum(w * gamma[:, :, s].sum(1))),
                'network': -float(np.sum(w * m.Dbar[s] * self.network[s] * cM[:, s])) / self.time_unit,
                'released by vacancy loops': float(np.sum(w * release)) / self.time_unit if s == m.iv else 0.0,
            }
            # the mobile species are at steady state: what remains leaves through the grain boundary
            # (and, for the clusters, by dissociation, which is not itemized)
            channels['grain boundary'] = -sum(channels.values())
            out[run.mobile_names[s]] = channels
        return out


def growth_strain(run):
    """Dose and average plastic distortion (growth strain) at every step, from F/F_0.txt."""
    rows = {}
    for line in (run.dir / 'F' / 'F_0.txt').read_text().split('\n'):
        parts = line.split()
        if len(parts) >= 12:
            rows[int(float(parts[0]))] = [float(x) for x in parts[:12]]  # a step written twice keeps its last record
    data = np.array([rows[k] for k in sorted(rows)])
    config = run.provenance['configuration']
    dose = data[:, 1] * run.provenance['derived']['time_unit_s'] * config['MATERIAL']['dose_rate_dpa_s']
    return dose, data[:, 3], data[:, 7], data[:, 11]


def dose_axis(run):
    doses = np.array(run.doses)
    return doses, doses > 0


def figures_volume_average(run, kin, plt, w):
    out = run.dir / 'volume_average'
    out.mkdir(exist_ok=True)
    doses, positive = dose_axis(run)
    x = doses[positive]
    idx = np.where(positive)[0]
    m = kin.model
    families = run.family_names
    interior = run.interior(0.8)
    shell = run.distance < 0.1 * run.distance.max()
    mean = lambda v, mask=None: float(np.sum(w * v)) if mask is None else float(np.sum(w[mask] * v[mask]) / np.sum(w[mask]))
    table = {'dose_dpa': list(x)}

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(out / name, dpi=150)
        plt.close(fig)

    # point defects
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for name in run.mobile_names:
        y = np.array([mean(run.quantity(name, i)) for i in idx]) / run.omega_SI * 1e-6
        ax.loglog(x, y, 'o-', label=name)
        table[name + '_cm-3'] = list(y)
    for k, f in enumerate(families[:2]):
        y = np.array([mean(run.quantity('C_' + f, i)) for i in idx]) / run.omega_SI * 1e-6
        ax.loglog(x, y, 's--', label='defects stored in %s loops' % ('<c>' if f == 'c' else '<a>1'))
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel(r'concentration [cm$^{-3}$]'); ax.set_title('point defects and stored content, volume means')
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'point_defects.png')

    # loop densities, sizes
    N = {f: np.array([mean(run.quantity('N_' + f, i)) for i in idx]) for f in families}
    Nint = {f: np.array([mean(run.quantity('N_' + f, i), interior) for i in idx]) for f in families}
    Nshell = {f: np.array([mean(run.quantity('N_' + f, i), shell) for i in idx]) for f in families}
    sizes = [kin.sizes(i) for i in idx]
    diameter, diameter_int = {}, {}
    for k, f in enumerate(families):
        # number-weighted mean diameter
        diameter[f] = np.array([np.sum(w * s[0][:, k] * s[4][:, k]) / max(np.sum(w * s[0][:, k]), 1e-300) for s in sizes])
        diameter_int[f] = np.array([np.sum((w * s[0][:, k] * s[4][:, k])[interior]) / max(np.sum((w * s[0][:, k])[interior]), 1e-300) for s in sizes])
        table['N_%s_m-3' % f] = list(N[f]); table['N_%s_interior_m-3' % f] = list(Nint[f]); table['d_%s_nm' % f] = list(diameter[f]); table['d_%s_interior_nm' % f] = list(diameter_int[f])
    Na = sum(N[f] for f in families if f != 'c')

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.loglog(x, N['c'], 'o-', color=FAMILY_COLORS['c'], label=r'$\langle c\rangle$ vacancy loops')
    ax.loglog(x, Na, 'o-', color=FAMILY_COLORS['a1'], label=r'$\langle a\rangle$ interstitial loops, three families')
    ax.axhspan(*MEASURED['N_a'], color=FAMILY_COLORS['a1'], alpha=0.12, label=r'measured $\langle a\rangle$ range (D1/M1 report)')
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel(r'number density [m$^{-3}$]'); ax.set_title('loop number density, volume mean')
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'loop_density.png')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, f, title in zip(axes, ('c', families[1] if len(families) > 1 else 'c'), (r'$\langle c\rangle$', r'$\langle a\rangle_1$')):
        ax.loglog(x, N[f], 'o-', color='k', label='grain')
        ax.loglog(x, Nint[f], 's--', color='#c62828', label='interior')
        ax.loglog(x, Nshell[f], '^:', color='#00838f', label='shell within 10 % of the half size of the boundary')
        ax.set_xlabel('dose (dpa)'); ax.set_ylabel(r'number density [m$^{-3}$]'); ax.set_title(title + ' loops by region')
        ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'loop_density_analysis.png')

    for name, with_measured in (('loop_sizes.png', False), ('loop_sizes_vs_experiment.png', True)):
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        ax.semilogx(x, diameter['c'], 'o-', color=FAMILY_COLORS['c'], label=r'$\langle c\rangle$, grain')
        ax.semilogx(x, diameter_int['c'], 'o--', color=FAMILY_COLORS['c'], alpha=0.6, label=r'$\langle c\rangle$, interior')
        if len(families) > 1:
            ax.semilogx(x, diameter[families[1]], 's-', color=FAMILY_COLORS['a1'], label=r'$\langle a\rangle$, grain')
            ax.semilogx(x, diameter_int[families[1]], 's--', color=FAMILY_COLORS['a1'], alpha=0.6, label=r'$\langle a\rangle$, interior')
        if with_measured:
            ax.axhspan(*MEASURED['d_c'], color=FAMILY_COLORS['c'], alpha=0.12, label=r'measured $\langle c\rangle$ range, 573-583 K (D1/M1 report)')
            ax.axhspan(*MEASURED['d_a'], color=FAMILY_COLORS['a1'], alpha=0.15, label=r'measured $\langle a\rangle$ saturation (D1/M1 report)')
        ax.set_xlabel('dose (dpa)'); ax.set_ylabel('number-weighted mean diameter [nm]'); ax.set_title('loop diameter')
        ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
        save(fig, name)

    # growth strain and its rate
    dose, e11, e22, e33 = growth_strain(run)
    ok = dose > 0
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for y, lab, col in ((e11, r'$\beta^P_{11}$  [2$\bar{1}\bar{1}$0]', '#c62828'), (e22, r'$\beta^P_{22}$  [01$\bar{1}$0]', '#2e7d32'), (e33, r'$\beta^P_{33}$  [0001]', '#1f4fbf')):
        ax.semilogx(dose[ok], 100.0 * y[ok], 'o-', ms=3, color=col, label=lab)
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel('growth strain [%]'); ax.set_title('irradiation growth: average plastic distortion of the grain')
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'irradiation_growth.png')
    table_growth = {'dose_dpa': list(dose), 'betaP_11': list(e11), 'betaP_22': list(e22), 'betaP_33': list(e33)}

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    if ok.sum() > 2:
        dm = 0.5 * (dose[ok][1:] + dose[ok][:-1])
        for y, lab, col in ((e11, 'a axis', '#c62828'), (e33, 'c axis', '#1f4fbf')):
            rate = np.diff(y[ok]) / np.diff(dose[ok])
            ax.loglog(dm[rate != 0], np.abs(rate[rate != 0]), 'o-', ms=3, color=col, label='|d strain / d dose|, ' + lab)
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel('strain rate [1/dpa]'); ax.set_title('growth strain rate')
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'strain_rates.png')

    # hardening, dispersed-barrier estimate
    mu = m.mu_SI
    hard = {}
    for k, f in enumerate(families):
        alpha = OBSTACLE_STRENGTH['c' if f == 'c' else 'a']
        hard[f] = np.array([TAYLOR_FACTOR * alpha * mu * m.b_SI * np.sqrt(max(np.sum(w * run.quantity('N_' + f, i) * s[4][:, k] * 1e-9), 0.0)) for i, s in zip(idx, sizes)]) / 1e6
    total = np.sqrt(sum(h ** 2 for h in hard.values()))
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.loglog(x, np.maximum(hard['c'], 1e-3), 'o-', color=FAMILY_COLORS['c'], label=r'$\langle c\rangle$ loops')
    ax.loglog(x, np.maximum(np.sqrt(sum(hard[f] ** 2 for f in families if f != 'c')), 1e-3), 's-', color=FAMILY_COLORS['a1'], label=r'$\langle a\rangle$ loops')
    ax.loglog(x, np.maximum(total, 1e-3), 'k-', lw=2, label='root sum of squares')
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel(r'$\Delta\sigma_y$ [MPa]')
    ax.set_title(r'hardening estimate, $\Delta\sigma = M\alpha\mu b\sqrt{N d}$  (M=%.2f, $\alpha$=%.2f)' % (TAYLOR_FACTOR, OBSTACLE_STRENGTH['a']), fontsize=9)
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'radiation_hardening.png')
    table['hardening_MPa'] = list(total)

    # stored defects and their shares
    stored = {f: np.array([mean(run.quantity('C_' + f, i)) for i in idx]) for f in families}
    mobile = {name: np.array([mean(run.quantity(name, i)) * abs(run.ms[j]) for i in idx]) for j, name in enumerate(run.mobile_names)}
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    vac = stored['c'] + sum(mobile[n] for j, n in enumerate(run.mobile_names) if run.ms[j] < 0)
    inter = sum(stored[f] for f in families if f != 'c') + sum(mobile[n] for j, n in enumerate(run.mobile_names) if run.ms[j] > 0)
    ax.loglog(x, vac, 'o-', color=FAMILY_COLORS['c'], label='vacancies: mobile and in loops')
    ax.loglog(x, inter, 's-', color=FAMILY_COLORS['a1'], label='interstitials: mobile and in loops')
    ax.loglog(x, x, 'k:', label='dose (all displaced atoms)')
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel('defects per atom'); ax.set_title('defects retained in the grain')
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'total_defect_evolution.png')

    for name, sign, title in (('vacancy_fractions.png', -1, 'vacancies'), ('interstitial_fractions.png', 1, 'interstitials')):
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        parts = {}
        for j, mn in enumerate(run.mobile_names):
            if np.sign(run.ms[j]) == sign:
                parts['mobile ' + mn] = mobile[mn]
        for k, f in enumerate(families):
            if m.pF[k] == sign:
                parts['in ' + f + ' loops'] = stored[f]
        whole = sum(parts.values())
        for lab, y in parts.items():
            ax.semilogx(x, y / np.maximum(whole, 1e-300), 'o-', ms=3, label=lab)
        ax.set_xlabel('dose (dpa)'); ax.set_ylabel('fraction of the retained ' + title); ax.set_ylim(-0.02, 1.02)
        ax.set_title('where the retained %s are' % title); ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
        save(fig, name)

    # fluxes to the sinks and balance of each species
    balances = [kin.balance(i, w) for i in idx]
    gammas = [kin.absorption(i) for i in idx]
    iv, ii = m.iv, list(m.ms).index(1.0)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    flux_v = np.array([-(b[run.mobile_names[iv]]['loops'] + b[run.mobile_names[iv]]['network']) for b in balances])
    flux_i = np.array([-sum(b[n]['loops'] + b[n]['network'] for j, n in enumerate(run.mobile_names) if run.ms[j] > 0) for b in balances])
    ax.loglog(x, flux_v, 'o-', color=FAMILY_COLORS['c'], label='vacancies to loops and network')
    ax.loglog(x, flux_i, 's-', color=FAMILY_COLORS['a1'], label='interstitials (i, 2i, 3i) to loops and network')
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel('defects per atom per second'); ax.set_title('fluxes to the sinks, volume means')
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'flux_evolution.png')
    table['flux_v_to_sinks'] = list(flux_v); table['flux_i_to_sinks'] = list(flux_i)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for k, f in enumerate(families[:2]):
        net = np.array([np.sum(w * (g[:, k, :] * (m.pF[k] * m.pM)[None, :]).sum(1)) for g in gammas])
        ax.semilogx(x, net, 'o-', color=FAMILY_COLORS[f], label='net growth of %s by absorption' % FAMILY_LABELS[f].split(',')[0])
        table['net_growth_%s' % f] = list(net)
    ax.axhline(0.0, color='0.5', lw=0.8)
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel('defects per atom per second'); ax.set_title('net flux to each loop type (like minus opposite defects)')
    ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    save(fig, 'net_flux_balance.png')

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    imbalance = flux_v - flux_i
    ax.semilogx(x, imbalance, 'ko-')
    ax.axhline(0.0, color='0.5', lw=0.8)
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel('defects per atom per second'); ax.set_title('vacancies minus interstitials arriving at the loops and the network')
    ax.grid(True, which='both', alpha=0.3)
    save(fig, 'defect_imbalance.png')
    table['defect_imbalance'] = list(imbalance)

    for species, name in ((run.mobile_names[iv], 'conservation_channels_v.png'), (run.mobile_names[ii], 'conservation_channels_i.png')):
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        for channel in balances[0][species]:
            y = np.array([b[species][channel] for b in balances])
            if np.any(y != 0):
                ax.loglog(x, np.abs(y), 'o-', ms=3, label='%s (%s)' % (channel, 'source' if y[-1] > 0 else 'loss'))
        ax.set_xlabel('dose (dpa)'); ax.set_ylabel('|rate|, defects per atom per second'); ax.set_title('balance of %s at steady state, volume means' % species)
        ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=7)
        save(fig, name)

    with open(out / 'additional_analysis.csv', 'w', newline='') as stream:
        writer = csv.writer(stream)
        keys = list(table)
        writer.writerow(keys)
        for r in range(len(x)):
            writer.writerow([repr(float(table[k][r])) for k in keys])
    with open(out / 'growth_strain.csv', 'w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(list(table_growth))
        for r in range(len(dose)):
            writer.writerow([repr(float(table_growth[k][r])) for k in table_growth])
    (out / 'provenance.md').write_text('Volume means of the run `%s`, computed by `simulations/postprocessing_extended.py` from its fields and material file.\n\n'
                                       '- weights: share of the volume nearest to each node (Monte Carlo)\n- interior: nodes farther than 80 %% of the half size from the boundary\n'
                                       '- hardening: dispersed barriers, M=%.2f, alpha=%.2f, root sum of squares over the families\n'
                                       '- measured ranges: those quoted in section 4.1 of the D1/M1 report\n'
                                       '- balance of the mobile species: the residual of production, reactions and sinks is what leaves through the grain boundary\n' % (run.dir.name, TAYLOR_FACTOR, OBSTACLE_STRENGTH['a']), encoding='utf-8')
    return balances, x


def boundary_flux_report(run, balances, x):
    lines = ['# Balance of the mobile species', '',
             'Volume means, in defects per atom per second. The mobile species are at steady state at each dose, so the sum of the '
             'channels is zero: the last column is what leaves through the grain boundary (for 2i and 3i it also holds their thermal dissociation, which is not itemized).', '']
    for species in balances[0]:
        channels = list(balances[0][species])
        lines += ['## %s' % species, '', '| dose [dpa] | ' + ' | '.join(channels) + ' | share lost to the boundary |', '|---:|' + '---:|' * (len(channels) + 1)]
        for dose, b in zip(x, balances):
            sources = sum(v for v in b[species].values() if v > 0)
            share = -b[species]['grain boundary'] / sources if sources > 0 else float('nan')
            lines.append('| %g | ' % dose + ' | '.join('%.3e' % b[species][c] for c in channels) + ' | %.1f %% |' % (100.0 * share))
        lines.append('')
    (run.dir / 'boundary_flux.md').write_text('\n'.join(lines), encoding='utf-8')


def figure_denuded_zone(run, kin, plt):
    last = len(run.doses) - 1
    n, c, size, R, d = kin.sizes(last)
    edges = np.logspace(np.log10(max(run.distance[run.distance > 0].min(), 0.3)), np.log10(run.distance.max()), 22)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for k, f in enumerate(run.family_names[:2]):
        N = run.quantity('N_' + f, last)
        xs, ys, ds = [], [], []
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (run.distance >= a) & (run.distance < b)
            if sel.any():
                xs.append(np.sqrt(a * b)); ys.append(N[sel].mean()); ds.append(np.sum(N[sel] * d[sel, k]) / max(N[sel].sum(), 1e-300))
        axes[0].loglog(xs, ys, 'o-', color=FAMILY_COLORS[f], label=FAMILY_LABELS[f].split(',')[0])
        axes[1].semilogx(xs, ds, 'o-', color=FAMILY_COLORS[f], label=FAMILY_LABELS[f].split(',')[0])
    axes[0].set_ylabel(r'number density [m$^{-3}$]'); axes[1].set_ylabel('mean diameter [nm]')
    for ax in axes:
        ax.set_xlabel('distance to the nearest face [nm]'); ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
    fig.suptitle('loops near the grain boundary at %g dpa' % run.doses[last], fontsize=10)
    fig.tight_layout()
    fig.savefig(run.dir / 'gb' / 'denuded_zone.png', dpi=150)
    plt.close(fig)


def figures_size_spectrum(run, kin, plt, w):
    out = run.dir / 'size_spectrum'
    out.mkdir(exist_ok=True)
    last = len(run.doses) - 1
    n, c, size, R, d = kin.sizes(last)
    regions = {'domain': np.ones(len(w), bool), 'interior': run.interior(0.5), 'shell': run.distance < 0.1 * run.distance.max()}
    styles = {'domain': ('k', '-'), 'interior': ('#c62828', '--'), 'shell': ('#00838f', ':')}
    grid = np.logspace(-0.5, 2.6, 300)
    bandwidth = 0.12  # in ln d: each node holds one size, and the kernel only smooths the sum over the nodes
    families = run.family_names[:2]

    def spectrum(k, mask):
        ok = mask & (d[:, k] > 0)
        if not ok.any():
            return np.zeros_like(grid)
        N = run.quantity('N_' + run.family_names[k], last)[ok] * w[ok] / w[mask].sum()
        z = (np.log(grid)[:, None] - np.log(d[ok, k])[None, :]) / bandwidth
        return (np.exp(-0.5 * z * z) / (bandwidth * np.sqrt(2.0 * np.pi)) * N[None, :]).sum(1)

    def panel(ax, k, names):
        for name in names:
            y = spectrum(k, regions[name])
            ax.loglog(grid, np.where(y > 0, y, np.nan), color=styles[name][0], ls=styles[name][1], lw=2, label=name)
        ax.set_xlabel('loop diameter $d$ [nm]'); ax.set_title('%s at %g dpa' % (FAMILY_LABELS[run.family_names[k]].split(',')[0], run.doses[last]))
        ax.grid(True, which='both', alpha=0.3); ax.legend(frameon=False, fontsize=8)
        top = np.nanmax([np.nanmax(spectrum(k, regions[nm])) for nm in names] + [1e-300])
        if top > 0:
            ax.set_ylim(top * 1e-5, top * 3)

    for name in list(regions) + ['regions']:
        fig, axes = plt.subplots(1, len(families), figsize=(5.5 * len(families), 4.0), squeeze=False)
        for j, f in enumerate(families):
            panel(axes[0, j], run.family_names.index(f), list(regions) if name == 'regions' else [name])
        axes[0, 0].set_ylabel(r'd$N$/dln$\,d$  [m$^{-3}$]')
        fig.tight_layout()
        fig.savefig(out / ('spectrum_%s.png' % name), dpi=150)
        plt.close(fig)


def sample_loops(run, kin, index, w, seed=0, maximum=4000):
    """Discrete loops drawn from the fields: for each family, as many loops as its field holds in the grain,
    placed in proportion to the local number of loops, each with the size of its node."""
    rng = np.random.default_rng(seed + index)
    n, c, size, R, d = kin.sizes(index)
    volume = run.volume_b3() * run.b_SI ** 3  # [m^3]
    spacing = (volume / len(run.nodes)) ** (1.0 / 3.0) / run.b_SI
    loops = []
    counts = {}
    for k, f in enumerate(run.family_names):
        lam = run.quantity('N_' + f, index) * w * volume
        lam = np.where(d[:, k] > 0, lam, 0.0)
        expected = lam.sum()
        number = int(round(expected))
        counts[f] = {'expected': float(expected), 'drawn': number}
        if number == 0:
            continue
        keep = min(number, maximum)
        nodes = rng.choice(len(lam), size=keep, p=lam / expected)
        centers = run.nodes[nodes] + 0.5 * spacing * (rng.random((keep, 3)) - 0.5)
        outside = run.depth(centers) < 0.0
        centers[outside] = run.nodes[nodes][outside]  # a loop stays in the grain
        centers = centers * run.b_SI * 1e9
        normal = kin.model.burgers[:, k] / kin.model.bMag[k]
        for center, node in zip(centers, nodes):
            loops.append({'family': f, 'x_nm': center[0], 'y_nm': center[1], 'z_nm': center[2], 'radius_nm': 0.5 * d[node, k],
                          'nx': normal[0], 'ny': normal[1], 'nz': normal[2], 'defects': size[node, k]})
        counts[f]['shown'] = keep
    return loops, counts


def ring(loop, points=40):
    n = np.array([loop['nx'], loop['ny'], loop['nz']])
    a = np.cross(n, [0.0, 0.0, 1.0] if abs(n[2]) < 0.9 else [1.0, 0.0, 0.0])
    a /= np.linalg.norm(a)
    b = np.cross(n, a)
    t = np.linspace(0.0, 2.0 * np.pi, points)
    return np.array([loop['x_nm'], loop['y_nm'], loop['z_nm']])[None, :] + loop['radius_nm'] * (np.cos(t)[:, None] * a + np.sin(t)[:, None] * b)


def draw_box(ax, lo, hi):
    for a in range(3):
        o = [i for i in range(3) if i != a]
        for u in (lo[o[0]], hi[o[0]]):
            for v in (lo[o[1]], hi[o[1]]):
                p0, p1 = np.array(lo, float), np.array(hi, float)
                p0[o[0]] = p1[o[0]] = u
                p0[o[1]] = p1[o[1]] = v
                ax.plot(*zip(p0, p1), color='0.35', lw=0.8)


def figure_loops_3d(run, loops, families, path, title, plt):
    fig = plt.figure(figsize=(6.4, 5.4))
    ax = fig.add_subplot(111, projection='3d')
    lo, hi = run.lo * run.b_SI * 1e9, run.hi * run.b_SI * 1e9
    for p0, p1 in run.outline():
        ax.plot(*zip(p0, p1), color='0.35', lw=0.8)
    shown = [l for l in loops if l['family'] in families]
    for loop in shown:
        r = ring(loop)
        ax.plot(r[:, 0], r[:, 1], r[:, 2], color=FAMILY_COLORS[loop['family']], lw=float(np.clip(0.05 * loop['radius_nm'], 0.6, 1.8)))
    ax.set_axis_off(); ax.set_box_aspect(tuple(hi - lo)); ax.view_init(elev=22, azim=-58)
    ax.set_title('%s\n%d loops' % (title, len(shown)), fontsize=10)
    handles = [plt.Line2D([], [], color=FAMILY_COLORS[f], lw=2, label=FAMILY_LABELS[f]) for f in families]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc='lower left')
    fig.savefig(path, dpi=130)
    plt.close(fig)


ZONE_AXES = {'0001': (2, (0, 1), r'$B \parallel [0001]$: $\langle c\rangle$ face-on, $\langle a\rangle$ edge-on'),
             '0110': (1, (0, 2), r'$B \parallel [01\bar{1}0]$'),
             '2110': (0, (1, 2), r'$B \parallel [2\bar{1}\bar{1}0]$')}


def tem_panel(ax, run, loops, zone, thickness, dose, rng):
    axis, plane, title = ZONE_AXES[zone]
    lo, hi = run.lo * run.b_SI * 1e9, run.hi * run.b_SI * 1e9
    mid = 0.5 * (lo[axis] + hi[axis])
    foil = [l for l in loops if abs((l['x_nm'], l['y_nm'], l['z_nm'])[axis] - mid) <= 0.5 * thickness]
    ax.imshow(0.72 + 0.03 * rng.standard_normal((80, 80)), cmap='gray', vmin=0, vmax=1, extent=(lo[plane[0]], hi[plane[0]], lo[plane[1]], hi[plane[1]]), interpolation='bilinear')
    count = {}
    for loop in foil:
        r = ring(loop)
        ax.fill(r[:, plane[0]], r[:, plane[1]], color='k', alpha=0.10, lw=0)
        ax.plot(r[:, plane[0]], r[:, plane[1]], color='k', lw=1.2, alpha=0.8)
        ax.plot((loop['x_nm'], loop['y_nm'], loop['z_nm'])[plane[0]], (loop['x_nm'], loop['y_nm'], loop['z_nm'])[plane[1]], '.', ms=3, color=FAMILY_COLORS[loop['family']])
        count[loop['family']] = count.get(loop['family'], 0) + 1
    ax.set_xlim(lo[plane[0]], hi[plane[0]]); ax.set_ylim(lo[plane[1]], hi[plane[1]]); ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title('projected loops, ' + title, fontsize=8)
    ax.text(0.02, 0.98, '%g dpa\n%.0f nm foil\n%d loops' % (dose, thickness, len(foil)), transform=ax.transAxes, va='top', fontsize=7, bbox=dict(facecolor='white', alpha=0.85, lw=0.3))
    bar = 0.25 * (hi[plane[0]] - lo[plane[0]])
    x0, y0 = lo[plane[0]] + 0.05 * (hi[plane[0]] - lo[plane[0]]), lo[plane[1]] + 0.06 * (hi[plane[1]] - lo[plane[1]])
    ax.plot([x0, x0 + bar], [y0, y0], 'k-', lw=3)
    ax.text(x0 + 0.5 * bar, y0 + 0.02 * (hi[plane[1]] - lo[plane[1]]), '%.0f nm' % bar, ha='center', fontsize=7)
    handles = [plt_line(f, count.get(f, 0)) for f in run.family_names]
    ax.legend(handles=handles, frameon=True, fontsize=6, loc='upper right')


def plt_line(family, number):
    from matplotlib.lines import Line2D
    return Line2D([], [], color=FAMILY_COLORS[family], marker='.', ls='', label='%s (%d)' % (FAMILY_LABELS[family].split(',')[0], number))


def figures_discrete(run, kin, plt, w, doses=None):
    out = run.dir / 'discrete_loops'
    out.mkdir(exist_ok=True)
    (run.dir / 'tem_slices').mkdir(exist_ok=True)
    (run.dir / '3d').mkdir(exist_ok=True)
    candidates = [i for i, dose in enumerate(run.doses) if dose > 0]
    if doses is None:
        chosen = candidates if len(candidates) <= 4 else [candidates[j] for j in np.unique(np.round(np.linspace(0, len(candidates) - 1, 4)).astype(int))]
    else:
        chosen = [i for i in candidates if any(abs(run.doses[i] - x) <= 1e-9 * max(x, 1e-30) for x in doses)]
    manifest = {'note': 'One size per family at each node: the spread of sizes is that of the field over the grain. Overlapping loops are not merged.', 'doses': {}}
    samples = {}
    for i in chosen:
        loops, counts = sample_loops(run, kin, i, w)
        samples[i] = loops
        tag = format_dose(run.doses[i])
        with open(out / ('loops_%sdpa.csv' % tag), 'w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=['family', 'x_nm', 'y_nm', 'z_nm', 'radius_nm', 'nx', 'ny', 'nz', 'defects'])
            writer.writeheader()
            writer.writerows(loops)
        volume = run.volume_b3() * run.b_SI ** 3
        entry = {}
        for f in run.family_names:
            radii = np.array([l['radius_nm'] for l in loops if l['family'] == f])
            spacing = (volume / max(counts[f]['drawn'], 1)) ** (1.0 / 3.0) * 1e9
            entry[f] = dict(counts[f], mean_radius_nm=float(radii.mean()) if radii.size else 0.0, spacing_nm=float(spacing),
                            diameter_over_spacing=float(2.0 * radii.mean() / spacing) if radii.size else 0.0)
        manifest['doses']['%g' % run.doses[i]] = entry
        figure_loops_3d(run, loops, run.family_names, out / ('loops_%sdpa.png' % tag), 'discrete loops at %g dpa' % run.doses[i], plt)
        for f in run.family_names:
            figure_loops_3d(run, loops, [f], out / ('loops_%s_%sdpa.png' % (f, tag)), '%s at %g dpa' % (FAMILY_LABELS[f], run.doses[i]), plt)
        for f in run.family_names[:2]:
            figure_loops_3d(run, loops, [f], run.dir / '3d' / ('loops_%s_%sdpa.png' % (f, tag)), '%s at %g dpa' % (FAMILY_LABELS[f], run.doses[i]), plt)
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1), encoding='utf-8')

    rng = np.random.default_rng(1)
    thickness = min(100.0, 0.5 * float((run.hi - run.lo)[0] * run.b_SI * 1e9))
    for zone in ZONE_AXES:
        fig, axes = plt.subplots(1, len(chosen), figsize=(5.2 * len(chosen), 5.4), squeeze=False)
        for ax, i in zip(axes[0], chosen):
            tem_panel(ax, run, samples[i], zone, thickness, run.doses[i], rng)
        fig.tight_layout()
        fig.savefig(run.dir / 'tem_slices' / ('tem_B_%s_montage.png' % zone), dpi=130)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(5.6, 5.8))
        tem_panel(ax, run, samples[chosen[-1]], zone, thickness, run.doses[chosen[-1]], rng)
        fig.tight_layout()
        fig.savefig(run.dir / 'tem_slices' / ('tem_B_%s.png' % zone), dpi=130)
        plt.close(fig)
    return manifest


def movies(run, plt, names=None):
    """One animated image per field: the mid-plane figures of 3d/, one frame per snapshot dose."""
    from PIL import Image
    out = run.dir / 'movies'
    out.mkdir(exist_ok=True)
    names = names or [run.mobile_names[0], run.mobile_names[1], 'N_' + run.family_names[0]] + (['N_' + run.family_names[1]] if run.nF > 1 else [])
    written = []
    for name in names:
        frames = [run.dir / '3d' / ('%s_%sdpa.png' % (name, format_dose(d))) for d in run.doses if d > 0]
        frames = [Image.open(f).convert('P', palette=Image.ADAPTIVE) for f in frames if f.is_file()]
        if len(frames) > 1:
            frames[0].save(out / ('%s.gif' % name), save_all=True, append_images=frames[1:], duration=900, loop=0)
            written.append(name)
    return written


def figure_loop_populations(run, kin, plt, w, dose, magnification=None, shown=28):
    """The loop families of the grain at one dose, as in Figure 27 of the D1/M1 report: the content of each family
    on the basal mid-plane, and some of its loops drawn in their habit plane with magnified radii."""
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from scipy.interpolate import LinearNDInterpolator
    index = int(np.argmin([abs(d - dose) for d in run.doses]))
    magnification = magnification or {'c': 1.5, 'a': 7.0}
    loops, counts = sample_loops(run, kin, index, w, seed=27)
    nm = run.b_SI * 1e9
    lo, hi = run.lo * nm, run.hi * nm
    points = run.nodes * nm
    resolution = 110
    grid = np.linspace(0.0, 1.0, resolution)
    A, B = np.meshgrid(grid, grid)
    P = np.zeros((resolution, resolution, 3))
    P[..., 0] = lo[0] + A * (hi[0] - lo[0]); P[..., 1] = lo[1] + B * (hi[1] - lo[1]); P[..., 2] = 0.5 * (lo[2] + hi[2])
    outside = run.depth(P.reshape(-1, 3) / nm).reshape(resolution, resolution) < 0.0
    order = [f for f in run.family_names if f != 'c'] + ['c']
    fig = plt.figure(figsize=(13, 11))
    rng = np.random.default_rng(27)
    for panel, family in enumerate(order):
        ax = fig.add_subplot(2, 2, panel + 1, projection='3d')
        values = run.quantity('C_' + family, index)
        positive = values[values > 0]
        if family == 'c':
            norm = colors.Normalize(vmin=float(positive.min()), vmax=float(positive.max()))
            V = LinearNDInterpolator(points, values)(P.reshape(-1, 3)).reshape(resolution, resolution)
            scale = 'lin'
        else:
            norm = colors.LogNorm(vmin=max(float(positive.min()), float(positive.max()) * 1e-2), vmax=float(positive.max()))
            V = np.exp(LinearNDInterpolator(points, np.log(np.maximum(values, norm.vmin)))(P.reshape(-1, 3)).reshape(resolution, resolution))
            scale = 'log'
        face = cm.jet(norm(np.nan_to_num(V, nan=norm.vmin)))
        face[outside | np.isnan(V), 3] = 0.0
        ax.plot_surface(P[..., 0], P[..., 1], P[..., 2], facecolors=face, rstride=1, cstride=1, shade=False, antialiased=False, linewidth=0)
        for p0, p1 in run.outline():
            ax.plot(*zip(p0, p1), color='0.4', lw=0.8)
        mine = [l for l in loops if l['family'] == family]
        factor = magnification['c' if family == 'c' else 'a']
        if len(mine) > shown:
            mine = [mine[i] for i in rng.choice(len(mine), shown, replace=False)]
        discs = []
        for loop in mine:
            big = dict(loop, radius_nm=factor * loop['radius_nm'])
            discs.append(ring(big, 24))
        ax.add_collection3d(Poly3DCollection(discs, facecolor=FAMILY_COLORS[family], edgecolor='k', linewidths=0.3, alpha=0.9))
        ax.set_axis_off(); ax.set_box_aspect(tuple(hi - lo)); ax.view_init(elev=20, azim=-60)
        name = r'$\langle c\rangle$' if family == 'c' else r'$\langle a\rangle_%s$' % family[1:]
        ax.set_title('%s loop population at %g dpa (platelet radii x%g, not to scale)\n%d of %d loops drawn' % (name, run.doses[index], factor, len(mine), counts[family]['drawn']), fontsize=9)
        bar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.jet), ax=ax, shrink=0.5, pad=0.02)
        bar.ax.set_title(scale, fontsize=8)
        bar.set_label('content of the family, defects per atom', fontsize=8)
    fig.tight_layout()
    (run.dir / 'figures').mkdir(exist_ok=True)
    path = run.dir / 'figures' / ('loop_populations_%sdpa.png' % format_dose(run.doses[index]))
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def gather_extended(run, plt, samples=2000000, progress=print, discrete_doses=None, population_dose=None):
    """Writes the extended figures. Returns the list of what was written, for the report."""
    w = run.weights(samples)
    kin = Kinetics(run)
    written = []
    balances, x = figures_volume_average(run, kin, plt, w)
    written.append('volume_average/')
    boundary_flux_report(run, balances, x)
    written.append('boundary_flux.md')
    figure_denuded_zone(run, kin, plt)
    figures_size_spectrum(run, kin, plt, w)
    written += ['gb/denuded_zone.png', 'size_spectrum/']
    progress('volume averages, balance and size spectra written')
    manifest = figures_discrete(run, kin, plt, w, discrete_doses)
    written += ['discrete_loops/', 'tem_slices/']
    progress('discrete loops and projected slices written')
    figure_loop_populations(run, kin, plt, w, population_dose if population_dose else run.doses[-1])
    written.append('figures/loop_populations')
    try:
        if movies(run, plt):
            written.append('movies/')
    except ImportError:
        progress('Pillow is not installed: no movies')
    return written, manifest
