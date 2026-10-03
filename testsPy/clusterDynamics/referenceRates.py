# Reference implementation of the rate equations of the immobile clusters
# (docs/reports/GW_Phase4_D1M1.pdf, sections 2.3 and 2.4), written independently of the C++ code
# and integrated with SciPy. checkCases.py compares it with the fields computed by MoDELib.
#
# Units: lengths in b, times in b/cs, as in MoDELib. n and c are per atom.

import re
import numpy as np
from scipy.integrate import solve_ivp

kB_eV = 8.617333262e-5  # [eV/K]
eV = 1.602176634e-19    # [J]


def value(fileName, key):
    for line in open(fileName):
        if re.match(r'^' + re.escape(key) + r' *=', line):
            return [float(x) for x in line.split('=')[1].split(';')[0].split()]
    raise RuntimeError(key + ' not found in ' + fileName)


def matrix(fileName, key, rows):
    text = open(fileName).read()
    block = re.search(r'^' + re.escape(key) + r' *=([^;]*);', text, re.M).group(1)
    return np.array([float(x) for x in block.split()]).reshape(rows, -1)


class Model:

    def __init__(self, materialFile, T):
        v = lambda key: np.array(value(materialFile, key))
        self.T = T
        self.kT = kB_eV * T  # [eV]
        self.b_SI = v('b_SI')[0]
        self.mu_SI = v('mu0_SI')[0]
        self.cs = np.sqrt(self.mu_SI / v('rho_SI')[0])
        time = self.b_SI / self.cs  # [s]
        self.kT_code = self.kT * eV / (self.mu_SI * self.b_SI**3)  # in units of mu*b^3
        self.omega = v('atomicVolume_SI')[0] / self.b_SI**3
        self.ms = v('mobileSpeciesVector')
        self.pM = np.sign(self.ms)
        nM = len(self.ms)
        Em = matrix(materialFile, 'mobileSpeciesEnergyMigration_eV', nM)
        D0 = matrix(materialFile, 'mobileSpeciesD0_SI', nM)
        D = D0 * np.exp(-Em / self.kT) / (self.b_SI * self.cs)  # components 11 12 13 22 23 33
        self.Dbar = np.array([np.linalg.det(np.array([[d[0], d[1], d[2]], [d[1], d[3], d[4]], [d[2], d[4], d[5]]]))**(1.0 / 3.0) for d in D])
        p = (D[:, 5] / D[:, 0])**(1.0 / 6.0)
        self.G = v('doseRate_dpaPerSec')[0] * time * v('mobileSpeciesSurvivingEfficiency')[0]
        f = v('mobileSpeciesCascadeFractions')
        self.epsV = max(1.0 - f[self.pM < 0].sum(), 0.0)
        self.epsI = max(1.0 - f[self.pM > 0].sum(), 0.0)
        self.Eb = v('mobileSpeciesBindingEnergy_eV')
        self.iv = list(self.ms).index(-1.0)

        # families
        self.pF = np.sign(v('immobileSpeciesVector'))
        nF = len(self.pF)
        lattice = np.array([[1.0, 0.5, 0.0], [0.0, np.sqrt(3.0) / 2.0, 0.0], [0.0, 0.0, np.sqrt(8.0 / 3.0)]])  # HEX, columns a1 a2 c
        latticeBurgers = matrix(materialFile, 'immobileSpeciesBurgers', 3)
        self.burgers = lattice @ latticeBurgers
        self.bMag = np.linalg.norm(self.burgers, axis=0)
        basal = (np.abs(latticeBurgers[2]) > 0) & (np.abs(latticeBurgers[0]) == 0) & (np.abs(latticeBurgers[1]) == 0)
        Z0 = matrix(materialFile, 'discreteDislocationBias', 2)  # rows: vacancy loops, interstitial loops
        self.Z = np.array([[Z0[0 if self.pF[k] < 0 else 1, m] * (p[m] if basal[k] else 0.5 * (p[m] + p[m]**-2)) for m in range(nM)] for k in range(nF)])
        try:
            self.Z = self.Z * v('loopSinkScale')[:, None]  # fitted factor on the sink strength of each family
        except RuntimeError:
            pass
        self.nmin, self.nmax, self.nNuc = v('nmin'), v('nmax'), v('nNuc')
        self.w0, self.n_s = v('w0')[0], v('n_s')[0]
        self.tau = v('tau0_vLoop_SI')[0] * np.exp(v('Ea_vLoop_eV')[0] / self.kT) / time
        self.cLL, self.cLN = v('cLL'), v('cLN')
        self.kLL, self.kLN = v('kappaLL')[0], v('kappaLN')[0]
        self.rhoN = v('rhoNetwork_SI')[0] * self.b_SI**2
        self.rmin = v('r_min') / self.b_SI
        try:
            self.embryosOnly = v('dissolveEmbryosOnly')[0] > 0
        except RuntimeError:
            self.embryosOnly = False
        # the two forms of the DisloCluster code, for comparison
        try:
            self.perReaction = v('loopNucleationPerReaction')[0] > 0
        except RuntimeError:
            self.perReaction = False
        try:
            self.climbBurgers = v('coalescenceClimbBurgers')[0] > 0
        except RuntimeError:
            self.climbBurgers = True

        # reactions between mobile species whose product is not mobile
        r = np.where(np.abs(self.ms) == 1, (3.0 * self.omega / 4.0 / np.pi)**(1.0 / 3.0), np.sqrt(np.abs(self.ms) * self.omega / np.pi))
        self.channels = []
        try:
            clustering = v('loopClusteringNucleation')[0] > 0
        except RuntimeError:
            clustering = True
        for a, b, prefactor in matrix(materialFile, 'reactionPrefactorMap', nM * (nM + 1) // 2):
            ia, ib = list(self.ms).index(a), list(self.ms).index(b)
            if clustering and prefactor > 0 and (a + b > self.ms.max() or a + b < self.ms.min()):
                K = prefactor * 4.0 * np.pi * (r[ia] + r[ib]) * (self.Dbar[ia] + self.Dbar[ib]) / self.omega
                self.channels.append((ia, ib, K, 1 if a + b > 0 else 0))

    def sigmoid(self, k, m):
        m0 = 0.5 * (self.nmin[k] + self.nmax[k]) - self.n_s
        return 1.0 / (1.0 + np.exp(-(m - m0) / (self.w0 * (self.nmax[k] - self.nmin[k]))))

    def loopRadius(self, k, m):
        return np.sqrt(m * self.omega / (np.pi * self.bMag[k]))

    def radius(self, k, m):  # Eq. (29)
        S = self.sigmoid(k, m)
        return (1.0 - S) * (m * self.omega / np.sqrt(8.0))**(1.0 / 3.0) + S * self.loopRadius(k, m)

    def chi(self, stress=None):  # Eq. (38). stress in units of mu
        nF = len(self.pF)
        w = np.ones(nF)
        if stress is not None:
            for k in range(nF):
                S = self.sigmoid(k, self.nNuc[k])
                bHat = self.burgers[:, k] / self.bMag[k]
                Omega = (S * np.outer(bHat, bHat) + (1.0 - S) * np.eye(3) / 3.0) * self.nNuc[k] * self.omega
                w[k] = np.exp(np.sum(stress * Omega) / self.kT_code)
        out = np.zeros(nF)
        for sign in (-1, 1):
            sel = self.pF == sign
            out[sel] = w[sel] / w[sel].sum()
        return out

    def rates(self, y, cM, chi):  # Eqs. (35) and (40)
        nF = len(self.pF)
        n, c = np.maximum(y[:nF], 0.0), np.maximum(y[nF:], 0.0)
        cM = np.maximum(cM, 0.0)
        C = cM / np.abs(self.ms)
        sNuc, events = [0.0, 0.0], [0.0, 0.0]
        for ia, ib, K, sign in self.channels:
            loss = K * C[ia] * C[ib]
            sNuc[sign] += abs(self.ms[ia]) * loss if ia == ib else (abs(self.ms[ia]) + abs(self.ms[ib])) * loss
            events[sign] += loss
        nDot, cDot = np.zeros(nF), np.zeros(nF)
        for k in range(nF):
            vac = self.pF[k] < 0
            source = chi[k] * (self.G * (self.epsV if vac else self.epsI) + sNuc[0 if vac else 1])
            nDot[k] = chi[k] * (self.G * (self.epsV if vac else self.epsI) / self.nNuc[k] + events[0 if vac else 1]) if self.perReaction else source / self.nNuc[k]
            cDot[k] = source
            if n[k] > 1e-30 and c[k] > 1e-30:
                m = c[k] / n[k]
                R = self.radius(k, m)
                x = (R / self.rmin[k] - 1.0) / 0.3 if self.rmin[k] > 0 else 1.0
                gate = 0.0 if x <= 0 else (1.0 if x >= 1 else x * x * (3.0 - 2.0 * x))
                gain = 0.0
                for mi in range(len(self.ms)):
                    like = self.pF[k] * self.pM[mi] > 0
                    Gamma = (1.0 if like else gate) * self.Dbar[mi] * 2.0 * np.pi * self.Z[k, mi] * R * n[k] / self.omega * cM[mi]  # Eq. (42)
                    cDot[k] += Gamma if like else -Gamma
                    gain += Gamma if like else 0.0
                if vac:
                    nDot[k] -= n[k] / self.tau  # Eq. (39)
                    cDot[k] -= (min(self.nNuc[k] * n[k], c[k]) if self.embryosOnly else c[k]) / self.tau  # Eq. (51), or the content of the embryos
                    cDot[k] -= self.Dbar[self.iv] * 2.0 * np.pi * self.Z[k, self.iv] * R * n[k] / self.omega * np.exp(-self.Eb[self.iv] / self.kT)  # Eq. (54)
                r = self.loopRadius(k, m)
                N = n[k] / self.omega
                vAbs = 0.5 * r * gain / c[k] * (1.0 if self.climbBurgers else self.bMag[k])
                rLL, rLN = min(r, N**(-1.0 / 3.0)), min(r, self.rhoN**-0.5)
                LL = self.cLL[k] * vAbs * N**(1.0 / 3.0) * (1.0 - np.exp(-self.kLL * 4.0 / 3.0 * np.pi * rLL**3 * N))   # Eqs. (58), (59)
                LN = self.cLN[k] * vAbs * np.sqrt(self.rhoN) * (1.0 - np.exp(-self.kLN * np.pi * rLN**2 * self.rhoN))
                nDot[k] -= (LL + LN) * n[k]  # Eq. (36)
                cDot[k] -= LN * c[k]         # Eq. (41)
        return np.concatenate([nDot, cDot])

    def integrate(self, y0, cM, chi, dt):
        # time scaled by dt, as numbers of order one suit the integrator
        solution = solve_ivp(lambda s, y: self.rates(y, cM, chi) * dt, (0.0, 1.0), y0, method='BDF', rtol=1e-10, atol=1e-28)
        if not solution.success:
            raise RuntimeError(solution.message)
        return solution.y[:, -1]
