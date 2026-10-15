# Checks the cases run by runCases.sh against the expected values.
#
#   python3 checkCases.py
#
# Rate equations of the D1/M1 report (CVODE): the fields are compared with referenceRates.py.
# The fields are read from the text files evl_<n>.txt, which hold six significant digits.

import os, re, sys
import numpy as np

def value(fileName, key):
    for line in open(fileName):
        if re.match(r'^' + re.escape(key) + r' *=', line):
            return [float(x) for x in line.split('=')[1].split(';')[0].split()]
    raise RuntimeError(key + ' not found in ' + fileName)

def fields(case, step):
    lines = open(case + '/evl/evl_%d.txt' % step).read().split('\n')
    nNodes = int(lines[9])  # the cluster-dynamics block is the last of the file
    rows = [l for l in lines[10:] if l.strip()]
    return np.array([[float(x) for x in r.split()] for r in rows[-nNodes:]])

results = []
def check(name, ok, detail):
    results.append(ok)
    print(('PASS  ' if ok else 'FAIL  ') + name + ': ' + detail)

have = lambda case: os.path.isfile(case + '/evl/evl_0.txt')
mSize, nF = 4, 4

if have('cvode'):
    import referenceRates
    case = 'cvode'
    mat = case + '/inputFiles/Zr_CD4opt.txt'
    T = value(case + '/inputFiles/polycrystal.txt', 'absoluteTemperature')[0]
    model = referenceRates.Model(mat, T)
    dt = value(case + '/inputFiles/DefectiveCrystal.txt', 'dtMax')[0]

    def compare(case, step, model, chi, nodes):
        # the fields of evl_<step> are those before the update of that step
        f0, f1 = fields(case, step), fields(case, step + 1)
        worst = 0.0
        for node in nodes:
            y0 = np.concatenate([f0[node, mSize:mSize + nF] * model.omega, f0[node, mSize + nF:]])
            y1 = model.integrate(y0, f0[node, :mSize], chi, dt)
            code = np.concatenate([f1[node, mSize:mSize + nF] * model.omega, f1[node, mSize + nF:]])
            worst = max(worst, np.abs(code / y1 - 1.0).max())
        return worst

    nNodes = len(fields(case, 0))
    sample = list(range(0, nNodes, 400))
    for step in (0, 2):
        error = compare(case, step, model, model.chi(), sample)
        check('cvode, step %d against the reference equations' % step, error < 1e-4, 'largest relative difference %.1e over %d nodes and %d fields' % (error, len(sample), 2 * nF))

    f = fields(case, 3)
    interior = f[:, 0] > 0.9 * f[:, 0].max()
    nC = f[interior, mSize] * model.omega
    expected = model.G * model.epsV / model.nNuc[0] * model.tau
    # production times lifetime is the density without coalescence, which can only lower it
    check('cvode, saturation of the vacancy-loop density', nC.max() <= expected * (1.0 + 1e-6) and nC.min() > 0.9 * expected,
          'n = %.5e per atom at the interior nodes; production times lifetime G*eps_vL/n_nuc*tau = %.5e is its upper bound' % (nC.mean(), expected))
    symmetry = max(np.abs(f[:, mSize + 2] / f[:, mSize + 1] - 1.0).max(), np.abs(f[:, mSize + 3] / f[:, mSize + 1] - 1.0).max())
    check('cvode, the three <a> families without stress', symmetry < 1e-9, 'largest relative difference of the number densities %.1e' % symmetry)
    check('cvode, no negative value', f.min() >= 0.0, 'smallest value of the fields %.1e' % f.min())
    b3 = model.b_SI**3
    print('      after %.1f dpa, interior: N_c = %.3e m^-3, N_a = %.3e m^-3 per family, defects per loop %.0f (c) and %.0f (a)' % (3 * dt * model.G, f[interior, mSize].mean() / b3, f[interior, mSize + 1].mean() / b3, (f[interior, mSize + nF] / f[interior, mSize] / model.omega).mean(), (f[interior, mSize + nF + 1] / f[interior, mSize + 1] / model.omega).mean()))

    if have('cvodeRestart'):
        a, b = fields('cvode', 3), fields('cvodeRestart', 1)
        scale = np.abs(a).max(0)
        error = (np.abs(a - b).max(0) / scale).max()
        check('cvodeRestart, a new study started from a configuration file', error < 1e-8, 'largest difference from the uninterrupted run %.1e, relative to the largest value of each field' % error)

    if have('cvodeStress'):
        case = 'cvodeStress'
        stress = np.zeros((3, 3))
        stress[0, 0] = value(case + '/inputFiles/ElasticDeformation.txt', 'ExternalStress0')[0]
        chi = model.chi(stress)
        error = compare(case, 0, model, chi, sample)
        check('cvodeStress, against the reference equations', error < 1e-4, 'nucleation fractions of the families %s, largest relative difference %.1e' % (np.round(chi, 4), error))

print('%d checks, %d failed' % (len(results), results.count(False)))
sys.exit(results.count(False))
