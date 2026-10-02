# Checks the cases run by runCases.sh against the expected values.
#
#   python3 checkCases.py [folder of a reference run of the case "reg"]
#
# Nucleation: the sources of number density and of content are compared with their analytical values.
# Conversion to discrete loops: the totals printed by the code are checked (number, stored defects, sink strength).
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

def log(case):
    return re.sub(r'\x1b\[[0-9;]*m', '', open(case + '/DDomp.log').read())

def conversion(case):
    # per family: field loops, drawn, coalesced, inserted, field defects, inserted defects, field sink, inserted sink
    text = log(case)
    num = r'([-+0-9.eE]+)'
    loops = re.findall(r'loops: +field ' + num + ', drawn ' + num + ', after coalescence ' + num + ', inserted ' + num, text)
    defects = re.findall(r'defects: +field ' + num + ', inserted ' + num + ', left in the field ' + num, text)
    sinks = re.findall(r'sink strength: field ' + num + ', inserted ' + num, text)
    added = int(re.findall(r'(\d+) loops added to the DislocationNetwork', text)[0])
    return np.array(loops, float), np.array(defects, float), np.array(sinks, float), added

results = []
def check(name, ok, detail):
    results.append(ok)
    print(('PASS  ' if ok else 'FAIL  ') + name + ': ' + detail)

have = lambda case: os.path.isfile(case + '/evl/evl_0.txt')
mSize, nF = 4, 4

if have('nuc') or have('clus'):
    mat = ('nuc' if have('nuc') else 'clus') + '/inputFiles/Zr4_Fitted.txt'
    b_SI = value(mat, 'b_SI')[0]
    cs = np.sqrt(value(mat, 'mu0_SI')[0] / value(mat, 'rho_SI')[0])  # shear wave speed [m/s]
    omega = np.sqrt(2.0) / 2.0  # atomic volume of the HEX lattice, in b^3
    G0 = value(mat, 'doseRate_dpaPerSec')[0] * b_SI / cs  # dose rate per unit of time b/cs
    eta = value(mat, 'mobileSpeciesSurvivingEfficiency')[0]
    ms = np.array(value(mat, 'mobileSpeciesVector'))

if have('nuc') and have('nuc0'):
    mat = 'nuc/inputFiles/Zr4_Fitted.txt'
    eps = np.array(value(mat, 'loopCascadeFractions'))
    nNuc = np.array(value(mat, 'nNuc'))
    dt = value('nuc/inputFiles/DefectiveCrystal.txt', 'dtMax')[0]
    f0, f1, z0, z1 = fields('nuc', 0), fields('nuc', 1), fields('nuc0', 0), fields('nuc0', 1)
    dN = f1[:, mSize:mSize + nF] - f0[:, mSize:mSize + nF]
    expected = G0 * eta * eps / nNuc / omega * dt
    error = np.abs(dN / expected - 1.0).max()
    check('cascade nucleation, number density', error < 2e-3, 'largest nodal error %.1e of G*eps/nNuc/Omega*dt = %s per b^3' % (error, expected))
    dc = f1[:, mSize + nF:] - z1[:, mSize + nF:]  # same mobile field and same absorption in the two cases
    expected = G0 * eta * eps * dt
    error = np.abs(dc / expected - 1.0).max()
    check('cascade nucleation, content', error < 5e-3, 'largest nodal error %.1e of G*eps*dt = %s' % (error, expected))
    check('no nucleation without the keys', np.abs(z1[:, mSize:mSize + nF] - z0[:, mSize:mSize + nF]).max() == 0.0, 'number densities of case nuc0 unchanged')

if have('clus') and have('clusRef'):
    dt = value('clus/inputFiles/DefectiveCrystal.txt', 'dtMax')[0]
    channels = re.findall(r'clustering channel (\S+) \+ (\S+), rate coefficient \(in 1/s\) (\S+)', log('clus'))
    f0, f1, r1 = fields('clus', 0), fields('clus', 1), fields('clusRef', 1)
    C = f0[:, :mSize] / np.abs(ms)  # cluster concentrations: the fields hold |n|*C
    number = np.zeros(len(f0))
    content = np.zeros(len(f0))
    for a, b, K in channels:
        ia, ib = list(ms).index(float(a)), list(ms).index(float(b))
        loss = float(K) * b_SI / cs * C[:, ia] * C[:, ib]
        number += 0.5 * loss if ia == ib else loss
        content += abs(ms[ia]) * loss if ia == ib else (abs(ms[ia]) + abs(ms[ib])) * loss
    share = 1.0 / 3.0  # three interstitial families, equal shares
    dN = (f1[:, mSize + 1] - f0[:, mSize + 1]).mean()
    dc = (f1[:, mSize + nF + 1] - r1[:, mSize + nF + 1]).mean()
    ratioN = dN / (share * number.mean() / omega * dt)
    ratioC = dc / (share * content.mean() * dt)
    check('clustering nucleation, number density', abs(ratioN - 1.0) < 0.1, 'channels %s, mean source/expected = %.3f' % ([c[0] + '+' + c[1] for c in channels], ratioN))
    check('clustering nucleation, content', abs(ratioC - 1.0) < 0.1, 'mean source/expected = %.3f, defects per new cluster %.2f' % (ratioC, dc / (dN * omega)))
    check('clustering nucleation, vacancy family', np.abs(f1[:, mSize] - f0[:, mSize]).max() == 0.0, 'no interstitial reaction feeds the <c> family')

for case, lumping in (('discAll', 1.0), ('disc', 100.0), ('discMin', 100.0)):
    if have(case):
        loops, defects, sinks, added = conversion(case)
        header = open(case + '/evl/evl_%d.txt' % (len([f for f in os.listdir(case + '/evl') if f.startswith('evl_')]) - 1)).read().split('\n')[:10]
        families = range(1, nF) if case == 'discMin' else range(nF)
        error = max(abs(defects[k, 1] / defects[k, 0] - 1.0) for k in families)
        check(case + ', stored defects', error < 1e-9, 'inserted/field - 1 = %.1e' % error)
        check(case + ', number of loops', all(abs(loops[k, 1] - loops[k, 0] / lumping) <= 0.5 for k in range(nF)), 'drawn %s for %s clusters, %s after coalescence' % (loops[:, 1].astype(int), np.round(loops[:, 0], 2), loops[:, 2].astype(int)))
        ratio = np.array([sinks[k, 1] / sinks[k, 0] for k in families])
        check(case + ', sink strength', np.abs(ratio * np.sqrt(lumping) - 1.0).max() < 1e-2, 'inserted/field = %s, expected %.4f' % (np.round(ratio, 5), 1.0 / np.sqrt(lumping)))
        if case != 'discAll':  # the last configuration of discAll is written before the conversion
            check(case + ', loops in the network', int(header[1]) == added == int(loops[:, 3].sum()), '%d loops in the last configuration' % int(header[1]))
        if case == 'discMin':
            check(case + ', minimumLoopSize', loops[0, 3] == 0 and defects[0, 1] == 0.0, 'the <c> family stays in the fields')
        if case == 'disc':
            steps = len(re.findall(r'^runID', log(case), re.M))
            check(case + ', climb steps after the conversion', 'simulation steps completed' in log(case), '%d steps completed' % steps)
        if case != 'discAll':
            F = np.loadtxt(case + '/F/F_0.txt', ndmin=2)
            print('      betaP_11 before and after the conversion: %.4e -> %.4e;  betaP_33: %.4e -> %.4e' % (F[0, 3], F[1, 3], F[0, 11], F[1, 11]))

if have('reg'):
    f0, f5 = fields('reg', 0), fields('reg', 5)
    check('reg, number densities', np.abs(f5[:, mSize:mSize + nF] - f0[:, mSize:mSize + nF]).max() == 0.0, 'constant without nucleation')
    if len(sys.argv) > 1:
        F, Fref = np.loadtxt('reg/F/F_0.txt', ndmin=2), np.loadtxt(sys.argv[1] + '/F/F_0.txt', ndmin=2)
        n = min(len(F), len(Fref))
        error = [np.abs(F[:n, j] / Fref[:n, j] - 1.0).max() for j in (3, 7, 11)]
        check('reg, growth strain against the reference', max(error) < 1e-3, 'largest relative difference of betaP_11, betaP_22, betaP_33: %s' % np.array(error))

print('%d checks, %d failed' % (len(results), results.count(False)))
sys.exit(results.count(False))
