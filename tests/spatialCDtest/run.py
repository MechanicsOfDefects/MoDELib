"""Run the unmodified spatial-CD physical setup and retain an auditable report."""
import argparse
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RTOL = 1.e-4
ATOL = 1.e-10


def scalar(path, name):
    match = re.search(r'^\s*' + re.escape(name) + r'\s*=\s*([^;]+);', path.read_text(), re.M)
    if match is None:
        raise ValueError(f'Missing {name} in {path}')
    return float(match.group(1))


def compare(actual, reference):
    if actual.shape != reference.shape or actual.shape != (120, 7):
        raise ValueError(f'Expected a complete 120-step history, got {actual.shape}')
    if not np.isfinite(actual).all():
        raise ValueError('Nonfinite output')
    if not np.allclose(actual[:, :3], reference[:, :3], rtol=1.e-12, atol=1.e-12):
        raise ValueError('Time/dose sampling changed')
    error = np.abs(actual[:, 3:] - reference[:, 3:])
    limits = ATOL + RTOL * np.abs(reference[:, 3:])
    if np.any(error > limits):
        step, component = np.unravel_index(np.argmax(error / limits), error.shape)
        raise ValueError(f'Strain/swelling regression at step {step}, component {component}: '
                         f'{actual[step, component+3]} versus {reference[step, component+3]}')
    if not np.allclose(actual[:, 3], actual[:, 4], rtol=1.e-9, atol=1.e-12):
        raise ValueError('In-plane strain symmetry lost')
    return dict(max_absolute_error=error.max(axis=0).tolist(), rtol=RTOL, atol=ATOL)


def history(report):
    values = np.atleast_2d(np.loadtxt(report / 'F/F_0.txt'))
    labels = (report / 'F/F_labels.txt').read_text().splitlines()
    if not np.array_equal(values[:, 0], np.arange(120)):
        raise ValueError('Missing, repeated, or unexpected run IDs')
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite simulation output')
    material = report / 'inputFiles/Zr4_Fitted.txt'
    t = values[:, labels.index('time [b/cs]')]
    conversion = scalar(material, 'b_SI') / np.sqrt(scalar(material, 'mu0_SI') / scalar(material, 'rho_SI'))
    dose = t * conversion * scalar(material, 'doseRate_dpaPerSec')
    strain = np.column_stack([values[:, labels.index(f'betaP_{ij}')] for ij in ['11', '22', '33']])
    for ij in ['12', '13', '23', '21', '31', '32']:
        if np.max(np.abs(values[:, labels.index(f'betaP_{ij}')])) > 1.e-10:
            raise ValueError('Unexpected shear plastic distortion')
    result = np.column_stack([t, dose, dose + .5, strain, strain.mean(axis=1)])
    np.savetxt(report / 'growth-swelling-data.csv', result, delimiter=',',
               header='time_b_over_cs,dose_dpa,plotted_dose_dpa,beta11,beta22,beta33,trace_beta_over_3', comments='')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--reports', type=Path, default=HERE / 'reports')
    args = parser.parse_args()
    if importlib.util.find_spec("matplotlib") is None:
        parser.error("Matplotlib is required to generate the test plots.")
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    report = args.reports.resolve() / stamp
    report.mkdir(parents=True)
    print(f'Report: {report}', flush=True)
    shutil.copytree(HERE / 'inputFiles', report / 'inputFiles')
    shutil.copy2(ROOT / 'Library/Meshes/unitCube_15K.msh', report / 'inputFiles/unitCube_15K.msh')
    for name in ['experimental-data.csv', 'reference.csv', 'stored-tutorial-result.pdf']:
        shutil.copy2(HERE / name, report / name)
    for name in ['F', 'config', '.cache']:
        (report / name).mkdir()
    settings = {
        'physics': 'ClusterDynamics', 'useFEM': True, 'temperature_K': 553,
        'steps': 120, 'dt_b_over_cs': 3.4e18, 'output_frequency': 1,
        'mesh': 'unitCube_15K.msh', 'boxScaling': [3093, 3093, 3093],
        'active_periodic_faces': [], 'discrete_microstructure': [],
        'experimental_series': {'EXPIa': 'Iodide purity, a', 'EXPIc': 'Iodide purity, c',
                                'EXPZa': 'Zone-refined purity, a', 'EXPZc': 'Zone-refined purity, c'},
        'experimental_strain_unit': '1e-4', 'plot_dose_offset_dpa': .5,
        'plot_convention': 'Original script: beta components in percent; swelling trace(beta)/3 unscaled.',
        'reference': 'Successful WeakProblem tutorial run from 2026-10-19, not a fit to experiment.',
        'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in sorted((report / 'inputFiles').iterdir())},
    }
    (report / 'settings.json').write_text(json.dumps(settings, indent=2) + '\n')
    env = dict(os.environ, OMP_NUM_THREADS='12', MPLBACKEND='Agg',
               MPLCONFIGDIR=str(report / '.cache'), XDG_CACHE_HOME=str(report / '.cache'))
    started = time.monotonic()
    result = {'passed': False, 'report': str(report)}
    try:
        with (report / 'simulation.log').open('w') as log:
            subprocess.run([str(args.executable.resolve()), str(report)], cwd=report,
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1650)
        actual = history(report)
        # Keep plots even when a numerical regression is detected.
        with (report / 'plot.log').open('w') as log:
            subprocess.run([sys.executable, str(HERE / 'plotGrowthSwellingDensity.py'), str(report)],
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=60)
        result.update(compare(actual, np.loadtxt(report / 'reference.csv', delimiter=',', skiprows=1)))
        result.update(passed=True, outputs=len(actual), final_values=actual[-1].tolist())
    except Exception as error:
        result['error'] = str(error)
        print(f'FAIL: {error}\nInspect {report}', file=sys.stderr)
    finally:
        result['elapsed_seconds'] = time.monotonic() - started
        (report / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(('PASS' if result['passed'] else 'FAIL') + f': spatialCDtest; report={report}', flush=True)
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
