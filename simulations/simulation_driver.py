"""Logic behind run_simulation.ipynb: build MoDELib, stage a case, run it, and gather its output.

The case is spatially resolved cluster dynamics of irradiated zirconium in a single grain whose
surface is a grain boundary: the mobile species are solved at steady state on the finite-element
mesh, and the immobile loop populations are integrated at every node by CVODE.

Everything is resolved from the location of this file. On Windows the executables are built and
run through WSL; on Linux and macOS they are built and run directly.
"""

import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import time
from pathlib import Path

SIMULATIONS_DIR = Path(__file__).resolve().parent
ROOT = SIMULATIONS_DIR.parent
LIBRARY = ROOT / 'Library'
OUTPUT_DIR = SIMULATIONS_DIR / 'output'
GEOMETRY_DIR = SIMULATIONS_DIR / 'geometry'
REFERENCE_DIR = SIMULATIONS_DIR / 'reference'
USE_WSL = os.name == 'nt'

PRESETS = {
    # the cube of 200 nm at the doses of the ENNDS reference run of the same geometry
    'reference': {'CASE': {'geometry': 'cube_200nm'},
                  'MARCH': {'doses': [1e-4, 1e-3, 1e-2, 0.1, 1.0, 2.0, 5.0, 10.0], 'steps_per_interval': 3}},
    # the hexagonal crystal of the D1/M1 report (its Figure 27), with the dose steps of the DisloCluster run: 1 dpa
    'hex400': {'CASE': {'geometry': 'hex_400nm'},
               'MARCH': {'doses': [1.0, 6.0, 11.0, 16.0, 21.0, 26.0], 'step_dpa': 1.0},
               # the runs of the report nucleated the loops from the cascades only (its section 5, item 1)
               'MATERIAL': {'overrides': {'loopClusteringNucleation': '0'}},
               'OUTPUT': {'population_dose': 6.0}},
    # the hexagonal crystal with the two terms of the DisloCluster code that differ from the report:
    # one cluster per clustering reaction, and a climb speed of coalescence without the Burgers-vector magnitude.
    # It checks the numerical accuracy of MoDELib against the DisloCluster run, term for term
    'hex400dc': {'CASE': {'geometry': 'hex_400nm', 'tag': 'hex400_DisloCluster_terms'},
                 'MARCH': {'doses': [1.0, 6.0, 11.0, 16.0, 21.0, 26.0], 'step_dpa': 1.0},
                 'MATERIAL': {'overrides': {'loopClusteringNucleation': '1', 'loopNucleationPerReaction': '1', 'coalescenceClimbBurgers': '0'}},
                 'OUTPUT': {'population_dose': 6.0}},
    # one interval: a check that the build and the workflow are in order
    'smoke': {'CASE': {'geometry': 'cube_200nm'}, 'MARCH': {'doses': [1e-4], 'steps_per_interval': 3}},
}

DEFAULTS = {
    'CASE': {'geometry': 'cube_200nm', 'tag': None},
    'MATERIAL': {'material_file': 'Zr_CD4opt.txt', 'temperature_K': 573.0, 'dose_rate_dpa_s': 1.0e-7, 'overrides': {}},
    'MARCH': {'doses': None, 'steps_per_interval': 3, 'step_dpa': None, 'immobile_integrator': 'cvode'},
    'RUN': {'threads': None, 'keep_all_configurations': False},
    'OUTPUT': {'mesh_cut_figures': True, 'figure_doses': None, 'voronoi_samples': 2000000, 'population_dose': None},
}


# ---------------------------------------------------------------------------------------------
# shell
def shell_path(path):
    """The path as the shell that runs the executables sees it."""
    path = Path(path).resolve()
    if USE_WSL:
        return '/mnt/' + path.drive[0].lower() + path.as_posix()[2:]
    return path.as_posix()


def shell(command, log=None, check=True):
    """Runs a bash command (in WSL on Windows). Returns its output; raises if it fails."""
    args = ['wsl.exe', '-e', 'bash', '-lc', command] if USE_WSL else ['bash', '-lc', command]
    result = subprocess.run(args, capture_output=True, text=True, encoding='utf-8', errors='replace')
    output = result.stdout + result.stderr
    if log is not None:
        Path(log).write_text(output, encoding='utf-8')
    if check and result.returncode != 0:
        raise RuntimeError('command failed (%d): %s\n%s' % (result.returncode, command, output[-2000:]))
    return output


# ---------------------------------------------------------------------------------------------
# build
def find_sundials():
    """The install prefix of SUNDIALS, or None when it is expected in a system location."""
    candidates = [os.environ.get('SUNDIALS_ROOT'), '$HOME/sundials-install', str(ROOT.parent / 'Libraries' / 'sundials-7.1.1')]
    for candidate in candidates:
        if candidate:
            prefix = candidate if candidate.startswith('$') else shell_path(candidate)
            if 'yes' in shell('test -d %s/lib/cmake/sundials && echo yes || true' % prefix, check=False):
                return prefix
    return None


def build(force=False, build_dir='build', jobs=None):
    """Configures and builds the command-line tools. Returns the executables and whether CVODE is in the build."""
    build_path = ROOT / build_dir
    if force and build_path.exists():
        shutil.rmtree(build_path)
    prefix = find_sundials()
    configure = 'cmake -S %s -B %s -DUSE_PYBIND11=OFF -DBUILD_DDQT=OFF' % (shell_path(ROOT), shell_path(build_path))
    if prefix:
        configure += ' -DCMAKE_PREFIX_PATH=%s' % prefix
    build_path.mkdir(exist_ok=True)
    log = shell(configure, log=build_path / 'configure.log')
    jobs = jobs or os.cpu_count() or 2
    shell('cmake --build %s -j%d' % (shell_path(build_path), jobs), log=build_path / 'build.log')
    executables = {
        'DDomp': build_path / 'tools' / 'DDomp' / 'DDomp',
        'microstructureGenerator': build_path / 'tools' / 'MicrostructureGenerator' / 'microstructureGenerator',
    }
    for name, path in executables.items():
        if not path.is_file():
            raise RuntimeError('%s was not built; see %s' % (name, build_path / 'build.log'))
    return {'executables': executables, 'sundials': 'SUNDIALS_VERSION' in log or '_MODEL_SUNDIALS_' in (build_path / 'CMakeCache.txt').read_text(errors='replace'),
            'sundials_prefix': prefix, 'build_dir': build_path}


# ---------------------------------------------------------------------------------------------
# configuration
def configuration(preset='reference', **sections):
    """The full configuration: the defaults, the preset, then the sections given (CASE, MATERIAL, MARCH, RUN, OUTPUT)."""
    if preset not in PRESETS:
        raise KeyError('unknown preset %r; choose from %s' % (preset, sorted(PRESETS)))
    config = json.loads(json.dumps(DEFAULTS))
    for name, values in PRESETS[preset].items():
        config[name].update(values)
    config['PRESET'] = preset
    for name, values in sections.items():
        if name not in DEFAULTS:
            raise KeyError('unknown section %r; the sections are %s' % (name, sorted(DEFAULTS)))
        for key, value in (values or {}).items():
            if key not in DEFAULTS[name]:
                raise KeyError('unknown key %r in %s; the keys are %s' % (key, name, sorted(DEFAULTS[name])))
            config[name][key] = value
    return config


def material_value(text, key):
    match = re.search(r'^' + re.escape(key) + r' *=([^;]*);', text, re.M)
    if not match:
        raise KeyError('%s not found in the material file' % key)
    return [float(x) for x in match.group(1).split()]


def set_variable(text, key, value):
    """Replaces the value of key=...; in the text of an input file. Raises if the key is absent."""
    pattern = re.compile(r'^(' + re.escape(key) + r' *=)[^;]*;', re.M)
    if not pattern.search(text):
        raise KeyError('variable %s not found' % key)
    return pattern.sub(lambda m: m.group(1) + str(value) + ';', text, count=1)


def validate(config, built=None):
    """Checks the configuration before anything runs. Returns the derived quantities."""
    geometry = GEOMETRY_DIR / config['CASE']['geometry']
    case = json.loads((geometry / 'case.json').read_text())
    mesh = next(geometry.glob('*.msh'))
    material = LIBRARY / 'Materials' / config['MATERIAL']['material_file']
    if not material.is_file():
        raise FileNotFoundError(material)
    text = material.read_text(encoding='utf-8')
    for key, value in config['MATERIAL']['overrides'].items():
        text = set_variable(text, key, value)
    doses = [float(d) for d in config['MARCH']['doses']]
    if not doses or any(b <= a for a, b in zip([0.0] + doses, doses)):
        raise ValueError('MARCH doses must increase from zero: %s' % doses)
    if config['MARCH']['step_dpa']:
        # steps of a given dose: each interval must hold a whole number of them
        steps = [(b - a) / float(config['MARCH']['step_dpa']) for a, b in zip([0.0] + doses, doses)]
        if any(abs(s - round(s)) > 1e-9 or round(s) < 1 for s in steps):
            raise ValueError('MARCH step_dpa=%s does not divide the dose intervals %s' % (config['MARCH']['step_dpa'], doses))
        steps = [int(round(s)) for s in steps]
    else:
        if int(config['MARCH']['steps_per_interval']) < 1:
            raise ValueError('MARCH steps_per_interval must be at least 1')
        steps = [int(config['MARCH']['steps_per_interval'])] * len(doses)
    if config['MARCH']['immobile_integrator'] not in ('cvode', 'euler'):
        raise ValueError('MARCH immobile_integrator is cvode or euler')
    if built is not None and config['MARCH']['immobile_integrator'] == 'cvode' and not built['sundials']:
        raise RuntimeError('immobile_integrator=cvode needs a build with SUNDIALS, which was not found')
    if case['geometry']['type'] not in ('cubic', 'hexagonal'):
        raise ValueError('the geometry type is cubic or hexagonal')
    # the mesh file spans the unit box: the extents of the bounding box scale it
    extents_nm = case['geometry'].get('extents_nm') or [case['geometry']['size_nm'] * a for a in (case['geometry'].get('aspect') or [1.0, 1.0, 1.0])]
    b_SI = material_value(text, 'b_SI')[0]
    cs = (material_value(text, 'mu0_SI')[0] / material_value(text, 'rho_SI')[0]) ** 0.5
    try:
        omega_SI = material_value(text, 'atomicVolume_SI')[0]
    except KeyError:
        omega_SI = b_SI ** 3 * 2.0 ** 0.5 / 2.0  # HEX lattice
    return {'geometry_dir': geometry, 'case': case, 'mesh': mesh, 'material_text': text, 'doses': doses,
            'b_SI': b_SI, 'cs_SI': cs, 'time_unit_s': b_SI / cs, 'omega_SI': omega_SI,
            'size_b': extents_nm[0] * 1e-9 / b_SI, 'extents_b': [e * 1e-9 / b_SI for e in extents_nm],
            'geometry_type': case['geometry']['type'], 'steps': steps,
            'mobile_species': [int(x) for x in material_value(text, 'mobileSpeciesVector')],
            'families': len(material_value(text, 'immobileSpeciesVector'))}


# ---------------------------------------------------------------------------------------------
# run directory and provenance
def git_state():
    try:
        commit = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(['git', '-C', str(ROOT), 'status', '--porcelain', '--untracked-files=no'], capture_output=True, text=True).stdout.strip())
        return commit, dirty
    except OSError:
        return '', False


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def new_run_directory(config):
    commit, dirty = git_state()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S')
    tag = config['CASE']['tag'] or '%s_%s' % (config['CASE']['geometry'], config['PRESET'])
    name = '%s_%s%s_%s' % (stamp, commit[:7] or 'nogit', '-dirty' if dirty else '', tag)
    run_dir = OUTPUT_DIR / name
    run_dir.mkdir(parents=True)
    return run_dir


def write_text_lf(path, text):
    """Input files are written with LF line endings: the mesh reader does not accept CRLF."""
    with open(path, 'w', encoding='utf-8', newline='\n') as stream:
        stream.write(text.replace('\r\n', '\n'))


def stage(config, derived, run_dir):
    """Writes the input files of the run. Returns their paths."""
    inputs = run_dir / 'inputFiles'
    for folder in (inputs, run_dir / 'evl', run_dir / 'F', run_dir / 'logs'):
        folder.mkdir(exist_ok=True)
    material_name = config['MATERIAL']['material_file']
    text = set_variable(derived['material_text'], 'doseRate_dpaPerSec', '%.15e' % config['MATERIAL']['dose_rate_dpa_s'])
    write_text_lf(inputs / material_name, text)
    write_text_lf(inputs / derived['mesh'].name, derived['mesh'].read_text(encoding='utf-8'))

    dc = (LIBRARY / 'DefectiveCrystal' / 'DefectiveCrystal.txt').read_text(encoding='utf-8')
    for key, value in (('physics', 'ClusterDynamics'), ('useFEM', 1), ('maxResolveSteps', 0), ('outputFrequency', 1), ('startAtTimeStep', -1)):
        dc = set_variable(dc, key, value)
    write_text_lf(inputs / 'DefectiveCrystal.txt', dc)

    cd = (LIBRARY / 'ClusterDynamics' / 'ClusterDynamics.txt').read_text(encoding='utf-8')
    cd = set_variable(cd, 'useClusterDynamicsFEM', 1)
    cd = set_variable(cd, 'immobileIntegrator', config['MARCH']['immobile_integrator'])
    write_text_lf(inputs / 'ClusterDynamics.txt', cd)

    extents = derived['extents_b']
    polycrystal = ('materialFile=%s;\nabsoluteTemperature=%.6f; # [K] simulation temperature\nmeshFile=%s; # mesh file\n'
                   'F=%.14f 0 0\n  0 %.14f 0\n  0 0 %.14f; # mesh deformation gradient, x = F*(X-X0)\n'
                   'X0=0 0 0; # mesh shift\nperiodicFaceIDs=; # every face is a grain boundary\n'
                   'C2G1 =\n1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0;\n'
                   % (material_name, config['MATERIAL']['temperature_K'], derived['mesh'].name, extents[0], extents[1], extents[2]))
    write_text_lf(inputs / 'polycrystal.txt', polycrystal)
    write_text_lf(inputs / 'initialMicrostructure.txt', '')
    (run_dir / 'config.json').write_text(json.dumps(config, indent=1), encoding='utf-8')
    return sorted(p for p in inputs.iterdir() if p.is_file())


def provenance(run_dir, config, derived, built, record):
    commit, dirty = git_state()
    data = {
        'application': 'MoDELib: spatially resolved cluster dynamics of irradiated zirconium',
        'run': run_dir.name,
        'status': record.get('status'),
        'started_utc': record.get('started'), 'finished_utc': record.get('finished'), 'wall_s': record.get('wall_s'),
        'modelib_commit': commit, 'dirty': dirty,
        'machine': {'system': platform.platform(), 'processor': platform.machine(), 'cpus': os.cpu_count(), 'python': platform.python_version(),
                    'runs_through_WSL': USE_WSL},
        'build': {'sundials': built['sundials'], 'sundials_prefix': built['sundials_prefix']},
        'configuration': config,
        'derived': {k: derived[k] for k in ('doses', 'steps', 'b_SI', 'cs_SI', 'time_unit_s', 'omega_SI', 'size_b', 'extents_b', 'geometry_type', 'mobile_species', 'families')},
        'input_files': {p.name: sha256(p) for p in sorted((run_dir / 'inputFiles').iterdir()) if p.is_file()},
        'executables': {name: sha256(path) for name, path in built['executables'].items()},
        'intervals': record.get('intervals', []),
    }
    (run_dir / 'provenance.json').write_text(json.dumps(data, indent=1), encoding='utf-8')
    lines = ['# Provenance: %s' % run_dir.name, '',
             '- **Status:** %s' % data['status'],
             '- **Started / finished (UTC):** %s / %s (%.1f s)' % (data['started_utc'], data['finished_utc'], data['wall_s'] or 0.0),
             '- **MoDELib:** commit `%s`%s' % (commit, ' (working tree modified)' if dirty else ''),
             '- **Machine:** %s, %s, %s CPUs; Python %s%s' % (platform.platform(), platform.machine(), os.cpu_count(), platform.python_version(), '; executables run through WSL' if USE_WSL else ''),
             '- **SUNDIALS in the build:** %s' % ('yes' if built['sundials'] else 'no'), '',
             '## Configuration', '', '```json', json.dumps(config, indent=1), '```', '',
             '## Input files (SHA-256)', '']
    lines += ['- `%s` %s' % item for item in data['input_files'].items()]
    lines += ['', '## Executables (SHA-256)', ''] + ['- `%s` %s' % item for item in data['executables'].items()]
    (run_dir / 'provenance.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return data


# ---------------------------------------------------------------------------------------------
# run
def run(config, derived, built, run_dir, progress=print):
    """Runs the march: one call of DDomp per dose interval, each continuing from the configuration
    file written by the previous one. Returns the record of the run."""
    record = {'status': 'running', 'started': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), 'intervals': []}
    t_start = time.time()
    threads = config['RUN']['threads'] or os.cpu_count() or 2
    environment = 'export OMP_NUM_THREADS=%d; cd %s; ' % (threads, shell_path(run_dir))
    executables = {name: shell_path(path) for name, path in built['executables'].items()}
    shell(environment + executables['microstructureGenerator'] + ' .', log=run_dir / 'logs' / 'microstructureGenerator.log')

    dc_file = run_dir / 'inputFiles' / 'DefectiveCrystal.txt'
    dose_unit = config['MATERIAL']['dose_rate_dpa_s'] * derived['time_unit_s']  # dpa per unit of time
    snapshots = {0: 0.0}
    done = 0
    previous = 0.0
    try:
        for index, dose in enumerate(derived['doses']):
            t0 = time.time()
            steps = derived['steps'][index]
            dt = (dose - previous) / steps / dose_unit
            # the configuration of step n holds the fields before the update of that step: one more step writes the fields at the end of the interval
            text = set_variable(set_variable(dc_file.read_text(encoding='utf-8'), 'Nsteps', done + steps + 1), 'dtMax', '%.15e' % dt)
            write_text_lf(dc_file, text)
            log = run_dir / 'logs' / ('DDomp_%02d.log' % index)
            output = shell(environment + executables['DDomp'] + ' .', log=log, check=False)
            done += steps
            if 'simulation steps completed' not in output or not (run_dir / 'evl' / ('evl_%d.txt' % done)).is_file():
                raise RuntimeError('DDomp did not complete the interval ending at %g dpa; see %s' % (dose, log))
            snapshots[done] = dose
            record['intervals'].append({'dose_from': previous, 'dose_to': dose, 'steps': steps, 'dt': dt, 'wall_s': time.time() - t0,
                                        'configuration': 'evl_%d.txt' % done})
            progress('interval %d of %d: %g -> %g dpa in %.0f s' % (index + 1, len(derived['doses']), previous, dose, time.time() - t0))
            previous = dose
        record['status'] = 'completed'
    except Exception as error:
        record['status'] = 'failed: %s' % error
        raise
    finally:
        record['finished'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
        record['wall_s'] = time.time() - t_start
        record['snapshots'] = {str(k): v for k, v in snapshots.items()}
        (run_dir / 'run_record.json').write_text(json.dumps(record, indent=1), encoding='utf-8')
        if not config['RUN']['keep_all_configurations']:
            keep = {'evl_%d.txt' % k for k in snapshots} | {'cdNodes.txt'}
            for path in (run_dir / 'evl').iterdir():
                if path.name not in keep:
                    path.unlink()
        provenance(run_dir, config, derived, built, record)
    return record


def find_runs():
    """The run directories, oldest first."""
    return sorted(p for p in OUTPUT_DIR.iterdir() if p.is_dir() and (p / 'provenance.json').is_file())


def latest_run():
    runs = find_runs()
    return runs[-1] if runs else None
