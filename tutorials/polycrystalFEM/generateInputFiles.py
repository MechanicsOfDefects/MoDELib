"""Generate an empty XY-periodic polycrystal with FEM strain-rate loading in z."""
import argparse
import os
from pathlib import Path
import sys

repository = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent,
                    help="Simulation folder (defaults to this tutorial folder).")
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)
os.chdir(args.output_dir.resolve())
sys.path.insert(0, str(repository / "python"))
from modlibUtils import *


# Create folder structure
folders=['config','F','inputFiles']
for x in folders:
    if not os.path.exists(x):
        os.makedirs(x)

# Make a local copy of DefectiveCrystal parameters file and modify that copy if necessary
DCfile='DefectiveCrystal.txt'
DCfileTemplate=str(repository)+'/Library/DefectiveCrystal/'+DCfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DCfileTemplate,'inputFiles/'+DCfile)
setInputVariable('inputFiles/'+DCfile,'physics','DislocationDynamics ElasticDeformation')
setInputVariable('inputFiles/'+DCfile,'useFEM','1')
setInputVariable('inputFiles/'+DCfile,'Nsteps','100')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'maxResolveSteps','0')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'dtMax','1e4')
setInputVariable('inputFiles/'+DCfile,'outputFrequency','1')  # output frequency

setInputVariable('inputFiles/'+DCfile,'startAtTimeStep','0')
p = Path('inputFiles') / DCfile
p.write_text('\n'.join('startAtTimeStep=0; # start from the empty frame 0' if line.startswith('startAtTimeStep') else line for line in p.read_text().splitlines()) + '\n')
setInputVector('inputFiles/'+DCfile,'periodicImageSize',np.array([0,0]),'image counts for the active periodic directions')

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate=str(repository)+'/Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','fixed') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dxMax','10') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','1') # don't filter velocity if noise is enabled
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1.0') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'Lmin','25')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','150')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'glideSolverType','Galerkin')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','none')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'crossSlipModel','0')  # crossSlipModel

# Make a local copy of material file, and modify that copy if necessary
materialFile='W.txt';
materialFileTemplate=str(repository)+'/Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
#setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','full<111>{110}')
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate=str(repository)+'/Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVariable('inputFiles/'+elasticDeformatinoFile, 'useElasticDeformationFEM', '1')
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.zeros(6),'applied stress')
# Voigt order: 11, 22, 33, 12, 23, 13. Strain rate is per simulation time b/cs.
strainRate33 = 1e-8
setInputVector('inputFiles/'+elasticDeformatinoFile, 'ExternalStrainRate',
               np.array([0.0, 0.0, strainRate33, 0.0, 0.0, 0.0]), 'applied strain rate in the global 33 direction')
setInputVector('inputFiles/'+elasticDeformatinoFile, 'stiffnessRatio',
               np.array([1e20, 0.0, 1e20, 0.0, 0.0, 0.0]), 'strain control in 11 and 33; zero stress in 22 and shear components')


# Create polycrystal.txt using local material file
meshFile='polycrystalCube_10_periodicXY.msh';
meshFileTemplate=str(repository)+'/Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
# Grain 1: crystal and global coordinates aligned.
pf.grain1globalX1 = np.array([1.0, 0.0, 0.0])
pf.grain1globalX3 = np.array([0.0, 0.0, 1.0])
pf.boxScaling=np.array([4e-6,4e-6,4e-6])/b_SI # scale of each mesh coordinate in Burgers vector units
pf.X0=np.array([0,0,0]) # Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
# Explicit example orientations for all 10 grains. These rational rotations
# are chosen for lattice-compatible periodic translations; they do not use the
# random orientations stored in the source mesh's ElsetOrientations section.
rx = np.array([[1., 0., 0.], [0., 0.8, -0.6], [0., 0.6, 0.8]])
ry = np.array([[0.8, 0., 0.6], [0., 1., 0.], [-0.6, 0., 0.8]])
rz = np.array([[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]])
grainOrientations = [np.eye(3), rx, rx.T, ry, ry.T, rz, rz.T,
                     rx @ ry, ry @ rz, rz @ rx]
# The source mesh has unit volume and unit X/Y repeats, with an irregular
# outer boundary. Scale to approximately (4 micrometers)^3, not a 4 um
# bounding cube. Round repeats to multiples of 25 for all orientations.
pf.compute()
repeatCounts = 25 * np.maximum(1, np.round(pf.boxScaling / (25 * pf.boxEdgesLatticeLengths)))
pf.boxScaling = repeatCounts * pf.boxEdgesLatticeLengths
pf.write('inputFiles')
with open('inputFiles/polycrystal.txt', 'a') as polycrystalFile:
    for grainID, orientation in enumerate(grainOrientations[1:], start=2):
        rows = [' '.join(f'{value:.17g}' for value in row) for row in orientation]
        polycrystalFile.write(f'C2G{grainID}=' + '\n     '.join(rows)
                             + '; # crystal-to-global rotation\n')

# Enable DD without inserting any defects. The generator writes an empty frame 0.
with open('inputFiles/initialMicrostructure.txt', 'w') as initialMicrostructureFile:
    initialMicrostructureFile.write('# No initial defects; no microstructureFile entries.\n')

lengthSI = pf.F[0, 0] * b_SI
print(f'XY repeat length and height: {lengthSI * 1e6:.6f} um')
print(f'Mesh volume: {(lengthSI * 1e6)**3:.6f} um^3')
print('Strain increment per step: 1e-4; last saved strain: 0.0099 at frame 99.')
print('Inputs ready: ' + str(Path.cwd()))
print('Run microstructureGenerator, then DDomp, using this simulation folder.')
