"""Generate a polycrystal with random prismatic loops and 33-direction strain loading."""
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
from countGrains import countGrains


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
setInputVariable('inputFiles/'+DCfile,'useFEM','0')
setInputVariable('inputFiles/'+DCfile,'Nsteps','10000')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'maxResolveSteps','0')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'dtMax','1e25')
setInputVariable('inputFiles/'+DCfile,'outputFrequency','1')  # output frequency

setInputVariable('inputFiles/'+DCfile,'startAtTimeStep','-1')
setInputVector('inputFiles/'+DCfile,'periodicImageSize',np.array([0,0,0]),'image counts for the active periodic directions')

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate=str(repository)+'/Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
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
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.zeros(6),'applied stress')
# Voigt order: 11, 22, 33, 12, 23, 13. Strain rate is per simulation time b/cs.
strainRate33 = 1e-8
setInputVector('inputFiles/'+elasticDeformatinoFile, 'ExternalStrainRate',
               np.array([0.0, 0.0, strainRate33, 0.0, 0.0, 0.0]), 'applied strain rate in the global 33 direction')
setInputVector('inputFiles/'+elasticDeformatinoFile, 'stiffnessRatio',
               np.array([0.0, 0.0, 1e20, 0.0, 0.0, 0.0]), 'strain control in 33; zero-stress control in the other components')


# Create polycrystal.txt using local material file
meshFile='poly30_periodicXYZ.msh';
meshFileTemplate=str(repository)+'/Library/Meshes/'+meshFile;
numberOfGrains = countGrains(meshFileTemplate)
print(f'Mesh contains {numberOfGrains} grains')
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
# Grain 1: crystal and global coordinates aligned.
pf.grain1globalX1 = np.array([1.0, 0.0, 0.0])
pf.grain1globalX3 = np.array([0.0, 0.0, 1.0])
pf.boxScaling=np.array([1e-6,1e-6,1e-6])/b_SI # scale of each mesh coordinate in Burgers vector units
pf.X0=np.array([0,0,0]) # Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
# Explicit example orientations for all grains detected in the mesh. These rational rotations
# are chosen for lattice-compatible periodic translations; they do not use the
# random orientations stored in the source mesh's ElsetOrientations section.
rx = np.array([[1., 0., 0.], [0., 0.8, -0.6], [0., 0.6, 0.8]])
ry = np.array([[0.8, 0., 0.6], [0., 1., 0.], [-0.6, 0., 0.8]])
rz = np.array([[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]])
grainOrientations = [np.eye(3), rx, rx.T, ry, ry.T, rz, rz.T,
                     rx @ ry, ry @ rz, rz @ rx]
# Preserve the original first ten orientations, then add distinct products
# of two elementary rotations. Their denominators still divide 25.
elementaryRotations = [rx, rx.T, ry, ry.T, rz, rz.T]
for first in elementaryRotations:
    for second in elementaryRotations:
        orientation = first @ second
        if not any(
                np.allclose(orientation, existing, rtol=0, atol=1e-12)
                for existing in grainOrientations):
            grainOrientations.append(orientation)
# Reuse this finite lattice-compatible orientation pool for larger meshes.
# The first 30 orientations are unchanged from the original tutorial.
orientationPool = grainOrientations
grainOrientations = [orientationPool[index % len(orientationPool)]
                     for index in range(numberOfGrains)]
# Products have denominators up to 25. Round the axial BCC repeat counts
# to multiples of 25 so all grain orientations admit the periodic translations.
pf.compute()
repeatCounts = 25 * np.maximum(1, np.round(pf.boxScaling / (25 * pf.boxEdgesLatticeLengths)))
meshScale = 2.0  # Double linear dimensions while preserving lattice-compatible repeats.
pf.boxScaling = meshScale * repeatCounts * pf.boxEdgesLatticeLengths
pf.write('inputFiles')
with open('inputFiles/polycrystal.txt', 'a') as polycrystalFile:
    for grainID, orientation in enumerate(grainOrientations[1:], start=2):
        rows = [' '.join(f'{value:.17g}' for value in row) for row in orientation]
        polycrystalFile.write(f'C2G{grainID}=' + '\n     '.join(rows)
                             + '; # crystal-to-global rotation\n')

# Random prismatic loops throughout all grains. Density is total dislocation
# line length per unit volume; radii are in meters.
microstructureName = 'prismaticLoopsDensity.txt'
microstructureFile = 'inputFiles/' + microstructureName
shutil.copy2(repository / 'Library/Microstructures/prismaticLoopsDensity.txt',
             microstructureFile)
setInputVariable(microstructureFile, 'targetDensity_SI', '1e13')  # m^-2
setInputVariable(microstructureFile, 'radiusDistributionMean_SI', '1e-7')  # 100 nm
setInputVariable(microstructureFile, 'radiusDistributionStd_SI', '0')  # fixed radius
setInputVector(microstructureFile, 'allowedGrainIDs', np.array([-1]), 'random placement throughout all grains')
setInputVariable(microstructureFile, 'allowedSlipSystemIDs', '-1')  # all slip systems
with open('inputFiles/initialMicrostructure.txt', 'w') as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile=' + microstructureName + ';\n')

# Generate the populated frame 0 with microstructureGenerator before running DDomp.
print('Inputs ready: ' + str(Path.cwd()))
print('Run microstructureGenerator, then DDomp, using this simulation folder.')
