"""Generate the uniform-load example with one rotated ellipsoidal inclusion."""
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
folders=['evl','F','inputFiles']
for x in folders:
    if not os.path.exists(x):
        os.makedirs(x)

# Make a local copy of DefectiveCrystal parameters file and modify that copy if necessary
DCfile='DefectiveCrystal.txt'
DCfileTemplate=str(repository)+'/Library/DefectiveCrystal/'+DCfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DCfileTemplate,'inputFiles/'+DCfile)
setInputVariable('inputFiles/'+DCfile,'physics','DislocationDynamics ElasticDeformation InclusionMicrostructure')
setInputVariable('inputFiles/'+DCfile,'useFEM','0')
setInputVariable('inputFiles/'+DCfile,'Nsteps','10000')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'maxResolveSteps','0')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'dtMax','1e25')
setInputVariable('inputFiles/'+DCfile,'outputFrequency','10')  # output frequency

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
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.01,0.0,0.0,0.0]),'applied stress')

# Create polycrystal.txt using local material file
meshFile='unitCube24.msh';
meshFileTemplate=str(repository)+'/Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
# Keep crystal and global coordinates aligned: C2G = identity.
pf.grain1globalX1 = np.array([1.0, 0.0, 0.0])
pf.grain1globalX3 = np.array([0.0, 0.0, 1.0])
pf.boxScaling=np.array([1e-6,1e-6,1e-6])/b_SI # length of box edges in Burgers vector units
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')
microstructureFiles = []

microstructureFileID=1;
microstructureName='prismaticLoopsDensity.txt';
microstructureFileTemplate=str(repository)+'/Library/Microstructures/'+microstructureName;
microstructureFileLocal='inputFiles/microstructureFile'+str(microstructureFileID)+'.txt'
microstructureFiles.append('microstructureFile'+str(microstructureFileID)+'.txt')
print("\033[1;32mCreating  microstructureFile"+str(microstructureFileID)+"\033[0m")
shutil.copy2(microstructureFileTemplate,microstructureFileLocal) # target filename is /dst/dir/file.ext
setInputVariable(microstructureFileLocal,'targetDensity_SI','1e12')
setInputVariable(microstructureFileLocal,'radiusDistributionMean_SI','1e-07')
setInputVariable(microstructureFileLocal,'radiusDistributionStd_SI','0e-07')
setInputVector(microstructureFileLocal,'allowedGrainIDs',np.array([-1]),'set of grain IDs where loops are allowed. Use -1 for all grains')
setInputVector(microstructureFileLocal,'allowedSlipSystemIDs',np.array([-1]),'set of slip system IDs whose Burgers vector are allowed to be the prism axis. Use -1 for all slip systems')

# Arithmetic semiaxis means/std in meters: full diameter 40 nm, thickness 3 nm.
axesLognormalDistributionMean_SI = np.array([20e-9, 20e-9, 1.5e-9])
axesLognormalDistributionStd_SI = np.zeros(3) # Circular bases, fixed dimensions.
principalEigenDistortion = np.array([1e-2, -5e-3, 2e-2])
# Individual precipitate example (disabled; uncomment to use it instead of density).
# crystalAxes = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
# crystalEigenDistortion = crystalAxes @ np.diag(principalEigenDistortion) @ crystalAxes.T
# inclusionFile = 'inputFiles/ellipsoidalInclusionIndividual.txt'
# shutil.copy2(repository / 'Library/Microstructures/ellipsoidalInclusionIndividual.txt',
#              inclusionFile)
# setInputVector(inclusionFile, 'center', np.zeros(3), '[b] center in simulation coordinates')
# setInputMatrix(inclusionFile, 'crystalAxes_SI', crystalAxes @ np.diag(axesLognormalDistributionMean_SI))
# setInputVector(inclusionFile, 'crystalEigenDistortion', crystalEigenDistortion.flatten(),
#                'dimensionless crystal-frame tensor: 11 12 13 21 22 23 31 32 33')
# setInputVariable(inclusionFile, 'velocityReductionFactor', '1.0')
# setInputVariable(inclusionFile, 'phaseID', '0')
# setInputVariable(inclusionFile, 'allowOverlap', '0')
# setInputVariable(inclusionFile, 'allowOutside', '0')
# microstructureFiles.append('ellipsoidalInclusionIndividual.txt')

# Three EllipsoidalInclusionDensitySpecifications, one for each {100} platelet base.
# Diameter = 40 nm and full thickness = 3 nm (oblate ellipsoids).
# Zero size spread keeps both in-plane semiaxes exactly equal, hence circular bases.
totalTargetDensity = 20e19 # [m^-3], split equally among the three families
patternVector_SI = np.zeros(3) # No positional plane pattern.
# Columns: first in-plane axis, second in-plane axis, platelet normal.
plateletFamilies = [
    ('100', np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])),
    ('010', np.array([[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]])),
    ('001', np.eye(3)),
]
for family, crystalAxes in plateletFamilies:
    crystalEigenDistortion = crystalAxes @ np.diag(principalEigenDistortion) @ crystalAxes.T
    inclusionName = 'ellipsoidalInclusionsDensity_' + family + '.txt'
    inclusionFile = 'inputFiles/' + inclusionName
    shutil.copy2(repository / 'Library/Microstructures/ellipsoidalInclusionsDensity.txt',
                 inclusionFile)
    setInputVariable(inclusionFile, 'targetDensity', str(totalTargetDensity / len(plateletFamilies)))
    setInputMatrix(inclusionFile, 'crystalAxes', crystalAxes)
    setInputVector(inclusionFile, 'axesLognormalDistributionMean_SI', axesLognormalDistributionMean_SI,
                   '[m] arithmetic mean semiaxis lengths')
    setInputVector(inclusionFile, 'axesLognormalDistributionStd_SI', axesLognormalDistributionStd_SI,
                   'zero spread: circular 40 nm diameter, 3 nm thickness')
    setInputVector(inclusionFile, 'patternVector_SI', patternVector_SI, 'zero: no positional pattern')
    setInputVector(inclusionFile, 'crystalEigenDistortion', crystalEigenDistortion.flatten(),
                   'dimensionless crystal-frame tensor: 11 12 13 21 22 23 31 32 33')
    setInputVariable(inclusionFile, 'velocityReductionFactor', '1.0')
    setInputVariable(inclusionFile, 'phaseID', '0')
    setInputVariable(inclusionFile, 'allowOverlap', '0')
    setInputVariable(inclusionFile, 'allowOutside', '0')
    setInputVariable(inclusionFile, 'maxAttempts', '10000')
    microstructureFiles.append(inclusionName)

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    for micrFile in microstructureFiles:
        initialMicrostructureFile.write('microstructureFile='+micrFile+';\n')
