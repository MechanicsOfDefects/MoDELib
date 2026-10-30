import argparse
import os
from pathlib import Path
import sys

repository = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description="Generate the annealing tutorial inputs.")
parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)
os.chdir(args.output_dir.resolve())
sys.path.insert(0, str(repository / "python"))
from modlibUtils import *
from countGrains import countGrains

# Annealing of vacancy and interstitial Frank loops by climb, coupled to the diffusion of single vacancies.
# The material file used here (Zr_BMD19.txt) has ONE mobile species and no immobile species, while the
# numbers of cluster-dynamics species are fixed when MoDELib is compiled (defaults: 4 mobile, 8 immobile).
# Build MoDELib for this tutorial with:
#   cmake -S . -B build -DMODELIB_CD_MSIZE=1 -DMODELIB_CD_ISIZE=0

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
setInputVariable('inputFiles/'+DCfile,'physics','DislocationDynamics ClusterDynamics ElasticDeformation')
setInputVariable('inputFiles/'+DCfile,'useFEM','1')
setInputVariable('inputFiles/'+DCfile,'Nsteps','2000')  # number of simulation steps
# Converge DD climb and CD concentration at fixed time and geometry.
setInputVariable('inputFiles/'+DCfile,'solvedQuantityMaxIterations','100')
setInputVariable('inputFiles/'+DCfile,'solvedQuantityRelativeTolerance','1e-3')
setInputVariable('inputFiles/'+DCfile,'solvedQuantityRelaxation','1.0')
setInputVariable('inputFiles/'+DCfile,'dtMax','1e35')
setInputVariable('inputFiles/'+DCfile,'outputFrequency','1')  # output frequency

# Make a local copy of material file, and modify that copy if necessary
materialFile='Zr_BMD19.txt';
materialFileTemplate=str(repository)+'/Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
#setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','full<111>{110}')
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate=str(repository)+'/Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dxMax','10') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1.0') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'Lmin','25')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','150')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'absoluteAreaThreshold','100')  # loop collapse area (in Burgers vector units ^2). Should be about dx^2
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','0')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'glideSolverType','none')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','Galerkin')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'crossSlipModel','0')  # crossSlipModel
setInputVariable('inputFiles/'+DDfile,'outputLoopLength','10')  # crossSlipModel

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate=str(repository)+'/Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVariable('inputFiles/'+elasticDeformatinoFile,'useElasticDeformationFEM','1')
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')

# Make a local copy of ClusterDynamics file, and modify that copy if necessary
clusterDynamicsFile='ClusterDynamics.txt';
clusterDynamicsFileTemplate=str(repository)+'/Library/ClusterDynamics/'+clusterDynamicsFile;
print("\033[1;32mCreating  clusterDynamicsFile\033[0m")
shutil.copy2(clusterDynamicsFileTemplate,'inputFiles/'+clusterDynamicsFile)
setInputVariable('inputFiles/'+clusterDynamicsFile,'useClusterDynamicsFEM','1')

# Create polycrystal.txt using local material file
meshFile='unitCube8188.msh';
#meshFile='poly10_periodicXYZ.msh';
meshFileTemplate=str(repository)+'/Library/Meshes/'+meshFile;
numberOfGrains = countGrains(meshFileTemplate)
print(f'Mesh contains {numberOfGrains} grains')
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=700;
pf.meshFile=meshFile
pf.boxScaling=np.array([1e-6,1e-6,1e-6])/b_SI # length of box edges in Burgers vector units
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
if numberOfGrains>1:
    pf.periodicFaceIDs=np.array([-1])
else:
    pf.periodicFaceIDs=np.array([])
pf.write('inputFiles')
# Keep the original grain-1 orientation and assign example rotations to the
# additional grains. These do not use the mesh's ElsetOrientations section.
rx = np.array([[1., 0., 0.], [0., 0.8, -0.6], [0., 0.6, 0.8]])
ry = np.array([[0.8, 0., 0.6], [0., 1., 0.], [-0.6, 0., 0.8]])
rz = np.array([[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]])
orientationPool = [np.eye(3), rx, rx.T, ry, ry.T, rz, rz.T,
                   rx @ ry, ry @ rz, rz @ rx]
with open('inputFiles/polycrystal.txt', 'a') as polycrystalFile:
    for grainID in range(2, numberOfGrains + 1):
        orientation = orientationPool[(grainID - 1) % len(orientationPool)] @ pf.C2G
        rows = [' '.join(f'{value:.17g}' for value in row) for row in orientation]
        polycrystalFile.write(f'C2G{grainID}=' + '\n     '.join(rows)
                             + '; # crystal-to-global rotation\n')
microstructureFiles = []

microstructureFileID=1;
microstructureName='frankLoopsDensity.txt';
microstructureFileTemplate=str(repository)+'/Library/Microstructures/'+microstructureName;
microstructureFileLocal='inputFiles/microstructureFile'+str(microstructureFileID)+'.txt'
microstructureFiles.append('microstructureFile'+str(microstructureFileID)+'.txt')
print("\033[1;32mCreating  microstructureFile"+str(microstructureFileID)+"\033[0m")
shutil.copy2(microstructureFileTemplate,microstructureFileLocal) # target filename is /dst/dir/file.ext
setInputVariable(microstructureFileLocal,'targetDensity_SI','1e13')
setInputVariable(microstructureFileLocal,'radiusDistributionMean_SI','5e-08')
setInputVariable(microstructureFileLocal,'radiusDistributionStd_SI','1e-09')
setInputVariable(microstructureFileLocal,'areVacancyLoops','1')
setInputVector(microstructureFileLocal,'allowedGrainIDs',np.array([-1]),'set of grain IDs where loops are allowed. Use -1 for all grains')
setInputVector(microstructureFileLocal,'allowedPlaneIDs',np.array([0]),'# set of plane IDs containing the loops. Use -1 for random plane')

microstructureFileID=2;
microstructureName='frankLoopsDensity.txt';
microstructureFileTemplate=str(repository)+'/Library/Microstructures/'+microstructureName;
microstructureFileLocal='inputFiles/microstructureFile'+str(microstructureFileID)+'.txt'
microstructureFiles.append('microstructureFile'+str(microstructureFileID)+'.txt')
print("\033[1;32mCreating  microstructureFile"+str(microstructureFileID)+"\033[0m")
shutil.copy2(microstructureFileTemplate,microstructureFileLocal) # target filename is /dst/dir/file.ext
setInputVariable(microstructureFileLocal,'targetDensity_SI','1e13')
setInputVariable(microstructureFileLocal,'radiusDistributionMean_SI','5e-08')
setInputVariable(microstructureFileLocal,'radiusDistributionStd_SI','1e-09')
setInputVariable(microstructureFileLocal,'areVacancyLoops','0')
setInputVector(microstructureFileLocal,'allowedGrainIDs',np.array([-1]),'set of grain IDs where loops are allowed. Use -1 for all grains')
setInputVector(microstructureFileLocal,'allowedPlaneIDs',np.array([1,2,3]),'# set of plane IDs containing the loops. Use -1 for random plane')

#microstructureFileID=3;
#microstructureName='frankLoopsIndividual.txt';
#microstructureFileTemplate=str(repository)+'/Library/Microstructures/'+microstructureName;
#microstructureFileLocal='inputFiles/microstructureFile'+str(microstructureFileID)+'.txt'
#microstructureFiles.append('microstructureFile'+str(microstructureFileID)+'.txt')
#print("\033[1;32mCreating  microstructureFile"+str(microstructureFileID)+"\033[0m")
#shutil.copy2(microstructureFileTemplate,microstructureFileLocal) # target filename is /dst/dir/file.ext
#setInputVector(microstructureFileLocal,'planeIDs',np.array([0,0,-1]),'set of grain IDs where loops are allowed. Use -1 for all grains')
#setInputVector(microstructureFileLocal,'loopRadii_SI',np.array([5e-08,5e-08,10e-08]),'set of grain IDs where loops are allowed. Use -1 for all grains')
#setInputVector(microstructureFileLocal,'loopSides',np.array([12,12,12]),'set of grain IDs where loops are allowed. Use -1 for all grains')
#setInputMatrix(microstructureFileLocal,'loopCenters',np.array([[0,0,-100],[0,0,100],[100,100,100]]))
#setInputVector(microstructureFileLocal,'isVacancyLoop',np.array([1,0,1]),'set of grain IDs where loops are allowed. Use -1 for all grains')


print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    for micrFile in microstructureFiles:
        initialMicrostructureFile.write('microstructureFile='+micrFile+';\n')
