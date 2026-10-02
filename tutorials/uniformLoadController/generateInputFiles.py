import sys
sys.path.append("../../python/")
from modlibUtils import *

# Create folder structure
folders=['evl','F','inputFiles']
for x in folders:
    if not os.path.exists(x):
        os.makedirs(x)

# Make a local copy of DefectiveCrystal parameters file and modify that copy if necessary
DCfile='DefectiveCrystal.txt'
DCfileTemplate='../../Library/DefectiveCrystal/'+DCfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DCfileTemplate,'inputFiles/'+DCfile)
setInputVariable('inputFiles/'+DCfile,'physics','DislocationDynamics ElasticDeformation')
setInputVariable('inputFiles/'+DCfile,'useFEM','0')
setInputVariable('inputFiles/'+DCfile,'Nsteps','10000')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'maxResolveSteps','0')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'dtMax','1e25')
setInputVariable('inputFiles/'+DCfile,'outputFrequency','10')  # output frequency

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate='../../Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dtMax','1e25')
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
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
#setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','full<111>{110}')
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.01,0.0,0.0,0.0]),'applied stress')

# Create polycrystal.txt using local material file
meshFile='unitCube24.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
#pf.alignToSlipSystem0=1
#pf.grain1globalX1=np.array([1,1,0])
#pf.grain1globalX3=np.array([0,0,1])
#pf.boxEdges=np.array([[1,1,0],[0,1,0],[0,0,1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([1e-6,1e-6,1e-6])/b_SI # length of box edges in Burgers vector units
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')
microstructureFiles = []

microstructureFileID=1;
microstructureName='prismaticLoopsDensity.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureName;
microstructureFileLocal='inputFiles/microstructureFile'+str(microstructureFileID)+'.txt'
microstructureFiles.append('microstructureFile'+str(microstructureFileID)+'.txt')
print("\033[1;32mCreating  microstructureFile"+str(microstructureFileID)+"\033[0m")
shutil.copy2(microstructureFileTemplate,microstructureFileLocal) # target filename is /dst/dir/file.ext
setInputVariable(microstructureFileLocal,'targetDensity_SI','1e13')
setInputVariable(microstructureFileLocal,'radiusDistributionMean_SI','1e-07')
setInputVariable(microstructureFileLocal,'radiusDistributionStd_SI','0e-07')
setInputVector(microstructureFileLocal,'allowedGrainIDs',np.array([-1]),'set of grain IDs where loops are allowed. Use -1 for all grains')
setInputVector(microstructureFileLocal,'allowedSlipSystemIDs',np.array([-1]),'set of slip system IDs whose Burgers vector are allowed to be the prism axis. Use -1 for all slip systems')

microstructureFileID=2;
microstructureName='frankLoopsDensity.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureName;
microstructureFileLocal='inputFiles/microstructureFile'+str(microstructureFileID)+'.txt'
microstructureFiles.append('microstructureFile'+str(microstructureFileID)+'.txt')
print("\033[1;32mCreating  microstructureFile"+str(microstructureFileID)+"\033[0m")
shutil.copy2(microstructureFileTemplate,microstructureFileLocal) # target filename is /dst/dir/file.ext
setInputVariable(microstructureFileLocal,'targetDensity_SI','1e13')
setInputVariable(microstructureFileLocal,'radiusDistributionMean_SI','5e-08')
setInputVariable(microstructureFileLocal,'radiusDistributionStd_SI','1e-09')
setInputVariable(microstructureFileLocal,'areVacancyLoops','0')
setInputVector(microstructureFileLocal,'allowedGrainIDs',np.array([-1]),'set of grain IDs where loops are allowed. Use -1 for all grains')
setInputVector(microstructureFileLocal,'allowedPlaneIDs',np.array([-1]),'# set of plane IDs containing the loops. Use -1 for random plane')

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    for micrFile in microstructureFiles:
        initialMicrostructureFile.write('microstructureFile='+micrFile+';\n')
