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
setInputVariable('inputFiles/'+DCfile,'outputFrequency','100')  # output frequency
setInputVector('inputFiles/'+DCfile,'periodicImageSize',np.array([2,2,2]),'n of images in each direction')
setInputVariable('inputFiles/'+DCfile,'EwaldLengthFactor','1.0')

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate='../../Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dxMax','1') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','1') # don't filter velocity if noise is enabled
setInputVariable('inputFiles/'+DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1.0') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'Lmin','5')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','20')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'computeElasticEnergyPerLength','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'glideSolverType','Galerkin')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','none')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'Nsteps','1000000')  # number of simulation steps
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','1')
setInputVariable('inputFiles/'+DDfile,'quadPerLength','1.0')

# Make a local copy of material file, and modify that copy if necessary
materialFile='W.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStressRate',np.array([0.0,0.0,0.0,0.0,0.0,1e-9]),'applied strain rate')
setInputVector('inputFiles/'+elasticDeformatinoFile,'stiffnessRatio',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'stiffness ratio')

# Create polycrystal.txt using local material file
meshFile='unitCube24.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
pf.alignToSlipSystem0=1;
pf.boxScaling=np.array([500e-9,100e-9,500e-9])/b_SI # length of box edges in Burgers vector units
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')
microstructureFiles = []

# make a local copy of microstructure file, and modify that copy if necessary
microstructureFileID=1;
microstructureName='periodicDipoleIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureName;
microstructureFileLocal='inputFiles/microstructureFile'+str(microstructureFileID)+'.txt'
microstructureFiles.append('microstructureFile'+str(microstructureFileID)+'.txt')
print("\033[1;32mCreating  microstructureFile"+str(microstructureFileID)+"\033[0m")
shutil.copy2(microstructureFileTemplate,microstructureFileLocal) # target filename is /dst/dir/file.ext
setInputVector(microstructureFileLocal,'slipSystemIDs',np.array([0,-1]),'slip system IDs for each dipole')
setInputVector(microstructureFileLocal,'exitFaceIDs',np.array([4,4]),'4 is for edge, 2 for screw')
setInputMatrix(microstructureFileLocal,'dipoleCenters',np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]]))
setInputVector(microstructureFileLocal,'nodesPerLine',np.array([10,10]),'number of extra nodes on each dipole')
setInputVector(microstructureFileLocal,'dipoleHeights',np.array([120,1200]),'height of each dipole, in number of planes')
setInputVector(microstructureFileLocal,'glideSteps',np.array([1.0,30.0]),'step of each dipole in the glide plane')

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    for micrFile in microstructureFiles:
        initialMicrostructureFile.write('microstructureFile='+micrFile+';\n')
