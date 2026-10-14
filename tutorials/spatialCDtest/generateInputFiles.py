import sys
sys.path.append("../../python/")
from modlibUtils import *

# Cluster-dynamics-only input generator (current MoDELib input layout, modeled on tutorials/annealing)

# Create folder structure
folders=['evl','F','inputFiles']
for x in folders:
    if not os.path.exists(x):
        os.makedirs(x)

# Make a local copy of DefectiveCrystal parameters file and modify that copy if necessary
DCfile='DefectiveCrystal.txt'
DCfileTemplate='../../Library/DefectiveCrystal/'+DCfile
print("\033[1;32mCreating  DCfile\033[0m")
shutil.copy2(DCfileTemplate,'inputFiles/'+DCfile)
setInputVariable('inputFiles/'+DCfile,'physics','ClusterDynamics')  # add 'ElasticDeformation' and/or 'DislocationDynamics' to couple
setInputVariable('inputFiles/'+DCfile,'useFEM','1')  # required for spatial (FEM) cluster dynamics
setInputVariable('inputFiles/'+DCfile,'Nsteps','120')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'maxResolveSteps','0')  # no re-solve needed with cluster dynamics alone
setInputVariable('inputFiles/'+DCfile,'dtMax','3.4e18')  # time step (in b/cs); CD-only runs always use dtMax
setInputVariable('inputFiles/'+DCfile,'outputFrequency','1')  # output frequency

# Make a local copy of material file, and modify that copy if necessary
materialFile='Zr4_Fitted.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','<a>{basal} <a>{prismatic}')

# Make a local copy of ClusterDynamics file, and modify that copy if necessary
clusterDynamicsFile='ClusterDynamics.txt';
clusterDynamicsFileTemplate='../../Library/ClusterDynamics/'+clusterDynamicsFile;
print("\033[1;32mCreating  clusterDynamicsFile\033[0m")
shutil.copy2(clusterDynamicsFileTemplate,'inputFiles/'+clusterDynamicsFile)
setInputVariable('inputFiles/'+clusterDynamicsFile,'useClusterDynamicsFEM','1')  # 1=spatial FEM solver, 0=uniform controllers

# Create polycrystal.txt using local material file
meshFile='unitCube_15K.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=553;
pf.meshFile=meshFile
pf.grain1globalX1=np.array([1,0,0])     # global x1 axis. Overwritten if alignToSlipSystem0=true
pf.grain1globalX3=np.array([0,0,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.boxEdges=np.array([[1,0,0],[0,1,0],[0,0,1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([3093,3093,3093]) # must be a vector of integers
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([])
pf.write('inputFiles')

# No discrete dislocation microstructure for a cluster-dynamics-only run
print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('')
