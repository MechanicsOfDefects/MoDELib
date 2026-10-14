import sys
sys.path.append("../../python/")
from modlibUtils import *

# Cluster dynamics of irradiated zirconium with nucleation, dissolution and coalescence of the loops.
# The immobile clusters are integrated implicitly by CVODE: MoDELib must be built with SUNDIALS.
# Model and parameters: docs/reports/GW_Phase4_D1M1.pdf, material file Zr_CD4opt.txt (fitted at 573 K).
# The species counts are the defaults of the build (4 mobile, 8 immobile).

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
setInputVariable('inputFiles/'+DCfile,'physics','ClusterDynamics')
setInputVariable('inputFiles/'+DCfile,'useFEM','1')  # required for spatial (FEM) cluster dynamics
setInputVariable('inputFiles/'+DCfile,'Nsteps','30')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'maxResolveSteps','0')  # no re-solve needed with cluster dynamics alone
setInputVariable('inputFiles/'+DCfile,'dtMax','6.958689643537561e18')  # time step (in b/cs): 0.1 dpa at 1e-7 dpa/s. The mobile species are solved once per step
setInputVariable('inputFiles/'+DCfile,'outputFrequency','1')  # output frequency

# Make a local copy of material file, and modify that copy if necessary
materialFile='Zr_CD4opt.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)

# Make a local copy of ClusterDynamics file, and modify that copy if necessary
clusterDynamicsFile='ClusterDynamics.txt';
clusterDynamicsFileTemplate='../../Library/ClusterDynamics/'+clusterDynamicsFile;
print("\033[1;32mCreating  clusterDynamicsFile\033[0m")
shutil.copy2(clusterDynamicsFileTemplate,'inputFiles/'+clusterDynamicsFile)
setInputVariable('inputFiles/'+clusterDynamicsFile,'useClusterDynamicsFEM','1')  # 1=spatial FEM solver, 0=uniform controllers
setInputVariable('inputFiles/'+clusterDynamicsFile,'immobileIntegrator','cvode')  # implicit integration of the immobile species

# Create polycrystal.txt using local material file
meshFile='unitCube_15K.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=573;
pf.meshFile=meshFile
pf.grain1globalX1=np.array([1,0,0])     # global x1 axis. Overwritten if alignToSlipSystem0=true
pf.grain1globalX3=np.array([0,0,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.boxEdges=np.array([[1,0,0],[0,1,0],[0,0,1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([3093,3093,3093]) # must be a vector of integers: 1 micron
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([])
pf.write('inputFiles')

# No discrete dislocation microstructure for a cluster-dynamics-only run
print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('')

# To continue this study in a new folder, copy the last configuration evl/evl_<n>.txt of this one
# to evl/evl_0.txt of the new folder after running microstructureGenerator there: the fields of
# the mobile and immobile species are read from it.
