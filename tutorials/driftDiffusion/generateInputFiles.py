import sys
sys.path.append("../../python/")
from modlibUtils import *

# Drift-diffusion (diffusion-bias) test: climb of a single prismatic <a> interstitial loop in Zr
# coupled to the spatial cluster dynamics of single vacancies.
# Ported from the DiffusionBias tutorials (drift, drift_E0, coGrowth, which shared these inputs).
#
# Requires MoDELib built with a single mobile species and no immobile species:
#   cmake .. -DMODELIB_CD_MSIZE=1 -DMODELIB_CD_ISIZE=0
#
# Drift-diffusion inputs in the material file (Zr_driftDiffusion.txt):
#   mobileSpeciesRelRelaxVol         sets the interaction energy E^e=-tr(sigma)*msVector*msRelRelaxVol*omega/3
#                                    (set it to 0 to switch off the stress-induced drift)
#   mobileSpeciesMigrationRelVolume  migration activation volumes for the stress-dependent diffusion tensor
#   lumpingFactor                    scaling of the lumped climb stiffness

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
setInputVariable('inputFiles/'+DCfile,'physics','DislocationDynamics ClusterDynamics ElasticDeformation')
setInputVariable('inputFiles/'+DCfile,'useFEM','1')  # required for spatial (FEM) cluster dynamics
setInputVariable('inputFiles/'+DCfile,'Nsteps','20')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'maxResolveSteps','0')
setInputVariable('inputFiles/'+DCfile,'dtMax','3.47935e30')  # max time step (in b/cs)
setInputVariable('inputFiles/'+DCfile,'outputFrequency','1')  # output frequency
setInputVector('inputFiles/'+DCfile,'periodicImageSize',np.array([1,1,1]),'n of images in each direction')
setInputVariable('inputFiles/'+DCfile,'EwaldLengthFactor','1')

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate='../../Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dxMax','1') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0')
setInputVariable('inputFiles/'+DDfile,'useStochasticForce','0')
setInputVariable('inputFiles/'+DDfile,'glideSolverType','none')  # the loop is sessile: climb only
setInputVariable('inputFiles/'+DDfile,'climbSolverType','Galerkin')
setInputVariable('inputFiles/'+DDfile,'quadPerLength','0.0001')  # quadrature points per unit length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1.0')
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','1')
setInputVariable('inputFiles/'+DDfile,'Lmin','450')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','1500')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')

# Make a local copy of material file, and modify that copy if necessary
materialFile='Zr_driftDiffusion.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# Make a local copy of ClusterDynamics file, and modify that copy if necessary
clusterDynamicsFile='ClusterDynamics.txt';
clusterDynamicsFileTemplate='../../Library/ClusterDynamics/'+clusterDynamicsFile;
print("\033[1;32mCreating  clusterDynamicsFile\033[0m")
shutil.copy2(clusterDynamicsFileTemplate,'inputFiles/'+clusterDynamicsFile)
setInputVariable('inputFiles/'+clusterDynamicsFile,'useClusterDynamicsFEM','1')  # 1=spatial FEM solver

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformationFile='ElasticDeformation.txt';
elasticDeformationFileTemplate='../../Library/ElasticDeformation/'+elasticDeformationFile;
print("\033[1;32mCreating  elasticDeformationFile\033[0m")
shutil.copy2(elasticDeformationFileTemplate,'inputFiles/'+elasticDeformationFile)
setInputVector('inputFiles/'+elasticDeformationFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')
setInputVector('inputFiles/'+elasticDeformationFile,'stiffnessRatio',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'stiffness ratio')

# Create polycrystal.txt using local material file
meshFile='unitCube_15K.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=625;
pf.meshFile=meshFile
pf.grain1globalX1=np.array([1,0,0])     # global x1 axis (C2G=identity, as in the DiffusionBias inputs)
pf.grain1globalX3=np.array([0,0,1])    # global x3 axis
pf.boxEdges=np.array([[1,0,0],[0,1,0],[0,0,1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([3093,3093,3093]) # 1 micron box, in Burgers vector units
pf.X0=np.array([0,0,0]) # Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([])
pf.write('inputFiles')

# Single prismatic <a> interstitial loop (100 Angstrom radius, 36 sides) at the center of the box.
# The loop is sessile, with normal parallel to its Burgers vector a1=[1 0 0] (HEX lattice, C2G=identity).
# Its plane x=1546.5 is a lattice plane perpendicular to a1 (plane spacing 1/2 b).
# The DiffusionBias tutorials generated this loop with the aLoop generator, which is not in this version of MoDELib.
microstructureFiles=[]
microstructureFileID=1;
microstructureFileLocal='inputFiles/microstructureFile'+str(microstructureFileID)+'.txt'
microstructureFiles.append('microstructureFile'+str(microstructureFileID)+'.txt')
loopCenter=np.array([1546.5,1546.5,1546.5])
loopRadius=100.0e-10/b_SI
loopSides=36
isVacancyLoop=0
theta=np.arange(loopSides)*2.0*np.pi/loopSides
loopPoints=np.zeros((loopSides,3))
loopPoints[:,0]=loopCenter[0]
loopPoints[:,1]=loopCenter[1]+loopRadius*np.cos(theta)
loopPoints[:,2]=loopCenter[2]+loopRadius*np.sin(theta)
loopNormal=np.array([1.0,0.0,0.0]) # right-handed normal of the point sequence above
loopBurgers=loopNormal if isVacancyLoop else -loopNormal # b.n>0 for vacancy loops, b.n<0 for interstitial loops
print("\033[1;32mCreating  microstructureFile"+str(microstructureFileID)+"\033[0m")
with open(microstructureFileLocal, "w") as mf:
    mf.write('type=PlanarLoop;\n')
    mf.write('style=individual;\n')
    mf.write('tag=aLoop;\n')
    mf.write('loopPoints=')
    for k in range(loopSides):
        mf.write(('' if k==0 else '           ')+'%.15e %.15e %.15e'%tuple(loopPoints[k])+(';\n' if k==loopSides-1 else '\n'))
    mf.write('burgers=%.15e %.15e %.15e; # Burgers vector of the loop\n'%tuple(loopBurgers))
    mf.write('normal=%.15e %.15e %.15e; # normal to the loop plane\n'%tuple(loopNormal))

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    for micrFile in microstructureFiles:
        initialMicrostructureFile.write('microstructureFile='+micrFile+';\n')
