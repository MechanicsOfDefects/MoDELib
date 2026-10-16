# Machine-learned mobility tutorial.
# Three prismatic loops on different slip systems in BCC Fe at 300 K, with the dislocation mobility given by a
# Gaussian-process regression model trained on MD (Library/Materials/Fe_320.txt, python/MobilityLaws/Fe320_GPR).
# The box is the MD simulation cell, loaded in compression along z at a constant strain rate.
import sys
sys.path.append("../../python/")
from modlibUtils import *
import numpy as np

# Create folder structure
folders=['evl','F','inputFiles']
for x in folders:
    if not os.path.exists(x):
        os.makedirs(x)

# Make a local copy of DefectiveCrystal parameters file and modify that copy if necessary
DCfile='DefectiveCrystal.txt'
DCfileTemplate='../../Library/DefectiveCrystal/'+DCfile
print("\033[1;32mCreating  DefectiveCrystalFile\033[0m")
shutil.copy2(DCfileTemplate,'inputFiles/'+DCfile)
setInputVariable('inputFiles/'+DCfile,'physics','DislocationDynamics ElasticDeformation')
setInputVariable('inputFiles/'+DCfile,'useFEM','0')
setInputVariable('inputFiles/'+DCfile,'Nsteps','80000')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'dtMax','0.5')  # max time step [b/cs]
setInputVariable('inputFiles/'+DCfile,'outputFrequency','100')  # output frequency
setInputVariable('inputFiles/'+DCfile,'outputBinary','0')
setInputVector('inputFiles/'+DCfile,'periodicImageSize',np.array([1,1,1]),'number of images for each periodic shift vector')
setInputVariable('inputFiles/'+DCfile,'EwaldLengthFactor','0')

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate='../../Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
# Time stepping
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dxMax','0.5') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0')
setInputVariable('inputFiles/'+DDfile,'velocityReductionFactor','0.75')
setInputVariable('inputFiles/'+DDfile,'useStochasticForce','0') # Langevin thermal noise disabled
setInputVariable('inputFiles/'+DDfile,'subcyclingBins','1 2 5 10 50 100') # subcycling, as in the LANL run (useSubCycling=1)
# Glide and climb solvers
setInputVariable('inputFiles/'+DDfile,'glideSolverType','Galerkin')
setInputVariable('inputFiles/'+DDfile,'climbSolverType','none')
# Dislocation fields and remeshing
setInputVariable('inputFiles/'+DDfile,'quadPerLength','0.1')  # number of quadrature points per unit length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'coreSize','2.0')
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','4.0') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','20')
setInputVariable('inputFiles/'+DDfile,'Lmin','5')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','25')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'relativeAreaThreshold','0.2')
# Junctions, cross slip and output
setInputVariable('inputFiles/'+DDfile,'maxJunctionIterations','1')
setInputVariable('inputFiles/'+DDfile,'crossSlipModel','1')
setInputVariable('inputFiles/'+DDfile,'crossSlipAngle_deg','30') # 30 instead of 10 so that near-screw dipole arms (20-30 deg) can cross slip and annihilate
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','0')

# Make a local copy of material file, and modify that copy if necessary
materialFile='Fe_320.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'mu0_SI','49e9')
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

t_dd2SI = 1.083093316853979e-13 # [s] time unit b/cs, value from the LANL run (it corresponds to mu0_SI=41e9)
StrainRate = 2.5e7*t_dd2SI # 2.5e7 1/s in code units

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformationFile='ElasticDeformation.txt';
elasticDeformationFileTemplate='../../Library/ElasticDeformation/'+elasticDeformationFile;
print("\033[1;32mCreating  elasticDeformationFile\033[0m")
shutil.copy2(elasticDeformationFileTemplate,'inputFiles/'+elasticDeformationFile)
setInputVector('inputFiles/'+elasticDeformationFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')
setInputVector('inputFiles/'+elasticDeformationFile,'ExternalStrainRate',np.array([0.0,0.0,-StrainRate,0.0,0.0,0.0]),'applied strain rate')
setInputVector('inputFiles/'+elasticDeformationFile,'stiffnessRatio',np.array([0.0,0.0,1e20,0.0,0.0,0.0]),'machine stiffness')

# Box dimensions of the MD simulation
XLength = 338.244807997994e-10/b_SI   # meters/b
YLength = 356.047166313678e-10/b_SI   # meters/b
ZLength = 1764.879600e-10/b_SI        # meters/b
print("XLength (in b): ", XLength, " YLength (in b): ", YLength, " ZLength (in b): ", ZLength)

# Create polycrystal.txt using local material file
meshFile='unitCube24.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
pf.grain1globalX1=np.array([3,1,0])
pf.grain1globalX3=np.array([0,0,1])
pf.boxEdges=np.array([[3,1,0],np.cross([0,0,1],[3,1,0]),[0,0,1]]) # i-th row is the direction of i-th box edge
pf.boxScaling=np.array([XLength,YLength,ZLength]) # length of box edges in Burgers vector units
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')

# Loop positions as fractions of the box (the box is centered at 0,0,0)
def frac_pos(L, f):
    return -L/2 + f*L

p30 = [frac_pos(XLength, 0.42), frac_pos(YLength, 0.55), frac_pos(ZLength, 0.13)]
p70 = [frac_pos(XLength, 0.63), frac_pos(YLength, 0.66), frac_pos(ZLength, 0.71)]

# Make a local copy of microstructure file, and modify that copy if necessary
microstructureFile1='prismaticLoopsIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile1;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1)
setInputVector('inputFiles/'+microstructureFile1,'slipSystemIDs',np.array([5,6,11]),'slip system IDs for each loop')
setInputVector('inputFiles/'+microstructureFile1,'loopRadii_SI',np.array([11e-9,11e-9,11e-9]),'radii of loops in SI')
setInputVector('inputFiles/'+microstructureFile1,'glideSteps',np.array([0,0,0]),'[b] distance from sessile loop (not important, used for debug)')
setInputMatrix('inputFiles/'+microstructureFile1,'loopCenters',np.array([[0,0,0], p30, p70]))

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
