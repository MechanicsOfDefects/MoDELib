"""Generate the defect-free bicrystal_48 tutorial for DDqt."""
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
setInputVariable('inputFiles/'+DCfile,'physics','DislocationDynamics ElasticDeformation')
setInputVariable('inputFiles/'+DCfile,'useFEM','0')
setInputVariable('inputFiles/'+DCfile,'Nsteps','0')  # visualization only; no simulation steps
setInputVariable('inputFiles/'+DCfile,'maxResolveSteps','0')  # number of simulation steps
setInputVariable('inputFiles/'+DCfile,'dtMax','1e25')
setInputVariable('inputFiles/'+DCfile,'outputFrequency','1')  # output frequency

setInputVariable('inputFiles/'+DCfile,'startAtTimeStep','0')

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

# Create polycrystal.txt using local material file
meshFile='bicrystal_48.msh';
meshFileTemplate=str(repository)+'/Library/Meshes/'+meshFile;
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
# Grain 2's rational rotation (cos=4/5, sin=3/5) requires the
# axial BCC repeat counts to be multiples of five for commensurate periods.
# Preserve the intended approximately 1-micrometer scale in both grains.
pf.compute()
repeatCounts = 5 * np.maximum(1, np.round(pf.boxScaling / (5 * pf.boxEdgesLatticeLengths)))
pf.boxScaling = repeatCounts * pf.boxEdgesLatticeLengths
pf.write('inputFiles')
# The mesh contains grains 1 and 2. Grain 2 is rotated 36.87 degrees about Y.
with open('inputFiles/polycrystal.txt', 'a') as polycrystalFile:
    polycrystalFile.write('C2G2=0.8 0 0.6\n     0 1 0\n     -0.6 0 0.8; # grain 2 crystal-to-global rotation\n')

# No dislocations, inclusions, or other defect microstructures.
microstructureFiles = []
print("Creating initialMicrostructureFile (empty)")
with open('inputFiles/initialMicrostructure.txt', 'w') as initialMicrostructureFile:
    for micrFile in microstructureFiles:
        initialMicrostructureFile.write('microstructureFile='+micrFile+';\n')

# Supply empty frame 0 so the mesh can be opened directly in DDqt.
# DDconfigIO has eleven block counts; DDauxIO has three.
Path('evl/evl_0.txt').write_text('0\n' * 11)
Path('evl/ddAux_0.txt').write_text('0\n' * 3)
print('Ready for DDqt: ' + str(Path.cwd()))
