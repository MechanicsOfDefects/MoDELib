import numpy as np
# Reader of the auxiliary files evl/ddAux_<n>.txt (DDauxIO::writeTxt).
# Three counts, then the blocks:
#   mesh nodes         nodeID  displacement(3)
#   quadrature points  sourceID sinkID qID r(3) j rl(3) stress(9) pkForce(3) stackingFaultForce(3)
#                      lineTensionForce(3) velocity(3) elasticEnergyPerLength coreEnergyPerLength cCD(mSize) cDD(mSize)
#   periodic patches   glidePlaneID  patchBoundaryID  shift(3)

class AUX:
    meshNodes=np.empty([0,4])
    gaussPoints=np.empty([0,0])
    sourceID=np.empty(0)
    sinkID=np.empty(0)
    qID=np.empty(0)
    r=np.empty([0,3])
    j=np.empty(0)
    rl=np.empty([0,3])
    stress=np.empty([0,3,3])
    pkForce=np.empty([0,3])
    stackingFaultForce=np.empty([0,3])
    lineTensionForce=np.empty([0,3])
    velocity=np.empty([0,3])
    elasticEnergyPerLength=np.empty(0)
    coreEnergyPerLength=np.empty(0)
    cCD=np.empty([0,0])  # cluster-dynamics concentrations of the mobile species at the quadrature point
    cDD=np.empty([0,0])  # concentrations in equilibrium with the dislocation
    periodicPatches=np.empty([0,5])

def readAUXtxt(filename):
    with open(filename+'.txt', "r") as auxFile:
        lines=auxFile.read().splitlines()
    numNodes=int(lines[0])
    numGPs=int(lines[1])
    numPGPP=int(lines[2])
    aux=AUX()
    k=3
    if numNodes>0:
        aux.meshNodes=np.loadtxt(lines[k:k+numNodes],ndmin=2)
    k+=numNodes
    if numGPs>0:
        data=np.loadtxt(lines[k:k+numGPs],ndmin=2)
        mSize=(data.shape[1]-33)//2  # number of mobile species of the build
        aux.gaussPoints=data
        aux.sourceID=data[:,0]
        aux.sinkID=data[:,1]
        aux.qID=data[:,2]
        aux.r=data[:,3:6]
        aux.j=data[:,6]
        aux.rl=data[:,7:10]
        aux.stress=data[:,10:19].reshape(-1,3,3)
        aux.pkForce=data[:,19:22]
        aux.stackingFaultForce=data[:,22:25]
        aux.lineTensionForce=data[:,25:28]
        aux.velocity=data[:,28:31]
        aux.elasticEnergyPerLength=data[:,31]
        aux.coreEnergyPerLength=data[:,32]
        aux.cCD=data[:,33:33+mSize]
        aux.cDD=data[:,33+mSize:33+2*mSize]
    k+=numGPs
    if numPGPP>0:
        aux.periodicPatches=np.loadtxt(lines[k:k+numPGPP],ndmin=2)
    return aux
