import numpy as np
# Reader of the configuration files evl/evl_<n>.txt and of evl/cdNodes.txt.
#
# evl_<n>.txt starts with one count per block, one per line, followed by the blocks in the same order
# (DDconfigIO::writeTxtStream). Files written since the ellipsoidal inclusions were added have eleven
# counts; older files have ten (no ellipsoidal block). Both are read:
#   network nodes      sID  P(3)  V(3)  climbVelocityScalar  velocityReduction  meshLocation
#   loops              sID  B(3)  N(3)  P(3)  grainID  loopType  loopLength(4)  slippedArea
#   loop links         loopID  sourceID  sinkID  hasNetworkLink  meshLocation
#   loop nodes         loopID  sID  P(3)  networkNodeID  periodicShift(3)  edgeID1  edgeID2
#   spherical inclusions, ellipsoidal inclusions (eleven-count files only), polyhedron inclusions,
#   polyhedron inclusion nodes, polyhedron inclusion edges
#   displacement       one row per finite-element node of ElasticDeformation
#   cluster dynamics   one row per finite-element node: the mobile contents (defects per atom), then the
#                      number densities of the immobile families (per unit volume, in b^-3), then their
#                      contents (defects per atom)

class EVL:
    # network nodes
    nodes=np.empty([0,3])            # positions
    nodesV=np.empty([0,3])           # velocities
    nodesClimbVScalar=np.empty([0])  # climb velocity
    # loops
    loopsID=np.empty([0])
    burgers=np.empty([0,3])
    normal=np.empty([0,3])
    position=np.empty([0,3])
    loopsArea=np.empty([0])          # slipped area
    # loop links and loop nodes
    loopLinks=np.empty([0,5])
    loopNodesLoopID=np.empty([0])
    loopNodesID=np.empty([0])
    loopNodesPos=np.empty([0,3])
    # other blocks, as read
    sphericalInclusions=np.empty([0,0])
    ellipsoidalInclusions=np.empty([0,0])
    polyhedronInclusions=np.empty([0,0])
    polyhedronInclusionNodes=np.empty([0,0])
    polyhedronInclusionEdges=np.empty([0,0])
    displacement=np.empty([0,0])
    clusterDynamics=np.empty([0,0])

def readBlock(lines,start,count):
    if count==0:
        return np.empty([0,0]),start
    return np.loadtxt(lines[start:start+count],ndmin=2),start+count

def readEVLtxt(filename):
    with open(filename+'.txt', "r") as evlFile:
        lines=evlFile.read().splitlines()
    nCounts=0 # the header: the leading lines that hold a single integer
    while nCounts<len(lines) and len(lines[nCounts].split())==1 and lines[nCounts].strip().lstrip('-').isdigit():
        nCounts+=1
    if nCounts==10: # file written before the ellipsoidal inclusions: no ellipsoidal block
        counts=[int(lines[j]) for j in range(5)]+[0]+[int(lines[j]) for j in range(5,10)]
    elif nCounts>=11:
        nCounts=11
        counts=[int(lines[j]) for j in range(11)]
    else:
        raise ValueError(filename+'.txt: unexpected header of %d counts'%nCounts)
    evl=EVL()
    k=nCounts
    data,k=readBlock(lines,k,counts[0])
    if counts[0]>0:
        evl.nodes=data[:,1:4]
        evl.nodesV=data[:,4:7]
        evl.nodesClimbVScalar=data[:,7]
    data,k=readBlock(lines,k,counts[1])
    if counts[1]>0:
        evl.loopsID=data[:,0]
        evl.burgers=data[:,1:4]
        evl.normal=data[:,4:7]
        evl.position=data[:,7:10]
        evl.loopsArea=data[:,-1]
    data,k=readBlock(lines,k,counts[2])
    if counts[2]>0:
        evl.loopLinks=data
    data,k=readBlock(lines,k,counts[3])
    if counts[3]>0:
        evl.loopNodesLoopID=data[:,0]
        evl.loopNodesID=data[:,1]
        evl.loopNodesPos=data[:,2:5]
    evl.sphericalInclusions,k=readBlock(lines,k,counts[4])
    evl.ellipsoidalInclusions,k=readBlock(lines,k,counts[5])
    evl.polyhedronInclusions,k=readBlock(lines,k,counts[6])
    evl.polyhedronInclusionNodes,k=readBlock(lines,k,counts[7])
    evl.polyhedronInclusionEdges,k=readBlock(lines,k,counts[8])
    evl.displacement,k=readBlock(lines,k,counts[9])
    evl.clusterDynamics,k=readBlock(lines,k,counts[10])
    return evl

# Positions of the cluster-dynamics nodes (rows of evl.clusterDynamics), in units of b
def readCDNodes(folder):
    data=np.loadtxt(folder+'/evl/cdNodes.txt',ndmin=2)
    return data[:,1:4]
