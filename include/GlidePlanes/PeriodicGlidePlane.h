/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlane_H_
#define model_PeriodicGlidePlane_H_


#include <memory>
#include <algorithm>
#include <limits>
#include <PolygonPatchIntersection.h>
#include <string>
#include <list>

#include <GlidePlane.h>
#include <GlidePlaneFactory.h>
//#include <PeriodicGlidePlaneFactory.h>
#include <TerminalColors.h>
#include <DiophantineSolver.h>
#include <SharedPtrFactories.h>
#include <GramSchmidt.h>


namespace model
{

template<int dim>
struct PeriodicGlidePlaneFactory;

template<int dim>
struct PeriodicPlanePatch;

template<int dim>
class PeriodicPlaneNode;

template<int dim>
struct PeriodicPatchBoundary;

template<int dim>
class PeriodicGlidePlane;

//    template<int dim>
//    class PeriodicPlaneEdge;

/**********************************************************************/
/**********************************************************************/
template<int dim>
struct PeriodicPlaneEdge
{
    
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    const PeriodicPlanePatch<dim>* const patch;
    const std::shared_ptr<PeriodicPlaneNode<dim>> source;
    const std::shared_ptr<PeriodicPlaneNode<dim>>   sink;
    const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection;
    const short int edgeID;
    const VectorDim deltaShift;
    
    PeriodicPlaneEdge<dim>* next;
    PeriodicPlaneEdge<dim>* prev;
    PeriodicPlaneEdge<dim>* twin;
    
    //static VectorDim getDeltaShift(const PeriodicPlanePatch<dim>* const patch,const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection);
    
    /**********************************************************************/
    PeriodicPlaneEdge(const PeriodicPlanePatch<dim>* const patch_in,
                      const std::shared_ptr<PeriodicPlaneNode<dim>>& source_in,
                      const std::shared_ptr<PeriodicPlaneNode<dim>>& sink_in,
                      const std::shared_ptr<const MeshBoundarySegment<dim>> meshIntersection_in,
                      const short int&)
    ;
    ~PeriodicPlaneEdge();
    
    std::string tag() const;
};

/**********************************************************************/
/**********************************************************************/
template<int dim>
struct NodalConnectivity
{
    typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
    
    PeriodicPlaneEdgeType* inEdge;
    PeriodicPlaneEdgeType* outEdge;
    
    NodalConnectivity() :
    /* init */ inEdge(nullptr)
    /* init */,outEdge(nullptr)
    {
        
    }
    
};

/**********************************************************************/
/**********************************************************************/
template<int dim>
class PeriodicPlaneNode : public StaticID<PeriodicPlaneNode<dim>>
/*                     */,public Eigen::Matrix<double,dim-1,1>
{
    
public:
    
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    typedef PeriodicPlaneNode<dim> PeriodicPlaneNodeType;
    typedef PeriodicPlaneEdge<dim> PeriodicPlaneEdgeType;
    typedef NodalConnectivity<dim> NodalConnectivityType;
    typedef std::map<size_t,NodalConnectivity<dim>> NodalConnectivityContainerType;
    typedef std::set<PeriodicPlaneEdgeType*> InOutEdgeContainerType;
    typedef VectorLowerDim KeyType;
    typedef CompareVectorsByComponent<double,dim-1,float> CompareType;
    
private:
    
    NodalConnectivityContainerType _patchConnectivities;
    NodalConnectivityContainerType _neighborConnectivities;
    
    
public:
    
    PeriodicPlaneNode(PeriodicPatchBoundary<dim>* const,const VectorLowerDim& pos);
    void addLink(PeriodicPlaneEdgeType* const link);
    void removeLink(PeriodicPlaneEdgeType* const link);
    const NodalConnectivityContainerType& patchConnectivities() const;
    NodalConnectivityContainerType& patchConnectivities();
    const NodalConnectivityContainerType& neighborConnectivities() const;
    NodalConnectivityContainerType& neighborConnectivities();
    InOutEdgeContainerType inEdges() const;
    InOutEdgeContainerType outEdges() const;
};

/**********************************************************************/
/**********************************************************************/
template<int dim>
struct PeriodicPlanePatch : public StaticID<PeriodicPlanePatch<dim>>
/*                      */, private std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>>
{
    
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef VectorDim KeyType;
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    typedef CompareVectorsByComponent<double,dim,float> CompareType;
    typedef std::vector<std::shared_ptr<PeriodicPlaneEdge<dim>>> PeriodicPlaneEdgeContainerType;
    
    PeriodicPatchBoundary<dim>* const patchBoundary;
    const VectorDim shift;
    const std::shared_ptr<GlidePlane<dim>> glidePlane;
    
    PeriodicPlanePatch(PeriodicPatchBoundary<dim>* const,const VectorDim&);
    ~PeriodicPlanePatch();
    static bool isRightHandedBoundary(const BoundingMeshSegments<dim>&, const Plane<dim>&);
    void addMeshIntersections(const BoundingMeshSegments<dim>&);
    const PeriodicPlaneEdgeContainerType& edges() const;
    std::shared_ptr<PeriodicPlaneEdge<dim>> edgeOnFaces(const std::set<const Polyhedron<dim,dim-1>*>& faces) const;
    int contains(const VectorLowerDim& test) const;
};





template<int dim>
struct PeriodicPatchBoundary : public StaticID<PeriodicPatchBoundary<dim>>
/*                          */,public KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlaneNode<dim>,typename PeriodicPlaneNode<dim>::CompareType>
/*                          */,private std::set<const PeriodicPlaneEdge<dim>*>
/*                          */,public KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlanePatch<dim>,typename PeriodicPlanePatch<dim>::CompareType>

{
    typedef Eigen::Matrix<double,dim,dim> MatrixDim;
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    typedef KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlaneNode<dim>,typename PeriodicPlaneNode<dim>::CompareType> NodeCointainerType;
    typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
    typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
    typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
    typedef KeyConstructableWeakPtrFactory<PeriodicPatchBoundary<dim>,PeriodicPlanePatch<dim>,typename PeriodicPlanePatch<dim>::CompareType> PatchContainerType;
    
    BoundariesContainerType _outerBoundaries;
    BoundariesContainerType _innerBoundaries;
    GlidePlaneFactory<dim>& glidePlaneFactory;
    const std::shared_ptr<GlidePlane<dim>> referencePlane;
    //        const MatrixDim L2G;
    
    
    //        static MatrixDim getL2G(VectorDim z);
    PeriodicPatchBoundary(GlidePlaneFactory<dim>&, const std::shared_ptr<GlidePlane<dim>>& referencePlane_in);
    //        GlidePlaneKey<dim> getGlidePlaneKey(const VectorDim& shift);
    //        std::shared_ptr<GlidePlane<dim>> getGlidePlane(const VectorDim& shift);
    void createNewBoundary(const PeriodicPlaneEdge<dim>* currentEdge,UntwinnedEdgeContainerType& untwinnedCopy);
    void updateBoundaries();
    BoundariesContainerType& outerBoundaries();
    const BoundariesContainerType& outerBoundaries() const;
    BoundariesContainerType& innerBoundaries();
    const BoundariesContainerType& innerBoundaries() const;
    const UntwinnedEdgeContainerType& untwinnedEdges() const;
    UntwinnedEdgeContainerType& untwinnedEdges();
    void addUntwinnedEdge(const PeriodicPlaneEdge<dim>* link);
    void removeUntwinnedEdge(const PeriodicPlaneEdge<dim>* link);
    NodeCointainerType& nodes();
    const NodeCointainerType& nodes() const;
    int isInsideOuterBoundary(const VectorLowerDim& test);
    bool isCompact() const;
    static VectorDim rightHandedNormal(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
    double rightHandedArea(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
    bool isRightHandedBoundary(const std::vector<const PeriodicPlaneEdge<dim>*>& bnd);
    VectorLowerDim getLocalPosition(const VectorDim& point,const VectorDim& shift) const;
    VectorDim getGlobalPosition(const VectorLowerDim& point) const;
    std::shared_ptr<PeriodicPlaneNode<dim> > getSharedNode(const VectorDim& pointDim,const VectorDim& shift);
    std::shared_ptr<PeriodicPlanePatch<dim>> getPatch(const VectorDim& shift);
    GlidePlaneKey<dim> getGlidePlaneKey(const VectorDim& shift);
    std::shared_ptr<GlidePlane<dim>> getGlidePlane(const VectorDim& shift);
    const PatchContainerType& patches() const;
    PatchContainerType& patches();
    
    
};




/**********************************************************************/
/**********************************************************************/
// ,
template<int dim>
class PeriodicGlidePlane : public PeriodicPatchBoundary<dim>
//    /*                      */,public KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>,PeriodicPlanePatch<dim>> // container of patches
{
    
    
public:
    
    typedef Eigen::Matrix<double,dim,dim> MatrixDim;
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef Eigen::Matrix<long int, dim, 1> VectorDimI;
    typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
    //        typedef KeyConstructableSharedPtrFactory<PeriodicGlidePlane<dim>,PeriodicPlanePatch<dim>> PatchContainerType;
    typedef std::set<const PeriodicPlaneEdge<dim>*> UntwinnedEdgeContainerType;
    typedef std::vector<const PeriodicPlaneEdge<dim>*> BoundaryContainerType;
    typedef std::vector<BoundaryContainerType>    BoundariesContainerType;
    typedef PeriodicGlidePlaneFactory<dim> PeriodicGlidePlaneFactoryType;
    typedef GlidePlane<dim> GlidePlaneType;
    typedef typename GlidePlaneType::KeyType KeyType;
    
    GlidePlaneFactory<dim>& glidePlaneFactory;
    PeriodicGlidePlaneFactoryType& periodicGlidePlaneFactory;
    const std::shared_ptr<GlidePlane<dim>> referencePlane;
    
    
    PeriodicGlidePlane(PeriodicGlidePlaneFactoryType* const pgpf,const KeyType& key_in);
    void print();
    //        VectorDim getGlidePlaneShiftfromReferencePlane(const GlidePlane<dim> *gp) const;
    
    template<typename T>
    std::vector<PolygonPatchIntersectionRecord<dim,T>> polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints,const bool& intersectInternalNodes=false);
    template<typename T>
    std::vector<PolygonPatchIntersectionRecord<dim,T>> polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim,const T* const>>& polyPoints,const bool& intersectInternalNodes=false);
    
    
    std::pair<VectorDim,bool> findPatch(const VectorLowerDim&,const VectorDim&);
    
    std::vector<std::shared_ptr<PeriodicPlanePatch<dim>>> filledPatches(const std::vector<VectorDim>& patchShifts);
    
};



template<int dim>
struct PeriodicPlanePatchIO
{
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    size_t glidePlaneID;
    size_t patchBoundaryID;
    VectorDim shift;
    
    PeriodicPlanePatchIO();
    PeriodicPlanePatchIO(const PeriodicPlanePatch<dim>& patch);
    PeriodicPlanePatchIO(std::stringstream& ss);
    //        template <class T>
    //        friend T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds);
    
};

template <int dim,class T>
T& operator << (T& os, const PeriodicPlanePatchIO<dim>& ds)
{
    os  << ds.glidePlaneID<<"\t"
    /**/<< ds.patchBoundaryID<<"\t"
    /**/<< std::setprecision(15)<<std::scientific<<ds.shift.transpose();
    return os;
}

template<int dim>
template<typename T>
std::vector<PolygonPatchIntersectionRecord<dim,T>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorDim,const T* const>>& polyPoints, const bool& intersectInternalNodes)
{
    std::vector<std::pair<VectorLowerDim,const T* const>> lowerPolyPoints;
    for(const auto& point : polyPoints)
    {
        assert(this->referencePlane->contains(point.first) && "reference plane does not cointain point");
        lowerPolyPoints.push_back(std::make_pair(referencePlane->localPosition(point.first-VectorDim::Zero()),point.second));
    }
    return polygonPatchIntersection(lowerPolyPoints,intersectInternalNodes);
}

template<int dim>
template<typename T>
std::vector<PolygonPatchIntersectionRecord<dim,T>> PeriodicGlidePlane<dim>::polygonPatchIntersection(const std::vector<std::pair<VectorLowerDim, const T* const>>& polyPoints, const bool& intersectInternalNodes)
{
    using Edge = PeriodicPlaneEdge<dim>;
    using Patch = PeriodicPlanePatch<dim>;
    using EdgeIDs = std::pair<short int, short int>;
    using EdgeMap = std::map<VectorDim, std::set<EdgeIDs>, CompareVectorsByComponent<double,dim,float>>;
    using Record = PolygonPatchIntersectionRecord<dim,T>;
    using Metadata = PolygonPatchMetadata<T>;
    using Arc = PolygonPatchArcKind;
    std::vector<Record> result;
    if(polyPoints.size() < 3)
    {
        return result;
    }
    // Always intersect geometry: equal endpoint patches do not imply that the
    // intervening polygon edge stays inside a non-convex periodic domain.
    (void)intersectInternalNodes;
    const auto cross = [](const VectorLowerDim& a, const VectorLowerDim& b)
    {
        return a[0]*b[1]-a[1]*b[0];
    };
    double scale = 1.0;
    double signedArea = 0.0;
    for(size_t i=0; i<polyPoints.size(); ++i)
    {
        const VectorLowerDim a = polyPoints[i].first-polyPoints.front().first;
        const VectorLowerDim b = polyPoints[(i+1)%polyPoints.size()].first-polyPoints.front().first;
        scale = std::max(scale,a.norm());
        signedArea += cross(a,b);
    }
    const double tolerance = 1.e-9*scale;
    // Junction contraction can leave distinct node identities at coincident
    // positions. Use the same edge tolerance as the walk below: fewer than
    // three retained edges cannot enclose a slipped area, even when roundoff
    // leaves a nonzero signed area for the original nearly collapsed triangle.
    size_t retainedEdges=0;
    for(size_t i=0; i<polyPoints.size(); ++i)
    {
        if((polyPoints[(i+1)%polyPoints.size()].first-polyPoints[i].first).norm()>tolerance)
        {
            ++retainedEdges;
        }
    }
    if(retainedEdges<3)
    {
        return result;
    }
    // Zero algebraic area is not a collapse proof (opposite lobes cancel).
    // Only discard a contour confined to a line. For cancellation, a local
    // oriented turn supplies the boundary orientation convention; successor
    // selection still uses intersection incidence and signed winding layers.
    VectorLowerDim longest=VectorLowerDim::Zero();
    for(const auto& vertex : polyPoints)
    {
        const VectorLowerDim d=vertex.first-polyPoints.front().first;
        if(d.squaredNorm()>longest.squaredNorm())
        {
            longest=d;
        }
    }
    double orientation=0.0;
    for(const auto& vertex : polyPoints)
    {
        const double turn=cross(longest,vertex.first-polyPoints.front().first);
        if(std::abs(turn)>tolerance*longest.norm())
        {
            orientation=turn;
            break;
        }
    }
    if(orientation==0.0)
    {
        return result;
    }
    const bool forward = (std::abs(signedArea)>tolerance*tolerance ? signedArea : orientation)>0.0;
    PeriodicPatchBoundary<dim> boundary(glidePlaneFactory,referencePlane);
    std::vector<std::shared_ptr<Patch>> heldPatches;
    const auto getPatch = [&](const VectorDim& shift)
    {
        auto patch = boundary.getPatch(shift);
        heldPatches.push_back(patch);
        return patch;
    };
    const auto traversable = [](const Edge* edge)
    {
        for(const auto* face : edge->meshIntersection->faces)
        {
            if(!face->periodicFacePair.second ||
               face->regionIDs != face->periodicFacePair.second->regionIDs)
            {
                return false;
            }
        }
        return edge->deltaShift.squaredNorm() > 0.0;
    };
    const auto crossSeam = [&](const Edge* edge)
    {
        auto patch = getPatch(edge->patch->shift+edge->deltaShift);
        if(!edge->twin)
        {
            throw std::runtime_error("Periodic polygon walk: same-grain edge has no twin after creating its neighbor.");
        }
        return patch;
    };
    const auto edgeIDs = [&](const Patch* patch, const VectorLowerDim& point)
    {
        EdgeIDs ids(-1,-1);
        for(const auto& edge : patch->edges())
        {
            const VectorLowerDim d = *edge->sink-*edge->source;
            const double t = d.dot(point-*edge->source)/d.squaredNorm();
            if(t>=-tolerance/d.norm() && t<=1.0+tolerance/d.norm() &&
               std::abs(cross(d,point-*edge->source))<=tolerance*d.norm())
            {
                if(ids.first<0)
                {
                    ids.first=edge->edgeID;
                }
                else
                {
                    ids.second=edge->edgeID;
                }
            }
        }
        return ids;
    };
    const auto originalSegment = [&](size_t index)
    {
        const size_t next=(index+1)%polyPoints.size();
        return PolygonOriginalSegment<T>{index,next,polyPoints[index].second,polyPoints[next].second};
    };
    const auto metadata = [&](const VectorLowerDim& point, std::optional<size_t> segment,
                              std::optional<size_t> vertex, Arc incoming,
                              std::optional<size_t> incomingSegment, bool seam)
    {
        Metadata value;
        value.originalVertex=vertex;
        value.incomingArc=incoming;
        value.crossesPeriodicSeam=seam;
        if(segment)
        {
            const auto origin=originalSegment(*segment);
            const VectorLowerDim d=polyPoints[origin.sinkIndex].first-polyPoints[origin.sourceIndex].first;
            const double t=d.squaredNorm()>0.0 ? d.dot(point-polyPoints[origin.sourceIndex].first)/d.squaredNorm() : 0.0;
            value.polygonLocation=PolygonSegmentLocation<T>{origin,std::clamp(t,0.0,1.0)};
        }
        if(incomingSegment)
        {
            value.incomingOriginalSegment=originalSegment(*incomingSegment);
        }
        return value;
    };
    const auto classify = [](Metadata& value)
    {
        if(value.originalVertex)
        {
            value.kind=PolygonPatchNodeKind::OriginalVertex;
        }
        else if(value.onPhysicalBoundary)
        {
            value.kind=value.polygonLocation ? PolygonPatchNodeKind::PhysicalBoundaryIntersection
                                            : PolygonPatchNodeKind::PhysicalBoundaryCorner;
        }
        else
        {
            value.kind=PolygonPatchNodeKind::PeriodicCrossing;
        }
    };
    const auto append = [&](const VectorLowerDim& point, const Patch* patch, const T* original, Metadata value)
    {
        const EdgeIDs geometricIDs=edgeIDs(patch,point);
        value.atPatchCorner=geometricIDs.second>=0;
        for(short id : {geometricIDs.first,geometricIDs.second})
        {
            if(id>=0 && !traversable(patch->edges()[id].get()))
            {
                value.onPhysicalBoundary=true;
            }
        }
        const EdgeIDs ids = original ? EdgeIDs(-1,-1) : geometricIDs;
        const auto samePort=[&](const Metadata& previous)
        {
            if((!previous.originalVertex && !previous.polygonLocation) ||
               (!value.originalVertex && !value.polygonLocation))
            {
                return true; // A corner annotation adds no distinct polygon port.
            }
            if(previous.originalVertex && value.originalVertex)
            {
                return previous.originalVertex==value.originalVertex;
            }
            return previous.polygonLocation && value.polygonLocation &&
                previous.polygonLocation->segment.sourceIndex==value.polygonLocation->segment.sourceIndex &&
                std::abs(previous.polygonLocation->parameter-value.polygonLocation->parameter)<=1.e-10;
        };
        if(!result.empty() && samePort(std::get<6>(result.back())) &&
           (std::get<0>(result.back())-point).norm()<=tolerance &&
           (std::get<1>(result.back())-patch->shift).norm()<=tolerance)
        {
            auto merged=std::get<6>(result.back());
            merged.crossesPeriodicSeam=merged.crossesPeriodicSeam || value.crossesPeriodicSeam;
            merged.onPhysicalBoundary=merged.onPhysicalBoundary || value.onPhysicalBoundary;
            merged.atPatchCorner=merged.atPatchCorner || value.atPatchCorner;
            if(value.polygonLocation)
            {
                merged.polygonLocation=value.polygonLocation;
            }
            if(value.originalVertex)
            {
                merged.originalVertex=value.originalVertex;
            }
            classify(merged);
            if(original)
            {
                result.pop_back();
                result.emplace_back(point,patch->shift,ids,EdgeMap(),0,original,merged);
            }
            else
            {
                std::get<6>(result.back())=merged;
            }
            return;
        }
        classify(value);
        result.emplace_back(point,patch->shift,ids,EdgeMap(),0,original,value);
    };
    // Geometry stage: split the directed polygon at patch intersections. Each
    // interval retains its original segment identity, including coincident arcs.
    const auto central=getPatch(VectorDim::Zero());
    VectorLowerDim lower=polyPoints.front().first;
    VectorLowerDim upper=lower;
    for(const auto& vertex : polyPoints)
    {
        lower=lower.cwiseMin(vertex.first);
        upper=upper.cwiseMax(vertex.first);
    }
    for(const auto& edge : central->edges())
    {
        lower=lower.cwiseMin(*edge->source);
        upper=upper.cwiseMax(*edge->source);
    }
    std::vector<std::shared_ptr<Patch>> patches{central};
    std::set<const Patch*> discovered{central.get()};
    for(size_t i=0; i<patches.size(); ++i)
    {
        for(const auto& edge : patches[i]->edges())
        {
            if(traversable(edge.get()))
            {
                const auto neighbor=crossSeam(edge.get());
                if(discovered.insert(neighbor.get()).second)
                {
                    VectorLowerDim lo=*neighbor->edges().front()->source;
                    VectorLowerDim hi=lo;
                    for(const auto& other : neighbor->edges())
                    {
                        lo=lo.cwiseMin(*other->source);
                        hi=hi.cwiseMax(*other->source);
                    }
                    if((lo.array()<=upper.array()+tolerance).all() &&
                       (hi.array()>=lower.array()-tolerance).all())
                    {
                        patches.push_back(neighbor);
                    }
                }
            }
        }
        if(patches.size()>100000)
        {
            throw std::runtime_error("Periodic polygon graph exceeded its patch limit.");
        }
    }
    struct Piece
    {
        size_t segment;
        double begin,end;
        const Patch* patch;
        VectorLowerDim source,sink;
    };
    std::vector<Piece> pieces;
    for(size_t segment=0; segment<polyPoints.size(); ++segment)
    {
        const VectorLowerDim a=polyPoints[segment].first;
        const VectorLowerDim d=polyPoints[(segment+1)%polyPoints.size()].first-a;
        if(d.norm()<=tolerance)
        {
            continue;
        }
        const double paramTolerance=tolerance/d.norm();
        std::vector<Piece> intervals;
        for(const auto& patch : patches)
        {
            double enter=0.0;
            double leave=1.0;
            bool outside=false;
            for(const auto& edge : patch->edges())
            {
                const VectorLowerDim direction=*edge->sink-*edge->source;
                const double f0=cross(direction,a-*edge->source);
                const double f1=cross(direction,a+d-*edge->source);
                // Same-grain periodic faces partition adjacent patches. A
                // geometric expansion here creates overlapping intervals at an
                // oblique seam/physical corner, possibly assigning a cut to the
                // wrong patch with no physical boundary port. Allow only
                // roundoff at these seams. Physical faces retain the geometric
                // contact tolerance so stationary boundary contours are stable.
                const double coordinateScale=std::max({1.0,a.norm(),(a+d).norm(),
                                                       edge->source->norm(),edge->sink->norm()});
                const double distanceTolerance=traversable(edge.get()) ?
                    64*std::numeric_limits<double>::epsilon()*coordinateScale : tolerance;
                const double eps=distanceTolerance*direction.norm();
                if(f0 < -eps && f1 < -eps)
                {
                    outside=true;
                    break;
                }
                if(f0 < -eps)
                {
                    enter=std::max(enter,f0/(f0-f1));
                }
                else if(f1 < -eps)
                {
                    leave=std::min(leave,f0/(f0-f1));
                }
            }
            if(!outside && leave-enter>paramTolerance)
            {
                intervals.push_back({segment,enter,leave,patch.get(),a+enter*d,a+leave*d});
            }
        }
        std::stable_sort(intervals.begin(),intervals.end(),[](const Piece& a,const Piece& b)
        {
            return a.begin<b.begin;
        });
        // A segment on a periodic seam belongs to both closed patches. Retain
        // one interval, preferring the incoming patch to avoid an artificial hop.
        for(size_t i=0; i<intervals.size();)
        {
            size_t end=i+1;
            while(end<intervals.size() && std::abs(intervals[end].begin-intervals[i].begin)<=paramTolerance &&
                  std::abs(intervals[end].end-intervals[i].end)<=paramTolerance)
            {
                ++end;
            }
            size_t chosen=i;
            if(!pieces.empty())
            {
                for(size_t j=i; j<end; ++j)
                {
                    if(intervals[j].patch==pieces.back().patch)
                    {
                        chosen=j;
                        break;
                    }
                }
            }
            pieces.push_back(intervals[chosen]);
            i=end;
        }
    }
    const auto nextPhysicalInDirection=[&](const Edge* edge,bool direction)
    {
        const Edge* next=direction ? edge->next : edge->prev;
        std::set<const Edge*> corners;
        while(traversable(next))
        {
            if(!corners.insert(next).second)
            {
                throw std::runtime_error("Periodic polygon graph has a cyclic seam corner.");
            }
            crossSeam(next);
            next=direction ? next->twin->next : next->twin->prev;
        }
        return next;
    };
    if(pieces.empty())
    {
        // No retained polygon arc does not imply an empty slipped area: the
        // polygon may surround an entire bounded grain cross-section. With no
        // crossing arcs, winding is constant in each convex patch interior.
        for(const auto& patch : patches)
        {
            VectorLowerDim interior=VectorLowerDim::Zero();
            for(const auto& edge : patch->edges())
            {
                interior+=*edge->source;
            }
            interior/=patch->edges().size();
            int winding=0;
            for(size_t i=0; i<polyPoints.size(); ++i)
            {
                const VectorLowerDim a=polyPoints[i].first-interior;
                const VectorLowerDim b=polyPoints[(i+1)%polyPoints.size()].first-interior;
                if(a.y()<=0 && b.y()>0 && cross(a,b)>0)
                {
                    ++winding;
                }
                else if(a.y()>0 && b.y()<=0 && cross(a,b)<0)
                {
                    --winding;
                }
            }
            if(winding==0)
            {
                continue;
            }
            const Edge* seed=nullptr;
            for(const auto& edge : patch->edges())
            {
                if(!traversable(edge.get()))
                {
                    seed=edge.get();
                    break;
                }
            }
            if(!seed)
            {
                throw std::runtime_error("Enclosed periodic patch has no physical boundary seed.");
            }
            const bool direction=winding>0;
            for(int layer=0; layer<std::abs(winding); ++layer)
            {
                const Edge* edge=seed;
                std::set<const Edge*> visited;
                do
                {
                    if(!visited.insert(edge).second || visited.size()>100000)
                    {
                        throw std::runtime_error("Enclosed physical boundary does not form a closed cycle.");
                    }
                    const VectorLowerDim corner=direction ? *edge->sink : *edge->source;
                    auto enclosed=metadata(corner,std::nullopt,std::nullopt,Arc::PhysicalBoundary,std::nullopt,false);
                    enclosed.enclosureLayer=layer;
                    append(corner,edge->patch,nullptr,enclosed);
                    const Edge* next=nextPhysicalInDirection(edge,direction);
                    if(next->patch!=edge->patch)
                    {
                        std::get<6>(result.back()).crossesPeriodicSeam=true;
                        enclosed.incomingArc=Arc::PeriodicSeam;
                        enclosed.crossesPeriodicSeam=true;
                        append(corner,next->patch,nullptr,enclosed);
                    }
                    edge=next;
                }
                while(edge!=seed);
            }
            // Retain one connected contour, as in the intersection-driven walk.
            return result;
        }
        return result; // The polygon is disjoint from every discovered patch.
    }
    const size_t none=std::numeric_limits<size_t>::max();
    const auto continuous=[&](size_t i,size_t j)
    {
        const auto& a=pieces[i];
        const auto& b=pieces[j];
        if((a.sink-b.source).norm()>tolerance)
        {
            return false;
        }
        if(a.segment==b.segment && b.begin>=a.end-1.e-10)
        {
            return true;
        }
        // Clipping can move an endpoint a roundoff-sized distance along its
        // original segment. Use the same physical tolerance as the skipped
        // short edges below, rather than requiring exact parameters 0 and 1.
        // Connectivity still follows the original ordered polygon path.
        if((a.sink-polyPoints[(a.segment+1)%polyPoints.size()].first).norm()>tolerance ||
           (b.source-polyPoints[b.segment].first).norm()>tolerance)
        {
            return false;
        }
        for(size_t k=(a.segment+1)%polyPoints.size(); k!=b.segment; k=(k+1)%polyPoints.size())
        {
            if((polyPoints[(k+1)%polyPoints.size()].first-polyPoints[k].first).norm()>tolerance)
            {
                return false;
            }
        }
        return true;
    };
    const auto nextPhysical=[&](const Edge* edge)
    {
        return nextPhysicalInDirection(edge,forward);
    };
    const auto boundaryPort=[&](const Patch* patch,const VectorLowerDim& point)
    {
        const EdgeIDs ids=edgeIDs(patch,point);
        const Edge* found=nullptr;
        for(short id : {ids.first,ids.second})
        {
            if(id>=0 && !traversable(patch->edges()[id].get()))
            {
                found=patch->edges()[id].get();
                const VectorLowerDim end=forward ? *found->sink : *found->source;
                if((point-end).norm()<=tolerance)
                {
                    found=nextPhysical(found);
                }
                break;
            }
        }
        if(!found)
        {
            throw std::runtime_error("Periodic polygon cut has no physical boundary port.");
        }
        return found;
    };
    struct Port
    {
        size_t piece;
        bool exit;
        const Edge* edge;
        double parameter;
        VectorLowerDim point;
        size_t polygonNext;
    };
    std::vector<Port> ports;
    std::vector<size_t> successor(pieces.size(),none);
    std::vector<size_t> exitPort(pieces.size(),none);
    std::vector<size_t> entryPort(pieces.size(),none);
    std::vector<bool> boundaryWrap(pieces.size(),false);
    std::vector<bool> boundaryReverse(pieces.size(),false);
    // Original polygon links define the port permutation. An outside chain
    // joins each exit to its original entry; it is never inferred by proximity.
    for(size_t i=0; i<pieces.size(); ++i)
    {
        const size_t j=(i+1)%pieces.size();
        if(continuous(i,j))
        {
            successor[i]=j;
        }
        else
        {
            const size_t out=ports.size();
            for(bool exit : {true,false})
            {
                const size_t piece=exit ? i : j;
                const VectorLowerDim point=exit ? pieces[piece].sink : pieces[piece].source;
                const Edge* edge=boundaryPort(pieces[piece].patch,point);
                const VectorLowerDim a=forward ? *edge->source : *edge->sink;
                const VectorLowerDim d=forward ? VectorLowerDim(*edge->sink-*edge->source) : VectorLowerDim(*edge->source-*edge->sink);
                ports.push_back({piece,exit,edge,d.dot(point-a)/d.squaredNorm(),point,ports.size()+1});
            }
            exitPort[i]=out;
            entryPort[j]=out+1;
        }
    }
    if(!ports.empty())
    {
        ports.back().polygonNext=0;
    }
    std::map<const Edge*,std::vector<size_t>> edgePorts;
    std::set<const Edge*> physicalEdges;
    for(const auto& patch : patches)
    {
        for(const auto& edge : patch->edges())
        {
            if(!traversable(edge.get()))
            {
                physicalEdges.insert(edge.get());
            }
        }
    }
    for(size_t i=0; i<ports.size(); ++i)
    {
        physicalEdges.insert(ports[i].edge);
        edgePorts[ports[i].edge].push_back(i);
    }
    // Edge incidence is the only boundary connectivity: next/prev, crossing
    // twins at periodic corners. Sorting is confined to intersection ordering.
    std::map<const Edge*,const Edge*> boundaryNext;
    std::map<const Edge*,const Edge*> boundaryPrevious;
    std::set<const Edge*> hasPredecessor;
    for(const auto* edge : physicalEdges)
    {
        const Edge* next=nextPhysical(edge);
        if(physicalEdges.count(next))
        {
            boundaryNext[edge]=next;
            boundaryPrevious[next]=edge;
            hasPredecessor.insert(next);
        }
        auto& list=edgePorts[edge];
        std::stable_sort(list.begin(),list.end(),[&](size_t a,size_t b)
        {
            return ports[a].parameter<ports[b].parameter;
        });
    }
    std::vector<const Edge*> seeds;
    for(const auto* edge : physicalEdges)
    {
        if(!hasPredecessor.count(edge))
        {
            seeds.push_back(edge);
        }
    }
    seeds.insert(seeds.end(),physicalEdges.begin(),physicalEdges.end());
    std::set<const Edge*> connected;
    for(const auto* seed : seeds)
    {
        if(connected.count(seed))
        {
            continue;
        }
        std::vector<std::vector<size_t>> buckets;
        const Edge* edge=seed;
        while(edge && connected.insert(edge).second)
        {
            const auto& list=edgePorts[edge];
            for(size_t i=0; i<list.size();)
            {
                size_t j=i+1;
                const double eps=tolerance/(*edge->sink-*edge->source).norm();
                while(j<list.size() && std::abs(ports[list[j]].parameter-ports[list[i]].parameter)<=eps)
                {
                    ++j;
                }
                buckets.emplace_back(list.begin()+i,list.begin()+j);
                i=j;
            }
            const auto next=boundaryNext.find(edge);
            edge=next==boundaryNext.end() ? nullptr : next->second;
        }
        if(buckets.empty())
        {
            continue;
        }
        // Count signed winding layers in boundary order. Only a closed edge
        // cycle may rotate its port list to a zero-winding starting point.
        int winding=0;
        int minimum=0;
        size_t first=0;
        for(size_t i=0; i<buckets.size(); ++i)
        {
            for(size_t port : buckets[i])
            {
                winding+=ports[port].exit ? 1 : -1;
            }
            if(winding<minimum)
            {
                minimum=winding;
                first=i+1;
            }
        }
        if(winding!=0)
        {
            throw std::runtime_error("Unbalanced physical-boundary intersection ports.");
        }
        bool zeroWindingCut=false;
        if(edge==seed)
        {
            // A closed boundary may meet a negative-winding lobe of a folded
            // polygon. Rotating to the minimum prefix count would add an
            // artificial full-grain winding. Find an actual zero-winding gap
            // using polygon intersections, then keep both signed layers.
            for(size_t k=0;k<buckets.size() && !zeroWindingCut;++k)
            {
                const auto& port=ports[buckets[k].front()];
                const auto& previous=ports[buckets[(k+buckets.size()-1)%buckets.size()].back()];
                // The gap before a corner port can lie entirely on preceding
                // edges. Follow boundary incidence to sample that gap; looking
                // only on the port's outgoing edge misses every such gap.
                std::vector<VectorLowerDim> samples;
                const Edge* gapEdge=previous.edge;
                VectorLowerDim gapStart=previous.point;
                bool wrap=gapEdge==port.edge && previous.parameter>=port.parameter;
                std::set<const Edge*> gapVisited;
                while(gapVisited.insert(gapEdge).second)
                {
                    const bool last=gapEdge==port.edge && !wrap;
                    const VectorLowerDim gapEnd=last ? port.point :
                        (forward ? *gapEdge->sink : *gapEdge->source);
                    if((gapEnd-gapStart).norm()>2*tolerance)
                    {
                        samples.push_back(0.5*(gapStart+gapEnd));
                    }
                    if(last)
                    {
                        break;
                    }
                    wrap=false;
                    gapEdge=boundaryNext.at(gapEdge);
                    gapStart=forward ? *gapEdge->source : *gapEdge->sink;
                    // Permit the final partial interval after wrapping once.
                    if(gapEdge==port.edge)
                    {
                        if((port.point-gapStart).norm()>2*tolerance)
                        {
                            samples.push_back(0.5*(gapStart+port.point));
                        }
                        break;
                    }
                }
                for(const auto& sample : samples)
                {
                    int sampleWinding=0;
                    bool onPolygon=false;
                    for(size_t j=0;j<polyPoints.size();++j)
                    {
                        const VectorLowerDim p=polyPoints[j].first-sample;
                        const VectorLowerDim q=polyPoints[(j+1)%polyPoints.size()].first-sample;
                        if(std::abs(cross(p,q))<=tolerance*(q-p).norm() && p.dot(q)<=tolerance*tolerance)
                        {
                            onPolygon=true;
                            break;
                        }
                        if(p.y()<=0 && q.y()>0 && cross(p,q)>0)
                        {
                            ++sampleWinding;
                        }
                        else if(p.y()>0 && q.y()<=0 && cross(p,q)<0)
                        {
                            --sampleWinding;
                        }
                    }
                    if(!onPolygon && sampleWinding==0)
                    {
                        std::rotate(buckets.begin(),buckets.begin()+k,buckets.end());
                        zeroWindingCut=true;
                        break;
                    }
                }
            }
        }
        if(edge!=seed || zeroWindingCut)
        {
            // An open chain has a fixed beginning: rotating its intersection
            // list would invent a connection between its two ends. Signed
            // winding layers may instead require a predecessor walk (an entry
            // followed by an exit), as in a self-intersecting moving loop.
            std::vector<size_t> exits;
            std::vector<std::set<size_t>> pendingEntries;
            const auto connectPorts=[&](size_t out,size_t in,bool reverse)
            {
                successor[ports[out].piece]=ports[in].piece;
                boundaryReverse[ports[out].piece]=reverse;
                if(edge==seed && ports[out].edge==ports[in].edge)
                {
                    const double delta=ports[in].parameter-ports[out].parameter;
                    const double eps=tolerance/(*ports[out].edge->sink-*ports[out].edge->source).norm();
                    boundaryWrap[ports[out].piece]=reverse ? delta>eps : delta<-eps;
                }
            };
            const auto chooseEntry=[&](size_t out,const std::set<size_t>& entries)
            {
                size_t in=ports[out].polygonNext;
                while(!entries.count(in) && in!=out)
                {
                    in=ports[in].polygonNext;
                }
                if(in==out)
                {
                    throw std::runtime_error("Disconnected polygon intersection ports.");
                }
                return in;
            };
            for(const auto& bucket : buckets)
            {
                std::vector<size_t> localExits;
                std::set<size_t> entries;
                for(size_t port : bucket)
                {
                    if(ports[port].exit)
                    {
                        localExits.push_back(port);
                    }
                    else
                    {
                        entries.insert(port);
                    }
                }
                // Resolve coincident ports by their original polygon
                // permutation before extending either signed winding layer.
                while(!localExits.empty() && !entries.empty())
                {
                    const size_t out=localExits.back();
                    localExits.pop_back();
                    const size_t in=chooseEntry(out,entries);
                    connectPorts(out,in,false);
                    entries.erase(in);
                }
                for(size_t out : localExits)
                {
                    if(pendingEntries.empty())
                    {
                        exits.push_back(out);
                    }
                    else
                    {
                        auto& candidates=pendingEntries.back();
                        const size_t in=chooseEntry(out,candidates);
                        connectPorts(out,in,true);
                        candidates.erase(in);
                        if(candidates.empty())
                        {
                            pendingEntries.pop_back();
                        }
                    }
                }
                while(!entries.empty())
                {
                    if(exits.empty())
                    {
                        pendingEntries.push_back(entries);
                        entries.clear();
                    }
                    else
                    {
                        const size_t out=exits.back();
                        exits.pop_back();
                        const size_t in=chooseEntry(out,entries);
                        connectPorts(out,in,false);
                        entries.erase(in);
                    }
                }
            }
            if(!exits.empty() || !pendingEntries.empty())
            {
                throw std::runtime_error("Unbalanced open-boundary intersection ports.");
            }
            continue;
        }
        std::vector<size_t> stack;
        for(size_t k=0; k<buckets.size(); ++k)
        {
            const auto& bucket=buckets[(first+k)%buckets.size()];
            std::set<size_t> entries;
            for(size_t port : bucket)
            {
                if(ports[port].exit)
                {
                    stack.push_back(port);
                }
                else
                {
                    entries.insert(port);
                }
            }
            while(!entries.empty())
            {
                if(stack.empty())
                {
                    throw std::runtime_error("Physical-boundary entry has no exit port.");
                }
                const size_t out=stack.back();
                stack.pop_back();
                // Coincident ports remain distinct. Follow the original port
                // permutation to resolve their sheet incidence, independent of
                // coordinates, node IDs, or the starting polygon vertex.
                size_t in=ports[out].polygonNext;
                while(!entries.count(in) && in!=out)
                {
                    in=ports[in].polygonNext;
                }
                if(in==out)
                {
                    throw std::runtime_error("Disconnected polygon intersection ports.");
                }
                successor[ports[out].piece]=ports[in].piece;
                // Intersection ordering records whether this boundary link
                // passes around the entire edge cycle before returning here.
                boundaryWrap[ports[out].piece]=ports[out].edge==ports[in].edge &&
                    ports[in].parameter<ports[out].parameter-tolerance/(*ports[out].edge->sink-*ports[out].edge->source).norm();
                entries.erase(in);
            }
        }
        if(!stack.empty())
        {
            throw std::runtime_error("Physical-boundary exit has no entry port.");
        }
    }
    // Walk the already connected graph. There is no geometric branch search
    // here: each polygon arc has exactly one successor, possibly via a face arc.
    size_t start=0;
    for(size_t i=0; i<pieces.size(); ++i)
    {
        if(pieces[i].begin==0.0)
        {
            start=i;
            break;
        }
    }
    const auto vertexAt=[&](const Piece& piece,bool sink)->std::optional<size_t>
    {
        if(sink && piece.end==1.0)
        {
            return (piece.segment+1)%polyPoints.size();
        }
        if(!sink && piece.begin==0.0)
        {
            return piece.segment;
        }
        return std::nullopt;
    };
    const auto appendPort=[&](const VectorLowerDim& point,const Patch* patch,size_t segment,
                        std::optional<size_t> vertex,Arc arc,bool seam)
    {
        append(point,patch,vertex ? polyPoints[*vertex].second : nullptr,
               metadata(point,segment,vertex,arc,arc==Arc::OriginalSegment ? std::optional<size_t>(segment) : std::nullopt,seam));
    };
    std::set<size_t> walked;
    size_t current=start;
    do
    {
        if(!walked.insert(current).second || successor[current]==none)
        {
            throw std::runtime_error("Invalid polygon intersection successor cycle.");
        }
        const auto& piece=pieces[current];
        const size_t next=successor[current];
        const auto& target=pieces[next];
        const bool boundaryArc=exitPort[current]!=none;
        appendPort(piece.sink,piece.patch,piece.segment,vertexAt(piece,true),Arc::OriginalSegment,false);
        if(boundaryArc)
        {
            const Edge* edge=ports[exitPort[current]].edge;
            const Edge* destination=ports[entryPort[next]].edge;
            const Patch* active=piece.patch;
            if(edge->patch!=active)
            {
                std::get<6>(result.back()).crossesPeriodicSeam=true;
                active=edge->patch;
                appendPort(piece.sink,active,piece.segment,vertexAt(piece,true),Arc::PeriodicSeam,true);
            }
            std::set<const Edge*> route;
            bool wrap=boundaryWrap[current];
            while(edge!=destination || wrap)
            {
                wrap=false;
                if(!route.insert(edge).second)
                {
                    throw std::runtime_error("Boundary arc cannot reach its assigned entry port.");
                }
                const VectorLowerDim corner=(forward!=boundaryReverse[current]) ? *edge->sink : *edge->source;
                append(corner,active,nullptr,metadata(corner,std::nullopt,std::nullopt,Arc::PhysicalBoundary,std::nullopt,false));
                const auto& adjacency=boundaryReverse[current] ? boundaryPrevious : boundaryNext;
                const auto following=adjacency.find(edge);
                if(following==adjacency.end())
                {
                    throw std::runtime_error("Boundary arc has no stored successor.");
                }
                edge=following->second;
                if(edge->patch!=active)
                {
                    std::get<6>(result.back()).crossesPeriodicSeam=true;
                    active=edge->patch;
                    append(corner,active,nullptr,metadata(corner,std::nullopt,std::nullopt,Arc::PeriodicSeam,std::nullopt,true));
                }
            }
            appendPort(target.source,active,target.segment,vertexAt(target,false),Arc::PhysicalBoundary,false);
            if(active!=target.patch)
            {
                std::get<6>(result.back()).crossesPeriodicSeam=true;
                appendPort(target.source,target.patch,target.segment,vertexAt(target,false),Arc::PeriodicSeam,true);
            }
        }
        else if(piece.patch!=target.patch)
        {
            std::get<6>(result.back()).crossesPeriodicSeam=true;
            appendPort(target.source,target.patch,target.segment,vertexAt(target,false),Arc::PeriodicSeam,true);
        }
        current=next;
    }
    while(current!=start);
    // Clipping can collapse a nondegenerate input polygon onto a physical
    // face. This is the COMPLETE walked contour, including every seam and
    // boundary bend, not merely its original anchors. With fewer than three
    // records its closed straight-edge path can only retrace a point or line.
    // Return the ordinary empty-contour result so reconciliation disconnects
    // only this loop's memberships, preserving shared network junctions.
    if(result.size()<3)
    {
        return {};
    }
    // Preserve the existing record contract: crossings between successive
    // original nodes carry the same edge map and a decreasing remaining count.
    size_t begin=0;
    for(size_t end=0; end<result.size(); ++end)
    {
        if(std::get<5>(result[end]) || end+1==result.size())
        {
            EdgeMap edges;
            size_t count=0;
            for(size_t i=begin; i<end; ++i)
            {
                edges[std::get<1>(result[i])].insert(std::get<2>(result[i]));
                ++count;
            }
            for(size_t i=begin; i<=end; ++i)
            {
                std::get<3>(result[i])=edges;
                if(i<end)
                {
                    std::get<4>(result[i])=--count;
                }
            }
            begin=end+1;
        }
    }
    return result;
}

template<typename T1,typename T2, typename T3,typename T4>
void printIntersections(const T1& ppi,const T2& usedPatches,const T3& lastPatch,int& removeMe, const T4& untwinnedEdges)
{
    std::ofstream pointsFile("points_"+std::to_string(removeMe)+".txt");
    for(const auto& tup : ppi)
    {
        //                        if(!node.second.expired())
        //                        {
        pointsFile<<"    "<<std::get<0>(tup).transpose()<<std::endl;
        //                        }
    }
    
    std::ofstream edgesFile("edges_"+std::to_string(removeMe)+".txt");
    for(const auto& patch : usedPatches)
    {
        for(const auto& edge : patch->edges())
        {
            if(untwinnedEdges.find(edge.get())!=untwinnedEdges.end())
            {
                edgesFile<<patch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<" 1"<<std::endl;
            }
            else
            {
                edgesFile<<patch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<" 0"<<std::endl;
            }
            
        }
    }
    
    std::ofstream lastPatchFile("lastPatch_"+std::to_string(removeMe)+".txt");
    for(const auto& edge : lastPatch->edges())
    {
        lastPatchFile<<lastPatch->sID<<" "<<edge->source->transpose()<<" "<<edge->sink->transpose()<<std::endl;
    }
    removeMe++;
}


}
#endif
