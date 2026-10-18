/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_PolygonPatchIntersection_H_
#define model_PolygonPatchIntersection_H_
#include <Eigen/Core>
#include <CompareVectorsByComponent.h>
#include <map>
#include <set>
#include <tuple>
#include <optional>

namespace model
{
    enum class PolygonPatchNodeKind
    {
        OriginalVertex,
        PeriodicCrossing,
        PhysicalBoundaryIntersection,
        PhysicalBoundaryCorner
    };

    enum class PolygonPatchArcKind
    {
        OriginalSegment,
        PhysicalBoundary,
        PeriodicSeam
    };

    template<typename T>
    struct PolygonOriginalSegment
    {
        // Indices refer to the supplied polygon ordering, including its closing
        // edge. Pointers are borrowed, just like legacy record field 5; they may
        // be null and do not extend the lifetime of the original objects.
        size_t sourceIndex;
        size_t sinkIndex;
        const T* source;
        const T* sink;
    };

    template<typename T>
    struct PolygonSegmentLocation
    {
        PolygonOriginalSegment<T> segment;
        double parameter; // x=(1-parameter)*source+parameter*sink in the input polygon
    };

    template<typename T>
    struct PolygonPatchMetadata
    {
        PolygonPatchNodeKind kind=PolygonPatchNodeKind::OriginalVertex;
        std::optional<size_t> originalVertex;
        // Full enclosure has no local segment port. Retain the winding sheet
        // so reconciliation can transfer complete pre-update contour ownership.
        std::optional<size_t> enclosureLayer;
        std::optional<PolygonSegmentLocation<T>> polygonLocation;
        // These flags are independent: a physical-boundary corner can also
        // cross a periodic seam. Physical includes cross-grain interfaces.
        bool onPhysicalBoundary=false;
        bool crossesPeriodicSeam=false;
        bool atPatchCorner=false;
        // Describes the directed arc from the preceding output record to this
        // record (cyclically). Seam arcs have zero length in the extended plane.
        PolygonPatchArcKind incomingArc=PolygonPatchArcKind::OriginalSegment;
        std::optional<PolygonOriginalSegment<T>> incomingOriginalSegment;
    };

    template<int dim,typename T>
    using PolygonPatchIntersectionRecord=std::tuple<
        Eigen::Matrix<double,dim-1,1>, Eigen::Matrix<double,dim,1>,
        std::pair<short int,short int>,
        std::map<Eigen::Matrix<double,dim,1>,std::set<std::pair<short int,short int>>,
                 CompareVectorsByComponent<double,dim,float>>,
        size_t,const T* const,PolygonPatchMetadata<T>>;
    // Fields 0..5 retain their existing layout. Field 6 supplies provenance;
    // boundary-walk arcs deliberately have no original segment attribution.
}
#endif
