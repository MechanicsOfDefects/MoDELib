/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_PlaneBoxIntersection_H_
#define model_PlaneBoxIntersection_H_
#include <FMM/BoxTree.h>

namespace model
{
    // Immutable cut of one tree level by an infinite plane. Currently 3D only.
    // Identities come from box entities, not rounded coordinates. The result
    // remains usable after the tree is rebuilt; revision identifies its snapshot.
    template<int dim>
    class PlaneBoxIntersection
    {
        static_assert(dim == 3);
      public:
        using Entity = BoxEntity<dim>;
        using Key = typename Entity::Key;
        using Vector = typename Entity::Vector;
        struct Node
        {
            std::shared_ptr<const Entity> support; // vertex or crossed edge
            Vector position;
        };
        using NodePtr = std::shared_ptr<const Node>;
        using SegmentKey = std::pair<Key, Key>;
        struct Segment
        {
            NodePtr first;
            NodePtr second;
        };
        using SegmentPtr = std::shared_ptr<const Segment>;
        struct Side
        {
            SegmentPtr segment;
            int orientation;
        };
        struct Polygon
        {
            std::shared_ptr<const Entity> box;
            std::vector<NodePtr> vertices; // right-handed about normal()
            std::vector<Side> boundary;
            double area;
        };
      private:
        Vector originValue;
        Vector normalValue;
        unsigned levelValue;
        std::uint64_t revisionValue;
        double tolerance;
        std::map<Key, NodePtr> nodesValue;
        std::map<SegmentKey, SegmentPtr> segmentsValue;
        std::map<Key, std::vector<NodePtr>> edgeCuts;
        std::map<Key, std::vector<SegmentPtr>> faceCutsValue;
        std::map<Key, Polygon> polygonsValue;
        NodePtr node(const std::shared_ptr<const Entity>&, const Vector&);
        SegmentPtr segment(NodePtr, NodePtr);
        const std::vector<NodePtr>& cutEdge(const std::shared_ptr<const Entity>&);
        const std::vector<SegmentPtr>& cutFace(const std::shared_ptr<const Entity>&);
      public:
        PlaneBoxIntersection(const BoxTree<dim>&, unsigned level, const Vector& origin, const Vector& normal);
        const Vector& origin() const;
        const Vector& normal() const;
        unsigned level() const;
        std::uint64_t revision() const;
        const std::map<Key, Polygon>& polygons() const;
        const std::map<Key, std::vector<SegmentPtr>>& faceCuts() const;
    };
    extern template class PlaneBoxIntersection<3>;
}
#endif
