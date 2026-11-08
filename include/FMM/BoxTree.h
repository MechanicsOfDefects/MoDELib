/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_BoxTree_H_
#define model_BoxTree_H_

#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>
#include <limits>
#include <algorithm>

namespace model
{
    template<int dim> class BoxTree;

    // An affine cell of order 0..dim. Identity uses exact dyadic coordinates,
    // never rounded physical positions. Boundary and refinement are distinct.
    template<int dim>
    class BoxEntity
    {
      public:
        using Vector = Eigen::Matrix<double, dim, 1>;
        using Key = std::array<std::uint32_t, 2 * dim>;
        struct Incidence
        {
            std::shared_ptr<BoxEntity> entity;
            int axis;
            int side;
            int orientation;
        };

      private:
        friend class BoxTree<dim>;
        Key keyValue;
        Vector lowerValue;
        Vector upperValue;
        int orderValue;
        double measureValue=1;
        std::vector<Incidence> boundaryValue;
        std::map<Key, std::weak_ptr<BoxEntity>> parentsValue;
        std::vector<std::shared_ptr<BoxEntity>> subdivisionsValue;
        std::weak_ptr<BoxEntity> refinementParentValue;
        BoxEntity(const Key&, const Vector&, const Vector&);

      public:
        BoxEntity(const BoxEntity&) = delete;
        BoxEntity& operator=(const BoxEntity&) = delete;
        const Key& key() const;
        int order() const;
        // Physical axis-aligned bounds (vertices have identical lower/upper).
        const Vector& lower() const;
        const Vector& upper() const;
        Vector center() const;
        double measure() const;
        const std::vector<Incidence>& boundary() const;
        const std::map<Key, std::weak_ptr<BoxEntity>>& parents() const;
        const std::vector<std::shared_ptr<BoxEntity>>& subdivisions() const;
        std::shared_ptr<BoxEntity> refinementParent() const;
        const Incidence& face(int axis, int side) const;
    };

    // Uniform dyadic hierarchy for the initial loop-FMM geometry layer.
    // Physical cells may be lattice-aligned parallelepipeds. Periodic images
    // and wrapping remain explicit responsibilities of the caller.
    // A rebuild replaces the hierarchy atomically and increments revision().
    template<int dim>
    class BoxTree
    {
        static_assert(dim >= 1 && dim <= 3);
      public:
        using Entity = BoxEntity<dim>;
        using Vector = typename Entity::Vector;
        using Key = typename Entity::Key;
        using Pointer = std::shared_ptr<Entity>;
        using Factory = std::map<Key, Pointer>;
        static constexpr unsigned maxDepth = 5;
        static constexpr std::uint32_t latticeSize = 1u << 20;

      private:
        Vector lowerValue;
        Vector upperValue;
        Vector originValue;
        Eigen::Matrix<double,dim,dim> basisValue, inverseValue;
        double diameterValue;
        struct State
        {
            Factory factory;
            std::vector<std::vector<Pointer>> levels;
        };
        State state;
        std::uint64_t revisionValue = 0;
        Pointer get(State&, const Key&) const;
        void subdivide(State&, const Pointer&) const;

      public:
        BoxTree(const Vector& lower, const Vector& upper, unsigned depth = 0);
        using Matrix = Eigen::Matrix<double,dim,dim>;
        // Physical origin and cell edge vectors. lower()/upper() remain AABB bounds.
        struct LatticeTag {};
        BoxTree(const Vector& origin, const Matrix& basis, unsigned depth, LatticeTag);
        const Matrix& basis() const { return basisValue; }
        const Matrix& inverseBasis() const { return inverseValue; }
        Vector fractional(const Vector& x) const { return inverseValue*(x-originValue); }
        Vector physical(const Vector& q) const { return originValue+basisValue*q; }
        Vector cellCenter(unsigned level, std::size_t lex) const
        {
            const unsigned n=1u<<level;
            Vector q;
            for(int c=0;c<dim;++c) { q[c]=(lex%n+0.5)/n; lex/=n; }
            return physical(q);
        }
        double cellDiameter(unsigned level) const { return diameterValue/(1u<<level); }
        double minimumFaceSpacing(unsigned level) const
        {
            double result=std::numeric_limits<double>::max();
            for(int c=0;c<dim;++c) { result=std::min(result,1.0/inverseValue.row(c).norm()); }
            return result/(1u<<level);
        }
        bool admissible(const Vector& displacement,unsigned level) const
        { return displacement.norm()>1.01*cellDiameter(level); }
        // Translation-invariant stencil shared by direct leaf interactions.
        std::vector<Eigen::Matrix<int,dim,1>> nearOffsets(unsigned level) const;
        bool contains(const Vector& x,double tolerance=0) const
        {
            if(!x.allFinite()) { return false; }
            const Vector q=fractional(x);
            for(int c=0;c<dim;++c)
            {
                const double eps=tolerance*inverseValue.row(c).norm();
                if(q[c]<-eps || q[c]>1+eps) { return false; }
            }
            return true;
        }
        BoxTree(const BoxTree&) = delete;
        BoxTree& operator=(const BoxTree&) = delete;
        const Vector& lower() const;
        const Vector& upper() const;
        unsigned depth() const;
        std::uint64_t revision() const;
        const Factory& entities() const;
        const std::vector<Pointer>& boxes(unsigned level) const;
        const Pointer& root() const;
        void setDepth(unsigned depth);
        // Same-level topological adjacency; outer faces return nullptr.
        Pointer neighbor(const Pointer& box, int axis, int side) const;
        // Half-open cells, with the root's upper boundary included in the last cell.
        Pointer locate(const Vector& point, unsigned level) const;
        // Unique edges of a displayed level (no overdraw on shared faces).
        std::vector<Pointer> edges(unsigned level) const;
    };

    extern template class BoxEntity<1>;
    extern template class BoxEntity<2>;
    extern template class BoxEntity<3>;
    extern template class BoxTree<1>;
    extern template class BoxTree<2>;
    extern template class BoxTree<3>;
}
#endif
