/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_Polyhedron_H_
#define model_Polyhedron_H_
#include <Eigen/Dense>
#include <StaticID.h>
#include <SharedPtrFactories.h>
#include <map>
#include <set>
#include <memory>
#include <vector>
#include <type_traits>

namespace model
{
    template <short int dim, short int order> class Simplex;
    template <int dim> struct Plane;
    template <int dim> class PolyhedralMesh;
    template <int dim> class SimplicialMesh;
    template <int dim, int order> struct Polyhedron;
    template <int dim, int order> struct PolyhedronFactory;
    template <int dim, int order>
    using PolyhedronKey = std::conditional_t<(order >= dim - 1), size_t, std::vector<size_t>>;

    template <int dim, int order, bool enabled = (order > 0)> struct PolyhedralChildren
    {
        using Child = Polyhedron<dim, order - 1>;
        using Children = std::map<PolyhedronKey<dim, order - 1>, std::shared_ptr<Child>>;
        Children children;
    };
    template <int dim, int order> struct PolyhedralChildren<dim, order, false>
    {
    };
    template <int dim, int order, bool enabled = (order < dim)> struct PolyhedralParents
    {
        using Parent = Polyhedron<dim, order + 1>;
        std::map<PolyhedronKey<dim, order + 1>, std::weak_ptr<Parent>> parents;
    };
    template <int dim, int order> struct PolyhedralParents<dim, order, false>
    {
    };

    // An oriented edge of a polygonal 2-cell, including regions in 2D.
    template <int dim> struct PolyhedralFaceIncidence
    {
        std::weak_ptr<Polyhedron<dim, 2>> face;
        std::shared_ptr<Polyhedron<dim, 1>> edge;
        size_t index;
        bool reversed;
        std::shared_ptr<Polyhedron<dim, 0>> source() const;
        std::shared_ptr<Polyhedron<dim, 0>> sink() const;
        Eigen::Matrix<double, dim, 1> tangent() const;
        Eigen::Matrix<double, dim, 1> outNormal(size_t regionID) const
            requires(dim == 3);
        PolyhedralFaceIncidence next() const;
        PolyhedralFaceIncidence prev() const;
        PolyhedralFaceIncidence otherFace(size_t regionID) const
            requires(dim == 3);
    };
    template <int dim, int order> struct PolyhedralGeometry
    {
    };
    template <int dim> struct PolyhedralGeometry<dim, 0>
    {
        Eigen::Matrix<double, dim, 1> position = Eigen::Matrix<double, dim, 1>::Zero();
    };
    template <int dim> struct PolyhedralGeometry<dim, 1>
    {
        std::map<PolyhedronKey<dim, 2>, size_t> faceIncidences;
    };
    template <int dim> struct PolyhedralGeometry<dim, 2>
    {
        std::vector<PolyhedralFaceIncidence<dim>> perimeter;
    };
    template <int dim, int order, bool face = (order == dim - 1)> struct PolyhedralFaceData
    {
        template <typename Key> explicit PolyhedralFaceData(const Key &)
        {
        }
    };
    template <int dim, int order> struct PolyhedralFaceData<dim, order, true>
    {
        using VectorDim = Eigen::Matrix<double, dim, 1>;
        const size_t sID;
        explicit PolyhedralFaceData(size_t id) : sID(id)
        {
        }
        std::pair<int, int> regionIDs{-1, -1};
        VectorDim n = VectorDim::Zero();
        VectorDim c = VectorDim::Zero();
        std::vector<const Simplex<dim, 0> *> hull;
        std::pair<VectorDim, const Polyhedron<dim, order> *> periodicFacePair{VectorDim::Zero(), nullptr};
    };
    template <int dim, int order, bool region = (order == dim)> struct PolyhedralRegionData
    {
        template <typename Key> explicit PolyhedralRegionData(const Key &)
        {
        }
    };
    template <int dim, int order> struct PolyhedralRegionData<dim, order, true>
    {
        const size_t regionID;
        explicit PolyhedralRegionData(size_t id) : regionID(id)
        {
        }
    };

    template <int dim, int order>
    struct Polyhedron : PolyhedralChildren<dim, order>,
                        PolyhedralParents<dim, order>,
                        PolyhedralGeometry<dim, order>,
                        PolyhedralFaceData<dim, order>,
                        PolyhedralRegionData<dim, order>,
                        std::set<const Simplex<dim, order> *>
    {
        static_assert(order >= 0 && order <= dim);
        using KeyType = PolyhedronKey<dim, order>;
        using CompareType = std::less<KeyType>;
        using VectorDim = Eigen::Matrix<double, dim, 1>;
        using SimplexType = Simplex<dim, order>;
        using SimplexContainerType = std::set<const SimplexType *>;
        using MeshFacesContainerType = std::map<size_t, std::shared_ptr<Polyhedron<dim, dim - 1>>>;
        const KeyType key;
        PolyhedronFactory<dim, order> &cellFactory;
        Polyhedron(PolyhedronFactory<dim, order> &factory, const KeyType &id)
            : PolyhedralFaceData<dim, order>(id), PolyhedralRegionData<dim, order>(id), key(id), cellFactory(factory)
        {
        }
        const SimplexContainerType &simplices() const
        {
            return *this;
        }
        SimplexContainerType &simplices()
        {
            return *this;
        }
        const MeshFacesContainerType &faces() const
            requires(order == dim)
        {
            return this->children;
        }
        MeshFacesContainerType &faces()
            requires(order == dim)
        {
            return this->children;
        }
        bool isExternal() const
            requires(order == dim - 1)
        {
            return this->parents.size() == 1;
        }
        void initialize(const SimplexType *simplex)
            requires(order == dim - 1);
        const SimplexContainerType &internalSimplices() const
            requires(order == dim - 1)
        {
            return *this;
        }
        const std::vector<const Simplex<dim, 0> *> &convexHull() const
            requires(order == dim - 1);
        void finalize()
            requires(order == dim - 1);
        VectorDim outNormal() const
            requires(order == dim - 1);
        VectorDim outNormal(const int &regionID) const
            requires(order == dim - 1);
        const VectorDim &center() const
            requires(order == dim - 1);
        Plane<dim> asPlane() const
            requires(order == dim - 1);
        VectorDim outNormal(const std::set<size_t> &ids) const
            requires(order == dim);
        void update()
            requires(order == dim);

      private:
        void buildFaces()
            requires(order == dim);
        void buildSingleFace(const Simplex<dim, dim - 1> *start, std::shared_ptr<Polyhedron<dim, dim - 1>> face,
                             std::map<std::vector<size_t>, const Simplex<dim, dim - 1> *> &remaining)
            requires(order == dim);
    };
    template <int dim, int order>
    struct PolyhedronFactory : KeyConstructableSharedPtrFactory<PolyhedronFactory<dim, order>, Polyhedron<dim, order>>
    {
        PolyhedralMesh<dim> *owner = nullptr;
    };
} // namespace model
#endif
