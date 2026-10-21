/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_PolyhedralMesh_H_
#define model_PolyhedralMesh_H_
#include <Polyhedron.h>
#include <tuple>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <utility>
namespace model
{
    template <int dim> class SimplicialMesh;
    template <int dim, typename Sequence> struct PolyhedralFactories;
    template <int dim, size_t... orders> struct PolyhedralFactories<dim, std::index_sequence<orders...>>
    {
        using Type = std::tuple<PolyhedronFactory<dim, int(orders)>...>;
    };
    template <int dim> class PolyhedralMesh
    {
        static_assert(dim >= 1);
        typename PolyhedralFactories<dim, std::make_index_sequence<dim + 1>>::Type factories;
        template <int order>
        void buildBoundary(const std::shared_ptr<Polyhedron<dim, order>> &cell,
                           const std::vector<std::vector<size_t>> &simplices,
                           const std::map<size_t, Eigen::Matrix<double, dim, 1>> &positions);
        template <int order> void validateLevel() const;
        template <int order> void attachSimplex(const Simplex<dim, order> *simplex);
        void buildPerimeter(const std::shared_ptr<Polyhedron<dim, 2>> &cell);
        void orientMeshFaces();

      public:
        using Vector = Eigen::Matrix<double, dim, 1>;
        struct FaceInput
        {
            size_t id;
            std::vector<size_t> regions;
            std::vector<std::array<size_t, dim>> simplices;
        };
        PolyhedralMesh()
        {
            std::apply([this](auto &...factory) { ((factory.owner = this), ...); }, factories);
        }
        PolyhedralMesh(const PolyhedralMesh &) = delete;
        PolyhedralMesh &operator=(const PolyhedralMesh &) = delete;
        PolyhedralMesh(const std::map<size_t, Vector> &positions, const std::vector<FaceInput> &faces);
        template <int order> const PolyhedronFactory<dim, order> &factory() const
        {
            return std::get<order>(factories);
        }
        std::shared_ptr<Polyhedron<dim, dim>> getRegion(size_t id)
        {
            return std::get<dim>(factories).get(id);
        }
        std::shared_ptr<Polyhedron<dim, dim - 1>> createFace(const Simplex<dim, dim - 1> *simplex);
        const Polyhedron<dim, dim> &region(size_t id) const;
        const Polyhedron<dim, dim - 1> &face(size_t id) const;
        void buildTopology();
        void validate() const;
    };
} // namespace model
#include <PolyhedralMeshImplementation.h>
#endif
