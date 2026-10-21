#include <PolyhedralMesh.h>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace
{
    void require(bool condition, const char *message)
    {
        if (!condition)
        {
            throw std::runtime_error(message);
        }
    }
    void interval()
    {
        using Mesh = model::PolyhedralMesh<1>;
        std::map<size_t, Mesh::Vector> positions;
        for (size_t k = 0; k < 3; ++k)
        {
            positions[k][0] = double(k);
        }
        Mesh mesh(positions, {{10, {0}, {{{0}}}}, {11, {0, 1}, {{{1}}}}, {12, {1}, {{{2}}}}});
        mesh.validate();
        require(mesh.face(10).isExternal() && !mesh.face(11).isExternal(), "1D external face classification incorrect");
        require(mesh.face(11).parents.size() == 2, "1D shared vertex lost its regions");
        require(mesh.region(0).children.at({11}) == mesh.region(1).children.at({11}), "1D vertex not shared");
        require(mesh.face(12).position[0] == 2, "1D face position incorrect");
    }
    void polygons()
    {
        using Mesh = model::PolyhedralMesh<2>;
        static_assert(
            std::is_same_v<decltype(std::declval<const Mesh &>().face(0)), const model::Polyhedron<2, 1> &>);
        const std::map<size_t, Mesh::Vector> positions{{0, {0, 0}}, {1, {1, 0}}, {2, {2, 0}},  {3, {0, 1}},
                                                       {4, {1, 1}}, {5, {2, 1}}, {6, {1, 0.5}}};
        const std::vector<Mesh::FaceInput> faces{
            {10, {0}, {{{0, 1}}}}, {11, {1}, {{{1, 2}}}}, {12, {0}, {{{0, 3}}}}, {13, {0, 1}, {{{1, 6}}, {{6, 4}}}},
            {14, {1}, {{{2, 5}}}}, {15, {0}, {{{3, 4}}}}, {16, {1}, {{{4, 5}}}}};
        std::weak_ptr<model::Polyhedron<2, 0>> weak;
        {
            Mesh mesh(positions, faces);
            mesh.validate();
            require(mesh.factory<2>().size() == 2 && mesh.factory<1>().size() == 7 &&
                        mesh.template factory<0>().size() == 6,
                    "2D factory counts incorrect");
            require(mesh.face(10).isExternal() && !mesh.face(13).isExternal(),
                    "2D external face classification incorrect");
            require(mesh.face(13).children.count({1}) && mesh.face(13).children.count({4}),
                    "2D face endpoints incorrect");
            require(mesh.region(0).children.at({13}) == mesh.region(1).children.at({13}), "2D face not shared");
            require(mesh.template factory<0>().at({1})->parents.size() == 3, "Triple edge incidence lost");
            for (const auto &region : mesh.factory<2>())
            {
                require(region.second->perimeter.size() == 4, "2D region perimeter missing");
                for (const auto &incidence : region.second->perimeter)
                {
                    require(incidence.sink() == incidence.next().source(), "2D perimeter disconnected");
                }
            }
            weak = mesh.template factory<0>().at({1});
        }
        require(weak.expired(), "2D ownership cycle");
        auto broken = faces;
        broken.pop_back();
        bool rejected = false;
        try
        {
            Mesh mesh(positions, broken);
        }
        catch (const std::runtime_error &)
        {
            rejected = true;
        }
        require(rejected, "Open polygon accepted");
    }
    template <int dim> void simplexBoundary()
    {
        using Mesh = model::PolyhedralMesh<dim>;
        std::map<size_t, typename Mesh::Vector> positions;
        positions[0] = Mesh::Vector::Zero();
        for (int k = 0; k < dim; ++k)
        {
            positions[k + 1] = Mesh::Vector::Unit(k);
        }
        std::vector<typename Mesh::FaceInput> faces;
        for (int excluded = 0; excluded <= dim; ++excluded)
        {
            typename Mesh::FaceInput face;
            face.id = 100 + excluded;
            face.regions = {0};
            std::array<size_t, dim> simplex;
            size_t index = 0;
            for (int k = 0; k <= dim; ++k)
            {
                if (k != excluded)
                {
                    simplex[index++] = k;
                }
            }
            face.simplices.push_back(simplex);
            faces.push_back(face);
        }
        Mesh mesh(positions, faces);
        mesh.validate();
        require(mesh.template factory<0>().size() == dim + 1, "Simplex vertices not shared");
        require(mesh.region(0).children.size() == dim + 1, "Simplex facets incorrect");
        if constexpr (dim >= 3)
        {
            for (const auto &face : mesh.template factory<2>())
            {
                for (const auto &incidence : face.second->perimeter)
                {
                    require(incidence.sink() == incidence.next().source(), "Polygon connectivity incorrect");
                    require(incidence.next().prev().edge == incidence.edge, "Polygon next/prev incorrect");
                }
            }
        }
    }
} // namespace
int main()
{
    interval();
    polygons();
    simplexBoundary<3>();
    simplexBoundary<4>();
    simplexBoundary<5>();
    std::cout << "Dimension-generic polyhedral topology passed (1D through 5D)\n";
}
