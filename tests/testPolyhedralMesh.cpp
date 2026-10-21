#include <PolyhedralMesh.h>
#include <MeshModule.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace
{
    using V = Eigen::Vector3d;
    void require(bool ok, const char *message)
    {
        if (!ok)
        {
            throw std::runtime_error(message);
        }
    }
    double sectionArea(const model::BoundingMeshSegments<3> &section, const model::Plane<3> &plane)
    {
        double area = 0;
        for (size_t k = 0; k < section.size(); ++k)
        {
            const auto &arc = *section[k];
            require(arc.P1 == section[(k + 1) % section.size()]->P0, "Section is not exactly connected");
            require(std::abs((arc.P0 - plane.P).dot(plane.unitNormal)) < 1.e-6, "Section vertex left plane");
            area += 0.5 * (arc.P0 - plane.P).cross(arc.P1 - plane.P).dot(plane.unitNormal);
        }
        require(area >= -1.e-10, "Section orientation is reversed");
        return area;
    }
    // Independent reference: sum planar cuts of the volume tetrahedra.
    // This geometric ordering is deliberately confined to the test oracle.
    double referenceArea(const model::Polyhedron<3,3> &region, const model::Plane<3> &plane)
    {
        double area = 0;
        for (const auto *tet : region.simplices())
        {
            const auto vertices = tet->vertexPositionMatrix();
            std::vector<Eigen::Vector2d> points;
            auto add = [&](const V &p)
            {
                const Eigen::Vector2d x = plane.localPosition(p);
                for (const auto &existing : points)
                {
                    if ((x - existing).norm() < 1.e-9)
                    {
                        return;
                    }
                }
                points.push_back(x);
            };
            for (int i = 0; i < 4; ++i)
            {
                const double a = (vertices.col(i) - plane.P).dot(plane.unitNormal);
                if (std::abs(a) < 1.e-10)
                {
                    add(vertices.col(i));
                }
                for (int j = i + 1; j < 4; ++j)
                {
                    const double b = (vertices.col(j) - plane.P).dot(plane.unitNormal);
                    if (a * b < 0 && std::abs(a) > 1.e-10 && std::abs(b) > 1.e-10)
                    {
                        add(vertices.col(i) + a / (a - b) * (vertices.col(j) - vertices.col(i)));
                    }
                }
            }
            if (points.size() < 3)
            {
                continue;
            }
            Eigen::Vector2d center = Eigen::Vector2d::Zero();
            for (const auto &p : points)
            {
                center += p / double(points.size());
            }
            std::sort(points.begin(), points.end(),
                      [&](const auto &a, const auto &b)
                      {
                          return std::atan2(a.y() - center.y(), a.x() - center.x()) <
                                 std::atan2(b.y() - center.y(), b.x() - center.x());
                      });
            for (size_t k = 0; k < points.size(); ++k)
            {
                const auto a = (points[k] - center).eval();
                const auto b = (points[(k + 1) % points.size()] - center).eval();
                area += 0.5 * (a.x() * b.y() - a.y() * b.x());
            }
        }
        return area;
    }
    void exercise(const std::string &path, const Eigen::Matrix3d &rotation, const V &shift)
    {
        std::weak_ptr<model::Polyhedron<3, 0>> weakVertex;
        std::weak_ptr<model::Polyhedron<3, 2>> weakFace;
        {
            model::SimplicialMesh<3> mesh(path, rotation, shift, {-1});
            const auto &topology = mesh.polyhedralMesh();
            topology.validate();
            weakVertex = topology.factory<0>().begin()->second;
            weakFace = topology.factory<2>().begin()->second;
            require(topology.factory<3>().size() == mesh.regions().size(), "Missing regions");
            for (const auto& region : mesh.regions())
            {
                require(region.second.get() == &topology.region(region.first), "Parallel region objects remain");
                require(&region.second->faces() == &region.second->children, "Faces are not the region children");
                for (const auto* simplex : region.second->simplices())
                {
                    require(simplex->region == region.second, "Simplex references a different region");
                }
                for (const auto& face : region.second->faces())
                {
                    require(face.second.get() == &topology.face(face.first), "Parallel face objects remain");
                    require(face.second->parents.at(region.first).lock() == region.second, "Face parent is not its region");
                    require(!face.second->simplices().empty(), "Face lost its boundary simplices");
                }
            }
            for (const auto& edge : topology.factory<1>())
            {
                require(edge.second->simplices().size() == 1, "Edge does not store its mesh simplex");
            }
            for (const auto& vertex : topology.factory<0>())
            {
                require(vertex.second->simplices().size() == 1, "Vertex does not store its mesh simplex");
            }
            bool sharedJunction = false;
            for (const auto &edge : topology.factory<1>())
            {
                sharedJunction = sharedJunction || edge.second->parents.size() > 2;
            }
            require(sharedJunction, "Fixture did not exercise an edge with more than two parent faces");
            for (const auto &f : topology.factory<2>())
            {
                const auto &face = *f.second;
                require((face.parents.size() == 1) == face.isExternal(), "External classification differs");
                for (const auto &incidence : face.perimeter)
                {
                    require(topology.factory<1>().at(incidence.edge->key) == incidence.edge, "Duplicate edge object");
                    for (const auto &vertex : incidence.edge->children)
                    {
                        require(topology.factory<0>().at(vertex.first) == vertex.second, "Duplicate vertex object");
                        require(vertex.second->parents.at(incidence.edge->key).lock() == incidence.edge,
                                "Missing edge parent");
                    }
                    for (const auto &parent : face.parents)
                    {
                        const auto other = incidence.otherFace(parent.first);
                        require(other.otherFace(parent.first).face.lock() == f.second, "Face switch not reciprocal");
                        require(incidence.next().prev().source() == incidence.source(), "Broken perimeter links");
                    }
                }
            }
            for (const auto &r : mesh.regions())
            {
                V center = V::Zero();
                for (const auto *tet : r.second->simplices())
                {
                    center += tet->center() / double(r.second->simplices().size());
                }
                for (const V &localNormal : {V(1, 0.317, -0.219), V(0.731, 1, 0.277)})
                {
                    const model::Plane<3> plane(center, rotation * localNormal);
                    const model::BoundingMeshSegments<3> section(mesh, r.first, plane);
                    const double expected = referenceArea(*r.second, plane);
                    require(std::abs(sectionArea(section, plane) - expected) < 1.e-7 * std::max(1.0, expected),
                            "Grain contour area differs from tetrahedral cuts");
                    for (const auto &segment : section)
                    {
                        for (const auto *face : segment->faces)
                        {
                            require(r.second->faces().count(face->sID), "Section includes another grain's face");
                        }
                    }
                    const model::Plane<3> opposite(center, -rotation * localNormal);
                    require(std::abs(sectionArea(model::BoundingMeshSegments<3>(mesh, r.first, opposite), opposite) -
                                     expected) < 1.e-7 * std::max(1.0, expected),
                            "Reversed normal changed contour area");
                }
                const auto &face = *topology.region(r.first).children.begin()->second;
                const auto &incidence = face.perimeter.front();
                const V vertex = incidence.source()->position;
                const model::Plane<3> vertexPlane(vertex, rotation * V(0.319, 0.711, 1.0));
                const double expected = referenceArea(*r.second, vertexPlane);
                require(std::abs(sectionArea(model::BoundingMeshSegments<3>(mesh, r.first, vertexPlane), vertexPlane) -
                                 expected) < 1.e-7 * std::max(1.0, expected),
                        "Vertex-hit contour changed area");
                const V normal = incidence.tangent().cross(rotation * V(0.419, 0.833, 0.231));
                const model::Plane<3> edgePlane(vertex, normal);
                const auto edgeSection = model::BoundingMeshSegments<3>(mesh, r.first, edgePlane);
                sectionArea(edgeSection, edgePlane);
                if (!edgeSection.empty())
                {
                    bool shared = false;
                    for (const auto &segment : edgeSection)
                    {
                        shared = shared || segment->faces.size() == 2;
                    }
                    require(shared, "Coincident edge lost its two supporting faces");
                }
            }
        }
        require(weakVertex.expired() && weakFace.expired(), "Polyhedral topology contains an ownership cycle");
    }
} // namespace
int main(int argc, char **argv)
{
    try
    {
        require(argc == 2, "Expected repository path");
        const std::string meshes = std::string(argv[1]) + "/Library/Meshes/";
        const Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.71, V(1, 2, 3).normalized()).toRotationMatrix();
        for (const std::string name : {"bicrystal_48.msh", "polycrystalCube_10_periodicX.msh",
                                       "polycrystalCube_10_periodicXY.msh", "poly30_periodicXYZ.msh"})
        {
            exercise(meshes + name, Eigen::Matrix3d::Identity(), V::Zero());
            exercise(meshes + name, rotation, V(10, -5, 3));
        }
        // A coincident face, edge tangency and vertex tangency of the box.
        model::SimplicialMesh<3> mesh(meshes + "bicrystal_48.msh", Eigen::Matrix3d::Identity(), V::Zero(), {});
        for (const auto &r : mesh.regions())
        {
            V lo = V::Constant(std::numeric_limits<double>::max());
            V hi = -lo;
            for (const auto *tet : r.second->simplices())
            {
                for (const auto *vertex : tet->vertices())
                {
                    lo = lo.cwiseMin(vertex->P0);
                    hi = hi.cwiseMax(vertex->P0);
                }
            }
            const model::Plane<3> facePlane(lo, V::UnitZ());
            const auto section = model::BoundingMeshSegments<3>(mesh, r.first, facePlane);
            require(std::abs(sectionArea(section, facePlane) - (hi.x() - lo.x()) * (hi.y() - lo.y())) < 1.e-8,
                    "Coincident face area incorrect");
            require(model::BoundingMeshSegments<3>(mesh, r.first, model::Plane<3>(lo, V::Ones())).empty(),
                    "Vertex tangency is not empty");
            require(model::BoundingMeshSegments<3>(mesh, r.first, model::Plane<3>(lo, V(1, 1, 0))).empty(),
                    "Edge tangency is not empty");
        }
        std::weak_ptr<model::Polyhedron<3,3>> oldRegion = mesh.regions().begin()->second;
        std::weak_ptr<model::Polyhedron<3, 0>> oldVertex = mesh.polyhedralMesh().factory<0>().begin()->second;
        mesh.readMesh(meshes + "bicrystal_48.msh", rotation, V(1, 2, 3), {});
        require(oldRegion.expired(), "Reload retained the previous region objects");
        require(oldVertex.expired(), "Reload retained the previous topology");
        mesh.polyhedralMesh().validate();
        std::cout << "PASS factory sharing, reciprocal incidences, grain contours, degeneracies and ownership\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
