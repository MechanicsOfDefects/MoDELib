/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_SurfaceMeshGeometry_H_
#define model_SurfaceMeshGeometry_H_
#include <SimplicialMesh.h>
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace model
{
    // A conforming refinement of the actual boundary simplices, not their convex hull.
    // Coordinates remain in simulation length units [b]. Never deform the reference mesh.
    struct SurfaceMeshGeometry
    {
        using Vector = Eigen::Vector3d;
        using Triangle = std::array<size_t, 3>;
        std::vector<Vector> referencePoints;
        std::vector<const Simplex<3, 3> *> parents;
        std::vector<Triangle> triangles;
        std::vector<const Simplex<3, 2> *> boundaryFaces;
        std::vector<Vector> samplePoints;
        std::vector<Triangle> walls;

        struct Trace
        {
            Vector point;
            Vector normal;
            std::set<const Simplex<3, 2> *> faces;
        };

        static const Simplex<3,3>* samplingParent(const Simplex<3,2>* face)
        {
            const Simplex<3,3>* result=nullptr;
            for(const auto& entry : face->parents())
            {
                if(!result || entry.second->region->regionID<result->region->regionID)
                {
                    result=entry.second;
                }
            }
            return result;
        }

        // Default to the exterior; explicit selection can include internal faces.
        static SurfaceMeshGeometry build(const SimplicialMesh<3>& mesh,double maxArea)
        {
            std::set<size_t> selected;
            for(const auto& region : mesh.regions())
            {
                for(const auto& face : region.second->faces())
                {
                    if(face.second->isExternal())
                    {
                        selected.insert(face.second->sID);
                    }
                }
            }
            return build(mesh,maxArea,selected);
        }

        static SurfaceMeshGeometry build(const SimplicialMesh<3> &mesh, double maxArea,
                                         const std::set<size_t> &selectedFaceIDs)
        {
            if(!std::isfinite(maxArea) || maxArea <= 0.0)
            {
                throw std::invalid_argument("Maximum triangle area must be finite and positive");
            }
            std::set<const Simplex<3, 2> *> selectedSimplices;
            for(const auto &region : mesh.regions())
            {
                for(const auto &face : region.second->faces())
                {
                    if(selectedFaceIDs.count(face.second->sID))
                    {
                        selectedSimplices.insert(face.second->simplices().begin(),face.second->simplices().end());
                    }
                }
            }
            std::vector<const Simplex<3, 2> *> faces(selectedSimplices.begin(),selectedSimplices.end());
            double largestArea=0.0;
            for(const auto* face : faces)
            {
                largestArea=std::max(largestArea,face->vol0);
            }
            if(faces.empty())
            {
                return {};
            }
            // A common subdivision count makes all shared edge samples identical.
            const double subdivisions = std::max(1.0, std::ceil(std::sqrt(largestArea / maxArea)));
            if(!std::isfinite(subdivisions) || subdivisions * subdivisions * faces.size() > 500000)
            {
                throw std::invalid_argument("Surface would exceed 500000 triangles; increase mesh size");
            }
            const size_t n = static_cast<size_t>(subdivisions);
            SurfaceMeshGeometry result;
            using Key = std::vector<std::pair<size_t, size_t>>;
            std::map<Key, size_t> pointIDs;
            for(const auto *face : faces)
            {
                const auto vertices = face->vertices();
                std::map<std::pair<size_t, size_t>, size_t> local;
                for(size_t i = 0; i <= n; ++i)
                {
                    for(size_t j = 0; j <= n - i; ++j)
                    {
                        const std::array<size_t, 3> weights{n - i - j, i, j};
                        Key key;
                        Vector position = Vector::Zero();
                        for(size_t k = 0; k < 3; ++k)
                        {
                            if(weights[k] != 0)
                            {
                                key.emplace_back(face->xID[k], weights[k]);
                                position += (static_cast<double>(weights[k]) / n) * vertices[k]->P0;
                            }
                        }
                        const auto inserted = pointIDs.emplace(key, result.referencePoints.size());
                        if(inserted.second)
                        {
                            result.referencePoints.push_back(position);
                            result.parents.push_back(samplingParent(face));
                        }
                        local.emplace(std::make_pair(i, j), inserted.first->second);
                    }
                }
                const Vector normal = face->outNormal(samplingParent(face)->region->regionID);
                auto append = [&](size_t a, size_t b, size_t c)
                {
                    if((result.referencePoints[b] - result.referencePoints[a])
                           .cross(result.referencePoints[c] - result.referencePoints[a])
                           .dot(normal) < 0)
                    {
                        std::swap(b, c);
                    }
                    result.triangles.push_back({a, b, c});
                    result.boundaryFaces.push_back(face);
                };
                for(size_t i = 0; i < n; ++i)
                {
                    for(size_t j = 0; j < n - i; ++j)
                    {
                        append(local.at({i, j}), local.at({i + 1, j}), local.at({i, j + 1}));
                        if(i + j + 1 < n)
                        {
                            append(local.at({i + 1, j}), local.at({i + 1, j + 1}), local.at({i, j + 1}));
                        }
                    }
                }
            }
            result.samplePoints = result.referencePoints;
            return result;
        }

        // Cut each refined boundary triangle into convex regions. Coincident trace
        // vertices deliberately have separate IDs and one-sided field samples.
        SurfaceMeshGeometry split(const std::vector<Trace> &inputTraces) const
        {
            if(inputTraces.empty())
            {
                return *this;
            }
            double scale = 1.0;
            for(const auto &p : referencePoints)
            {
                scale = std::max(scale, p.norm());
            }
            const double tolerance = scale * 1.e-10;
            std::vector<Trace> traces;
            for(const auto &trace : inputTraces)
            {
                bool merged = false;
                for(auto &existing : traces)
                {
                    if(existing.normal.cross(trace.normal).norm() < 1.e-10 &&
                       std::abs(existing.normal.dot(trace.point - existing.point)) < tolerance)
                    {
                        existing.faces.insert(trace.faces.begin(), trace.faces.end());
                        merged = true;
                        break;
                    }
                }
                if(!merged)
                {
                    traces.push_back(trace);
                }
            }
            SurfaceMeshGeometry result;
            using Polygon = std::vector<Vector>;
            using EdgeKey = std::array<long long, 10>;
            struct Edge
            {
                size_t a;
                size_t b;
                bool positive;
            };
            std::map<EdgeKey, Edge> seams;
            for(size_t t = 0; t < triangles.size(); ++t)
            {
                const auto *face = boundaryFaces[t];
                const auto &tri = triangles[t];
                std::vector<size_t> active;
                for(size_t k = 0; k < traces.size(); ++k)
                {
                    if(traces[k].faces.count(face))
                    {
                        active.push_back(k);
                    }
                }
                std::vector<Polygon> polygons{
                    {referencePoints[tri[0]], referencePoints[tri[1]], referencePoints[tri[2]]}};
                for(const auto k : active)
                {
                    const auto &trace = traces[k];
                    std::vector<Polygon> next;
                    for(const auto &polygon : polygons)
                    {
                        double lo = 0.0;
                        double hi = 0.0;
                        for(const auto &p : polygon)
                        {
                            const double d = trace.normal.dot(p - trace.point);
                            lo = std::min(lo, d);
                            hi = std::max(hi, d);
                        }
                        if(lo >= -tolerance || hi <= tolerance)
                        {
                            next.push_back(polygon);
                            continue;
                        }
                        Polygon positive;
                        Polygon negative;
                        for(size_t j = 0; j < polygon.size(); ++j)
                        {
                            const Vector &a = polygon[j];
                            const Vector &b = polygon[(j + 1) % polygon.size()];
                            double da = trace.normal.dot(a - trace.point);
                            double db = trace.normal.dot(b - trace.point);
                            if(std::abs(da) <= tolerance)
                            {
                                da = 0.0;
                            }
                            if(std::abs(db) <= tolerance)
                            {
                                db = 0.0;
                            }
                            if(da >= 0.0)
                            {
                                positive.push_back(a);
                            }
                            if(da <= 0.0)
                            {
                                negative.push_back(a);
                            }
                            if(da * db < 0.0)
                            {
                                const Vector intersection = a + da / (da - db) * (b - a);
                                positive.push_back(intersection);
                                negative.push_back(intersection);
                            }
                        }
                        next.push_back(std::move(positive));
                        next.push_back(std::move(negative));
                    }
                    polygons = std::move(next);
                    if(polygons.size() + result.triangles.size() > 500000)
                    {
                        throw std::runtime_error("Too many slip-trace cuts; increase surface mesh size");
                    }
                }
                for(const auto &polygon : polygons)
                {
                    Vector center = Vector::Zero();
                    for(const auto &p : polygon)
                    {
                        center += p;
                    }
                    center /= polygon.size();
                    const size_t first = result.referencePoints.size();
                    for(const auto &p : polygon)
                    {
                        double fraction = 0.0;
                        for(const auto k : active)
                        {
                            if(std::abs(traces[k].normal.dot(p - traces[k].point)) <= tolerance)
                            {
                                const double distance = std::abs(traces[k].normal.dot(center - p));
                                // Stay outside the relative FLT_EPSILON dead band used
                                // by the dislocation solid-angle evaluation.
                                if(distance > tolerance)
                                {
                                    fraction = std::max(fraction, std::min(0.25, 2.e-6 * scale / distance));
                                }
                            }
                        }
                        result.referencePoints.push_back(p);
                        // A convex combination stays in the original boundary simplex
                        // and its FEM parent, including at intersecting slip traces.
                        result.samplePoints.push_back(Vector(p + fraction * (center - p)));
                        result.parents.push_back(samplingParent(face));
                    }
                    for(size_t j = 1; j + 1 < polygon.size(); ++j)
                    {
                        result.triangles.push_back({first, first + j, first + j + 1});
                        result.boundaryFaces.push_back(face);
                    }
                    for(size_t j = 0; j < polygon.size(); ++j)
                    {
                        const size_t next = (j + 1) % polygon.size();
                        const Vector &a = polygon[j];
                        const Vector &b = polygon[next];
                        for(const auto k : active)
                        {
                            const auto &trace = traces[k];
                            const double side = trace.normal.dot(center - trace.point);
                            if(std::abs(side) <= tolerance || std::abs(trace.normal.dot(a - trace.point)) > tolerance ||
                               std::abs(trace.normal.dot(b - trace.point)) > tolerance)
                            {
                                continue;
                            }
                            std::array<long long, 3> qa, qb;
                            for(size_t d = 0; d < 3; ++d)
                            {
                                qa[d] = std::llround(a[d] / tolerance);
                                qb[d] = std::llround(b[d] / tolerance);
                            }
                            if(qb < qa)
                            {
                                std::swap(qa, qb);
                            }
                            const Vector normal = face->outNormal(samplingParent(face)->region->regionID);
                            const EdgeKey key{static_cast<long long>(k),
                                              qa[0],
                                              qa[1],
                                              qa[2],
                                              qb[0],
                                              qb[1],
                                              qb[2],
                                              std::llround(normal[0] * 1.e8),
                                              std::llround(normal[1] * 1.e8),
                                              std::llround(normal[2] * 1.e8)};
                            const Edge edge{first + j, first + next, side > 0.0};
                            const auto inserted = seams.emplace(key, edge);
                            if(!inserted.second)
                            {
                                const auto other = inserted.first->second;
                                if(other.positive != edge.positive)
                                {
                                    // Adjacent outward polygons traverse their common edge
                                    // oppositely. Reverse those edges on the connecting wall.
                                    result.walls.push_back({other.b, other.a, edge.b});
                                    result.walls.push_back({other.b, edge.b, edge.a});
                                }
                                seams.erase(inserted.first);
                            }
                        }
                    }
                    if(result.triangles.size() + result.walls.size() > 500000)
                    {
                        throw std::runtime_error("Surface with slip walls exceeds 500000 triangles");
                    }
                }
            }
            return result;
        }

        std::vector<Vector> displaced(const std::vector<Vector> &displacement, double amplification) const
        {
            if(displacement.size() != referencePoints.size() || !std::isfinite(amplification) || amplification < 0)
            {
                throw std::invalid_argument("Invalid surface displacement or amplification");
            }
            std::vector<Vector> points;
            points.reserve(referencePoints.size());
            for(size_t k = 0; k < referencePoints.size(); ++k)
            {
                const Vector point = referencePoints[k] + amplification * displacement[k];
                if(!point.allFinite())
                {
                    throw std::runtime_error("Nonfinite surface displacement; check the field or amplification");
                }
                points.push_back(point);
            }
            return points;
        }
    };
} // namespace model
#endif
