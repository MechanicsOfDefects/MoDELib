/* Dimension-generic template definitions for PolyhedralMesh. */
#ifndef model_PolyhedralMeshImplementation_H_
#define model_PolyhedralMeshImplementation_H_
namespace model
{
    // Only normal orientation and the corresponding grain checks are specialized.
    template <> void PolyhedralMesh<3>::orientMeshFaces();
    template <> void PolyhedralMesh<3>::validate() const;
    template <int dim> std::shared_ptr<Polyhedron<dim, 0>> PolyhedralFaceIncidence<dim>::source() const
    {
        return reversed ? edge->children.rbegin()->second : edge->children.begin()->second;
    }
    template <int dim> std::shared_ptr<Polyhedron<dim, 0>> PolyhedralFaceIncidence<dim>::sink() const
    {
        return reversed ? edge->children.begin()->second : edge->children.rbegin()->second;
    }
    template <int dim> Eigen::Matrix<double, dim, 1> PolyhedralFaceIncidence<dim>::tangent() const
    {
        return (sink()->position - source()->position).normalized();
    }
    template <int dim> PolyhedralFaceIncidence<dim> PolyhedralFaceIncidence<dim>::next() const
    {
        const auto owner = face.lock();
        if (!owner)
        {
            throw std::runtime_error("Expired polyhedral face incidence");
        }
        return owner->perimeter.at((index + 1) % owner->perimeter.size());
    }
    template <int dim> PolyhedralFaceIncidence<dim> PolyhedralFaceIncidence<dim>::prev() const
    {
        const auto owner = face.lock();
        if (!owner)
        {
            throw std::runtime_error("Expired polyhedral face incidence");
        }
        return owner->perimeter.at((index + owner->perimeter.size() - 1) % owner->perimeter.size());
    }
    template <int dim> const Polyhedron<dim, dim> &PolyhedralMesh<dim>::region(size_t regionID) const
    {
        return *factory<dim>().at({regionID});
    }
    template <int dim> const Polyhedron<dim, dim - 1> &PolyhedralMesh<dim>::face(size_t faceID) const
    {
        return *factory<dim - 1>().at({faceID});
    }

    template <int dim>
    template <int order>
    void PolyhedralMesh<dim>::buildBoundary(const std::shared_ptr<Polyhedron<dim, order>> &cell,
                                            const std::vector<std::vector<size_t>> &simplices,
                                            const std::map<size_t, Vector> &positions)
    {
        if constexpr (order == 0)
        {
            if (simplices.size() != 1 || simplices.front().size() != 1)
            {
                throw std::runtime_error("A polyhedral vertex requires one mesh vertex");
            }
            cell->position = positions.at(simplices.front().front());
        }
        else
        {
            using Key = std::vector<size_t>;
            std::map<Key, size_t> counts;
            std::set<Key> unique;
            for (auto simplex : simplices)
            {
                std::sort(simplex.begin(), simplex.end());
                if (simplex.size() != order + 1 ||
                    std::adjacent_find(simplex.begin(), simplex.end()) != simplex.end() ||
                    !unique.insert(simplex).second)
                {
                    throw std::runtime_error("Invalid or duplicate boundary simplex");
                }
                for (const auto id : simplex)
                {
                    positions.at(id);
                }
                for (size_t k = 0; k < simplex.size(); ++k)
                {
                    auto child = simplex;
                    child.erase(child.begin() + k);
                    ++counts[child];
                }
            }
            for (const auto &entry : counts)
            {
                if (entry.second > 2)
                {
                    throw std::runtime_error("Non-manifold boundary simplex");
                }
                if (entry.second == 1)
                {
                    const auto child = std::get<order - 1>(factories).get(entry.first);
                    cell->children.emplace(child->key, child);
                    child->parents.emplace(cell->key, cell);
                    buildBoundary<order - 1>(child, {entry.first}, positions);
                }
            }
            if constexpr (order == 1)
            {
                if (cell->children.size() != 2)
                {
                    throw std::runtime_error("A polyhedral edge requires two endpoints");
                }
                // Reject disconnected cycles hidden beside a valid open chain.
                std::map<size_t, std::set<size_t>> adjacency;
                for (const auto &simplex : simplices)
                {
                    adjacency[simplex[0]].insert(simplex[1]);
                    adjacency[simplex[1]].insert(simplex[0]);
                }
                std::set<size_t> visited;
                std::vector<size_t> pending{simplices.front().front()};
                while (!pending.empty())
                {
                    const auto id = pending.back();
                    pending.pop_back();
                    if (visited.insert(id).second)
                    {
                        if (adjacency.at(id).size() > 2)
                        {
                            throw std::runtime_error("Branched polyhedral edge");
                        }
                        for (const auto neighbor : adjacency.at(id))
                        {
                            pending.push_back(neighbor);
                        }
                    }
                }
                if (visited.size() != adjacency.size())
                {
                    throw std::runtime_error("Disconnected polyhedral edge");
                }
            }
            if constexpr (order == 2)
            {
                buildPerimeter(cell);
            }
        }
    }

    template <int dim> void PolyhedralMesh<dim>::buildPerimeter(const std::shared_ptr<Polyhedron<dim, 2>> &cell)
    {
        using Key = PolyhedronKey<dim, 0>;
        if (cell->perimeter.empty())
        {
            std::map<Key, std::vector<std::shared_ptr<Polyhedron<dim, 1>>>> adjacency;
            for (const auto &entry : cell->children)
            {
                for (const auto &vertex : entry.second->children)
                {
                    adjacency[vertex.first].push_back(entry.second);
                }
            }
            for (const auto &entry : adjacency)
            {
                if (entry.second.size() != 2)
                {
                    throw std::runtime_error("Open or branched polygon perimeter");
                }
            }
            if (cell->children.size() < 3)
            {
                throw std::runtime_error("Polygon requires at least three edges");
            }
            const auto first = cell->children.begin()->second;
            auto edge = first;
            auto vertex = first->children.begin()->second;
            const auto start = vertex;
            std::set<PolyhedronKey<dim, 1>> visited;
            do
            {
                if (!visited.insert(edge->key).second)
                {
                    throw std::runtime_error("Polygon perimeter repeats an edge");
                }
                cell->perimeter.push_back(
                    {cell, edge, cell->perimeter.size(), vertex != edge->children.begin()->second});
                edge->faceIncidences.emplace(cell->key, cell->perimeter.size() - 1);
                vertex = cell->perimeter.back().sink();
                const auto &neighbors = adjacency.at(vertex->key);
                edge = neighbors[0] == edge ? neighbors[1] : neighbors[0];
            } while (edge != first);
            if (vertex != start || visited.size() != cell->children.size())
            {
                throw std::runtime_error("Polygon requires one closed perimeter");
            }
            if constexpr (dim == 2)
            {
                double area = 0;
                const auto origin = start->position;
                for (const auto &incidence : cell->perimeter)
                {
                    const Vector a = incidence.source()->position - origin;
                    const Vector b = incidence.sink()->position - origin;
                    area += a[0] * b[1] - a[1] * b[0];
                }
                if (area < 0)
                {
                    std::reverse(cell->perimeter.begin(), cell->perimeter.end());
                    for (auto &incidence : cell->perimeter)
                    {
                        incidence.reversed = !incidence.reversed;
                    }
                }
                for (size_t k = 0; k < cell->perimeter.size(); ++k)
                {
                    cell->perimeter[k].index = k;
                    cell->perimeter[k].edge->faceIncidences.at(cell->key) = k;
                }
            }
        }
    }
    template <int dim> void PolyhedralMesh<dim>::orientMeshFaces()
    {
        // Higher-dimensional orientation is not needed for parent/child topology.
    }

    template <int dim>
    PolyhedralMesh<dim>::PolyhedralMesh(const std::map<size_t, Vector> &positions, const std::vector<FaceInput> &faces)
        : PolyhedralMesh()
    {
        std::set<size_t> ids;
        for (const auto &input : faces)
        {
            if (!ids.insert(input.id).second || input.regions.empty() || input.regions.size() > 2 ||
                (input.regions.size() == 2 && input.regions[0] == input.regions[1]) || input.simplices.empty())
            {
                throw std::runtime_error("Invalid polyhedral face input");
            }
            const auto face = std::get<dim - 1>(factories).get({input.id});
            std::vector<std::vector<size_t>> simplices;
            for (const auto &simplex : input.simplices)
            {
                simplices.emplace_back(simplex.begin(), simplex.end());
            }
            buildBoundary<dim - 1>(face, simplices, positions);
            for (const auto id : input.regions)
            {
                const auto region = std::get<dim>(factories).get({id});
                region->children.emplace(face->key, face);
                face->parents.emplace(region->key, region);
            }
        }
        if constexpr (dim == 2)
        {
            for (const auto &region : factory<2>())
            {
                buildPerimeter(region.second);
            }
        }
        validateLevel<dim>();
    }

    template <int dim>
    std::shared_ptr<Polyhedron<dim, dim - 1>> PolyhedralMesh<dim>::createFace(const Simplex<dim, dim - 1> *simplex)
    {
        const StaticID<Polyhedron<dim, dim - 1>> identifier;
        const size_t id = identifier.sID;
        const auto face = std::get<dim - 1>(factories).get(id);
        face->initialize(simplex);
        for (const auto regionID : simplex->regionIDs())
        {
            const auto region = getRegion(regionID);
            face->parents.emplace(regionID, region);
            region->children.emplace(id, face);
        }
        return face;
    }
    template <int dim> template <int order> void PolyhedralMesh<dim>::attachSimplex(const Simplex<dim, order> *simplex)
    {
        if constexpr (order < dim - 1)
        {
            const std::vector<size_t> key(simplex->xID.begin(), simplex->xID.end());
            const auto found = std::get<order>(factories).find(key);
            if (found != std::get<order>(factories).end())
            {
                found->second->simplices().insert(simplex);
            }
        }
        if constexpr (order > 0)
        {
            for (const auto &child : simplex->children())
            {
                attachSimplex<order - 1>(child.get());
            }
        }
    }
    template <int dim> void PolyhedralMesh<dim>::buildTopology()
    {
        for (const auto &entry : factory<dim - 1>())
        {
            const auto &face = entry.second;
            std::map<size_t, Vector> positions;
            std::vector<std::vector<size_t>> simplices;
            for (const auto *simplex : face->simplices())
            {
                simplices.emplace_back(simplex->xID.begin(), simplex->xID.end());
                for (const auto *vertex : simplex->vertices())
                {
                    positions.emplace(vertex->xID[0], vertex->P0);
                }
            }
            buildBoundary<dim - 1>(face, simplices, positions);
            for (const auto *simplex : face->simplices())
            {
                attachSimplex<dim - 1>(simplex);
            }
        }
        if constexpr (dim == 2)
        {
            for (const auto &region : factory<dim>())
            {
                buildPerimeter(region.second);
            }
        }
        orientMeshFaces();
        validate();
    }

    template <int dim> template <int order> void PolyhedralMesh<dim>::validateLevel() const
    {
        for (const auto &entry : factory<order>())
        {
            const auto &cell = *entry.second;
            if constexpr (order > 0)
            {
                if (cell.children.size() < size_t(order + 1))
                {
                    throw std::runtime_error("Polyhedral cell has too few boundary cells");
                }
                for (const auto &child : cell.children)
                {
                    if (factory<order - 1>().at(child.first) != child.second ||
                        child.second->parents.at(cell.key).lock() != entry.second)
                    {
                        throw std::runtime_error("Nonreciprocal polyhedral child incidence");
                    }
                }
                if constexpr (order > 1)
                {
                    std::map<PolyhedronKey<dim, order - 2>, size_t> counts;
                    for (const auto &child : cell.children)
                    {
                        for (const auto &grandchild : child.second->children)
                        {
                            ++counts[grandchild.first];
                        }
                    }
                    for (const auto &count : counts)
                    {
                        if (count.second != 2)
                        {
                            throw std::runtime_error("Open or non-manifold polyhedral cell boundary");
                        }
                    }
                }
            }
            if constexpr (order < dim)
            {
                if (cell.parents.empty() || (order == dim - 1 && cell.parents.size() > 2))
                {
                    throw std::runtime_error("Invalid polyhedral parent count");
                }
                for (const auto &parent : cell.parents)
                {
                    const auto owner = parent.second.lock();
                    if (!owner || factory<order + 1>().at(parent.first) != owner ||
                        owner->children.at(cell.key) != entry.second)
                    {
                        throw std::runtime_error("Nonreciprocal polyhedral parent incidence");
                    }
                }
            }
        }
        if constexpr (order > 0)
        {
            validateLevel<order - 1>();
        }
    }
    template <int dim> void PolyhedralMesh<dim>::validate() const
    {
        validateLevel<dim>();
    }
} // namespace model
#endif
