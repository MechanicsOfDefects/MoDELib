#include <FMM/PlaneBoxIntersection.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace model
{
    template<int dim>
    typename PlaneBoxIntersection<dim>::NodePtr PlaneBoxIntersection<dim>::node(
        const std::shared_ptr<const Entity>& support, const Vector& position)
    {
        auto found = nodesValue.find(support->key());
        if(found == nodesValue.end())
        {
            found = nodesValue.emplace(support->key(), std::make_shared<Node>(Node{support, position})).first;
        }
        return found->second;
    }

    template<int dim>
    typename PlaneBoxIntersection<dim>::SegmentPtr PlaneBoxIntersection<dim>::segment(NodePtr a, NodePtr b)
    {
        if(b->support->key() < a->support->key())
        {
            std::swap(a, b);
        }
        const SegmentKey key(a->support->key(), b->support->key());
        auto found = segmentsValue.find(key);
        if(found == segmentsValue.end())
        {
            found = segmentsValue.emplace(key, std::make_shared<Segment>(Segment{a, b})).first;
        }
        return found->second;
    }

    template<int dim>
    const std::vector<typename PlaneBoxIntersection<dim>::NodePtr>&
    PlaneBoxIntersection<dim>::cutEdge(const std::shared_ptr<const Entity>& edge)
    {
        auto [found, inserted] = edgeCuts.emplace(edge->key(), std::vector<NodePtr>{});
        if(!inserted)
        {
            return found->second;
        }
        const auto a = edge->boundary().at(0).entity;
        const auto b = edge->boundary().at(1).entity;
        const double da = normalValue.dot(a->lower() - originValue);
        const double db = normalValue.dot(b->lower() - originValue);
        if(std::abs(da) <= tolerance)
        {
            found->second.push_back(node(a, a->lower()));
        }
        if(std::abs(db) <= tolerance)
        {
            found->second.push_back(node(b, b->lower()));
        }
        if(found->second.empty() && ((da < 0 && db > 0) || (da > 0 && db < 0)))
        {
            found->second.push_back(node(edge, a->lower() + (da / (da - db)) * (b->lower() - a->lower())));
        }
        return found->second;
    }

    template<int dim>
    const std::vector<typename PlaneBoxIntersection<dim>::SegmentPtr>&
    PlaneBoxIntersection<dim>::cutFace(const std::shared_ptr<const Entity>& face)
    {
        auto [found, inserted] = faceCutsValue.emplace(face->key(), std::vector<SegmentPtr>{});
        if(!inserted)
        {
            return found->second;
        }
        std::map<Key, NodePtr> points;
        bool coplanar = true;
        for(const auto& incidence : face->boundary())
        {
            const auto& cut = cutEdge(incidence.entity);
            coplanar = coplanar && cut.size() == 2;
            for(const auto& point : cut)
            {
                points.emplace(point->support->key(), point);
            }
        }
        if(coplanar)
        {
            for(const auto& incidence : face->boundary())
            {
                const auto& cut = cutEdge(incidence.entity);
                found->second.push_back(segment(cut[0], cut[1]));
            }
        }
        else if(points.size() == 2)
        {
            found->second.push_back(segment(points.begin()->second, points.rbegin()->second));
        }
        else if(points.size() > 2)
        {
            throw std::runtime_error("Inconsistent plane/box face intersection");
        }
        return found->second;
    }

    template<int dim>
    PlaneBoxIntersection<dim>::PlaneBoxIntersection(const BoxTree<dim>& tree, unsigned level,
        const Vector& origin, const Vector& normal)
        : originValue(origin), normalValue(normal), levelValue(level), revisionValue(tree.revision())
    {
        if(!origin.allFinite() || !normal.allFinite() || normal.stableNorm() == 0)
        {
            throw std::invalid_argument("Plane/box cut needs a finite origin and nonzero finite normal");
        }
        normalValue /= normal.stableNorm();
        const double scale = std::max({tree.lower().cwiseAbs().maxCoeff(), tree.upper().cwiseAbs().maxCoeff(),
                                      origin.cwiseAbs().maxCoeff(), (tree.upper() - tree.lower()).norm()});
        tolerance = 64 * std::numeric_limits<double>::epsilon() * scale;
        for(const auto& box : tree.boxes(level))
        {
            const Vector halfSize = (box->upper() - box->lower()) * 0.5;
            if(std::abs(normalValue.dot(box->center() - originValue)) > normalValue.cwiseAbs().dot(halfSize) + tolerance)
            {
                continue;
            }
            // Half-open ownership: a plane coincident with an internal upper
            // face belongs to the neighbor on its positive-coordinate side.
            bool owned = true;
            for(int axis = 0; axis < dim; ++axis)
            {
                const auto& face = box->face(axis, 1).entity;
                bool coplanar = true;
                for(const auto& edge : face->boundary())
                {
                    coplanar = coplanar && cutEdge(edge.entity).size() == 2;
                }
                if(coplanar && box->key()[dim + axis] != BoxTree<dim>::latticeSize)
                {
                    owned = false;
                }
            }
            if(!owned)
            {
                continue;
            }
            std::map<SegmentKey, SegmentPtr> sides;
            for(const auto& face : box->boundary())
            {
                for(const auto& side : cutFace(face.entity))
                {
                    sides.emplace(SegmentKey(side->first->support->key(), side->second->support->key()), side);
                }
            }
            if(sides.size() < 3) // point/edge tangencies have no area
            {
                continue;
            }
            std::map<Key, std::vector<SegmentPtr>> adjacency;
            for(const auto& [key, side] : sides)
            {
                adjacency[key.first].push_back(side);
                adjacency[key.second].push_back(side);
            }
            for(const auto& [key, links] : adjacency)
            {
                if(links.size() != 2)
                {
                    throw std::runtime_error("Plane/box cut is not a closed topological contour");
                }
            }
            Polygon polygon{box, {}, {}, 0};
            const auto start = sides.begin()->second->first;
            auto current = start;
            SegmentPtr previous;
            do
            {
                polygon.vertices.push_back(current);
                const auto& links = adjacency.at(current->support->key());
                const auto next = links[0] == previous ? links[1] : links[0];
                const int sign = next->first == current ? 1 : -1;
                polygon.boundary.push_back(Side{next, sign});
                current = sign == 1 ? next->second : next->first;
                previous = next;
                if(polygon.boundary.size() > sides.size())
                {
                    throw std::runtime_error("Plane/box contour failed to close");
                }
            } while(current != start);
            if(polygon.boundary.size() != sides.size())
            {
                throw std::runtime_error("Plane/box contour has disconnected components");
            }
            Vector areaVector = Vector::Zero();
            for(size_t i = 1; i + 1 < polygon.vertices.size(); ++i)
            {
                areaVector += (polygon.vertices[i]->position - start->position).cross(
                    polygon.vertices[i + 1]->position - start->position) * 0.5;
            }
            const double signedArea = normalValue.dot(areaVector);
            if(std::abs(signedArea) <= tolerance * tolerance)
            {
                continue;
            }
            if(signedArea < 0)
            {
                // Reverse the walk, retaining its starting node.
                std::reverse(polygon.vertices.begin() + 1, polygon.vertices.end());
                std::reverse(polygon.boundary.begin(), polygon.boundary.end());
                for(auto& side : polygon.boundary)
                {
                    side.orientation *= -1;
                }
            }
            polygon.area = std::abs(signedArea);
            polygonsValue.emplace(box->key(), std::move(polygon));
        }
    }
    template<int dim> const typename PlaneBoxIntersection<dim>::Vector& PlaneBoxIntersection<dim>::origin() const
    {
        return originValue;
    }
    template<int dim> const typename PlaneBoxIntersection<dim>::Vector& PlaneBoxIntersection<dim>::normal() const
    {
        return normalValue;
    }
    template<int dim> unsigned PlaneBoxIntersection<dim>::level() const
    {
        return levelValue;
    }
    template<int dim> std::uint64_t PlaneBoxIntersection<dim>::revision() const
    {
        return revisionValue;
    }
    template<int dim> const std::map<typename PlaneBoxIntersection<dim>::Key, typename PlaneBoxIntersection<dim>::Polygon>& PlaneBoxIntersection<dim>::polygons() const
    {
        return polygonsValue;
    }
    template<int dim> const std::map<typename PlaneBoxIntersection<dim>::Key, std::vector<typename PlaneBoxIntersection<dim>::SegmentPtr>>& PlaneBoxIntersection<dim>::faceCuts() const
    {
        return faceCutsValue;
    }
    template class PlaneBoxIntersection<3>;
}
