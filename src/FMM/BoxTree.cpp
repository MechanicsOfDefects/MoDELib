#include <FMM/BoxTree.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace model
{
    template<int dim>
    BoxEntity<dim>::BoxEntity(const Key& key, const Vector& lower, const Vector& upper)
        : keyValue(key), lowerValue(lower), upperValue(upper), orderValue(0)
    {
        for(int axis = 0; axis < dim; ++axis)
        {
            orderValue += key[axis] != key[dim + axis];
        }
    }

    template<int dim> const typename BoxEntity<dim>::Key& BoxEntity<dim>::key() const
    {
        return keyValue;
    }
    template<int dim> int BoxEntity<dim>::order() const
    {
        return orderValue;
    }
    template<int dim> const typename BoxEntity<dim>::Vector& BoxEntity<dim>::lower() const
    {
        return lowerValue;
    }
    template<int dim> const typename BoxEntity<dim>::Vector& BoxEntity<dim>::upper() const
    {
        return upperValue;
    }
    template<int dim> typename BoxEntity<dim>::Vector BoxEntity<dim>::center() const
    {
        return lowerValue + 0.5 * (upperValue - lowerValue);
    }
    template<int dim> double BoxEntity<dim>::measure() const
    {
        double result = 1.0;
        for(int axis = 0; axis < dim; ++axis)
        {
            if(keyValue[axis] != keyValue[dim + axis])
            {
                result *= upperValue[axis] - lowerValue[axis];
            }
        }
        return result;
    }
    template<int dim> const std::vector<typename BoxEntity<dim>::Incidence>& BoxEntity<dim>::boundary() const
    {
        return boundaryValue;
    }
    template<int dim> const std::map<typename BoxEntity<dim>::Key, std::weak_ptr<BoxEntity<dim>>>&
    BoxEntity<dim>::parents() const
    {
        return parentsValue;
    }
    template<int dim> const std::vector<std::shared_ptr<BoxEntity<dim>>>& BoxEntity<dim>::subdivisions() const
    {
        return subdivisionsValue;
    }
    template<int dim> std::shared_ptr<BoxEntity<dim>> BoxEntity<dim>::refinementParent() const
    {
        return refinementParentValue.lock();
    }
    template<int dim> const typename BoxEntity<dim>::Incidence& BoxEntity<dim>::face(int axis, int side) const
    {
        for(const auto& incidence : boundaryValue)
        {
            if(incidence.axis == axis && incidence.side == side)
            {
                return incidence;
            }
        }
        throw std::out_of_range("Box entity has no boundary incidence for this axis and side");
    }

    template<int dim>
    BoxTree<dim>::BoxTree(const Vector& lower, const Vector& upper, unsigned depth)
        : lowerValue(lower), upperValue(upper)
    {
        const Vector span = upper - lower;
        if(!lower.allFinite() || !upper.allFinite() || !span.allFinite() || (span.array() <= 0.0).any())
        {
            throw std::invalid_argument("BoxTree bounds must be finite and strictly increasing");
        }
        setDepth(depth);
    }
    template<int dim> const typename BoxTree<dim>::Vector& BoxTree<dim>::lower() const
    {
        return lowerValue;
    }
    template<int dim> const typename BoxTree<dim>::Vector& BoxTree<dim>::upper() const
    {
        return upperValue;
    }
    template<int dim> unsigned BoxTree<dim>::depth() const
    {
        return static_cast<unsigned>(state.levels.size() - 1);
    }
    template<int dim> std::uint64_t BoxTree<dim>::revision() const
    {
        return revisionValue;
    }
    template<int dim> const typename BoxTree<dim>::Factory& BoxTree<dim>::entities() const
    {
        return state.factory;
    }
    template<int dim> const std::vector<typename BoxTree<dim>::Pointer>& BoxTree<dim>::boxes(unsigned level) const
    {
        return state.levels.at(level);
    }
    template<int dim> const typename BoxTree<dim>::Pointer& BoxTree<dim>::root() const
    {
        return state.levels.front().front();
    }

    template<int dim>
    typename BoxTree<dim>::Pointer BoxTree<dim>::get(State& target, const Key& key) const
    {
        const auto found = target.factory.find(key);
        if(found != target.factory.end())
        {
            return found->second;
        }
        Vector lower;
        Vector upper;
        for(int axis = 0; axis < dim; ++axis)
        {
            lower[axis] = lowerValue[axis] + (upperValue[axis] - lowerValue[axis]) *
                (static_cast<double>(key[axis]) / latticeSize);
            upper[axis] = lowerValue[axis] + (upperValue[axis] - lowerValue[axis]) *
                (static_cast<double>(key[dim + axis]) / latticeSize);
        }
        Pointer entity(new Entity(key, lower, upper));
        target.factory.emplace(key, entity);
        int rank = 0;
        for(int axis = 0; axis < dim; ++axis)
        {
            if(key[axis] != key[dim + axis])
            {
                for(int side = 0; side < 2; ++side)
                {
                    Key childKey = key;
                    childKey[axis] = childKey[dim + axis] = key[side * dim + axis];
                    auto child = get(target, childKey);
                    const int orientation = (rank % 2 == 0 ? 1 : -1) * (side == 1 ? 1 : -1);
                    entity->boundaryValue.push_back({child, axis, side, orientation});
                    child->parentsValue.emplace(key, entity);
                }
                ++rank;
            }
        }
        return entity;
    }

    template<int dim>
    void BoxTree<dim>::subdivide(State& target, const Pointer& entity) const
    {
        if(entity->order() == 0 || !entity->subdivisionsValue.empty())
        {
            return;
        }
        std::vector<int> axes;
        for(int axis = 0; axis < dim; ++axis)
        {
            if(entity->key()[axis] != entity->key()[dim + axis])
            {
                axes.push_back(axis);
            }
        }
        for(unsigned mask = 0; mask < (1u << axes.size()); ++mask)
        {
            Key key = entity->key();
            for(size_t rank = 0; rank < axes.size(); ++rank)
            {
                const int axis = axes[rank];
                const auto midpoint = (key[axis] + key[dim + axis]) / 2;
                if(mask & (1u << rank))
                {
                    key[axis] = midpoint;
                }
                else
                {
                    key[dim + axis] = midpoint;
                }
            }
            auto child = get(target, key);
            child->refinementParentValue = entity;
            entity->subdivisionsValue.push_back(child);
        }
        // The same factory identifies each coarse-face subdivision with the
        // corresponding fine-box face; edges and vertices are shared as well.
        for(const auto& incidence : entity->boundary())
        {
            subdivide(target, incidence.entity);
        }
    }

    template<int dim> void BoxTree<dim>::setDepth(unsigned requested)
    {
        if(requested > maxDepth)
        {
            throw std::invalid_argument("BoxTree uniform depth exceeds the supported memory bound");
        }
        if(!state.levels.empty() && depth() == requested)
        {
            return;
        }
        // Validate representability before allocating or changing the live tree.
        for(int axis = 0; axis < dim; ++axis)
        {
            const double step = (upperValue[axis] - lowerValue[axis]) / (1u << requested);
            if(lowerValue[axis] + step == lowerValue[axis] || upperValue[axis] - step == upperValue[axis])
            {
                throw std::invalid_argument("BoxTree depth is not resolvable at these physical coordinates");
            }
        }
        State next;
        next.levels.resize(requested + 1);
        Key rootKey{};
        for(int axis = 0; axis < dim; ++axis)
        {
            rootKey[dim + axis] = latticeSize;
        }
        next.levels[0].push_back(get(next, rootKey));
        for(unsigned level = 0; level < requested; ++level)
        {
            for(const auto& box : next.levels[level])
            {
                subdivide(next, box);
                next.levels[level + 1].insert(next.levels[level + 1].end(),
                    box->subdivisions().begin(), box->subdivisions().end());
            }
        }
        state.factory.swap(next.factory);
        state.levels.swap(next.levels);
        ++revisionValue;
    }

    template<int dim>
    typename BoxTree<dim>::Pointer BoxTree<dim>::neighbor(const Pointer& box, int axis, int side) const
    {
        if(!box || box->order() != dim)
        {
            throw std::invalid_argument("BoxTree neighbor requires a full-dimensional box");
        }
        const auto own = state.factory.find(box->key());
        if(own == state.factory.end() || own->second != box)
        {
            throw std::invalid_argument("BoxTree neighbor received a foreign or stale box");
        }
        for(const auto& entry : box->face(axis, side).entity->parents())
        {
            auto parent = entry.second.lock();
            if(parent && parent != box && parent->order() == dim &&
               parent->key()[dim] - parent->key()[0] == box->key()[dim] - box->key()[0])
            {
                return parent;
            }
        }
        return nullptr;
    }

    template<int dim>
    typename BoxTree<dim>::Pointer BoxTree<dim>::locate(const Vector& point, unsigned level) const
    {
        if(level > depth())
        {
            throw std::out_of_range("BoxTree level is not built");
        }
        if(!point.allFinite() || (point.array() < lowerValue.array()).any() ||
            (point.array() > upperValue.array()).any())
        {
            return nullptr;
        }
        const unsigned count = 1u << level;
        const auto width = latticeSize / count;
        Key key{};
        for(int axis = 0; axis < dim; ++axis)
        {
            const auto index = std::min(count - 1, static_cast<unsigned>(std::floor(
                (point[axis] - lowerValue[axis]) / (upperValue[axis] - lowerValue[axis]) * count)));
            key[axis] = index * width;
            key[dim + axis] = (index + 1) * width;
        }
        return state.factory.at(key);
    }

    template<int dim> std::vector<typename BoxTree<dim>::Pointer> BoxTree<dim>::edges(unsigned level) const
    {
        std::map<Key, Pointer> unique;
        const std::function<void(const Pointer&)> collect = [&](const Pointer& entity)
        {
            if(entity->order() == 1)
            {
                unique.emplace(entity->key(), entity);
            }
            else
            {
                for(const auto& incidence : entity->boundary())
                {
                    collect(incidence.entity);
                }
            }
        };
        for(const auto& box : boxes(level))
        {
            collect(box);
        }
        std::vector<Pointer> result;
        result.reserve(unique.size());
        for(const auto& entry : unique)
        {
            result.push_back(entry.second);
        }
        return result;
    }

    template class BoxEntity<1>;
    template class BoxEntity<2>;
    template class BoxEntity<3>;
    template class BoxTree<1>;
    template class BoxTree<2>;
    template class BoxTree<3>;
}
