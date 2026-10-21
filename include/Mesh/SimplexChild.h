/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_SimplexChild_H_
#define model_SimplexChild_H_
#include <memory>
#include <map>
#include <set>
#include <SimplexTraits.h>
#include <BoundarySimplex.h>
namespace model
{
    template <short int dim, short int order> class SimplexChild
    {
      public:
        using SimplexType = typename SimplexTraits<dim, order>::SimplexType;
        using ParentSimplexType = typename SimplexTraits<dim, order + 1>::SimplexType;
        using ParentSimplexIDType = typename SimplexTraits<dim, order + 1>::SimplexIDType;
        using ParentContainerType = std::map<ParentSimplexIDType, std::shared_ptr<ParentSimplexType>>;
        using ParentLinks = std::map<ParentSimplexIDType, std::weak_ptr<ParentSimplexType>>;
        using SiblingsContainerType =
            std::map<typename SimplexTraits<dim, order>::SimplexIDType, std::shared_ptr<SimplexType>>;

      private:
        ParentLinks parentLinks_;

      public:
        void addParent(const std::shared_ptr<ParentSimplexType> &parent)
        {
            const auto found = parentLinks_.find(parent->xID);
            if (found != parentLinks_.end())
            {
                if (const auto existing = found->second.lock(); existing && existing != parent)
                {
                    throw std::runtime_error("Conflicting simplex parent identity");
                }
            }
            parentLinks_.insert_or_assign(parent->xID, parent);
        }
        const ParentLinks &parentLinks() const
        {
            return parentLinks_;
        }
        // A temporary locked view; expired entries never affect topology queries.
        ParentContainerType parents() const
        {
            ParentContainerType result;
            for (const auto &entry : parentLinks_)
            {
                if (const auto parent = entry.second.lock())
                {
                    result.emplace(entry.first, parent);
                }
            }
            return result;
        }
        SiblingsContainerType siblings() const
        {
            SiblingsContainerType result;
            for (const auto &parent : parents())
            {
                for (const auto &child : parent.second->children())
                {
                    result.emplace(child->xID, child);
                }
            }
            return result;
        }
        bool isBoundarySimplex() const
        {
            return BoundarySimplex<dim, dim - order>::isBoundarySimplex(*this);
        }
        bool isRegionBoundarySimplex() const
        {
            return BoundarySimplex<dim, dim - order>::isRegionBoundarySimplex(*this);
        }
        bool isInRegion(const size_t &id) const
        {
            return regionIDs().count(id) != 0;
        }
        std::set<size_t> regionIDs() const
        {
            std::set<size_t> result;
            for (const auto &parent : parents())
            {
                const auto ids = parent.second->regionIDs();
                result.insert(ids.begin(), ids.end());
            }
            return result;
        }
    };
} // namespace model
#endif
