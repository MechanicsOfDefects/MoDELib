/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_SimplexFactory_H_
#define model_SimplexFactory_H_
#include <WeakPtrFactories.h>
#include <SimplexTraits.h>
#include <tuple>
#include <stdexcept>
namespace model
{
    template <int dim> class SimplicialMesh;
    template <short int dim, short int order> class Simplex;

    template <short int dim, short int order>
    struct SimplexFactory : WeakPtrFactory<SimplexFactory<dim, order>, Simplex<dim, order>>
    {
        using Base = WeakPtrFactory<SimplexFactory<dim, order>, Simplex<dim, order>>;
        using KeyType = typename SimplexTraits<dim, order>::SimplexIDType;
        using Pointer = std::shared_ptr<Simplex<dim, order>>;
        SimplicialMesh<dim> *mesh = nullptr;

        SimplexFactory() = default;
        SimplexFactory(const SimplexFactory &) = delete;
        SimplexFactory &operator=(const SimplexFactory &) = delete;

      private:
        using Base::create;

      public:
        template <typename... Args> Pointer getOrCreate(const KeyType &input, const Args &...args)
        {
            const auto key = SimplexTraits<dim, order>::sortID(input);
            if (const auto existing = Base::get(key))
            {
                if constexpr (order == dim)
                {
                    static_assert(sizeof...(Args) == 1);
                    if (existing->region->regionID != size_t(std::get<0>(std::tie(args...))))
                    {
                        throw std::runtime_error("A volume simplex cannot belong to two regions");
                    }
                }
                return existing;
            }
            const auto simplex = Base::create(key, args...);
            if constexpr (order > 0)
            {
                // Construction is complete: weak parents can now reference its control block.
                for (const auto &child : simplex->children())
                {
                    child->addParent(simplex);
                }
            }
            return simplex;
        }
        auto liveSimplices() const
        {
            std::map<KeyType, Pointer> result;
            for (const auto &entry : *this)
            {
                if (const auto simplex = entry.second.lock())
                {
                    result.emplace(entry.first, simplex);
                }
            }
            return result;
        }
        void pruneExpired()
        {
            for (auto entry = this->begin(); entry != this->end();)
            {
                if (entry->second.expired())
                {
                    entry = this->erase(entry);
                }
                else
                {
                    ++entry;
                }
            }
        }
        auto faces(const KeyType &id)
            requires(order > 0)
        {
            typename SimplexTraits<dim, order>::BaseArrayType result;
            for (size_t k = 0; k < result.size(); ++k)
            {
                result[k] =
                    mesh->template simplexFactory<order - 1>().getOrCreate(SimplexTraits<dim, order>::faceID(id, k));
            }
            return result;
        }
    };
    template <int dim, typename Sequence> struct SimplexFactoryTuple;
    template <int dim, size_t... orders> struct SimplexFactoryTuple<dim, std::index_sequence<orders...>>
    {
        using Type = std::tuple<SimplexFactory<dim, short(orders)>...>;
    };
    template <int dim> class SimplexFactories
    {
        typename SimplexFactoryTuple<dim, std::make_index_sequence<dim + 1>>::Type factories;

      public:
        explicit SimplexFactories(SimplicialMesh<dim> *owner)
        {
            std::apply([owner](auto &...factory) { ((factory.mesh = owner), ...); }, factories);
        }
        template <short int order> auto &factory()
        {
            return std::get<order>(factories);
        }
        template <short int order> const auto &factory() const
        {
            return std::get<order>(factories);
        }
        void clear()
        {
            std::apply([](auto &...factory) { (factory.clear(), ...); }, factories);
        }
    };
} // namespace model
#endif
