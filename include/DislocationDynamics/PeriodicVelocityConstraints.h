#ifndef model_PeriodicVelocityConstraints_H_
#define model_PeriodicVelocityConstraints_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <limits>
#include <map>
#include <vector>
#include <stdexcept>

namespace model
{
    // Homogeneous scalar constraints. Each pivot expresses one velocity in
    // terms of the remaining velocities. Repeated constraints are redundant;
    // different endpoint pairs impose compatibility, rather than adding two
    // interpolants to the same row. The three Cartesian components share Z.
    class PeriodicVelocityConstraints
    {
        using Row=std::map<size_t,double>;
        struct Pivot
        {
            size_t node;
            Row row;
        };
        size_t count;
        std::vector<Pivot> pivots;
        std::map<size_t,size_t> pivotIndex;

    public:
        explicit PeriodicVelocityConstraints(size_t countIn) : count(countIn)
        {
        }

        void interpolate(size_t node,size_t a,size_t b,double wa,double wb)
        {
            if(node>=count || a>=count || b>=count || !std::isfinite(wa) || !std::isfinite(wb))
            {
                throw std::runtime_error("Invalid periodic velocity interpolation.");
            }
            Row row;
            row[node]+=1.0;
            row[a]-=wa;
            row[b]-=wb;
            const double tolerance=128*std::numeric_limits<double>::epsilon();
            // Substitute earlier pivots in their creation order. Their rows can
            // reference later pivots, but never earlier ones.
            std::map<size_t,size_t> pending;
            for(const auto& value : row)
            {
                const auto found=pivotIndex.find(value.first);
                if(found!=pivotIndex.end())
                {
                    pending.emplace(found->second,value.first);
                }
            }
            while(!pending.empty())
            {
                const auto item=*pending.begin();
                pending.erase(pending.begin());
                const double factor=row[item.second];
                row.erase(item.second);
                for(const auto& value : pivots[item.first].row)
                {
                    row[value.first]+=factor*value.second;
                    const auto found=pivotIndex.find(value.first);
                    if(found!=pivotIndex.end())
                    {
                        pending.emplace(found->second,value.first);
                    }
                }
            }
            for(auto it=row.begin();it!=row.end();)
            {
                if(std::abs(it->second)<=tolerance)
                {
                    it=row.erase(it);
                }
                else
                {
                    ++it;
                }
            }
            if(row.empty())
            {
                return;
            }
            auto pivot=std::max_element(row.begin(),row.end(),[](const auto& x,const auto& y)
            {
                return std::abs(x.second)<std::abs(y.second);
            });
            // Every interpolation preserves rigid translation. Cancellation
            // between almost identical relations must not create an artificial
            // absolute-velocity constraint. Restore the zero coefficient sum
            // before normalizing a small residual row.
            long double otherSum=0.0L;
            for(const auto& value : row)
            {
                if(value.first!=pivot->first)
                {
                    otherSum+=value.second;
                }
            }
            pivot->second=-static_cast<double>(otherSum);
            if(std::abs(pivot->second)<=tolerance)
            {
                row.erase(pivot);
                if(row.empty())
                {
                    return;
                }
            }
            pivot=std::max_element(row.begin(),row.end(),[](const auto& x,const auto& y)
            {
                return std::abs(x.second)<std::abs(y.second);
            });
            const auto preferred=row.find(node);
            if(preferred!=row.end() && std::abs(preferred->second)>=0.5*std::abs(pivot->second))
            {
                pivot=preferred;
            }
            const size_t id=pivot->first;
            const double coefficient=pivot->second;
            row.erase(pivot);
            for(auto& value : row)
            {
                value.second/=-coefficient;
            }
            pivotIndex.emplace(id,pivots.size());
            pivots.push_back({id,std::move(row)});
        }

        template<int dim,class Triplets>
        size_t append(Triplets& triplets) const
        {
            std::vector<Row> expressions(count);
            size_t freeCount=0;
            for(size_t i=0;i<count;++i)
            {
                if(!pivotIndex.count(i))
                {
                    expressions[i].emplace(freeCount++,1.0);
                }
            }
            for(auto it=pivots.rbegin();it!=pivots.rend();++it)
            {
                for(const auto& dependency : it->row)
                {
                    for(const auto& value : expressions[dependency.first])
                    {
                        expressions[it->node][value.first]+=dependency.second*value.second;
                    }
                }
            }
            for(size_t i=0;i<count;++i)
            {
                for(const auto& value : expressions[i])
                {
                    for(int d=0;d<dim;++d)
                    {
                        triplets.emplace_back(dim*i+d,dim*value.first+d,value.second);
                    }
                }
            }
            return freeCount;
        }
    };

    template<int dim,class Network,class Triplets>
    size_t periodicVelocityConstraints(const Network& network,Triplets& triplets)
    {
        PeriodicVelocityConstraints constraints(network.networkNodes().size());
        // networkID() walks the ordered container. Cache the identical ordering
        // once rather than repeating that linear walk for every occurrence.
        std::map<size_t,size_t> indices;
        for(const auto& entry : network.networkNodes())
        {
            indices.emplace(entry.second.lock()->sID,indices.size());
        }
        for(const auto& entry : network.networkNodes())
        {
            const auto node=entry.second.lock();
            const bool subdivision=!node->loopNodes().empty() &&
                std::all_of(node->loopNodes().begin(),node->loopNodes().end(),[](const auto* ln)
                {
                    return ln->isPeriodicSubdivision();
                });
            if(!subdivision)
            {
                continue;
            }
            for(const auto* ln : node->loopNodes())
            {
                const auto a=ln->periodicPrev();
                const auto b=ln->periodicNext();
                if(a && b)
                {
                    const double da=(ln->get_P()-a->get_P()).norm();
                    const double db=(ln->get_P()-b->get_P()).norm();
                    if(da+db<=std::numeric_limits<double>::epsilon())
                    {
                        throw std::runtime_error("Degenerate periodic velocity interpolation endpoints.");
                    }
                    constraints.interpolate(indices.at(node->sID),indices.at(a->networkNode->sID),
                                            indices.at(b->networkNode->sID),db/(da+db),da/(da+db));
                }
            }
        }
        return constraints.template append<dim>(triplets);
    }
}
#endif
