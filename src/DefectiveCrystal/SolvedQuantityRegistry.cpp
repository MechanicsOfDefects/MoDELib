#include <DefectiveCrystal/SolvedQuantityRegistry.h>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <set>
#include <stdexcept>

namespace model
{
    namespace
    {
        void updateContribution(const std::shared_ptr<SolvedQuantityContribution>& contribution,
                                const SolvedQuantityContext& context)
        {
            const auto start=std::chrono::steady_clock::now();
            std::cout<<"Solving "<<contribution->provider<<":"<<contribution->name<<" "<<std::flush;
            contribution->update(context);
            std::cout<<" ["<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" sec]"<<std::endl;
        }
    }

    void SolvedQuantityOptions::validate() const
    {
        if(maxIterations==0 || !std::isfinite(relativeTolerance) || relativeTolerance<0.0
           || !std::isfinite(relaxation) || relaxation<=0.0 || relaxation>1.0)
        {
            throw std::runtime_error("Invalid SolvedQuantity convergence options.");
        }
    }

    void SolvedQuantityGroup::solve(const SolvedQuantityOptions& options,double time)
    {
        iterations=0;
        error=0.0;
        if(!cyclic)
        {
            updateContribution(contributions.front(),{time,0});
            iterations=1;
            return;
        }
        std::vector<Eigen::VectorXd> saved;
        for(const auto& c : contributions)
        {
            if(!c->state || !c->restore || !c->absoluteTolerance)
            {
                throw std::runtime_error("Cyclic contribution lacks convergence/restore callbacks: "+c->name);
            }
            saved.push_back(c->state());
        }
        std::string worst;
        Eigen::Index worstComponent=-1;
        double worstDelta=0.0;
        double worstTolerance=0.0;
        try
        {
            for(size_t iteration=0;iteration<options.maxIterations;++iteration)
            {
                error=0.0;
                worst.clear();
                worstComponent=-1;
                worstDelta=0.0;
                worstTolerance=0.0;
                for(const auto& c : contributions)
                {
                    const Eigen::VectorXd before=c->state();
                    updateContribution(c,{time,iteration});
                    const Eigen::VectorXd after=c->state();
                    const Eigen::VectorXd absolute=c->absoluteTolerance();
                    if(before.size()!=after.size() || absolute.size()!=after.size()
                       || !before.allFinite() || !after.allFinite() || !absolute.allFinite()
                       || (absolute.array()<0.0).any())
                    {
                        throw std::runtime_error("Invalid convergence state for "+c->name);
                    }
                    for(Eigen::Index k=0;k<after.size();++k)
                    {
                        const double delta=std::abs(after(k)-before(k));
                        const double tolerance=absolute(k)+options.relativeTolerance*std::max(std::abs(after(k)),std::abs(before(k)));
                        const double e=tolerance>0.0? delta/tolerance : (delta==0.0? 0.0 : INFINITY);
                        if(worstComponent<0 || e>error)
                        {
                            error=e;
                            worst=c->provider+":"+c->name;
                            worstComponent=k;
                            worstDelta=delta;
                            worstTolerance=tolerance;
                        }
                    }
                    // Test the unrelaxed discrepancy; damping must not fake convergence.
                    if(options.relaxation!=1.0)
                    {
                        c->restore(before+options.relaxation*(after-before));
                    }
                }
                iterations=iteration+1;
                std::cout<<"Coupled quantities: iteration "<<iterations<<", normalized error "<<error;
                if(worstComponent>=0)
                {
                    std::cout<<", max at "<<worst<<", component "<<worstComponent
                             <<" (zero-based state index), |update|="<<worstDelta
                             <<", tolerance="<<worstTolerance;
                }
                std::cout<<std::endl;
                if(error<=1.0)
                {
                    return;
                }
            }
            throw std::runtime_error("SolvedQuantityGroup did not converge after "+std::to_string(iterations)
                                     +" iterations; worst contribution "+worst+", component "+std::to_string(worstComponent)
                                     +", normalized error "+std::to_string(error));
        }
        catch(...)
        {
            for(size_t k=0;k<contributions.size();++k)
            {
                contributions[k]->restore(saved[k]);
            }
            throw;
        }
    }

    void SolvedQuantityRegistry::clear()
    {
        operations.clear();
        quantities.clear();
        groups.clear();
        finalized=false;
    }

    void SolvedQuantityRegistry::add(SolvedQuantityContribution c)
    {
        const std::string key=c.provider+":"+c.name;
        if(!c.update || c.quantities.empty() || operations.count(key))
        {
            throw std::runtime_error("Invalid or duplicate solved quantity operation: "+key);
        }
        auto ptr=std::make_shared<SolvedQuantityContribution>(std::move(c));
        operations.emplace(key,ptr);
        for(const auto& name : ptr->quantities)
        {
            auto& quantity=quantities[name];
            quantity.name=name;
            quantity.contributions.push_back(ptr);
        }
        finalized=false;
    }

    void SolvedQuantityRegistry::finalize()
    {
        options.validate();
        // Validate declarations even when a dependency is currently inactive.
        for(const auto& entry : operations)
        {
            for(const auto& d : entry.second->dependencies)
            {
                bool found=false;
                const auto q=quantities.find(d.quantity);
                if(q!=quantities.end())
                {
                    for(const auto& p : q->second.contributions)
                    {
                        found=found || ((d.provider.empty()? p!=entry.second : p->provider==d.provider));
                    }
                }
                if(d.required && !found)
                {
                    throw std::runtime_error("Missing provider for "+d.quantity+" required by "+entry.first);
                }
            }
        }
        finalized=true;
    }

    void SolvedQuantityRegistry::solve(double time)
    {
        if(!finalized)
        {
            throw std::runtime_error("Call finalizeSolvedQuantities before solve.");
        }
        options.validate();
        std::vector<std::shared_ptr<SolvedQuantityContribution>> nodes;
        std::map<const SolvedQuantityContribution*,size_t> ids;
        for(const auto& entry : operations)
        {
            ids.emplace(entry.second.get(),nodes.size());
            nodes.push_back(entry.second);
        }
        std::vector<std::set<size_t>> edges(nodes.size());
        for(size_t i=0;i<nodes.size();++i)
        {
            for(const auto& d : nodes[i]->dependencies)
            {
                if(d.enabled && !d.enabled())
                {
                    continue;
                }
                const auto q=quantities.find(d.quantity);
                if(q!=quantities.end())
                {
                    for(const auto& p : q->second.contributions)
                    {
                        if(d.provider.empty()? p!=nodes[i] : p->provider==d.provider)
                        {
                            edges[i].insert(ids.at(p.get()));
                        }
                    }
                }
            }
        }
        // Tarjan on consumer->provider edges emits dependency groups first.
        // Rebuild at each solve: DD can switch between glide and climb steps.
        const auto previousGroups=std::move(groups);
        groups.clear();
        std::vector<int> index(nodes.size(),-1),low(nodes.size(),-1);
        std::vector<size_t> stack;
        std::vector<bool> onStack(nodes.size(),false);
        int next=0;
        std::function<void(size_t)> visit=[&](size_t v)
        {
            index[v]=low[v]=next++;
            stack.push_back(v);
            onStack[v]=true;
            for(const size_t w : edges[v])
            {
                if(index[w]<0)
                {
                    visit(w);
                    low[v]=std::min(low[v],low[w]);
                }
                else if(onStack[w])
                {
                    low[v]=std::min(low[v],index[w]);
                }
            }
            if(low[v]==index[v])
            {
                SolvedQuantityGroup group;
                size_t w;
                do
                {
                    w=stack.back();
                    stack.pop_back();
                    onStack[w]=false;
                    group.contributions.push_back(nodes[w]);
                }
                while(w!=v);
                group.cyclic=group.contributions.size()>1 || edges[v].count(v);
                std::sort(group.contributions.begin(),group.contributions.end(),[](const auto& a,const auto& b)
                {
                    return a->provider+":"+a->name < b->provider+":"+b->name;
                });
                groups.push_back(std::move(group));
            }
        };
        for(size_t v=0;v<nodes.size();++v)
        {
            if(index[v]<0)
            {
                visit(v);
            }
        }
        const bool sameSequence=previousGroups.size()==groups.size()
            && std::equal(groups.begin(),groups.end(),previousGroups.begin(),[](const auto& a,const auto& b)
            {
                return a.cyclic==b.cyclic && a.contributions==b.contributions;
            });
        if(!sameSequence)
        {
            std::cout<<"Solve sequence ("<<groups.size()<<" groups):"<<std::endl;
            for(size_t k=0;k<groups.size();++k)
            {
                const auto& group=groups[k];
                std::cout<<"  Group "<<k+1<<" ["<<(group.cyclic? "iterate to convergence" : "solve once")<<"]"<<std::endl;
                for(const auto& contribution : group.contributions)
                {
                    std::cout<<"    "<<contribution->provider<<":"<<contribution->name<<std::endl;
                }
            }
        }
        for(size_t k=0;k<groups.size();++k)
        {
            std::cout<<"Solving group "<<k+1<<"/"<<groups.size()<<std::endl;
            groups[k].solve(options,time);
        }
    }
}
