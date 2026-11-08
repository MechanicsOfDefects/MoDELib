#include <DefectiveCrystal/SolvedQuantityRegistry.h>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <sstream>
#include <iomanip>
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
        if(andersonDepth>20 || maxIterations==0 || !std::isfinite(relativeTolerance) || relativeTolerance<0.0
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
        // Accelerate the complete (possibly damped) sequential sweep, not each physics.
        Eigen::Index stateSize=0;
        for(const auto& state : saved) stateSize+=state.size();
        const auto gather=[&]() -> Eigen::VectorXd
        {
            Eigen::VectorXd result(stateSize);
            Eigen::Index offset=0;
            for(size_t i=0;i<contributions.size();++i)
            {
                const Eigen::VectorXd state=contributions[i]->state();
                if(state.size()!=saved[i].size()) throw std::runtime_error("Coupled state size changed during solve.");
                result.segment(offset,state.size())=state;
                offset+=state.size();
            }
            return result;
        };
        std::vector<Eigen::VectorXd> residualHistory,imageHistory;
        Eigen::VectorXd scales=Eigen::VectorXd::Ones(stateSize);
        double previousResidual=0.0;
        bool acceleratedLast=false;
        std::string worst;
        Eigen::Index worstComponent=-1;
        double worstDelta=0.0;
        double worstTolerance=0.0;
        try
        {
            for(size_t iteration=0;iteration<options.maxIterations;++iteration)
            {
                const Eigen::VectorXd sweepBefore=gather();
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
                    Eigen::Index controlling=-1,maxAbsoluteComponent=-1;
                    double contributionError=0.0,maxAbsoluteUpdate=0.0,controllingTolerance=0.0;
                    for(Eigen::Index k=0;k<after.size();++k)
                    {
                        const double delta=std::abs(after(k)-before(k));
                        const double tolerance=absolute(k)+options.relativeTolerance*std::max(std::abs(after(k)),std::abs(before(k)));
                        const double e=tolerance>0.0? delta/tolerance : (delta==0.0? 0.0 : INFINITY);
                        if(controlling<0 || e>contributionError)
                        {
                            controlling=k;contributionError=e;controllingTolerance=tolerance;
                        }
                        if(maxAbsoluteComponent<0 || delta>maxAbsoluteUpdate)
                        {
                            maxAbsoluteComponent=k;maxAbsoluteUpdate=delta;
                        }
                        if(worstComponent<0 || e>error)
                        {
                            error=e;
                            worst=c->provider+":"+c->name;
                            worstComponent=k;
                            worstDelta=delta;
                            worstTolerance=tolerance;
                        }
                    }
                    const Eigen::VectorXd update=after-before;
                    const double updateNorm=update.stableNorm();
                    const double beforeNorm=before.stableNorm(),afterNorm=after.stableNorm();
                    const double fieldNorm=std::max(beforeNorm,afterNorm);
                    const double relativeUpdate=fieldNorm>0.0?updateNorm/fieldNorm:(updateNorm==0.0?0.0:INFINITY);
                    std::ostringstream diagnostic;
                    diagnostic<<std::scientific<<std::setprecision(8)
                              <<"  Coupling diagnostic "<<c->provider<<":"<<c->name
                              <<" (unrelaxed): ||update||2="<<updateNorm
                              <<", ||before||2="<<beforeNorm<<", ||after||2="<<afterNorm
                              <<", relative L2 update="<<relativeUpdate
                              <<", max |update|="<<maxAbsoluteUpdate
                              <<" at component "<<maxAbsoluteComponent;
                    if(controlling>=0)
                    {
                        diagnostic<<", normalized max="<<contributionError
                                  <<" at component "<<controlling
                                  <<", before="<<before(controlling)<<", after="<<after(controlling)
                                  <<", signed update="<<update(controlling)
                                  <<", tolerance="<<controllingTolerance;
                    }
                    diagnostic<<", relaxation="<<options.relaxation;
                    std::cout<<diagnostic.str()<<std::endl;
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
                if(options.andersonDepth && stateSize && iteration+1<options.maxIterations)
                {
                    const Eigen::VectorXd image=gather();
                    const Eigen::VectorXd residual=image-sweepBefore;
                    // Fix one scale per contribution for this solve: units and large
                    // blocks must not dominate the least squares or change its metric.
                    if(iteration==0)
                    {
                        Eigen::Index offset=0;
                        for(size_t i=0;i<contributions.size();++i)
                        {
                            const Eigen::Index n=saved[i].size();
                            const double scale=std::max({sweepBefore.segment(offset,n).stableNorm(),
                                image.segment(offset,n).stableNorm(),
                                contributions[i]->absoluteTolerance().stableNorm()});
                            scales.segment(offset,n).setConstant(scale>0.0?scale:1.0);
                            offset+=n;
                        }
                    }
                    const Eigen::VectorXd scaled=residual.cwiseQuotient(scales);
                    const double norm=scaled.stableNorm();
                    const bool restart=acceleratedLast && norm>2.0*previousResidual;
                    acceleratedLast=false;
                    if(restart)
                    {
                        residualHistory.clear();imageHistory.clear();
                        std::cout<<"  Anderson acceleration restarted: sweep residual grew; using ordinary sweep."<<std::endl;
                    }
                    if(!restart && !residualHistory.empty())
                    {
                        const Eigen::Index m=residualHistory.size();
                        Eigen::MatrixXd differences(stateSize,m);
                        for(Eigen::Index j=0;j<m;++j) differences.col(j)=residualHistory[j]-scaled;
                        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(differences);
                        qr.setThreshold(1.e-10);
                        const Eigen::VectorXd alpha=qr.solve(-scaled);
                        const double last=1.0-alpha.sum();
                        Eigen::VectorXd candidate=last*image;
                        for(Eigen::Index j=0;j<m;++j) candidate+=alpha(j)*imageHistory[j];
                        // Bound extrapolation and reject degenerate/nonfinite fits.
                        if(qr.rank()>0 && alpha.allFinite() && candidate.allFinite()
                           && alpha.lpNorm<1>()+std::abs(last)<=100.0
                           && (scaled+differences*alpha).stableNorm()<norm)
                        {
                            Eigen::Index offset=0;
                            for(size_t i=0;i<contributions.size();++i)
                            {
                                contributions[i]->restore(candidate.segment(offset,saved[i].size()));
                                offset+=saved[i].size();
                            }
                            acceleratedLast=true;
                            std::cout<<"  Anderson acceleration applied after iteration "<<iterations
                                     <<" (history="<<m<<", rank="<<qr.rank()<<")."<<std::endl;
                        }
                        else
                        {
                            residualHistory.clear();imageHistory.clear();
                            std::cout<<"  Anderson acceleration skipped: unreliable extrapolation; history reset."<<std::endl;
                        }
                    }
                    residualHistory.push_back(scaled);imageHistory.push_back(image);
                    if(residualHistory.size()>options.andersonDepth)
                    {
                        residualHistory.erase(residualHistory.begin());imageHistory.erase(imageHistory.begin());
                    }
                    previousResidual=norm;
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
