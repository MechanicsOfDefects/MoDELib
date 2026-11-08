#ifndef model_GrainClimbPartition_H_
#define model_GrainClimbPartition_H_
#include <DislocationDynamics/ScreenedDislocationAssembly.h>
#include <map>
#include <set>
#include <stdexcept>
namespace model
{
    struct GrainClimbPartition
    {
        std::size_t grain=0;
        std::vector<std::size_t> nodes; // Local to network indexing.
        std::vector<bool> active;
        ClimbAssemblySnapshot snapshot;
        Eigen::VectorXd gather(const Eigen::VectorXd& value) const
        {
            Eigen::VectorXd result(nodes.size());
            for(std::size_t i=0;i<nodes.size();++i) { result[i]=value[nodes[i]]; }
            return result;
        }
    };
    inline std::vector<GrainClimbPartition> partitionClimbByGrain(const ClimbAssemblySnapshot& input)
    {
        if(input.speciesCount==0 || input.rhs.size()!=input.speciesCount)
        {
            throw std::invalid_argument("Invalid grain climb species/RHS dimensions");
        }
        for(const auto& rhs:input.rhs)
        {
            if(rhs.size()!=Eigen::Index(input.nodeCount) || !rhs.allFinite())
            {
                throw std::invalid_argument("Invalid grain climb RHS");
            }
        }
        std::map<std::size_t,GrainClimbPartition> groups;
        std::map<std::size_t,std::size_t> activeOwners;
        const auto observe=[&](const auto& entry,std::size_t grain)
        {
            for(int k=0;k<2;++k)
            {
                if(entry.nodes[k]>=input.nodeCount) { throw std::invalid_argument("Invalid grain climb node"); }
                if(entry.weights[k]!=0.0)
                {
                    const auto owner=activeOwners.emplace(entry.nodes[k],grain);
                    if(!owner.second && owner.first->second!=grain)
                    {
                        throw std::runtime_error("Independent grain climb requires no shared nonzero climb DOFs across grains");
                    }
                }
            }
        };
        for(const auto& target:input.targets)
        {
            observe(target,target.grain);
            groups[target.grain].snapshot.targets.push_back(target);
        }
        for(const auto& channel:input.channels)
        {
            for(const auto& source:channel.sources) { observe(source,channel.grain); }
            for(const auto& source:channel.segments) { observe(source,channel.grain); }
            groups[channel.grain].snapshot.channels.push_back(channel);
        }
        std::vector<GrainClimbPartition> result;
        for(auto& item:groups)
        {
            auto& group=item.second;
            if(group.snapshot.targets.empty()) { continue; }
            group.grain=item.first;
            std::set<std::size_t> ids;
            const auto add=[&](const auto& entry) { ids.insert(entry.nodes.begin(),entry.nodes.end()); };
            for(const auto& target:group.snapshot.targets) { add(target); }
            for(const auto& channel:group.snapshot.channels)
            {
                for(const auto& source:channel.sources) { add(source); }
                for(const auto& source:channel.segments) { add(source); }
            }
            group.nodes.assign(ids.begin(),ids.end());
            for(const auto node:group.nodes)
            {
                const auto owner=activeOwners.find(node);
                group.active.push_back(owner!=activeOwners.end() && owner->second==group.grain);
            }
            std::map<std::size_t,std::size_t> local;
            for(std::size_t i=0;i<group.nodes.size();++i) { local.emplace(group.nodes[i],i); }
            const auto remap=[&](auto& entry)
            {
                for(auto& node:entry.nodes) { node=local.at(node); }
            };
            for(auto& target:group.snapshot.targets) { remap(target); }
            for(auto& channel:group.snapshot.channels)
            {
                for(auto& source:channel.sources) { remap(source); }
                for(auto& source:channel.segments) { remap(source); }
            }
            group.snapshot.nodeCount=group.nodes.size();
            group.snapshot.speciesCount=input.speciesCount;
            for(const auto& rhs:input.rhs) { group.snapshot.rhs.push_back(group.gather(rhs)); }
            result.push_back(std::move(group));
        }
        return result;
    }
    // Restrict constraints, preserving any periodic/subdivision relation wholly
    // inside a grain. A shared reduced variable would couple the grain solves.
    inline std::vector<Eigen::SparseMatrix<double>> grainClimbConstraints(
        const Eigen::SparseMatrix<double>& global,const std::vector<GrainClimbPartition>& groups)
    {
        std::map<Eigen::Index,std::size_t> owners;
        std::vector<Eigen::SparseMatrix<double>> result;
        for(std::size_t g=0;g<groups.size();++g)
        {
            std::map<Eigen::Index,Eigen::Index> rows,cols;
            for(std::size_t i=0;i<groups[g].nodes.size();++i)
            {
                if(groups[g].active[i]) { rows.emplace(groups[g].nodes[i],i); }
            }
            std::vector<Eigen::Triplet<double>> entries;
            for(Eigen::Index c=0;c<global.outerSize();++c)
            {
                for(Eigen::SparseMatrix<double>::InnerIterator it(global,c);it;++it)
                {
                    const auto row=rows.find(it.row());
                    if(row==rows.end() || it.value()==0.0) { continue; }
                    const auto owner=owners.emplace(c,g);
                    if(!owner.second && owner.first->second!=g)
                    {
                        throw std::runtime_error("Climb constraint couples different grains; independent solves are invalid");
                    }
                    const auto column=cols.emplace(c,cols.size());
                    entries.emplace_back(row->second,column.first->second,it.value());
                }
            }
            Eigen::SparseMatrix<double> map(groups[g].nodes.size(),cols.size());
            map.setFromTriplets(entries.begin(),entries.end());
            result.push_back(std::move(map));
        }
        return result;
    }
}
#endif
