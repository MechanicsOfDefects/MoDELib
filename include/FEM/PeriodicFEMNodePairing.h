#ifndef model_PeriodicFEMNodePairing_H_
#define model_PeriodicFEMNodePairing_H_
#include <Eigen/Core>
#include <array>
#include <map>
#include <vector>
#include <limits>
#include <stdexcept>
#include <algorithm>
namespace model
{
    // Geometric relations only: no DOFs are eliminated and no field is tied.
    template<class Element>
    struct PeriodicFEMNodePairing
    {
        static constexpr int dim=Element::dim;
        using Vector=Eigen::Matrix<double,dim,1>;
        struct Pair
        {
            size_t source,sink;
            Vector shift;
        };
        std::vector<Pair> pairs; // unique undirected edges, source < sink
        std::vector<size_t> representatives;
        std::vector<Vector> shifts; // node position = representative + shift

        template<class FE>
        explicit PeriodicFEMNodePairing(const FE& fe)
        {
            const double tolerance=1.e-8*(fe.mesh.xMax()-fe.mesh.xMin()).norm()+
                64*std::numeric_limits<double>::epsilon()*std::max(fe.mesh.xMin().norm(),fe.mesh.xMax().norm());
            std::map<std::pair<size_t,size_t>,Vector> unique;
            for(const auto& entry : fe.mesh.periodicSimplexPairs())
            {
                const auto* source=entry.first;
                const auto& relation=entry.second;
                const auto* target=relation.partner;
                if(source->parents().size()!=1 || target->parents().size()!=1)
                {
                    throw std::runtime_error("Periodic FEM pairing requires external boundary simplices.");
                }
                const auto& a=fe.elements().at(source->parents().begin()->first);
                const auto& b=fe.elements().at(target->parents().begin()->first);
                const auto localIndices=[](const auto& element,const auto* face)
                {
                    std::array<size_t,dim> indices{};
                    for(size_t i=0;i<dim;++i)
                    {
                        bool found=false;
                        for(size_t j=0;j<dim+1;++j)
                        {
                            if(face->xID(i)==element.simplex.xID(j))
                            {
                                indices[i]=j;
                                found=true;
                                break;
                            }
                        }
                        if(!found)
                        {
                            throw std::runtime_error("Boundary vertex missing from parent element.");
                        }
                    }
                    return indices;
                };
                const auto ia=localIndices(a,source);
                const auto ib=localIndices(b,target);
                using Key=std::array<int,dim>;
                std::map<Key,size_t> targetNodes;
                for(size_t n=0;n<Element::nodesPerElement;++n)
                {
                    Key key{};
                    int sum=0;
                    for(size_t i=0;i<dim;++i)
                    {
                        key[i]=Element::baryStarsAndBars(n,ib[i]);
                        sum+=key[i];
                    }
                    if(sum==Element::order)
                    {
                        if(!targetNodes.emplace(key,b.node(n).gID).second)
                        {
                            throw std::runtime_error("Duplicate FEM boundary barycentric node.");
                        }
                    }
                }
                size_t matched=0;
                for(size_t n=0;n<Element::nodesPerElement;++n)
                {
                    Key key{};
                    int sum=0;
                    for(size_t i=0;i<dim;++i)
                    {
                        const int weight=Element::baryStarsAndBars(n,ia[i]);
                        key[relation.vertexPermutation[i]]=weight;
                        sum+=weight;
                    }
                    if(sum!=Element::order)
                    {
                        continue;
                    }
                    const auto found=targetNodes.find(key);
                    if(found==targetNodes.end())
                    {
                        throw std::runtime_error("Missing periodic higher-order FEM node.");
                    }
                    ++matched;
                    size_t u=a.node(n).gID,v=found->second;
                    Vector shift=relation.shift;
                    if((fe.node(u).P0+shift-fe.node(v).P0).norm()>tolerance || u==v)
                    {
                        throw std::runtime_error("Periodic FEM node translation mismatch.");
                    }
                    if(u>v)
                    {
                        std::swap(u,v);
                        shift=-shift;
                    }
                    const auto inserted=unique.emplace(std::make_pair(u,v),shift);
                    if(!inserted.second && (inserted.first->second-shift).norm()>tolerance)
                    {
                        throw std::runtime_error("Conflicting periodic FEM node translations.");
                    }
                }
                if(matched!=targetNodes.size())
                {
                    throw std::runtime_error("Incomplete periodic FEM face pairing.");
                }
            }
            std::vector<std::vector<std::pair<size_t,Vector>>> adjacency(fe.nodeSize());
            for(const auto& entry : unique)
            {
                const auto [u,v]=entry.first;
                pairs.push_back({u,v,entry.second});
                adjacency[u].emplace_back(v,entry.second);
                adjacency[v].emplace_back(u,-entry.second);
            }
            representatives.assign(fe.nodeSize(),fe.nodeSize());
            shifts.assign(fe.nodeSize(),Vector::Zero());
            for(size_t root=0;root<fe.nodeSize();++root)
            {
                if(representatives[root]!=fe.nodeSize())
                {
                    continue;
                }
                representatives[root]=root;
                std::vector<size_t> queue{root};
                for(size_t k=0;k<queue.size();++k)
                {
                    const size_t u=queue[k];
                    for(const auto& edge : adjacency[u])
                    {
                        const size_t v=edge.first;
                        const Vector offset=shifts[u]+edge.second;
                        if(representatives[v]==fe.nodeSize())
                        {
                            representatives[v]=root;
                            shifts[v]=offset;
                            queue.push_back(v);
                        }
                        else if(representatives[v]!=root || (shifts[v]-offset).norm()>tolerance)
                        {
                            throw std::runtime_error("Inconsistent periodic FEM edge/corner cycle.");
                        }
                    }
                }
            }
        }
        PeriodicFEMNodePairing()=default;
    };
}
#endif
