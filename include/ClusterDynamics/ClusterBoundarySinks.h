#ifndef model_ClusterBoundarySinks_H_
#define model_ClusterBoundarySinks_H_
#include <map>
#include <set>
#include <vector>
#include <Eigen/Core>
namespace model
{
    // Normals of physical sink faces incident on each higher-order FEM node.
    // Same-grain active periodic seams are not sinks; periodic GBs are.
    template<class FE>
    auto clusterBoundarySinks(const FE& fe)
    {
        std::map<size_t,std::vector<Eigen::Vector3d>> result;
        std::set<size_t> visited;
        using Element=typename FE::ElementType;
        for(const auto& region : fe.mesh.regions())
        {
            for(const auto& entry : region.second->faces())
            {
                const auto& face=*entry.second;
                if(!visited.insert(face.sID).second)
                {
                    continue;
                }
                const auto* partner=face.periodicFacePair.second;
                if(face.isExternal() && partner && face.regionIDs.first==partner->regionIDs.first)
                {
                    continue;
                }
                std::set<size_t> nodes;
                for(const auto* triangle : face.simplices())
                {
                    const auto& element=fe.elements().at(triangle->parents().begin()->first);
                    const size_t opposite=element.simplex.childOrder(triangle->xID);
                    for(size_t n=0;n<Element::nodesPerElement;++n)
                    {
                        if(Element::baryStarsAndBars(n,opposite)==0)
                        {
                            nodes.insert(element.node(n).gID);
                        }
                    }
                }
                for(const auto node : nodes)
                {
                    result[node].push_back(face.outNormal(face.regionIDs.first));
                }
            }
        }
        return result;
    }
}
#endif
