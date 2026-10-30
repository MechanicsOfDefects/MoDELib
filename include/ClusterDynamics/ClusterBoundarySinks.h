#ifndef model_ClusterBoundarySinks_H_
#define model_ClusterBoundarySinks_H_
#include <map>
#include <algorithm>
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
        std::map<size_t,std::map<size_t,Eigen::Vector3d>> classFaces;
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
                    const size_t faceID=partner?std::min(face.sID,partner->sID):face.sID;
                    classFaces[fe.periodicRepresentative(node)].emplace(faceID,face.outNormal(face.regionIDs.first));
                }
            }
        }
        // A periodic point sees the union of incident physical faces at all
        // copies. Deduplicate paired faces by identity, never by normal alone.
        for(const auto& node : fe.nodes())
        {
            const auto found=classFaces.find(fe.periodicRepresentative(node.gID));
            if(found!=classFaces.end())
            {
                for(const auto& face : found->second)
                {
                    result[node.gID].push_back(face.second);
                }
            }
        }
        return result;
    }
}
#endif
