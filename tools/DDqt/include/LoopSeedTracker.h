#ifndef model_LoopSeedTracker_H_
#define model_LoopSeedTracker_H_
#include <Eigen/Dense>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
namespace model
{
    // Positions are unwrapped loop coordinates; shift maps each outgoing edge to the RVE.
    struct LoopSeedNode
    {
        size_t id;
        Eigen::Vector3d position,shift;
        bool seedEdge=true;
    };
    class LoopSeedTracker
    {
        struct State
        {
            Eigen::Vector3d origin;
            size_t count;
        };
        std::map<size_t,std::map<size_t,State>> frames;
        double spacing=-1;
    public:
        std::vector<std::pair<size_t,Eigen::Vector3d>> seeds(size_t frame,size_t loop,
            const std::vector<LoopSeedNode>& nodes,const Eigen::Vector3d& normal,
            double requestedSpacing,double offset)
        {
            if(spacing!=requestedSpacing)
            {
                frames.clear();
                spacing=requestedSpacing;
            }
            if(nodes.empty())
            {
                return {};
            }
            size_t origin=0;
            for(size_t i=1;i<nodes.size();++i)
            {
                if(nodes[i].id<nodes[origin].id)
                {
                    origin=i;
                }
            }
            const State* previous=nullptr;
            auto frameIt=frames.find(frame);
            if(frameIt!=frames.end() && frameIt->second.count(loop))
            {
                previous=&frameIt->second.at(loop);
            }
            else
            {
                // Use the closest previously visited earlier frame. If seeking
                // backwards before all visited frames, use the earliest known one.
                auto it=frames.lower_bound(frame);
                while(it!=frames.begin())
                {
                    --it;
                    if(it->second.count(loop))
                    {
                        previous=&it->second.at(loop);
                        break;
                    }
                }
                if(!previous)
                {
                    for(const auto& item:frames)
                    {
                        if(item.second.count(loop))
                        {
                            previous=&item.second.at(loop);
                            break;
                        }
                    }
                }
            }
            if(previous)
            {
                double best=(nodes[origin].position-previous->origin).squaredNorm();
                for(size_t i=0;i<nodes.size();++i)
                {
                    const double distance=(nodes[i].position-previous->origin).squaredNorm();
                    if(distance<best || (distance==best && nodes[i].id<nodes[origin].id))
                    {
                        best=distance;
                        origin=i;
                    }
                }
            }
            double length=0;
            std::vector<double> lengths(nodes.size());
            for(size_t i=0;i<nodes.size();++i)
            {
                lengths[i]=(nodes[(i+1)%nodes.size()].position-nodes[i].position).norm();
                length+=lengths[i];
            }
            if(length<=0 || normal.norm()==0)
            {
                return {};
            }
            const size_t count=previous?previous->count:static_cast<size_t>(std::min(1.e6,std::max(1.0,std::ceil(length/spacing))));
            frames[frame][loop]={nodes[origin].position,count};
            std::vector<std::pair<size_t,Eigen::Vector3d>> result;
            size_t edge=0;
            double start=0;
            const Eigen::Vector3d u=normal.normalized();
            for(size_t seed=0;seed<count;++seed)
            {
                const double arc=(seed+0.5)*length/count;
                while(edge+1<nodes.size() && start+lengths[(origin+edge)%nodes.size()]<=arc)
                {
                    start+=lengths[(origin+edge)%nodes.size()];
                    ++edge;
                }
                const size_t i=(origin+edge)%nodes.size();
                if(!nodes[i].seedEdge || lengths[i]<=0)
                {
                    continue;
                }
                const Eigen::Vector3d tangent=(nodes[(i+1)%nodes.size()].position-nodes[i].position)/lengths[i];
                const Eigen::Vector3d v=tangent.cross(u).normalized();
                const Eigen::Vector3d center=nodes[i].position+(arc-start)*tangent+nodes[i].shift;
                result.emplace_back(4*seed,center+offset*u);
                result.emplace_back(4*seed+1,center+offset*v);
                result.emplace_back(4*seed+2,center-offset*u);
                result.emplace_back(4*seed+3,center-offset*v);
            }
            return result;
        }
    };
}
#endif
