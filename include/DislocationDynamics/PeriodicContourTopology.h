#ifndef model_PeriodicContourTopology_H_
#define model_PeriodicContourTopology_H_

#include <vector>
#include <cstddef>

namespace model
{
    // Inspect the committed unwrapped contour, before moving or contracting
    // endpoints. Coincident seam copies are one geometric vertex; both copies
    // must retain a genuine bend, including a reversal along the same line.
    template <class Nodes>
    void retainPeriodicContourBends(const Nodes& nodes,const double tolerance)
    {
        for(size_t i=0;i<nodes.size();++i)
        {
            if(!nodes[i]->isPeriodicSubdivision())
            {
                continue;
            }
            size_t previous=(i+nodes.size()-1)%nodes.size();
            size_t next=(i+1)%nodes.size();
            while(previous!=i && (nodes[previous]->get_P()-nodes[i]->get_P()).norm()<=tolerance)
            {
                previous=(previous+nodes.size()-1)%nodes.size();
            }
            while(next!=i && (nodes[next]->get_P()-nodes[i]->get_P()).norm()<=tolerance)
            {
                next=(next+1)%nodes.size();
            }
            const auto incoming=(nodes[i]->get_P()-nodes[previous]->get_P()).eval();
            const auto outgoing=(nodes[next]->get_P()-nodes[i]->get_P()).eval();
            if(previous!=i && next!=i &&
               (incoming.dot(outgoing)<=0 ||
                incoming.cross(outgoing).norm()>tolerance*(incoming.norm()+outgoing.norm())))
            {
                nodes[i]->periodicContourVertex=true;
            }
        }
    }

    // With two anchors, require both ordered paths to traverse the same
    // unwrapped segment monotonically in opposite directions. Area alone is
    // insufficient: a self-intersecting contour can have zero signed area.
    template <class Nodes>
    bool retracesContourSegment(const Nodes& nodes,const std::vector<size_t>& anchors,const double tolerance)
    {
        if(anchors.empty())
        {
            return false;
        }
        if(anchors.size()==1)
        {
            for(const auto& node : nodes)
            {
                if((node->get_P()-nodes[anchors.front()]->get_P()).norm()>tolerance)
                {
                    return false;
                }
            }
            return true;
        }
        if(anchors.size()!=2)
        {
            return false;
        }
        for(size_t path=0;path<2;++path)
        {
            const size_t begin=anchors[path];
            const size_t end=anchors[1-path];
            const auto origin=nodes[begin]->get_P();
            const auto chord=(nodes[end]->get_P()-origin).eval();
            const double length=chord.norm();
            double previous=0;
            for(size_t i=(begin+1)%nodes.size();i!=end;i=(i+1)%nodes.size())
            {
                const auto offset=(nodes[i]->get_P()-origin).eval();
                if(length<=tolerance)
                {
                    if(offset.norm()>tolerance)
                    {
                        return false;
                    }
                }
                else
                {
                    const double distance=offset.dot(chord)/length;
                    if(distance<previous-tolerance || distance>length+tolerance ||
                       (offset-distance*chord/length).norm()>tolerance)
                    {
                        return false;
                    }
                    previous=distance;
                }
            }
        }
        return true;
    }
}
#endif
