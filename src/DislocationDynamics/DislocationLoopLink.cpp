/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopLink_cpp_
#define model_DislocationLoopLink_cpp_

#include <DislocationLoopLink.h>

namespace model
{
    
    template <int dim>
    DislocationLoopLink<dim>::DislocationLoopLink(const std::shared_ptr<LoopNodeType>& so,
                                             const std::shared_ptr<LoopNodeType>& si,
                                             const std::shared_ptr<LoopType>& pL) :
    /* init */ LoopLink<LoopLinkType>(so,si,pL)
    {
        
        
    }
    
    template <int dim>
    bool DislocationLoopLink<dim>::hasNetworkLink() const
    {
        // return !((this->source->get_P()-this->sink->get_P()).squaredNorm()<FLT_EPSILON
        // /*    */ && this->source->periodicPlaneEdge
        // /*    */ &&   this->sink->periodicPlaneEdge
        // /*    */ && this->source->loop()->sID==this->sink->loop()->sID);
        
        if(this->source->networkNode==this->sink->networkNode)
        {
            return false;
        }
        
        if (   this->source->periodicPlaneEdge.first
            && this->sink->periodicPlaneEdge.first
            && this->source->periodicPlanePatch()!=this->sink->periodicPlanePatch())
        {
            return false;
        }
        
//        if (this->source->periodicPlaneEdge.first && this->sink->periodicPlaneEdge.first)
//        {
//            // bool temp((this->source->loop()->sID==this->sink->loop()->sID));
//            // if (temp)
//            // {
//            // std::cout<<"Checking Twin"<<std::endl;
//            // std::cout<<this->source->periodicPlaneEdge->twin->tag()<<"=>= " <<this->sink->periodicPlaneEdge->tag()<<std::endl;
//            // std::cout<<this->sink->periodicPlaneEdge->twin->tag()<<"=>= " <<this->source->periodicPlaneEdge->tag()<<std::endl;
//            // std::cout<<" Case A"<<std::flush;
//            if (this->source->periodicPlaneEdge.first->twin && this->sink->periodicPlaneEdge.first->twin)
//            {
//                // std::cout<<".1"<<std::endl;
//                // temp *= ((this->source->periodicPlaneEdge->twin == this->sink->periodicPlaneEdge.get()) &&
//                //          (this->sink->periodicPlaneEdge->twin == this->source->periodicPlaneEdge.get()));
//                const bool firstEdgeMatch((this->source->periodicPlaneEdge.first->twin == this->sink->periodicPlaneEdge.first.get()) &&
//                                          (this->sink->periodicPlaneEdge.first->twin == this->source->periodicPlaneEdge.first.get()));
//                if (firstEdgeMatch && this->source->periodicPlaneEdge.second && this->sink->periodicPlaneEdge.second)
//                {
//                    if (this->source->periodicPlaneEdge.second->twin && this->sink->periodicPlaneEdge.second->twin)
//                    {
//                        const bool secondEdgeMatch((this->source->periodicPlaneEdge.second->twin == this->sink->periodicPlaneEdge.second.get()) &&
//                                                   (this->sink->periodicPlaneEdge.second->twin == this->source->periodicPlaneEdge.second.get()));
//                        assert(secondEdgeMatch && "Second edge must match provided first edge match for network link creation");
//                    }
//                }
//
//                return !firstEdgeMatch;
//            }
//        }
        return true;
    }

    template <int dim>
    std::shared_ptr<PeriodicPlanePatch<dim>> DislocationLoopLink<dim>::periodicPlanePatch() const
    {
        
        return (this->source->periodicPlanePatch()==this->sink->periodicPlanePatch())? this->source->periodicPlanePatch() : std::shared_ptr<PeriodicPlanePatch<dim>>(nullptr);
        
    }

template <int dim>
CatmullRomSplineSegment<dim> DislocationLoopLink<dim>::spline() const
{
    const auto controlPoint=[&](const bool before)->VectorDim
    {
        const LoopNodeType* node=before ? this->source.get() : this->sink.get();
        bool backwards=before;
        VectorDim offset=VectorDim::Zero();
        std::set<std::pair<const LoopNodeType*,bool>> visited;
        while(visited.emplace(node,backwards).second)
        {
            const auto* adjacent=backwards ? node->prev.second : node->next.second;
            if(const auto* twin=adjacent->twin())
            {
                // Preserve the existing same-loop periodic continuation.
                const auto* neighbor=backwards ? twin->source->periodicPrev() : twin->sink->periodicNext();
                if(neighbor)
                {
                    return offset+neighbor->get_P()-node->periodicPlanePatch()->shift+neighbor->periodicPlanePatch()->shift;
                }
                break;
            }
            if(const auto* twin=adjacent->twinnedLink())
            {
                // A two-link zero-Burgers connection is artificial. Switch at
                // this endpoint, then follow the other loop away from it.
                // Match network endpoints rather than assuming equal loop
                // orientations: opposite Burgers vectors can also cancel.
                const bool atSource=twin->source->networkNode==node->networkNode;
                const auto* other=atSource ? twin->source.get() : twin->sink.get();
                offset+=node->get_P()-other->get_P();
                node=other;
                backwards=atSource;
                continue;
            }
            const auto* neighbor=backwards ? node->periodicPrev() : node->periodicNext();
            if(neighbor)
            {
                return offset+neighbor->get_P();
            }
            break;
        }
        // A closed chain of cancelling links has no unique continuation.
        // Retain the loop-local control point rather than circulating forever.
        const auto* neighbor=before ? this->source->periodicPrev() : this->sink->periodicNext();
        if(!neighbor)
        {
            throw std::runtime_error("Cannot determine spline control point.");
        }
        return neighbor->get_P();
    };
    return CatmullRomSplineSegment<dim>(controlPoint(true),this->source->get_P(),this->sink->get_P(),controlPoint(false));
}


    template <int dim>
    void DislocationLoopLink<dim>::initFromFile(const std::string& fileName)
    {
        verboseDislocationLoopLink=TextFileParser(fileName).readScalar<int>("verboseDislocationLoopLink",false);
    }
    
    template <int dim>
    const typename DislocationLoopLink<dim>::LoopLinkType * DislocationLoopLink<dim>::twinnedLink() const
    {
        if (hasNetworkLink())
        {
            if (this->networkLink()->loopLinks().size()==2 && this->networkLink()->hasZeroBurgers())
            {
                if (*this->networkLink()->loopLinks().begin()==this)
                {
                    return (*this->networkLink()->loopLinks().rbegin());
                }
                else
                {
                    return (*this->networkLink()->loopLinks().begin());
                }
                
            }
            else
            {
                return nullptr;
            }
            
        }
        else
        {
            return nullptr;
        }
        
    }

    template <int dim>
    int DislocationLoopLink<dim>::verboseDislocationLoopLink=0;

    template class DislocationLoopLink<3>;
    
}
#endif
