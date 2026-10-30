/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopNode_cpp_
#define model_DislocationLoopNode_cpp_

#include <DislocationLoopNode.h>

namespace model
{
    
    template <int dim>
    DislocationLoopNode<dim>::DislocationLoopNode(typename DislocationLoopNode<dim>::LoopNetworkType* const net,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim>::LoopType>& loop,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim>::NetworkNodeType>& networkNode,
                                                                           const VectorDim& P,
                                                                           const std::shared_ptr<PeriodicPlanePatch<dim>>& patch_in,
                                                                           const std::pair<const std::shared_ptr<PeriodicPlaneEdge<dim>>,const std::shared_ptr<PeriodicPlaneEdge<dim>>>& edge_in
//                                                                           const std::set<std::shared_ptr<PeriodicPlanePatch<dim>>>& auxPatch_in
                                                                           ) :
    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
    /* init */,SplineNodeType(P)
    /* init */,_periodicPlanePatch(patch_in)
    /* init */,periodicPlaneEdge(edge_in)
    {
        VerboseDislocationLoopNode(1,"Creating LoopNode "<<this->tag()<<std::endl;);
        if (periodicPlaneEdge.second!=nullptr)
        {
            if (periodicPlaneEdge.first == nullptr)
            {
                assert(false && "Inconsisten definition of periodic plane edges for the dislocation loop node");
            }
        }
        
    }

    template <int dim>
    DislocationLoopNode<dim>::DislocationLoopNode(typename DislocationLoopNode<dim>::LoopNetworkType* const net,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim>::LoopType>& loop,
                                                                           const std::shared_ptr<typename DislocationLoopNode<dim>::NetworkNodeType>& networkNode,
                                                                           const LoopLinkType* const loopLink) :
    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
    /* init */,SplineNodeType(loopLink->periodicPlanePatch()? networkNode->get_P()-loopLink->periodicPlanePatch()->shift : networkNode->get_P())
    /* init */,_periodicPlanePatch(loopLink->periodicPlanePatch())
    /* init */,periodicPlaneEdge(std::make_pair(nullptr,nullptr))
    {
        VerboseDislocationLoopNode(1,"Creating LoopNode without PeriodicPlaneEdge "<<this->tag()<<std::endl;);
        
    }
    
    template <int dim>
    bool DislocationLoopNode<dim>::isContractableTo(const LoopNodeType* const other) const
    {
        return LoopNode<LoopNodeType>::isContractableTo(other) && _periodicPlanePatch==other->periodicPlanePatch();
    }

    
    template <int dim>
    std::shared_ptr<DislocationLoopNode<dim>> DislocationLoopNode<dim>::clone(const std::shared_ptr<typename DislocationLoopNode<dim>::LoopType>& otherLoop,
                                                                                                                                const std::shared_ptr<typename DislocationLoopNode<dim>::NetworkNodeType>& otherNetworkNode) const
    {
        auto result=std::shared_ptr<DislocationLoopNode<dim>>(new DislocationLoopNode(this->p_network(),otherLoop,otherNetworkNode,this->get_P(),_periodicPlanePatch,periodicPlaneEdge));
        result->periodicContourVertex=periodicContourVertex;
        return result;
    }
    
    template <int dim>
    void DislocationLoopNode<dim>::initFromFile(const std::string& fileName)
    {
        verboseDislocationLoopNode=TextFileParser(fileName).readScalar<int>("verboseDislocationLoopNode",false);
    }

    template <int dim>
    const DislocationLoopNode<dim>* DislocationLoopNode<dim>::periodicPrev() const
    {
        if (this->prev.first)
        {
            auto currentPrev(this->prev.first);
            while (currentPrev->isPeriodicSubdivision()) //First is sufficient
            {
                if (currentPrev == this)
                {
                    return nullptr;
                }
                currentPrev = currentPrev->prev.first;
            }
            return currentPrev;
        }
        return nullptr;
    }
    
    template <int dim>
    const DislocationLoopNode<dim>* DislocationLoopNode<dim>::periodicNext() const
    {
        if (this->next.first)
        {
            auto currentNext(this->next.first);
            while (currentNext->isPeriodicSubdivision())
            {
                if (currentNext == this)
                {
                    return nullptr;
                }
                currentNext = currentNext->next.first;
            }
            return currentNext;
        }
        return nullptr;
    }
    
    template <int dim>
    std::vector<DislocationLoopNode<dim>*> DislocationLoopNode<dim>::boundaryPrev() const
    {
        assert(!this->periodicPlaneEdge.first);
        std::vector<DislocationLoopNode<dim>*> temp;
        if (this->prev.first)
        {
            auto currentPrev(this->prev.first);
            while (currentPrev->isPeriodicSubdivision())
            {
                temp.push_back(currentPrev);
                currentPrev = currentPrev->prev.first;
            }
        }
        return temp;
    }
    
    template <int dim>
    std::vector<DislocationLoopNode<dim>*> DislocationLoopNode<dim>::boundaryNext() const
    {
                assert(!this->periodicPlaneEdge.first);
        std::vector<DislocationLoopNode<dim>*> temp;
        if (this->next.first)
        {
            auto currentNext(this->next.first);
            while (currentNext->isPeriodicSubdivision())
            {
                temp.push_back(currentNext);
                currentNext = currentNext->next.first;
            }
        }

        return temp;
    }

    template <int dim>
    void DislocationLoopNode<dim>::set_P(const typename DislocationLoopNode<dim>::VectorLowerDim& newP)
    {
       VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->tag()<<" set_P (lowerDim)"<<std::endl;);
        assert(periodicPlaneEdge.first);

        VectorDim globalPosition(this->loop()->periodicGlidePlane->referencePlane->globalPosition(newP));
        VectorDim localPosition(VectorDim::Zero());

        if (periodicPlaneEdge.second==nullptr)
        {
            localPosition=periodicPlaneEdge.first->meshIntersection->snap(globalPosition + periodicPlaneEdge.first->patch->shift);
            // if (this->networkNode->meshFaces().size()==2)
            // {
            //     //This means that the node was created at a mesh face and not at the edge
            //     //Clear the confinement so that it can move away from the edge
            //     this->networkNode->meshFaces().clear();
            // }
        }
        else
        {
            SegmentSegmentDistance<3> ssd(periodicPlaneEdge.first->meshIntersection->P0,periodicPlaneEdge.first->meshIntersection->P1,
            periodicPlaneEdge.second->meshIntersection->P0,periodicPlaneEdge.second->meshIntersection->P1);
            assert(ssd.dMin<FLT_EPSILON &&  "The two periodic plane edges must intersect");
            assert((periodicPlaneEdge.first->patch->shift-periodicPlaneEdge.second->patch->shift).norm()<FLT_EPSILON && "Two patch shifts must match ");
            localPosition=0.5 * (ssd.x0 + ssd.x1);
        }
        globalPosition = localPosition - periodicPlaneEdge.first->patch->shift;

        // This global local is necessary so as to have the precision for the boundary node and eliminate any floating point error

        SplineNodeType::set_P(globalPosition);
        this->networkNode->set_P(localPosition);
        // this->networkNode->updateConfinement();
    }

//Yash's Version temporary

    // template <int dim>
    // void DislocationLoopNode<dim>::set_P(const typename DislocationLoopNode<dim>::VectorLowerDim& newP)
    // {
    //     VerboseDislocationLoopNode(2,"DislocationLoopNode "<<this->tag()<<" set_P (lowerDim)"<<std::endl;);
    //     assert(periodicPlaneEdge);

    //     //If we can get rid of the boundary junctions we do not need the snap operation
    //     const VectorDim globalRVEPosition(this->networkNode->snapToGlidePlanes(this->loop()->periodicGlidePlane->referencePlane->globalPosition(newP)+periodicPlaneEdge->patch->shift));
    //     SplineNodeType::set_P(globalRVEPosition-periodicPlaneEdge->patch->shift);
    //     // std::cout<<"this->networkNode "<<this->networkNode->sID<<this->networkNode->glidePlanes().size()<<std::endl;
    //     this->networkNode->set_P(globalRVEPosition);
    // }

    
    template <int dim>
    bool DislocationLoopNode<dim>::isPhysicalBoundaryNode() const
    {
        for(const auto& edge : {periodicPlaneEdge.first,periodicPlaneEdge.second})
        {
            if(edge)
            {
                for(const auto* face : edge->meshIntersection->faces)
                {
                    if(!face->periodicFacePair.second || face->regionIDs!=face->periodicFacePair.second->regionIDs)
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    template <int dim>
    bool DislocationLoopNode<dim>::isPeriodicSubdivision() const
    {
        return periodicPlaneEdge.first && !periodicContourVertex && !isPhysicalBoundaryNode();
    }

    template <int dim>
    void DislocationLoopNode<dim>::set_P(const typename DislocationLoopNode<dim>::VectorDim& newP)
    {
        if(!this->network().trialNodeMotion)
        {
            // Public motion moves every loop representation of the junction.
            const VectorDim mapped=newP+(_periodicPlanePatch ? _periodicPlanePatch->shift : VectorDim::Zero());
            this->networkNode->trySet_P(mapped);
            return;
        }
        if(!isPeriodicSubdivision())
        {
            // The trial point may be outside a physical face. Never assign it
            // to the network node before the retained contour has been found.
            SplineNodeType::set_P(this->loop()->glidePlane ? this->loop()->glidePlane->snapToPlane(newP) : newP);
        }
    }

    template <int dim>
    std::shared_ptr<PeriodicPlanePatch<dim>> DislocationLoopNode<dim>::periodicPlanePatch() const
    {
        return _periodicPlanePatch;
    }

    template <int dim>
    void DislocationLoopNode<dim>::updateGeometry()
    {

    }

    template <int dim>
    bool DislocationLoopNode<dim>::isGeometricallyRemovable(const double &Lmin,const double &Lmax,const double& relAreaTh)
    {
        const auto pPrev(periodicPrev());
        const auto pNext(periodicNext());
        if((pNext->get_P()-pPrev->get_P()).norm()<0.95*Lmax)
        {
            const double deltaArea(0.5 * (this->get_P() - pPrev->get_P()).cross(pNext->get_P() - this->get_P()).norm());

            // VectorDim vectorArea(VectorDim::Zero());
            // if (this->loop()->loopLinks().size())
            // {
            //     const VectorDim P0((*this->loop()->loopLinks().begin())->source->get_P());
            //     for (const auto &loopLink : this->loop()->loopLinks())
            //     {
            //         vectorArea += 0.5 * (loopLink->source->get_P() - P0).cross(loopLink->sink->get_P() - loopLink->source->get_P());
            //     }
            // }

            // const double loopArea(vectorArea.norm());
            const double loopArea(this->loop()->currentAreaVector().norm());
            const double loopAreaRate(this->loop()->currentAreaRate());
            // std::cout<<" For loop node "<<this->tag()<<" with "<<this->loop()->loopLinks().size()<<" loopLinks. loop area is "<<std::setprecision(15)<<loopArea<<std::endl;
            if (loopArea > std::pow(Lmin, 2))
            {
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " loopArea= " << loopArea << std::endl;);
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " deltaArea= " << deltaArea << std::endl;);
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " relAreaTh= " << relAreaTh << std::endl;);
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " removable? " << (deltaArea / loopArea < relAreaTh) << std::endl;);

                return deltaArea / loopArea < relAreaTh;
            }
            else if (loopArea > FLT_EPSILON && loopAreaRate < 0.0)
            {
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " removable (small loop)" << std::endl;);
                return true;
            }
            else
            {
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (zero area loop)" << std::endl;);
                return false;
            }
        }
        else
        {
            VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (short segment)" << std::endl;);
            return false;
        }
    }
    
    //Version which can classiy the boundary nodes as well...
    template <int dim>
    std::pair<bool,size_t> DislocationLoopNode<dim>::isRemovable(const double &Lmin,const double &Lmax, const double &relAreaTh)
    {
       //The pair will be populated as follows:
        //If the node is removable and it does not have a twin the pair will contain the same node
        //Otherwise it will contrain the twin (in the removal then we have to remove both the current node and twin)

        VerboseDislocationLoopNode(2, " Checking if DislocationLoopNode " << this->tag() << " is removable " << std::endl;);
        if (periodicPlaneEdge.first)
        { // a boundary node
            VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " onPlaneEdge" << std::endl;);
            const auto pPrev(periodicPrev());
            const auto pNext(periodicNext());
            if (pPrev && pNext)
            {
                VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (both periodicPrev and periodicNext exist)" << std::endl;);
                return std::make_pair(false,0);
            }
            else
            {
                const double loopArea(this->loop()->currentAreaVector().norm());
                if (loopArea < std::pow(Lmin, 2))
                {
                    VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << "  removable (on Edge via Geometric condition)" << std::endl;);
                    return std::make_pair(true,this->sID);
                }
                else
                {
                    VerboseDislocationLoopNode(2, "  DislocationLoopNode " << this->tag() << " NOT removable (on Edge)" << std::endl;);
                    return std::make_pair(false, 0);
                }
                
            }

        }
        else
        { // not a boundary node
            VerboseDislocationLoopNode(2, " Dislocation LoopNode Not on Boundary " << std::endl;);

            const auto pPrev(periodicPrev());
            const auto pNext(periodicNext());
            // if ((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax)
            // {
            if (pPrev && pNext)
            {
                VerboseDislocationLoopNode(2, "  Has pPrev and pNext" << std::endl;);

                if (this->networkNode->loopNodes().size() == 1)
                {
                    VerboseDislocationLoopNode(2, "  Has 1 loop nodes" << std::endl;);
                    if ((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax
//                    if (((pPrev->get_P() - pNext->get_P()).norm() < this->network().networkRemesher.Lmax || pPrev->periodicPlanePatch()!=pNext->periodicPlanePatch())
                        && isGeometricallyRemovable(Lmin, Lmax,relAreaTh))
                    {
                        return std::make_pair(true,this->sID);
                    }
                    else if (this->networkNode->isOnExternalBoundary())
                    {
                        //This means that an internal node is exactly on the the boundary based on the geometrically removable condition, the node can be removed
                        if (isGeometricallyRemovable(Lmin, Lmax, relAreaTh))
                        {
                            return std::make_pair(true,this->sID);
                        }
                        else
                        {
                            return std::make_pair(false,0);
                        }
                        
                    }
                    else
                    {
                        return std::make_pair(false,0);
                    }
                }
                else if (this->networkNode->loopNodes().size() >= 2)
                {
                    VerboseDislocationLoopNode(2, "  Has 2 loop nodes" << std::endl;);
                    const auto prevTwin(this->prev.second->twinnedLink());
                    const auto nextTwin(this->next.second->twinnedLink());

                    if (prevTwin && nextTwin)
                    {
                        VerboseDislocationLoopNode(2, "  Has prevTwin and nextTwin" << std::endl;);

                        if (prevTwin->loop == nextTwin->loop && (prevTwin->periodicPlanePatch() == nextTwin->periodicPlanePatch()))
                       {
                           VerboseDislocationLoopNode(2, "  Same patches" << std::endl;);
                            //Check based on neighbors of the networkNode
                            //Here the other node will be inserted
                            size_t prevtwinID((prevTwin->source->networkNode==this->networkNode)? prevTwin->source->sID : prevTwin->sink->sID);
                            size_t nexttwinID((nextTwin->source->networkNode==this->networkNode)? nextTwin->source->sID : nextTwin->sink->sID);
                            

                            if (prevtwinID==nexttwinID)
                            {
                                VerboseDislocationLoopNode(2, "  Same IDs" << std::endl;);
                                return std::make_pair(true,prevtwinID);
                            }
                            else
                            {
                                // The other loop may visit this network node at
                                // distinct occurrences. They are not a removable pair.
                                return std::make_pair(false,0);
                            }
                            
                        }
                        else
                        {
                            VerboseDislocationLoopNode(2, "  Different patches" << std::endl;);
                            return std::make_pair(false, 0);
                        }
                    }
                    else
                    {
                        VerboseDislocationLoopNode(2, "  No prevTwin or nextTwin" << std::endl;);
                        return std::make_pair(false,0);
                    }
                }
                else
                {
                    VerboseDislocationLoopNode(2, "this->networkNode->loopNodes().size()="<<this->networkNode->loopNodes().size()<< std::endl;);
                    return std::make_pair(false,0);
                }
            }
            else
            {
                VerboseDislocationLoopNode(2, "  No pPrev or pNext" << std::endl;);
                return std::make_pair(false,0);
            }

           
        }
    }

    template <int dim>
    bool DislocationLoopNode<dim>::isMovableTo(const VectorDim& X) const
    {
       VerboseDislocationLoopNode(4,"checking if DislocationLoopNode "<<this->sID<< " isMovable:"<<std::endl;);

        bool isMovable=true;
        
        if(periodicPlaneEdge.first)
        {
            isMovable=(isMovable && ((X-this->get_P()).squaredNorm()<FLT_EPSILON));
        }
        
        if(isMovable && this->loop()->glidePlane)
        {
            isMovable=(isMovable && this->loop()->glidePlane->contains(X));
        }
        else
        {
            isMovable=false;
        }
        
        return isMovable;
    }
    
    template <int dim>
    int DislocationLoopNode<dim>::verboseDislocationLoopNode=0;

    template class DislocationLoopNode<3>;
}
#endif
