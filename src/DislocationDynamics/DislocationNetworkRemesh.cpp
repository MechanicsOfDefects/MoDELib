/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetworkRemesh_cpp_
#define model_DislocationNetworkRemesh_cpp_

#include <DislocationNetworkRemesh.h>

#include <DislocationNetwork.h>

namespace model
{
    /**************************************************************************/
    template <int dim>
    DislocationNetworkRemesh<dim>::DislocationNetworkRemesh(DislocationNetwork<dim>& DN_in):
    /* init */ DN(DN_in)
    /* init */,Lmax(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("Lmax",true))
    /* init */,Lmin(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("Lmin",true))
    /* init */,absoluteAreaThreshold(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("absoluteAreaThreshold",true))
    /* init */,relativeAreaThreshold(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("relativeAreaThreshold",true))
    /* init */,remeshFrequency(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("remeshFrequency",true))
    {
        
        assert(Lmin<=Lmax);
        assert(Lmax>3.0*Lmin);
        assert(Lmin>=0.0);
        assert(relativeAreaThreshold>=0.0);
    }
    
    /**************************************************************************/
    template <int dim>
    void DislocationNetworkRemesh<dim>::remesh(const long int &runID)
    { /*! Performs remeshByContraction and then remeshByExpansion.
       * This order guarantees that 2-vertex NetworkComponents are expanded.
       */
        if (remeshFrequency)
        {
            if (!(runID % remeshFrequency))
            {
                remeshByRemoval();
                //                    remeshByContraction();
                remeshByExpansion();
                contract0chordSegments();
                remove0AreaLoopAcrossBnd();
                removeCollapsedLoops();
//                contractBoundaryNodes(); //Added by Yash
            }
        }
    }

template <int dim>
void DislocationNetworkRemesh<dim>::removeCollapsedLoops()
{
        std::cout << "        Remeshing network: removeCollapsedLoops... " << std::flush;
//        std::cout<<"Removing collapsed via climb... " << std::flush;
        const auto t0 = std::chrono::system_clock::now();
        std::deque<size_t> loopIDs;
        for(const auto& loop : DN.loops())
        {
            const auto current=loop.second.lock();
            const VectorDim area=current->currentAreaVector();
            const double magnitude=area.norm();
            // rightHandedUnitNormal still holds the last committed timestep.
            // Do not call updateGeometry here: it would overwrite that history.
            const bool flipped=magnitude>FLT_EPSILON &&
                current->rightHandedUnitNormal().dot(area/magnitude)<-FLT_EPSILON;
            if(current->isIsolated() && (flipped || magnitude<absoluteAreaThreshold))
            {
                loopIDs.push_back(loop.second.lock()->sID);
            }
        }
        for(const auto& loopID : loopIDs)
        {
            DN.deleteLoop(loopID);
        }
        std::cout<<"("<<loopIDs.size()<<" removed)"<<std::endl;
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }
    
    template <int dim>
    size_t DislocationNetworkRemesh<dim>::simplifyPhysicalBoundaryChains()
    {
        using LoopNode=typename DislocationNetwork<dim>::LoopNodeType;
        size_t removed=0;
        bool changed=true;
        while(changed)
        {
            changed=false;
            std::vector<size_t> candidates;
            for(const auto& entry : DN.networkNodes())
            {
                candidates.push_back(entry.first);
            }
            for(const size_t id : candidates)
            {
                const auto node=DN.networkNodes().get(id);
                if(!node || node->neighbors().size()!=2 || node->loopNodes().empty())
                {
                    continue;
                }
                std::vector<std::shared_ptr<LoopNode>> representations;
                std::set<size_t> loops;
                std::pair<size_t,size_t> endpoints;
                bool removable=true;
                for(const auto* ln : node->loopNodes())
                {
                    const auto* prev=ln->prev.first;
                    const auto* next=ln->next.first;
                    const auto& edge=ln->periodicPlaneEdge.first;
                    // Keep corners, seam crossings, and changes in loop incidence.
                    // Every representation must traverse the same physical edge in
                    // one patch and connect the same two network neighbors.
                    if(!prev || !next || !edge || ln->periodicPlaneEdge.second ||
                       !ln->isPhysicalBoundaryNode() || edge->meshIntersection->faces.size()!=1 ||
                       ln->periodicPlanePatch()!=prev->periodicPlanePatch() ||
                       ln->periodicPlanePatch()!=next->periodicPlanePatch() ||
                       (prev->periodicPlaneEdge.first!=edge && prev->periodicPlaneEdge.second!=edge) ||
                       (next->periodicPlaneEdge.first!=edge && next->periodicPlaneEdge.second!=edge) ||
                       ln->loop()->loopLinks().size()<=3 || !loops.insert(ln->loop()->sID).second)
                    {
                        removable=false;
                        break;
                    }
                    const auto pair=std::minmax(prev->networkNode->sID,next->networkNode->sID);
                    const std::pair<size_t,size_t> current(pair.first,pair.second);
                    if(current.first==current.second || current.first==id || current.second==id ||
                       (!representations.empty() && endpoints!=current))
                    {
                        removable=false;
                        break;
                    }
                    const VectorDim chord=next->get_P()-prev->get_P();
                    const double length2=chord.squaredNorm();
                    if(length2<=0.0)
                    {
                        removable=false;
                        break;
                    }
                    const VectorDim offset=ln->get_P()-prev->get_P();
                    const double t=offset.dot(chord)/length2;
                    const double tolerance=1.e-10*std::max(1.0,std::sqrt(length2));
                    if(t<=0.0 || t>=1.0 || (offset-t*chord).norm()>tolerance)
                    {
                        removable=false;
                        break;
                    }
                    endpoints=current;
                    representations.push_back(DN.loopNodes().get(ln->sID));
                }
                // Do not merge this chain into an existing, unrelated network edge.
                if(!removable || DN.networkLinks().get(endpoints))
                {
                    continue;
                }
                // Hold all representations alive until every loop has been reconnected.
                // The network link factory shares the replacement segment between them.
                for(const auto& ln : representations)
                {
                    DN.removeLoopNode(ln);
                }
                ++removed;
                changed=true;
            }
        }
        return removed;
    }

    /**************************************************************************/
    template <int dim>
    void DislocationNetworkRemesh<dim>::remeshByRemoval()
    {
        const auto t0 = std::chrono::system_clock::now();
        std::cout << "        Remeshing network: removing... " << std::flush;
        simplifyPhysicalBoundaryChains();
        DN.danglingBoundaryLoopNodes.clear();
        
        std::vector<std::shared_ptr<typename DislocationNetwork<dim>::LoopNodeType>> tempLoopNodes;
        for (const auto &loopNode : DN.loopNodes())
        {
            tempLoopNodes.push_back(loopNode.second.lock());
        }

        std::set<std::tuple<size_t,size_t,size_t>> bndLoopNodes; //First is the ID of the loop node, second is the periodicPrev and third is the periodic Next
        //Definition of bndLoopNode as a set of tuple is very important as remesh may remove the loop nodes hence changing the connectivity of the loopnodes at the boundary
        //IF that happens the corresponding boundary node should be classified as a dangling loop node
        //Also believe that this is only required here as nothing else can change the configuration of the loops near the boundary
        size_t Nremoved = 0;
        std::set<const typename DislocationNetwork<dim>::LoopNodeType*> removedLoopNodes;

        // std::deque<size_t> toBeRemoved;
        for (const auto &loopNode : DN.loopNodes())
        {
            //                std::cout<<"node "<<node.second->sID<<" "<<node.second->isSimpleBoundaryNode()<<" "<<node.second->isSimpleGrainBoundaryNode()<<std::endl;
//            if (!loopNode.second.lock()->networkNode->masterNode && removedLoopNodes.find(loopNode.second.lock().get())==removedLoopNodes.end())
            if (removedLoopNodes.find(loopNode.second.lock().get())==removedLoopNodes.end())
            {
                
                const auto isRemovableLoopNode(loopNode.second.lock()->isRemovable(Lmin, Lmax, relativeAreaThreshold));
                // std::cout<<" Trying to remove "<<loopNode.second.lock()->tag()<<" => "<<isRemovableLoopNode.first<<std::endl;
                if (isRemovableLoopNode.first)
                {
                    //Do this for all the loopNode of the networkNode 
                    if (isRemovableLoopNode.second==loopNode.second.lock()->sID)
                    {
                        //Only the current node is needed to be removed
                        if (!loopNode.second.lock()->periodicPlaneEdge.first) //No need to update the boundary nodes for the nodes at the boundary
                        {
                            for (const auto &bndNode : loopNode.second.lock()->boundaryPrev())
                            {
                                // std::cout << " For " << loopNode.second.lock()->tag() << " bndPrev is " << bndNode->tag() << std::endl;
                                bndLoopNodes.emplace(bndNode->sID,bndNode->periodicPrev()->sID,bndNode->periodicNext()->sID);
                            }
                            for (const auto &bndNode : loopNode.second.lock()->boundaryNext())
                            {
                                // std::cout << " For " << loopNode.second.lock()->tag() << " bndNext is " << bndNode->tag() << std::endl;
                                bndLoopNodes.emplace(bndNode->sID,bndNode->periodicPrev()->sID,bndNode->periodicNext()->sID);
                            }
                        }

                        // std::cout<<" Starting to remove the loop Node "<<std::endl;
//                        std::cout<<"Removing "<<loopNode.second.lock()->tag()<<std::endl;
                        DN.removeLoopNode(loopNode.second.lock()->sID);
                        // std::cout<<" Removed "<<loopNode.second.lock()->tag()<<std::endl;

                        removedLoopNodes.insert(loopNode.second.lock().get());
                        Nremoved++;
                    }
                    else
                    {
                        //Both the current and the twinned loopNode are needed to be removed
                        
                        //Collect corresponding to the twinned loopNode and current loopNode
                        for (const auto& loopN : loopNode.second.lock()->networkNode->loopNodes())
                        {
                            if (loopN->sID==loopNode.second.lock()->sID || loopN->sID==isRemovableLoopNode.second)
                            {
                                if (!loopN->periodicPlaneEdge.first)
                                {
                                    for (const auto &bndNode : loopN->boundaryPrev())
                                    {
                                        bndLoopNodes.emplace(bndNode->sID,bndNode->periodicPrev()->sID,bndNode->periodicNext()->sID);

                                    }
                                    for (const auto &bndNode : loopN->boundaryNext())
                                    {
                                        bndLoopNodes.emplace(bndNode->sID,bndNode->periodicPrev()->sID,bndNode->periodicNext()->sID);
                                    }
                                }
//                                std::cout<<"Removing LoopNode"<<loopN->tag()<<std::endl;
                                DN.removeLoopNode(loopN->sID);
                                removedLoopNodes.insert(loopN);
                                Nremoved++;
                            }
                        }
                    }
                }
            }
        }
        
        tempLoopNodes.clear(); //Very important to clear the loop node here
        std::cout << " (" << Nremoved << " LoopNodes removed) " << std::flush;
        Nremoved=0;
        for (const auto &nodeID : bndLoopNodes)
        {
            auto nodeIter(DN.loopNodes().find(std::get<0>(nodeID)));
            if (nodeIter != DN.loopNodes().end())
            {
                auto bndNode(nodeIter->second.lock());

                const auto pPrev(bndNode->periodicPrev());
                const auto pNext(bndNode->periodicNext());
                if (pPrev && pNext)
                {
                    if (pPrev->sID == std::get<1>(nodeID) && pNext->sID== std::get<2>(nodeID))
                    {
                        //If the periodic Prev and periodic next are the same as the remesh may change this
                        const auto pPrevLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pPrev->get_P()));
                        const auto pNextLocal(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(pNext->get_P()));
                        //

                        SegmentSegmentDistance<dim - 1> ssd1(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.first->source, *bndNode->periodicPlaneEdge.first->sink);
                        if (ssd1.dMin < FLT_EPSILON)
                        {
                            //                    VerboseDislocationLoopNode(3,"dMin= "<<ssd.dMin<<std::endl;);
                            if (bndNode->periodicPlaneEdge.second)
                            {
                                SegmentSegmentDistance<dim - 1> ssd2(pPrevLocal, pNextLocal, *bndNode->periodicPlaneEdge.second->source, *bndNode->periodicPlaneEdge.second->sink);
                                const VectorLowerDim ssd1Pos(0.5 * (ssd1.x0 + ssd1.x1));
                                const VectorLowerDim ssd2Pos(0.5 * (ssd2.x0 + ssd2.x1));
                                assert((ssd1Pos - ssd2Pos).norm() < FLT_EPSILON && "The two positions must match ");
                            }
                            bndNode->set_P(VectorLowerDim(0.5 * (ssd1.x0 + ssd1.x1)));
                        }
                        else
                        {
                            bndNode->set_P(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                            DN.danglingBoundaryLoopNodes.insert(bndNode.get());
                        }
                    }
                    else
                    {
                        bndNode->set_P(bndNode->loop()->periodicGlidePlane->referencePlane->localPosition(bndNode->get_P()));
                        DN.danglingBoundaryLoopNodes.insert(bndNode.get());
                    }
                    
                    
                }
                else
                {
                    //Sole boundary node as the survivor..Remove them
                    DN.removeLoopNode(bndNode->sID);
                    Nremoved++;
                }
                
                
            }
        }

        std::cout << " (" << Nremoved << " bndLoopNodes removed). Updating boundary nodes " << std::flush;

        if (DN.danglingBoundaryLoopNodes.size())
        {
            DN.updateBoundaryNodes();
        }

        std::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]." << defaultColor << std::endl;
    }
    
    /**************************************************************************/
    template <int dim>
    void DislocationNetworkRemesh<dim>::remeshByExpansion()
    {
        const auto t0 = std::chrono::system_clock::now();
        std::cout << "        Remeshing network: expanding... " << std::flush;
        std::set<std::pair<size_t, size_t>> toBeExpanded;
        for (const auto &linkIter : DN.networkLinks())
        {
            const auto sharedLink(linkIter.second.lock());
            
            if (!sharedLink->hasZeroBurgers() && !sharedLink->isSessile() && !sharedLink->isBoundarySegment() && !sharedLink->isGrainBoundarySegment())
            {
                const VectorDim chord(sharedLink->chord()); // this is sink->get_P() - source->get_P()
                const double chordLength(chord.norm());
                
//                if (sharedLink->source->loopNodes().size() > 1 && sharedLink->sink->loopNodes().size() > 1
//                    /*&& chord.dot(dv)>vTolexp*chordLength*dv.norm()*/
//                    && chordLength > 3.0 * Lmin)
//                { // also expands a straight line to generate glissile segment
//                    toBeExpanded.insert(std::make_pair(sharedLink->source->sID, sharedLink->sink->sID));
//                }
                
                // Expand segments shorter than Lmax
                if (chordLength > Lmax)
                {
                    toBeExpanded.insert(std::make_pair(sharedLink->source->sID, sharedLink->sink->sID));
                }
            }
        }
        
        // Call Network::expand
        unsigned int Nexpanded(0);
        const double expand_at(0.5);
        for (const auto &expIter : toBeExpanded)
        {
            //            const size_t i(expIter->first);
            //            const size_t j(expIter->second);
            const auto source(DN.networkNodes().get(expIter.first));
            const auto sink(DN.networkNodes().get(expIter.second));
            
            const auto Lij(DN.networkLinks().get(expIter));
            if (Lij)
            {
                VectorDim expandPoint(Lij->get_r(expand_at));
                
                if ((expandPoint - source->get_P()).squaredNorm() > FLT_EPSILON && (expandPoint - sink->get_P()).squaredNorm() > FLT_EPSILON)
                {
                    
                    const auto newNetNode(DN.networkNodes().create(expandPoint, 0.5 * (source->get_V() + sink->get_V()),0.5 * (source->climbVelocityScalar + sink->climbVelocityScalar) ,0.5 * (source->velocityReduction() + sink->velocityReduction())));
//                    std::cout<<"Expanding "<<Lij->tag()<<std::endl;

                    DN.expandNetworkLink(Lij, newNetNode);
                    Nexpanded++;
                }
            }
        }
        std::cout << " (" << Nexpanded << " expanded)" << std::flush;
        std::cout << magentaColor << std::setprecision(3) << std::scientific << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]." << defaultColor << std::endl;
    }
    
    /**************************************************************************/
    template <int dim>
    void DislocationNetworkRemesh<dim>::contractBoundaryNodes() //This function contracts boundary nodes if they are at the same position and share atleast one same neightbors
    {
        std::cout << "Contracting BoundaryNodes... " << std::flush;
        const auto t0 = std::chrono::system_clock::now();
        
        typedef std::pair<std::shared_ptr<typename DislocationNetwork<dim>::NetworkNodeType>, std::shared_ptr<typename DislocationNetwork<dim>::NetworkNodeType>> NetNodePairType;
        typedef std::set<NetNodePairType> BndSetType;
        
        std::map<NetNodePairType, BndSetType> bndNodesContractionMap;
        
        for (const auto &netNode : DN.networkNodes())
        {
            if (netNode.second.lock()->isBoundaryNode())
            {
                //                std::cout<<"NetNode is "<<netNode.second.lock()->tag()<<std::endl;
                
                std::set<std::shared_ptr<typename DislocationNetwork<dim>::NetworkNodeType>> internalNodes; // used to be regular ptrs
                for (const auto &loopNode : netNode.second.lock()->loopNodes())
                {
                    internalNodes.insert(loopNode->periodicPrev()->networkNode);
                    internalNodes.insert(loopNode->periodicNext()->networkNode);
                }
                
                if (internalNodes.size() == 2)
                {
                    
                    //                    std::cout<<"internal nodes are "<<(*internalNodes.begin())->tag()<<","<<(*internalNodes.rbegin())->tag()<<std::endl;
                    
                    const NetNodePairType key((*internalNodes.begin())->sID < (*internalNodes.rbegin())->sID ? std::make_pair(*internalNodes.begin(), *internalNodes.rbegin()) : std::make_pair(*internalNodes.rbegin(), *internalNodes.begin())); // Line that was causing problem
                    
                    if (bndNodesContractionMap.find(key) == bndNodesContractionMap.end())
                    { // key not found
                        std::map<std::set<const Polyhedron<dim,dim-1> *>, std::set<typename DislocationNetwork<dim>::NetworkNodeType *>> bndNodeMap;
                        
                        for (const auto &loopNode : (*internalNodes.begin())->loopNodes())
                        {
                            if (loopNode->periodicPrev()->networkNode == *internalNodes.rbegin())
                            {
                                for (const auto &bndNode : loopNode->boundaryPrev())
                                {
                                    //                                                        std::cout<<"adding bndNode "<<bndNode->networkNode->tag()<<std::endl;
                                    bndNodeMap[bndNode->networkNode->meshFaces()].insert(bndNode->networkNode.get());
                                }
                            }
                            else if (loopNode->periodicNext()->networkNode == *internalNodes.rbegin())
                            {
                                for (const auto &bndNode : loopNode->boundaryNext())
                                {
                                    //                                                        std::cout<<"adding bndNode "<<bndNode->networkNode->tag()<<std::endl;
                                    bndNodeMap[bndNode->networkNode->meshFaces()].insert(bndNode->networkNode.get());
                                }
                            }
                        }
                        
                        for (const auto &pair : bndNodeMap)
                        {
                            //                                                std::cout<<"Nodes of face"<<std::endl;
                            //                                                for(const auto& face : pair.first)
                            //                                                {
                            //                                                    std::cout<<face->outNormal().transpose()<<std::endl;
                            //                                                }
                            
                            //                                                std::cout<<"nodes size="<<pair.second.size()<<std::endl;
                            if (pair.second.size() > 1)
                            {
                                for (const auto &tempNode : pair.second)
                                {
                                    //                                                        std::cout<<"node size="<<tempNode->sID<<std::endl;
                                    if (tempNode != *pair.second.begin())
                                    {
                                        if ((tempNode->get_P() - (*pair.second.begin())->get_P()).norm() < FLT_EPSILON)
                                        {
                                            const NetNodePairType contractPair(tempNode->sID < (*pair.second.begin())->sID ? std::make_pair(tempNode, *pair.second.begin()) : std::make_pair(*pair.second.begin(), tempNode));
                                            //                                                                std::cout<<"contract pair="<<contractPair.first->sID<<" "<<contractPair.second->sID<<std::endl;
                                            bndNodesContractionMap[key].insert(contractPair);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //        std::cout<<"Finished selection"<<std::endl;
        
        for (const auto &pair : bndNodesContractionMap)
        {
            std::cout << "Contracting bnd nodes of intenal nodes " << pair.first.first->tag() << " " << pair.first.second->tag() << std::endl;
            for (const auto &bndPair : pair.second)
            {
                std::cout << "Contracting bnd nodes  " << bndPair.first->tag() << " " << bndPair.second->tag() << std::endl;
                DN.contract(bndPair.first, bndPair.second);
            }
        }
        std::cout << magentaColor << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]" << defaultColor << std::endl;
    }
    
    /**************************************************************************/
    template <int dim>
    bool DislocationNetworkRemesh<dim>::zeroChordCandidate(
        const std::shared_ptr<typename TraitsType::NetworkLinkType>& link) const
    {
        if(!link || link->hasZeroBurgers())
        {
            return false;
        }
        bool physical=false;
        for(const auto& node : {link->source,link->sink})
        {
            for(const auto* ln : node->loopNodes())
            {
                // Periodic subdivisions are rebuilt from loop topology, not
                // independently contracted in RVE coordinates.
                if(ln->isPeriodicSubdivision())
                {
                    return false;
                }
                physical=physical || ln->isPhysicalBoundaryNode();
            }
        }
        const VectorDim chord=link->sink->get_P()-link->source->get_P();
        if(chord.norm()<=FLT_EPSILON)
        {
            return true;
        }
        if(physical)
        {
            return false;
        }
        const auto sourcePlanes=link->source->glidePlanes();
        const auto sinkPlanes=link->sink->glidePlanes();
        return sourcePlanes.size()==dim-1 && sourcePlanes==sinkPlanes &&
            chord.squaredNorm()+chord.dot(link->sink->get_V()-link->source->get_V())*
                DN.ddBase.simulationParameters.dt<FLT_EPSILON;
    }

    template <int dim>
    size_t DislocationNetworkRemesh<dim>::contractZeroChordCandidates(
        const std::vector<std::pair<size_t,size_t>>& candidates)
    {
        size_t contracted=0;
        for(const auto& key : candidates)
        {
            const auto link=DN.networkLinks().get(key);
            // Earlier contractions can move endpoints and change confinement.
            if(zeroChordCandidate(link))
            {
                // DN.contract owns the complete boundary reconciliation and
                // shared-loop update. No per-occurrence boundary motion here.
                contracted+=DN.contract(link->source,link->sink);
            }
        }
        return contracted;
    }

    template <int dim>
    void DislocationNetworkRemesh<dim>::contract0chordSegments()
    {
        std::cout<<"        Contracting zero-chord segments... "<<std::flush;
        const auto t0=std::chrono::system_clock::now();
        std::set<std::pair<double,std::pair<size_t,size_t>>> ordered;
        for(const auto& entry : DN.networkLinks())
        {
            const auto link=entry.second.lock();
            if(zeroChordCandidate(link))
            {
                ordered.emplace((link->sink->get_P()-link->source->get_P()).norm(),link->key);
            }
        }
        std::vector<std::pair<size_t,size_t>> candidates;
        for(const auto& entry : ordered)
        {
            candidates.push_back(entry.second);
        }
        const auto contracted=contractZeroChordCandidates(candidates);
        std::cout<<"("<<contracted<<" contracted)"<<magentaColor<<" ["
                 <<std::chrono::duration<double>(std::chrono::system_clock::now()-t0).count()
                 <<" sec]"<<defaultColor<<std::endl;
    }

    template <int dim>
    void DislocationNetworkRemesh<dim>::remove0AreaLoopAcrossBnd()
    {
        std::cout << "        Removing zero-area loops across boundary... " << std::flush;
        const auto t0 = std::chrono::system_clock::now();
        size_t nRemoved(0);

        std::vector<size_t> loopToBeRemoved;
        for(const auto& entry : DN.loops())
        {
            const auto loop=entry.second.lock();
            if(loop->glidePlane && loop->currentAreaVector().norm()<FLT_EPSILON)
            {
                size_t anchors=0;
                size_t periodicSubdivisions=0;
                for(const auto* link : loop->linkSequence())
                {
                    if(link->source->isPeriodicSubdivision())
                    {
                        ++periodicSubdivisions;
                    }
                    else
                    {
                        // Physical surfaces and grain boundaries are retained
                        // contour vertices, including cross-grain periodic faces.
                        ++anchors;
                    }
                }
                if(anchors<3 && periodicSubdivisions>0)
                {
                    loopToBeRemoved.push_back(loop->sID);
                }
            }
        }
        // Remove only this loop's memberships. Other loops sharing network
        // nodes and segments keep their topology and Burgers contributions.
        for(const auto id : loopToBeRemoved)
        {
            DN.deleteLoop(id);
            ++nRemoved;
        }
        std::cout << "(" << nRemoved << " Loops Removed)" << std::flush;
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }
    
//    /**************************************************************************/
//    template <int dim>
//    double DislocationNetworkRemesh<dim>::minMeshSize(const SimplicialMesh<dim> &mesh)
//    {
//        return std::min(mesh.xMax(0) - mesh.xMin(0), std::min(mesh.xMax(1) - mesh.xMin(1), mesh.xMax(2) - mesh.xMin(2)));
//    }
    
    template class DislocationNetworkRemesh<3>;

    
} // namespace model
#endif

