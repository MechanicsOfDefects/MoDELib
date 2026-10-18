/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationJunctionFormation_cpp_
#define model_DislocationJunctionFormation_cpp_

#include <utility> // for std::pair
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <SegmentSegmentDistance.h>
#include <SweepPlane.h>
#include <TypeTraits.h>
#include <StressStraight.h>

#include <DislocationJunctionFormation.h>

#ifndef NDEBUG
#define VerboseJunctions(N,x) if(verboseJunctions>=N){std::cout<<x;}
#else
#define VerboseJunctions(N,x)
#endif


namespace model
{

    template <typename DislocationNetworkType>
    void DislocationJunctionFormation<DislocationNetworkType>::insertIntersection(std::deque<IntersectionTypeContainerType> &intersectionContainer,
                                                                                  const NetworkLinkType *const linkA,
                                                                                  const NetworkLinkType *const linkB,
                                                                                  const SegmentSegmentDistance<dim> &ssd,
                                                                                  const double &currentcCollisionTOL,
                                                                                  const bool &bndJunction, const bool &gbndJunction,const int event)
    {
        VerboseJunctions(2, "insertIntersection " << linkA->tag() << " [ " << linkA->source->get_P().transpose() << ", " << linkA->sink->get_P().transpose()
                         << " ]"
                         << "," << linkB->tag() << " [ " << linkB->source->get_P().transpose() << ", " << linkB->sink->get_P().transpose()
                         << " ]" << std::endl;);
        VerboseJunctions(2, "insertIntersection "
                         << "ssd.dMin=" << ssd.dMin << std::endl;);
        VerboseJunctions(2, "insertIntersection "
                         << "currentcCollisionTOL=" << currentcCollisionTOL << std::endl;);
        
        if (ssd.dMin < currentcCollisionTOL)
        {
            bool isValidJunction(bndJunction || gbndJunction);
            if (!isValidJunction && !linkA->isBoundarySegment() && !linkB->isBoundarySegment() && !linkA->isGrainBoundarySegment() && !linkB->isGrainBoundarySegment())
            { // Check force condition for internal segments
                
                const VectorDim chordA(linkA->sink->get_P() - linkA->source->get_P());
                const double LA(chordA.norm());
                const VectorDim chordB(linkB->sink->get_P() - linkB->source->get_P());
                const double LB(chordB.norm());
                
                if (LA > FLT_EPSILON && LB > FLT_EPSILON)
                {
                    if (ssd.dMin > FLT_EPSILON)
                    {
                        StressStraight<dim> stressA(DN.ddBase.poly,ssd.x0 - infiniteLineLength / LA * chordA,
                                                    ssd.x0 + infiniteLineLength / LA * chordA,
                                                    linkA->burgers(),DN.ddBase.EwaldLength);
                        
                        StressStraight<dim> stressB(DN.ddBase.poly,ssd.x1 - infiniteLineLength / LB * chordB,
                                                    ssd.x1 + infiniteLineLength / LB * chordB,
                                                    linkB->burgers(),DN.ddBase.EwaldLength);
                        
                        VerboseJunctions(3, "Non-intersecting pair" << std::endl;);
                        const VectorDim forceOnA = (stressB.stress(ssd.x0) * linkA->burgers()).cross(chordA);
                        const VectorDim forceOnB = (stressA.stress(ssd.x1) * linkB->burgers()).cross(chordB);
                        const VectorDim dxShift(ssd.x1 - ssd.x0);
                        if (forceOnA.dot(dxShift) > FLT_EPSILON && forceOnB.dot(dxShift) < -FLT_EPSILON)
                        {
                            VerboseJunctions(3, "attractive pair 1" << std::endl;);
                            isValidJunction = true; // for non-parallel lines this neglects the energy of rotation
                        }
                        else
                        {
                            VerboseJunctions(3, "non-attractive pair" << std::endl;);
                        }
                    }
                    else
                    {
                        VerboseJunctions(3, "Intersecting pair .. Determining via frank rule" << std::endl;);
                        isValidJunction=(linkA->burgers().dot(linkB->burgers())*linkA->chord().dot(linkB->chord())<=0.0);
                    }
                }
            }
            
            if (isValidJunction)
            {
    #ifdef _OPENMP
                intersectionContainer[omp_get_thread_num()].emplace_back(linkA->key,
                                                                         linkB->key,
                                                                         ssd,event);
    #else
                intersectionContainer[0].emplace_back(linkA->key,
                                                      linkB->key,
                                                      ssd,event);
    #endif
            }
        }
        else
        {
            VerboseJunctions(3, "dMin=" << ssd.dMin << ", collisionTol=" << currentcCollisionTOL << std::endl;);
        }
    }

    template <typename DislocationNetworkType>
    void DislocationJunctionFormation<DislocationNetworkType>::findIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                                                                                 const size_t& nThreads)

    {/*! @param[in]  avoidNodeIntersection
      *  Computes all the intersections between the edges of the DislocationNetwork
      */
        
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"		Finding collisions "<<std::flush;
        
        // Use SweepPlane to compute possible intersections
        SweepPlane<NetworkLinkType,dim> swp;
        for(const auto& link : DN.networkLinks())
        {
            const auto sharedLink(link.second.lock());
            
            if(   (!sharedLink->hasZeroBurgers() )
               || sharedLink->isBoundarySegment() || sharedLink->isGrainBoundarySegment() )
                
            {
                const double x0=sharedLink->source->get_P()(0);
                const double x1=sharedLink->sink->get_P()(0);
                swp.addSegment(std::min(x0,x1)-0.5*collisionTol,
                               std::max(x0,x1)+0.5*collisionTol,*sharedLink);
            }
        }
        swp.computeIntersectionPairs();
        std::cout<<"("<<swp.potentialIntersectionPairs().size()<<" sweep-line pairs) "<<defaultColor<<std::flush;
        
        
        std::deque<std::pair<const NetworkLinkType*,const NetworkLinkType*>> reducedIntersectionPairs;
        for(size_t k=0;k<swp.potentialIntersectionPairs().size();++k)
        {
            const auto& linkA(swp.potentialIntersectionPair(k).first);
            const auto& linkB(swp.potentialIntersectionPair(k).second);
            const VectorDim& sourceA(linkA->source->get_P());
            const VectorDim&   sinkA(linkA->sink->get_P());
            const VectorDim& sourceB(linkB->source->get_P());
            const VectorDim&   sinkB(linkB->sink->get_P());
                        
            // SweepPlane sweeps along x direction, so now check possible overlap in y and z direction
            bool yDontOverlap((std::min(sourceA(1),sinkA(1))>std::max(sourceB(1),sinkB(1))+collisionTol) || (std::max(sourceA(1),sinkA(1))<std::min(sourceB(1),sinkB(1))-collisionTol));
            bool zDontOverlap((std::min(sourceA(2),sinkA(2))>std::max(sourceB(2),sinkB(2))+collisionTol) || (std::max(sourceA(2),sinkA(2))<std::min(sourceB(2),sinkB(2))-collisionTol));
                        
            if(!yDontOverlap && !zDontOverlap)
            {
                reducedIntersectionPairs.emplace_back(linkA,linkB);
            }
        }
        
        std::cout<<" ("<<reducedIntersectionPairs.size()<<" reduced pairs) "<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        
        const auto t1= std::chrono::system_clock::now();
        std::cout<<"        Selecting junctions ("<<nThreads<<" threads): "<<std::flush;
        
        //! 2- loop over all links and determine their intersections
        // #ifdef _OPENMP
        // #pragma omp parallel for
        // #endif
        for(size_t k=0;k<reducedIntersectionPairs.size();++k)
        {
            
            const auto& linkA(reducedIntersectionPairs[k].first);
            const auto& linkB(reducedIntersectionPairs[k].second);
            
            findPairIntersections(intersectionContainer,linkA,linkB);
        }

        int nIntersections=0;
        for (const auto& intersectionByThreadContainer : intersectionContainer)
        {
            nIntersections+=intersectionByThreadContainer.size();
        }
        std::cout<<nIntersections<<" physical junctions "<<std::flush;
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;

    }

    template <typename DislocationNetworkType>
    void DislocationJunctionFormation<DislocationNetworkType>::findPairIntersections(
        std::deque<IntersectionTypeContainerType>& intersectionContainer,
        const NetworkLinkType* linkA,const NetworkLinkType* linkB)
    {
            VerboseJunctions(2,"checking intersection "<<linkA->tag()<<","<<linkB->tag()<<std::endl;);
            
            
            const bool linkAisBnd(linkA->isBoundarySegment());
            const bool linkBisBnd(linkB->isBoundarySegment());
            const bool linkAisGBnd(linkA->isGrainBoundarySegment());
            const bool linkBisGBnd(linkB->isGrainBoundarySegment());
            
            
            bool bndJunction=false;
            bool gbndJunction=false;
            const auto facesB=linkB->meshFaces();
            // Boundary merging requires a common physical face, not merely a
            // boundary label on one segment or parallel faces elsewhere.
            if(linkA->chord().norm()>FLT_EPSILON && linkB->chord().norm()>FLT_EPSILON &&
               linkA->chord().normalized().cross(linkB->chord().normalized()).norm()<FLT_EPSILON &&
               linkA->glidePlaneNormal().cross(linkB->glidePlaneNormal()).norm()<FLT_EPSILON)
            {
                for(const auto* face : linkA->meshFaces())
                {
                    if(facesB.count(face))
                    {
                        const bool grainBoundary=face->regionIDs.first!=face->regionIDs.second ||
                            (face->periodicFacePair.second && face->regionIDs!=face->periodicFacePair.second->regionIDs);
                        gbndJunction=gbndJunction || grainBoundary;
                        bndJunction=bndJunction || (face->isExternal() && !face->periodicFacePair.second);
                    }
                }
            }

            VerboseJunctions(2,"links are bnd: "<<linkAisBnd<<" "<<linkBisBnd<<std::endl;);
            VerboseJunctions(2,"links are grainBnd: "<<linkAisGBnd<<" "<<linkBisGBnd<<std::endl;);
            VerboseJunctions(2,"bndJunction: "<<bndJunction<<std::endl;);
            VerboseJunctions(2,"gbndJunction: "<<gbndJunction<<std::endl;);
            
            const bool intersectionIsSourceSource(linkA->source->sID==linkB->source->sID);
            const bool intersectionIsSourceSink(  linkA->source->sID==linkB->  sink->sID);
            const bool intersectionIsSinkSource(  linkA->  sink->sID==linkB->source->sID);
            const bool intersectionIsSinkSink(    linkA->  sink->sID==linkB->  sink->sID);
            
            double currentcCollisionTOL=collisionTol;
                        
            const auto pcPlanes(linkA->parallelAndCoincidentGlidePlanes(linkB->glidePlanes()));
            VerboseJunctions(2,"pcPlanes.size()= "<<pcPlanes.size()<<std::endl;);
            for(const auto& pair : pcPlanes)
            {
                if(pair.first==pair.second)
                {// a common coincident plane found
                    const double cosTheta=(intersectionIsSourceSource||intersectionIsSinkSink)? linkA->chord().normalized().dot(linkB->chord().normalized()) : ((intersectionIsSourceSink ||intersectionIsSinkSource)? -linkA->chord().normalized().dot(linkB->chord().normalized()) : -1.0);
                    if(cosTheta<0.7)
                    {// links are either disconnected, or connected to a node and forming a large angle at that node. Reduce tolerance
                        currentcCollisionTOL=FLT_EPSILON;
                    }
                    break;
                }
                else
                {// segments on parallel (non-coincident) planes
                    currentcCollisionTOL=FLT_EPSILON;
                }
            }
            
            if(intersectionIsSourceSource)
            {
                // intersect sink of A with link B
                SegmentSegmentDistance<dim> ssdA(linkA->sink->get_P(), // fake degenerete segment at sink of A
                                                 linkA->sink->get_P(),  // fake degenerete segment at sink of A
                                                 linkB->source->get_P(),
                                                 linkB->sink->get_P(),
                                                 1.0);
                insertIntersection(intersectionContainer,linkA,linkB,ssdA,currentcCollisionTOL,bndJunction,gbndJunction,1);
                
                // intersect sink of B with link A
                SegmentSegmentDistance<dim> ssdB(linkA->source->get_P(),
                                                 linkA->sink->get_P(),
                                                 linkB->sink->get_P(),    // fake degenerete segment at sink of B
                                                 linkB->sink->get_P(),    // fake degenerete segment at sink of B
                                                 1.0);
                insertIntersection(intersectionContainer,linkA,linkB,ssdB,currentcCollisionTOL,bndJunction,gbndJunction,2);
                
            }
            else if(intersectionIsSourceSink)
            {
                // intersect sink of A with link B
                SegmentSegmentDistance<dim> ssdA(linkA->sink->get_P(), // fake degenerete segment at sink of A
                                                 linkA->sink->get_P(),  // fake degenerete segment at sink of A
                                                 linkB->source->get_P(),
                                                 linkB->sink->get_P(),
                                                 1.0);
                insertIntersection(intersectionContainer,linkA,linkB,ssdA,currentcCollisionTOL,bndJunction,gbndJunction,1);
                
                // intersect source of B with link A
                SegmentSegmentDistance<dim> ssdB(linkA->source->get_P(),
                                                 linkA->sink->get_P(),
                                                 linkB->source->get_P(),    // fake degenerete segment at source of B
                                                 linkB->source->get_P(),    // fake degenerete segment at source of B
                                                 0.0);
                insertIntersection(intersectionContainer,linkA,linkB,ssdB,currentcCollisionTOL,bndJunction,gbndJunction,2);
            }
            else if(intersectionIsSinkSource)
            {
                
                // intersect source of A with link B
                SegmentSegmentDistance<dim> ssdA(linkA->source->get_P(), // fake degenerete segment at source of A
                                                 linkA->source->get_P(),  // fake degenerete segment at source of A
                                                 linkB->source->get_P(),
                                                 linkB->sink->get_P(),
                                                 0.0);
                insertIntersection(intersectionContainer,linkA,linkB,ssdA,currentcCollisionTOL,bndJunction,gbndJunction,1);
                
                
                // intersect sink of B with link A
                SegmentSegmentDistance<dim> ssdB(linkA->source->get_P(),
                                                 linkA->sink->get_P(),
                                                 linkB->sink->get_P(),    // fake degenerete segment at sink of B
                                                 linkB->sink->get_P(),    // fake degenerete segment at sink of B
                                                 1.0);
                insertIntersection(intersectionContainer,linkA,linkB,ssdB,currentcCollisionTOL,bndJunction,gbndJunction,2);
            }
            else if(intersectionIsSinkSink)
            {
                // intersect source of A with link B
                SegmentSegmentDistance<dim> ssdA(linkA->source->get_P(), // fake degenerete segment at source of A
                                                 linkA->source->get_P(),  // fake degenerete segment at source of A
                                                 linkB->source->get_P(),
                                                 linkB->sink->get_P(),
                                                 0.0);
                insertIntersection(intersectionContainer,linkA,linkB,ssdA,currentcCollisionTOL,bndJunction,gbndJunction,1);
                
                
                // intersect source of B with link A
                SegmentSegmentDistance<dim> ssdB(linkA->source->get_P(),
                                                 linkA->sink->get_P(),
                                                 linkB->source->get_P(),    // fake degenerete segment at source of B
                                                 linkB->source->get_P(),    // fake degenerete segment at source of B
                                                 0.0);
                insertIntersection(intersectionContainer,linkA,linkB,ssdB,currentcCollisionTOL,bndJunction,gbndJunction,2);
                
            }
            else
            {// segments don't share vertices
                
                SegmentSegmentDistance<dim> ssd(linkA->source->get_P(),
                                                linkA->sink->get_P(),
                                                linkB->source->get_P(),
                                                linkB->sink->get_P()
                                                );
                insertIntersection(intersectionContainer,linkA,linkB,ssd,currentcCollisionTOL,bndJunction,gbndJunction);
            }
    }

    //     Yash's version(Working)
    template <typename DislocationNetworkType>
    std::pair<std::shared_ptr<typename DislocationJunctionFormation<DislocationNetworkType>::NetworkNodeType>, bool> DislocationJunctionFormation<DislocationNetworkType>::junctionNode(const double &t,
                                                                                                                                                                                        const VectorDim &x,
                                                                                                                                                                                        const std::shared_ptr<NetworkLinkType> &L)
    {
        //            VerboseJunctions(4, "JunctionNode for segment: " << key.first << "->" << key.second << " @ " << t << std::endl;);
        
        auto Nclose = t < 0.5 ? L->source : L->sink;
        auto Nfar = t >= 0.5 ? L->source : L->sink;
        VerboseJunctions(4, "Nclose=" << Nclose->sID << " (GP Size)=> " << Nclose->glidePlanes().size() << std::endl;);
        VerboseJunctions(4, "Nfar=" << Nfar->sID << " (GP Size)=> " << Nfar->glidePlanes().size() << std::endl;);
        //            if(t>FLT_EPSILON && t<1.0-FLT_EPSILON && !Nclose->isBoundaryNode() && !Nfar->isBoundaryNode()) //This condition is the problem It needs to be enabled
        if (t > FLT_EPSILON && t < 1.0 - FLT_EPSILON)
        { // intersection point is not an end node
            if ((Nclose->get_P() - x).norm() > DN.networkRemesher.Lmin)
            { // far enough from Nclose
                VerboseJunctions(4, "JunctionNode case a" << std::endl;);
                
                const VectorDim snappedPosition(L->glidePlanes().size() > 0 ? L->snapToGlidePlanes(x) : x);
                const auto newNetNode(DN.networkNodes().create(snappedPosition, L->source->get_V() * (1.0 - t) + L->sink->get_V() * t,L->source->climbVelocityScalar*(1.0-t)+L->sink->climbVelocityScalar*t,L->source->velocityReduction() * (1.0 - t) + L->sink->velocityReduction() * t));
                DN.expandNetworkLink(L, newNetNode);
                return std::make_pair(newNetNode, true);
            }
            else
            { // close to from Nclose
                if (Nclose->isMovableTo(x) && !Nclose->isBoundaryNode())
                {
                    VerboseJunctions(4, "JunctionNode case b" << std::endl;);
                    return std::make_pair(Nclose, false);
                }
                else
                {
                    if ((Nfar->get_P() - x).norm() > DN.networkRemesher.Lmin)
                    {
                        VerboseJunctions(4, "JunctionNode case c" << std::endl;);
                        
                        const VectorDim snappedPosition(L->glidePlanes().size() > 0 ? L->snapToGlidePlanes(x) : x);
                        
                        const auto newNetNode(DN.networkNodes().create(snappedPosition, L->source->get_V() * (1.0 - t) + L->sink->get_V() * t, L->source->climbVelocityScalar*(1.0-t)+L->sink->climbVelocityScalar*t,L->source->velocityReduction() * (1.0 - t) + L->sink->velocityReduction() * t));
                        // std::cout << "x is " << x.transpose() << std::endl;
                        
                        // std::cout<<"Expanding network Links c "<<L->tag()<<" [ "<<L->source->get_P().transpose()<<" -> "<<L->sink->get_P().transpose()<<" ] GlidePlane size "
                        // <<" [ "<<L->source->glidePlanes().size()<<" --> "<<L->sink->glidePlanes().size() <<" ] "
                        // <<L->glidePlanes().size()<<" loopLinks size C"<<L->loopLinks().size()<<std::endl;
                        
                        DN.expandNetworkLink(L, newNetNode);
                        // std::cout<<"Expanded network Links"<<std::endl;
                        
                        return std::make_pair(newNetNode, true);
                    }
                    else
                    {
                        if (Nfar->isMovableTo(x) && !Nfar->isBoundaryNode())
                        {
                            VerboseJunctions(4, "JunctionNode case d" << std::endl;);
                            return std::make_pair(Nfar, false);
                        }
                        else
                        {
                            VerboseJunctions(4, "JunctionNode case e" << std::endl;);
                            // std::cout<<"x is "<<x.transpose()<<std::endl;
                            //   std::cout<<"Expanding network Links E "<<L->tag()<<" [ "<<L->source->get_P().transpose()<<" -> "<<L->sink->get_P().transpose()<<" ] GlidePlane size "
                            // <<" [ "<<L->source->glidePlanes().size()<<" --> "<<L->sink->glidePlanes().size() <<" ] "
                            // <<L->glidePlanes().size()<<" loopLinks size E"<<L->loopLinks().size()<<std::endl;
                            
                            const VectorDim snappedPosition(L->glidePlanes().size() > 0 ? L->snapToGlidePlanes(x) : x);
                            
                            const auto newNetNode(DN.networkNodes().create(snappedPosition, L->source->get_V() * (1.0 - t) + L->sink->get_V() * t, L->source->climbVelocityScalar*(1.0-t)+L->sink->climbVelocityScalar*t,L->source->velocityReduction() * (1.0 - t) + L->sink->velocityReduction() * t));
                            DN.expandNetworkLink(L, newNetNode);
                            return std::make_pair(newNetNode, true);
                        }
                    }
                }
            }
        }
        else
        {// junction node is an existing node
            bool physicalEndpoint=false;
            for(const auto* ln : Nclose->loopNodes())
            {
                physicalEndpoint=physicalEndpoint || ln->isPhysicalBoundaryNode();
            }
            if (!Nclose->isBoundaryNode() ||
                (physicalEndpoint && (Nclose->get_P()-x).norm()<=FLT_EPSILON))
            {
                VerboseJunctions(4, "JunctionNode case f" << std::endl;);
                return std::make_pair(Nclose, false);
            }
            else
            {
                VerboseJunctions(4, "JunctionNode case g (returning nullptr)" << std::endl;);
                return std::make_pair(nullptr, false);
            }
        }
    }

    //New Implemnetation with update of the boundary nodes
    template <typename DislocationNetworkType>
    size_t DislocationJunctionFormation<DislocationNetworkType>::contractJunctions(const std::deque<IntersectionTypeContainerType>& intersectionContainer)
    {
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"        : "<<std::flush;
        
        DN.danglingBoundaryLoopNodes.clear();
        size_t nContracted=0;
        size_t nDeferred=0;
        std::set<size_t> changedNodes;

        for (const auto& intersectionByThreadContainer : intersectionContainer)
        {
            for (const auto& intersection : intersectionByThreadContainer)
            {
                const KeyType& key1(std::get<0>(intersection));
                const KeyType& key2(std::get<1>(intersection));
                
                const auto L1(DN.networkLinks().get(key1));
                const auto L2(DN.networkLinks().get(key2));
                
                if(L1 && L2) // Links exist
                {
                    // Inspect whole incident contours before querying derived boundary
                    // topology. A shared node makes two otherwise distinct loops
                    // dependent. IDs also survive loop splitting during contraction.
                    std::set<size_t> footprint;
                    std::set<size_t> visitedLoops;
                    for(const auto& node : {L1->source,L1->sink,L2->source,L2->sink})
                    {
                        for(const auto* occurrence : node->loopNodes())
                        {
                            if(visitedLoops.insert(occurrence->loop()->sID).second)
                            {
                                for(const auto* link : occurrence->loop()->linkSequence())
                                {
                                    footprint.insert(link->source->networkNode->sID);
                                }
                            }
                        }
                    }
                    const bool dependent=std::any_of(footprint.begin(),footprint.end(),[&](const size_t id)
                    {
                        return changedNodes.count(id)!=0;
                    });
                    if(dependent)
                    {
                        ++nDeferred;
                        continue;
                    }
                    // Recheck only queued pairs, not the global sweep. Earlier
                    // contractions can move endpoints or change link membership
                    // without changing the two endpoint IDs.
                    std::deque<IntersectionTypeContainerType> refreshed(1);
                    findPairIntersections(refreshed,L1.get(),L2.get());
                    if(refreshed.front().empty())
                    {
                        continue;
                    }
                    // Keep the two endpoint-versus-segment events distinct.
                    // If connectivity changed, use the newly valid pair event.
                    auto refreshedEvent=refreshed.front().begin();
                    for(auto it=refreshed.front().begin();it!=refreshed.front().end();++it)
                    {
                        if(std::get<3>(*it)==std::get<3>(intersection))
                        {
                            refreshedEvent=it;
                            break;
                        }
                    }
                    const auto ssd=std::get<2>(*refreshedEvent);
                    const double t=ssd.t;
                    const double u=ssd.u;
                    VerboseJunctions(1,"forming junction "<<key1.first<<"->"<<key1.second<<", "
                                     /*                   */ <<key2.first<<"->"<<key2.second<<", "
                                     /*                   */ <<"dMin="<<ssd.dMin<<", "
                                     /*                   */ <<"@ ("<<t<<","<<u<<"), "<<std::flush);
                    VerboseJunctions(4, "JunctionNode for segment: " << key1.first << "->" << key1.second << " @ " << t << std::endl;);
                    std::pair<std::shared_ptr<NetworkNodeType>,bool> Ni(junctionNode(t,ssd.x0,L1));
                    VerboseJunctions(4, "JunctionNode for segment: " << key2.first << "->" << key2.second << " @ " << u << std::endl;);
                    std::pair<std::shared_ptr<NetworkNodeType>,bool> Nj(junctionNode(u,ssd.x1,L2));
                    
                    if (Ni.first && Nj.first)
                    {
                        
                        VerboseJunctions(1, "contracting " << Ni.first->sID << " " << Nj.first->sID << std::endl;);
                        
                        if (Ni.first->sID != Nj.first->sID)
                        {
                            const bool success = DN.nodeContractor.contract(Ni.first,Nj.first,false);
                            VerboseJunctions(1, "success= " << success << std::endl;);
                            
                            nContracted += success;
                            if (!success)
                            {
                                if (Ni.second)
                                {
                                    DN.removeNetworkNode(Ni.first->sID);
                                }
                                if (Nj.second)
                                {
                                    DN.removeNetworkNode(Nj.first->sID);
                                }
                            }
                            else
                            {
                                changedNodes.insert(footprint.begin(),footprint.end());
                                changedNodes.insert(Ni.first->sID);
                                changedNodes.insert(Nj.first->sID);
                            }
                        }
                    }
                    else
                    {
                        //Either Ni or Nj might have been created...
                        //Get rid of that
                        if (Ni.first)
                        {
                            if (Ni.second)
                            {
                                DN.removeNetworkNode(Ni.first->sID);
                            }
                        }
                        
                        if (Nj.first)
                        {
                            if (Nj.second)
                            {
                                DN.removeNetworkNode(Nj.first->sID);
                            }
                        }
                        
                    }
                    
                }
            }
        }
        
        if(nContracted || !DN.danglingBoundaryLoopNodes.empty())
        {
            DN.updateBoundaryNodes();
        }
        std::cout<<" ("<<nContracted<<" contracted, "<<nDeferred<<" deferred)"<<magentaColor
                 <<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()
                 <<" sec]"<<defaultColor<<std::endl;

        return nContracted;
    }

    //This version also compares the network links burgers across the boundary
    template <typename DislocationNetworkType>
    void DislocationJunctionFormation<DislocationNetworkType>::glissileJunctions(const double &dx)
    {
        const auto t0 = std::chrono::system_clock::now();
        std::cout << "        Forming Glissile Junctions: " << std::flush;
        
        std::deque<std::tuple<std::shared_ptr<NetworkNodeType>, std::shared_ptr<NetworkNodeType>, size_t, size_t>> glissDeq;
                
        for (const auto &linkIter : DN.networkLinks())
        {
            const auto link(linkIter.second.lock());
            
            //    if (link->isSessile() && link->loopLinks().size() > 1) // a junction
            if (   link->loopLinks().size() > 1 // a junction
                && link->chordLength() > dx
                && link->isSessile()
                )
            {
                                
                const VectorDim unitChord(link->chord() / link->chordLength());
                
                if (!link->isGrainBoundarySegment() && !link->isBoundarySegment())
                {
                    
                    VerboseJunctions(2, "glissele junction, segment " << link->tag() << std::endl;);
                    
                    for (const auto &gr : link->grains())
                    {
                        for (size_t s = 0; s < gr->slipSystems().size(); ++s)
                        {
                            const auto &slipSystem(gr->slipSystems().at(s));
                            if ((slipSystem->s.cartesian() - link->burgers()).norm() < FLT_EPSILON && fabs(slipSystem->n.cartesian().normalized().dot(unitChord)) < FLT_EPSILON)
                            {
                                VerboseJunctions(3, "glissDeq, emplacing" << std::endl;);
                                
                                glissDeq.emplace_back(link->source, link->sink, gr->grainID, s);
                            }
                        }
                    }
                }
            }
        }
        
        size_t formedJunctions=0;
        for(const auto& [source,sink,grainID,slipID] : glissDeq)
        {
            const auto link=DN.networkLinks().get(std::make_pair(source->sID,sink->sID));
            const auto& slip=DN.ddBase.poly.grain(grainID)->slipSystems().at(slipID);
            // Earlier seeds may already have consumed this periodic continuation.
            if(!link || link->hasZeroBurgers() || !link->isSessile() ||
               link->isBoundarySegment() || link->isGrainBoundarySegment() ||
               (link->burgers()-slip->s.cartesian()).norm()>FLT_EPSILON)
            {
                continue;
            }

            // Orient each contributing loop's walk from network source to sink.
            // Only artificial periodic subdivisions are traversed: physical
            // faces and original loop nodes terminate the continuation.
            std::deque<const LoopNodeType*> path;
            VectorDim originShift=VectorDim::Zero();
            bool valid=true;
            for(const auto* occurrence : link->loopLinks())
            {
                const bool forward=occurrence->source->networkNode==source;
                const auto* a=(forward? occurrence->source : occurrence->sink).get();
                const auto* b=(forward? occurrence->sink : occurrence->source).get();
                std::deque<const LoopNodeType*> candidate{a,b};
                std::set<const LoopNodeType*> visited{a,b};
                for(const bool left : {true,false})
                {
                    const auto* current=left? candidate.front() : candidate.back();
                    while(current->isPeriodicSubdivision())
                    {
                        const bool next=left? !forward : forward;
                        const auto* adjacent=next? current->next.first : current->prev.first;
                        if(!adjacent || !visited.insert(adjacent).second)
                        {
                            valid=false;
                            break;
                        }
                        if(left)
                        {
                            candidate.push_front(adjacent);
                        }
                        else
                        {
                            candidate.push_back(adjacent);
                        }
                        current=adjacent;
                    }
                    if(!valid)
                    {
                        break;
                    }
                }
                if(!valid)
                {
                    break;
                }
                const VectorDim shift=a->periodicPlanePatch()->shift;
                if(path.empty())
                {
                    path=candidate;
                    originShift=shift;
                }
                else
                {
                    // Sharing must agree through the whole continuation, not
                    // just at the original pair of network nodes. Do not choose
                    // an unrelated loop visiting a shared endpoint.
                    if(path.size()!=candidate.size())
                    {
                        valid=false;
                        break;
                    }
                    for(size_t k=0;k<path.size();++k)
                    {
                        if(path[k]->networkNode!=candidate[k]->networkNode ||
                           (path[k]->periodicPlanePatch()->shift-originShift-
                            candidate[k]->periodicPlanePatch()->shift+shift).norm()>FLT_EPSILON)
                        {
                            valid=false;
                            break;
                        }
                    }
                }
                if(!valid)
                {
                    break;
                }
            }
            if(!valid || path.size()<2)
            {
                continue;
            }

            const VectorDim midpoint=0.5*(source->get_P()+sink->get_P());
            const auto glidePlane=DN.ddBase.glidePlaneFactory.getFromKey(GlidePlaneKey<dim>(midpoint,slip->n));
            const auto periodicPlane=DN.ddBase.periodicGlidePlaneFactory.get(glidePlane->key);
            using EdgePair=std::pair<std::shared_ptr<PeriodicPlaneEdge<dim>>,std::shared_ptr<PeriodicPlaneEdge<dim>>>;
            struct SeedPoint
            {
                std::shared_ptr<NetworkNodeType> node;
                VectorDim position;
                std::shared_ptr<PeriodicPlanePatch<dim>> patch;
                EdgePair edges;
            };
            std::vector<SeedPoint> prepared;
            std::set<const NetworkLinkType*> physicalLinks;
            size_t midpointAfter=path.size();
            for(size_t k=0;k<path.size();++k)
            {
                const auto* node=path[k];
                const VectorDim shift=node->periodicPlanePatch()->shift-originShift;
                const VectorDim position=node->networkNode->get_P()-shift;
                if(std::abs((position-glidePlane->P).dot(glidePlane->unitNormal))>FLT_EPSILON)
                {
                    valid=false;
                    break;
                }
                const auto patch=periodicPlane->getPatch(shift);
                EdgePair edges;
                for(const auto& edge : patch->edges())
                {
                    // Geometry assigns metadata on the target plane; it never
                    // chooses the continuation. Corners retain both patch edges.
                    if(edge->meshIntersection->contains(node->networkNode->get_P()))
                    {
                        if(!edges.first)
                        {
                            edges.first=edge;
                        }
                        else if(!edges.second)
                        {
                            edges.second=edge;
                        }
                        else
                        {
                            valid=false;
                        }
                    }
                }
                if(node->isPeriodicSubdivision() && !edges.first)
                {
                    valid=false;
                }
                prepared.push_back({node->networkNode,position,patch,edges});
                if(k+1<path.size())
                {
                    const auto* next=path[k+1];
                    const auto* loopLink=node->next.first==next? node->next.second : node->prev.second;
                    const auto segment=loopLink->networkLink();
                    if(segment)
                    {
                        const VectorDim orientedB=(segment->source==node->networkNode? 1.0 : -1.0)*segment->burgers();
                        if(segment->isBoundarySegment() || segment->isGrainBoundarySegment() ||
                           (orientedB-link->burgers()).norm()>FLT_EPSILON ||
                           !physicalLinks.insert(segment.get()).second)
                        {
                            valid=false;
                        }
                        if(segment==link)
                        {
                            midpointAfter=k;
                        }
                    }
                }
                if(!valid)
                {
                    break;
                }
            }
            if(!valid || midpointAfter==path.size())
            {
                continue;
            }

            // Commit only after all contributors, translations and target-plane
            // metadata have been checked. The return side shares the original
            // network nodes, cancelling the junction on every periodic piece.
            const auto seed=DN.loops().create(slip->s.cartesian(),glidePlane);
            std::vector<std::shared_ptr<LoopNodeType>> nodes;
            for(auto iter=prepared.rbegin();iter!=prepared.rend();++iter)
            {
                nodes.push_back(DN.loopNodes().create(seed,iter->node,iter->position,iter->patch,iter->edges));
            }
            for(size_t k=0;k+1<prepared.size();++k)
            {
                const auto& point=prepared[k];
                if(k>0)
                {
                    const auto freeNode=DN.networkNodes().clone(point.node->sID);
                    nodes.push_back(DN.loopNodes().create(seed,freeNode,point.position,point.patch,point.edges));
                }
                if(k==midpointAfter)
                {
                    const auto freeNode=DN.networkNodes().create(midpoint);
                    nodes.push_back(DN.loopNodes().create(seed,freeNode,midpoint,periodicPlane->getPatch(VectorDim::Zero()),EdgePair()));
                }
            }
            DN.insertLoop(seed,nodes);
            ++formedJunctions;
        }
        std::cout << "(" << formedJunctions << " junctions)" << magentaColor << " [" << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]" << defaultColor << std::endl;
    }

    template <typename DislocationNetworkType>
    DislocationJunctionFormation<DislocationNetworkType>::DislocationJunctionFormation(DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    /* init */,maxJunctionIterations(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("maxJunctionIterations",true))
    /* init */,verboseJunctions(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("verboseJunctions",true))
    /* init */,infiniteLineLength(10000.0)
    {
        
    }

    template <typename DislocationNetworkType>
    void DislocationJunctionFormation<DislocationNetworkType>::initFromFile(const std::string& fileName)
    {
        collisionTol=TextFileParser(fileName).readScalar<double>("collisionTol",true);
    }

    template <typename DislocationNetworkType>
    void DislocationJunctionFormation<DislocationNetworkType>::formJunctions(const double& dx)
    {
    #ifdef _OPENMP
        // const size_t nThreads = omp_get_max_threads();
        const size_t nThreads = 1;
    #else
        const size_t nThreads = 1;
    #endif
        
        
        
        size_t nContracted=1;
        size_t iterations=0;
        while(nContracted && iterations<maxJunctionIterations)
        {
            std::deque<IntersectionTypeContainerType> intersectionContainer;
            intersectionContainer.resize(nThreads);
            findIntersections(intersectionContainer,nThreads);
            
            nContracted=contractJunctions(intersectionContainer);
            glissileJunctions(dx);
            iterations++;
        }
        
//        breakJunctions();
        
    }

    template <typename DislocationNetworkType>
    void DislocationJunctionFormation<DislocationNetworkType>::breakJunctions()
    {
        for (const auto &nodeIter : DN.networkNodes())
        {
            const auto netNode(nodeIter.second.lock()); // shared_ptr to DislocationNode
            
            //    if (link->isSessile() && link->loopLinks().size() > 1) // a junction
            if (   netNode->loopNodes().size()==2 // a junction
                && netNode->neighbors().size()==4
                )
            {
                std::cout<<"Detected 4-point junction at node "<<netNode->sID<<std::endl;
//                (*netNode->loopNodes().rbegin())->reset(netNode->clone());
            }
        }
    }


    template <typename DislocationNetworkType>
    double DislocationJunctionFormation<DislocationNetworkType>::collisionTol=10.0;

    template class DislocationJunctionFormation<DislocationNetwork<3>>;

}
#endif
