/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationCrossSlip_cpp_
#define model_DislocationCrossSlip_cpp_

#include <numbers>
#include <memory>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <DislocationCrossSlip.h>
#include <CrossSlipModels.h>
#include <OneSidedMotion.h>

namespace model
{
    template <typename DislocationNetworkType>
    std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>> DislocationCrossSlip<DislocationNetworkType>::getModel(const PolycrystallineMaterialBase& material,const DDtraitsIO& traitsIO)
    {
        const int crossSlipModel(TextFileParser(traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("crossSlipModel",true));
        switch (crossSlipModel)
        {
            case 0:
            {// no cross-slip
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                break;
            }
                
            case 1:
            {// deterministic cross-slip model based on max glide PK force
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(new AthermalCrossSlipModel<DislocationNetworkType>(traitsIO));
                break;
            }
                
            case 2:
            {
                if(material.crystalStructure=="HEX")
                {
                    return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(new EscaigCrossSlipModelHEX<DislocationNetworkType>(material,traitsIO));
                }
                else
                {
                    throw std::runtime_error("Unknown cross-slip model "+std::to_string(crossSlipModel)+" for crystal structure "+material.crystalStructure);
                    return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                }
                break;
            }
                
            default:
            {
                throw std::runtime_error("Unknown cross-slip model "+std::to_string(crossSlipModel));
                return std::shared_ptr<BaseCrossSlipModel<DislocationNetworkType>>(nullptr);
                break;
            }
        }
    }

    template <typename DislocationNetworkType>
    DislocationCrossSlip<DislocationNetworkType>::DislocationCrossSlip(DislocationNetworkType& DN_in) :
    /* init */ DN(DN_in)
    /* init */,verboseCrossSlip(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("verboseCrossSlip",true))
    {
    }

    template <typename DislocationNetworkType>
    void  DislocationCrossSlip<DislocationNetworkType>::findCrossSlipSegments(const double dt)
    {
        
        if(DN.crossSlipModel)
        {
            std::cout<<"CrossSlip: finding branches "<<std::flush;
            const auto t0= std::chrono::steady_clock::now();
            csNodes.clear();
            for(const auto& loopIter : DN.loops())
            {
                const auto loop(loopIter.second.lock());
                loop->crossSlipBranches(csNodes,dt);
            }
            std::cout<<csNodes.size()<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::steady_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
    }

    template <typename DislocationNetworkType>
    void DislocationCrossSlip<DislocationNetworkType>::execute()
    {
        if(!DN.crossSlipModel || !DN.glideSolver || DN.isClimbStep())
        {
            return;
        }
        auto phaseStart=std::chrono::steady_clock::now();
        const auto reportPhase=[&](const char* phase)
        {
            const auto now=std::chrono::steady_clock::now();
            std::cout<<"CrossSlip: "<<phase<<" ["<<std::chrono::duration<double>(now-phaseStart).count()<<" sec]"<<std::endl;
            phaseStart=std::chrono::steady_clock::now();
        };
        CrossSlipBranchContainerType alignedBranches;
        std::set<size_t> reservedNodes;
        for(const auto& branch : csNodes)
        {
            if(branch.first.size()<(branch.thermal? 2 : 3))
            {
                continue;
            }
            const auto loop=branch.first.back()->loop();
            if(!loop->glidePlane)
            {
                continue;
            }
            const auto& system=loop->grain.slipSystems().at(branch.second);
            VectorDim midpoint=VectorDim::Zero();
            double branchLength=0.0;
            for(size_t k=1;k<branch.first.size();++k)
            {
                const auto& a=branch.first[k-1]->get_P();
                const auto& b=branch.first[k]->get_P();
                const double length=(b-a).norm();
                midpoint+=0.5*length*(a+b);
                branchLength+=length;
            }
            if(branchLength<=FLT_EPSILON)
            {
                continue;
            }
            midpoint/=branchLength;
            const long spacing=std::max(1L,static_cast<long>(20.0/system->n.planeSpacing()));
            const long height=std::lround(double(system->n.closestPlaneIndexOfPoint(midpoint))/spacing)*spacing;
            const VectorDim planePoint=height*system->n.planeSpacing()*system->unitNormal;
            PlanePlaneIntersection<dim> intersection(loop->glidePlane->P,loop->glidePlane->unitNormal,
                                                     planePoint,system->unitNormal);
            // Prepare each NETWORK displacement once. No positions or links
            // change until every occurrence and every affected loop is checked.
            std::map<size_t,VectorDim> targets;
            bool valid=true;
            std::set<LoopType*> affected;
            for(const auto& node : branch.first)
            {
                const auto net=node->networkNode;
                const VectorDim delta=intersection.P+(node->get_P()-intersection.P).dot(intersection.d)*intersection.d-node->get_P();
                const VectorDim target=net->get_P()+delta;
                const auto prior=targets.find(net->sID);
                if(reservedNodes.count(net->sID) || (prior!=targets.end() && (prior->second-target).norm()>FLT_EPSILON))
                {
                    valid=false;
                    break;
                }
                typename NetworkNodeType::VectorOfNormalsType normals;
                for(const auto* plane : net->glidePlanes())
                {
                    normals.push_back(plane->unitNormal);
                }
                for(int row=0;row<DN.nodalVelocityConstraints.rows();++row)
                {
                    normals.push_back(DN.nodalVelocityConstraints.row(row));
                }
                if(!target.allFinite() || !net->isMovableTo(target) ||
                   (oneSidedMotion<dim>(delta,normals,net->physicalBoundaryNormals())-delta).norm()>FLT_EPSILON)
                {
                    valid=false;
                    break;
                }
                for(const auto* occurrence : net->loopNodes())
                {
                    affected.insert(occurrence->loop().get());
                }
                targets.emplace(net->sID,target);
            }
            for(const auto* affectedLoop : affected)
            {
                if(!valid)
                {
                    break;
                }
                std::vector<std::pair<VectorDim,const LoopNodeType* const>> polygon;
                for(const auto* link : affectedLoop->linkSequence())
                {
                    const auto& node=link->source;
                    if(!node->isPeriodicSubdivision())
                    {
                        VectorDim position=node->get_P();
                        const auto target=targets.find(node->networkNode->sID);
                        if(target!=targets.end())
                        {
                            position+=target->second-node->networkNode->get_P();
                        }
                        polygon.emplace_back(position,node.get());
                    }
                }
                const auto contour=affectedLoop->periodicGlidePlane->polygonPatchIntersection(polygon,true);
                std::vector<double> coverage(polygon.size(),0.0);
                for(size_t k=0;k<contour.size();++k)
                {
                    const auto& metadata=std::get<6>(contour[k]);
                    if(metadata.incomingOriginalSegment)
                    {
                        coverage.at(metadata.incomingOriginalSegment->sourceIndex)+=
                            (std::get<0>(contour[k])-std::get<0>(contour[(k+contour.size()-1)%contour.size()])).norm();
                    }
                }
                if(contour.empty())
                {
                    valid=false;
                }
                for(size_t k=0;k<polygon.size();++k)
                {
                    const double length=(polygon[(k+1)%polygon.size()].first-polygon[k].first).norm();
                    if(std::abs(length-coverage[k])>1.e-8*std::max(1.0,length))
                    {
                        valid=false;
                    }
                }
            }
            if(branch.thermal)
            {
                double alignedLength=0.0;
                for(size_t k=1;k<branch.first.size();++k)
                {
                    alignedLength+=std::abs((branch.first[k]->get_P()-branch.first[k-1]->get_P()).dot(intersection.d));
                }
                if(alignedLength+FLT_EPSILON<branch.thermal->length)
                {
                    valid=false;
                }
            }
            if(!valid)
            {
                continue;
            }
            for(const auto& [id,target] : targets)
            {
                const auto net=DN.networkNodes().get(id);
                const VectorDim delta=target-net->get_P();
                for(auto* occurrence : net->loopNodes())
                {
                    static_cast<typename LoopNodeType::SplineNodeType&>(*occurrence).set_P(occurrence->get_P()+delta);
                }
                // Network positions are committed only by the boundary walk,
                // which maps same-grain periodic exits back into the RVE.
                reservedNodes.insert(id);
            }
            alignedBranches.push_back(branch);
        }
        csNodes.clear();
        reportPhase("alignment validation and commit");
        std::cout<<"CrossSlip: "<<alignedBranches.size()<<" aligned branches"<<std::endl;
        if(alignedBranches.empty())
        {
            return;
        }
        DN.updateBoundaryNodes();
        reportPhase("boundary update");
        std::vector<std::shared_ptr<LoopType>> newLoops;
        for(auto& branch : alignedBranches)
        {
            // Reconciliation can replace a physical-face occurrence while
            // preserving its network node. Resolve live occurrences before
            // following prev/next connectivity.
            for(auto& anchor : branch.first)
            {
                const auto parent=anchor->loop();
                std::shared_ptr<LoopNodeType> live;
                for(const auto* link : parent->linkSequence())
                {
                    if(link->source->networkNode==anchor->networkNode &&
                       (link->source->get_P()-anchor->get_P()).norm()<FLT_EPSILON)
                    {
                        live=link->source;
                        break;
                    }
                }
                if(!live)
                {
                    throw std::runtime_error("Cross-slip alignment lost a branch anchor");
                }
                anchor=live;
            }
            const auto parent=branch.first.back()->loop();
            if(branch.thermal)
            {
                // Split only after successful transactional alignment. Walk live
                // loop connectivity, including periodic seam connectors; they
                // carry zero physical length and retain their original patches.
                const auto first=branch.first.front();
                const auto last=branch.first.back();
                double total=0.0;
                auto* cursor=first.get();
                size_t walked=0;
                while(cursor!=last.get())
                {
                    if(!cursor || !cursor->next.first || ++walked>parent->loopLinks().size())
                    {
                        throw std::runtime_error("Thermal branch lost connectivity");
                    }
                    total+=(cursor->next.first->get_P()-cursor->get_P()).norm();
                    cursor=cursor->next.first;
                }
                const double length=std::min(branch.thermal->length,total);
                const double center=std::clamp(branch.thermal->centerFraction*total,0.5*length,total-0.5*length);
                std::deque<std::shared_ptr<LoopNodeType>> nucleus;
                for(const double distance : {center-0.5*length,center,center+0.5*length})
                {
                    double offset=0.0;
                    cursor=first.get();
                    walked=0;
                    std::shared_ptr<LoopNodeType> selected;
                    while(cursor!=last.get())
                    {
                        if(!cursor || !cursor->next.first || ++walked>parent->loopLinks().size())
                        {
                            throw std::runtime_error("Thermal nucleus lost connectivity");
                        }
                        const auto* link=cursor->next.second;
                        const auto source=link->source;
                        const auto sink=link->sink;
                        const double chord=(sink->get_P()-source->get_P()).norm();
                        if(distance<=offset+FLT_EPSILON)
                        {
                            selected=source;
                            break;
                        }
                        if(chord>FLT_EPSILON && distance<offset+chord-FLT_EPSILON)
                        {
                            const double u=(distance-offset)/chord;
                            const auto a=source->networkNode;
                            const auto b=sink->networkNode;
                            const VectorDim position=(1.0-u)*a->get_P()+u*b->get_P();
                            const auto net=DN.networkNodes().create(position,(1.0-u)*a->get_V()+u*b->get_V(),
                                (1.0-u)*a->climbVelocityScalar+u*b->climbVelocityScalar,
                                (1.0-u)*a->velocityReduction()+u*b->velocityReduction());
                            selected=DN.loopNodes().create(parent,net,link);
                            DN.expandLoopLink(*link,selected);
                            break;
                        }
                        offset+=chord;
                        cursor=sink.get();
                    }
                    nucleus.push_back(selected? selected : last);
                }
                branch.first=std::move(nucleus);
            }
            const auto& system=parent->grain.slipSystems().at(branch.second);
            const auto anchor=branch.first.back();
            const VectorDim originShift=anchor->networkNode->get_P()-anchor->get_P();
            const auto plane=DN.ddBase.glidePlaneFactory.getFromKey(
                GlidePlaneKey<dim>(anchor->networkNode->get_P(),system->n));
            const auto seed=DN.loops().create(system->s.cartesian(),plane);
            std::deque<std::shared_ptr<LoopNodeType>> nodes;
            auto* current=anchor.get();
            size_t visited=0;
            while(true)
            {
                if(!current || ++visited>parent->loopLinks().size())
                {
                    throw std::runtime_error("Cross-slip branch lost its loop connectivity");
                }
                // Transfer the actual patch translation, not an alternating
                // edge counter. Physical faces never imply a patch transition.
                const VectorDim shift=current->periodicPlanePatch()->shift-originShift;
                const auto patch=seed->periodicGlidePlane->getPatch(shift);
                std::pair<std::shared_ptr<PeriodicPlaneEdge<dim>>,std::shared_ptr<PeriodicPlaneEdge<dim>>> edges;
                if(current->periodicPlaneEdge.first)
                {
                    edges.first=patch->edgeOnFaces(current->periodicPlaneEdge.first->meshIntersection->faces);
                    if(!edges.first)
                    {
                        throw std::runtime_error("Cross-slip patch missing first boundary face");
                    }
                }
                if(current->periodicPlaneEdge.second)
                {
                    edges.second=patch->edgeOnFaces(current->periodicPlaneEdge.second->meshIntersection->faces);
                    if(!edges.second)
                    {
                        throw std::runtime_error("Cross-slip patch missing second boundary face");
                    }
                }
                const VectorDim position=current->networkNode->get_P()-shift;
                nodes.emplace_back(DN.loopNodes().create(seed,current->networkNode,position,patch,edges));
                if(current==branch.first.front().get())
                {
                    break;
                }
                if(current!=anchor.get())
                {
                    nodes.emplace_front(DN.loopNodes().create(seed,DN.networkNodes().clone(current->networkNode->sID),position,patch,edges));
                }
                current=current->prev.first;
            }
            DN.insertLoop(seed,nodes);
            newLoops.push_back(seed);
        }
        reportPhase("seed insertion");
        // Seeds deliberately have zero area. Do not clip them before the
        // first physical motion; the ordinary moveNodes walk handles that.
        for(const auto& entry : DN.loops())
        {
            entry.second.lock()->updateGeometry();
        }
        reportPhase("loop geometry");
        // Recompute all fields: alignment also changed the parent/junctions.
        for(const auto& entry : DN.networkLinks())
        {
            entry.second.lock()->createQuadraturePoints(false);
        }
        reportPhase("quadrature creation");
        {
            // Topology is fixed here. Release all shared references before
            // returning to any operation that can change network connectivity.
            std::vector<std::shared_ptr<NetworkLinkType>> links;
            links.reserve(DN.networkLinks().size());
            for(const auto& entry : DN.networkLinks())
            {
                links.push_back(entry.second.lock());
            }
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(size_t k=0;k<links.size();++k)
            {
                links[k]->updateQuadraturePoints(false);
            }
        }
        reportPhase("field refresh");
        if(DN.glideSolver)
        {
            const Eigen::VectorXd velocities=DN.glideSolver->getNodeVelocities();
            if(velocities.size()!=static_cast<Eigen::Index>(dim*DN.networkNodes().size()))
            {
                throw std::runtime_error("Cross-slip velocity solve returned wrong size");
            }
            size_t k=0;
            for(const auto& entry : DN.networkNodes())
            {
                entry.second.lock()->set_V(velocities.template segment<dim>(dim*k),false);
                ++k;
            }
        }
        reportPhase("velocity solve and storage");
        for(const auto& entry : DN.loops())
        {
            entry.second.lock()->updateRates();
        }
        reportPhase("rate update");
        csNodes.clear();
    }

    template class DislocationCrossSlip<DislocationNetwork<3>>;
}
#endif
