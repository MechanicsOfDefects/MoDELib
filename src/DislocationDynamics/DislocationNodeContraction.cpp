/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeContraction_cpp_
#define model_DislocationNodeContraction_cpp_

#include <memory>
#include <tuple>
#include <Eigen/Dense>
#include <TypeTraits.h>
#include <GlidePlane.h>
//#include <BoundingMeshSegments.h>
#include <Grain.h>
#include <FiniteLineSegment.h>
#include <PlanesIntersection.h>
#include <DislocationNetwork.h>
#include <DislocationNodeContraction.h>
#include <OneSidedMotion.h>
#include <PeriodicContourTopology.h>

#ifndef NDEBUG
#define VerboseNodeContraction(N,x) if(verboseNodeContraction>=N){std::cout<<x;}
#else
#define VerboseNodeContraction(N,x)
#endif

namespace model
{

    template <int dim>
    DislocationNodeContraction<dim>::DislocationNodeContraction(DislocationNetwork<dim>& DN_in) :
    /* init */ DN(DN_in)
    /* init */,verboseNodeContraction(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("verboseNodeContraction",true))
    {
        
    }

    template <int dim>
    std::vector<std::set<std::shared_ptr<typename DislocationNodeContraction<dim>::NetworkNodeType>>> DislocationNodeContraction<dim>::contractBoundary(std::shared_ptr<NetworkNodeType> nA,
                                                                                                                                                                  std::shared_ptr<NetworkNodeType> nB)
    {
        using NodeSet=std::set<std::shared_ptr<NetworkNodeType>>;
        std::vector<NodeSet> groups;
        const auto collect=[&](const auto& chainA,const auto& chainB)
        {
            for(const auto* a : chainA)
            {
                for(const auto* b : chainB)
                {
                    // Only coincident subdivisions of the same periodic branch
                    // may be shared. Physical faces terminate these chains.
                    if(a->isPeriodicSubdivision() && b->isPeriodicSubdivision() &&
                       a->networkNode!=b->networkNode &&
                       a->networkNode->meshFaces()==b->networkNode->meshFaces() &&
                       (a->networkNode->get_P()-b->networkNode->get_P()).norm()<FLT_EPSILON)
                    {
                        NodeSet merged{a->networkNode,b->networkNode};
                        for(auto it=groups.begin();it!=groups.end();)
                        {
                            bool overlaps=false;
                            for(const auto& node : *it)
                            {
                                if(merged.count(node))
                                {
                                    overlaps=true;
                                    break;
                                }
                            }
                            if(overlaps)
                            {
                                merged.insert(it->begin(),it->end());
                                it=groups.erase(it);
                            }
                            else
                            {
                                ++it;
                            }
                        }
                        groups.push_back(std::move(merged));
                    }
                }
            }
        };
        for(const auto* a : nA->loopNodes())
        {
            if(a->periodicPlaneEdge.first)
            {
                continue;
            }
            for(const auto* b : nB->loopNodes())
            {
                if(b->periodicPlaneEdge.first)
                {
                    continue;
                }
                const auto ap=a->periodicPrev();
                const auto an=a->periodicNext();
                const auto bp=b->periodicPrev();
                const auto bn=b->periodicNext();
                if(ap && bp && ap->networkNode==bp->networkNode)
                {
                    collect(a->boundaryPrev(),b->boundaryPrev());
                }
                if(an && bn && an->networkNode==bn->networkNode)
                {
                    collect(a->boundaryNext(),b->boundaryNext());
                }
                if(ap && bn && ap->networkNode==bn->networkNode)
                {
                    collect(a->boundaryPrev(),b->boundaryNext());
                }
                if(an && bp && an->networkNode==bp->networkNode)
                {
                    collect(a->boundaryNext(),b->boundaryPrev());
                }
            }
        }
        return groups;
    }

    template <int dim>
    bool DislocationNodeContraction<dim>::canReach(const std::shared_ptr<NetworkNodeType>& node,const VectorDim& target) const
    {
        if(!target.allFinite() || !node->isMovableTo(target))
        {
            return false;
        }
        const VectorDim requested=target-node->get_P();
        typename NetworkNodeType::VectorOfNormalsType normals;
        for(const auto* plane : node->glidePlanes())
        {
            normals.push_back(plane->unitNormal);
        }
        for(int row=0;row<DN.nodalVelocityConstraints.rows();++row)
        {
            normals.push_back(DN.nodalVelocityConstraints.row(row));
        }
        const VectorDim actual=oneSidedMotion<dim>(requested,normals,node->physicalBoundaryNormals());
        if((actual-requested).norm()>FLT_EPSILON)
        {
            return false;
        }
        for(const auto* ln : node->loopNodes())
        {
            if(ln->periodicPlaneEdge.first && requested.norm()>FLT_EPSILON)
            {
                return false;
            }
            // Keep each representation in its current convex host-grain patch.
            // Contraction must not silently trim away one of its endpoints.
            if(!DN.ddBase.mesh.searchRegion(ln->loop()->grain.grainID,target).first)
            {
                return false;
            }
        }
        return true;
    }

    template <int dim>
    bool DislocationNodeContraction<dim>::commitContraction(std::shared_ptr<NetworkNodeType> a,
        std::shared_ptr<NetworkNodeType> b,const VectorDim& target,const bool reconcileBoundary)
    {
        if(!a || !b || a==b || a->loopNodes().empty() || b->loopNodes().empty() ||
           !b->isContractableTo(a) || !canReach(a,target) || !canReach(b,target))
        {
            return false;
        }
        // All ordinary rejection paths are above this point. Commit the exact
        // common point without running boundary reconciliation between moves.
        using LoopNodeType=typename DislocationNetwork<dim>::LoopNodeType;
        using NodePtr=std::shared_ptr<LoopNodeType>;
        struct Subdivision
        {
            NodePtr node,source,sink;
            double parameter;
        };
        std::vector<Subdivision> subdivisions;
        std::set<size_t> visitedLoops;
        const double tolerance=1.e-8*(DN.ddBase.mesh.xMax()-DN.ddBase.mesh.xMin()).norm();
        for(const auto& endpoint : {a,b})
        {
            for(const auto* occurrence : endpoint->loopNodes())
            {
                const auto loop=occurrence->loop();
                if(!visitedLoops.insert(loop->sID).second)
                {
                    continue;
                }
                std::vector<NodePtr> nodes;
                for(const auto* link : loop->linkSequence())
                {
                    nodes.push_back(link->source);
                }
                retainPeriodicContourBends(nodes,tolerance);
                std::vector<size_t> anchors;
                for(size_t i=0;i<nodes.size();++i)
                {
                    if(!nodes[i]->isPeriodicSubdivision())
                    {
                        anchors.push_back(i);
                    }
                }
                for(size_t k=0;k<anchors.size();++k)
                {
                    const size_t begin=anchors[k];
                    const size_t end=anchors[(k+1)%anchors.size()];
                    const VectorDim chord=nodes[end]->get_P()-nodes[begin]->get_P();
                    std::vector<Subdivision> edgeSubdivisions;
                    bool straight=true;
                    double previous=0;
                    const double length=chord.norm();
                    for(size_t i=(begin+1)%nodes.size();i!=end;i=(i+1)%nodes.size())
                    {
                        const VectorDim offset=nodes[i]->get_P()-nodes[begin]->get_P();
                        const double parameter=chord.squaredNorm()>0 ? offset.dot(chord)/chord.squaredNorm() : 0;
                        if((offset-parameter*chord).norm()>tolerance ||
                           (parameter-previous)*length < -tolerance || (parameter-1)*length>tolerance)
                        {
                            straight=false;
                        }
                        previous=parameter;
                        edgeSubdivisions.push_back({nodes[i],nodes[begin],nodes[end],parameter});
                    }
                    if(straight)
                    {
                        subdivisions.insert(subdivisions.end(),edgeSubdivisions.begin(),edgeSubdivisions.end());
                    }
                    else
                    {
                        // Local near-collinearity is not enough for a long
                        // chain: retain it if the complete path is not straight.
                        for(const auto& subdivision : edgeSubdivisions)
                        {
                            subdivision.node->periodicContourVertex=true;
                        }
                    }
                }
            }
        }
        for(const auto& node : {a,b})
        {
            for(auto* ln : node->loopNodes())
            {
                const VectorDim shift=ln->periodicPlanePatch() ? ln->periodicPlanePatch()->shift : VectorDim::Zero();
                static_cast<typename LoopNodeType::SplineNodeType&>(*ln).set_P(target-shift);
            }
            node->set_P(target);
        }
        // Transport only proven straight-edge subdivisions using their old
        // segment ancestry, before cutLoop changes the contour connectivity.
        // These are trial loop coordinates; boundary reconciliation recreates
        // the actual seam intersections and commits their network positions.
        for(const auto& subdivision : subdivisions)
        {
            const VectorDim position=(1-subdivision.parameter)*subdivision.source->get_P()+
                subdivision.parameter*subdivision.sink->get_P();
            static_cast<typename LoopNodeType::SplineNodeType&>(*subdivision.node).set_P(position);
        }
        subdivisions.clear();
        const auto groups=contractBoundary(a,b);
        if(!DN.contractNetworkNodes(a,b))
        {
            throw std::logic_error("Validated node contraction unexpectedly rejected.");
        }
        for(const auto& group : groups)
        {
            const auto survivor=*std::min_element(group.begin(),group.end(),[](const auto& x,const auto& y)
            {
                return x->sID<y->sID;
            });
            for(const auto& node : group)
            {
                if(node!=survivor && !node->loopNodes().empty() && !survivor->loopNodes().empty())
                {
                    DN.contractNetworkNodes(survivor,node);
                }
            }
        }
        if(reconcileBoundary)
        {
            DN.updateBoundaryNodes();
        }
        return true;
    }

    template <int dim>
    bool DislocationNodeContraction<dim>::contractSecondAndVirtual(std::shared_ptr<NetworkNodeType> nA,
        std::shared_ptr<NetworkNodeType> nB,const bool reconcileBoundary)
    {
        return commitContraction(nA,nB,nA->get_P(),reconcileBoundary);
    }

    template <int dim>
    bool DislocationNodeContraction<dim>::contractYoungest(std::shared_ptr<NetworkNodeType> nA,
                                                                              std::shared_ptr<NetworkNodeType> nB,const bool reconcileBoundary)
    {
        return nA->sID<nB->sID? contractSecondAndVirtual(nA,nB,reconcileBoundary) : contractSecondAndVirtual(nB,nA,reconcileBoundary);
        
    }

    template <int dim>
    bool DislocationNodeContraction<dim>::contractToPosition(std::shared_ptr<NetworkNodeType> nA,
                                                                                std::shared_ptr<NetworkNodeType> nB,
                                                                                const VectorDim& X,
                                                                                const double& maxRange,const bool reconcileBoundary)
    {
        
        if((nA->get_P()-X).norm()+(nB->get_P()-X).norm()>=maxRange)
        {
            return false;
        }
        return nA->sID<nB->sID ? commitContraction(nA,nB,X,reconcileBoundary) : commitContraction(nB,nA,X,reconcileBoundary);
    }


    template <int dim>
    bool DislocationNodeContraction<dim>::contract(std::shared_ptr<NetworkNodeType> nA,
                                                                      std::shared_ptr<NetworkNodeType> nB,const bool reconcileBoundary)
    {
        if(!nA || !nB || nA==nB || nA->loopNodes().empty() || nB->loopNodes().empty())
        {
            return false;
        }
        VerboseNodeContraction(1, "DislocationNodeContraction::contract " << nA->sID << " " << nB->sID << std::endl;);
        if (   !nA->isBoundaryNode()
            && !nA->isGrainBoundaryNode()
            && !nB->isBoundaryNode()
            && !nB->isGrainBoundaryNode()
            ) //Added for avoiding node contraction with boundary node
        {
            
            const bool nAisMovable = canReach(nA,nB->get_P());
            const bool nBisMovable = canReach(nB,nA->get_P());
            
            if (nAisMovable && nBisMovable)
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 1a" << std::endl;);
                return contractYoungest(nA, nB,reconcileBoundary);
            }
            else if (nAisMovable && !nBisMovable)
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 1b" << std::endl;);
                return contractSecondAndVirtual(nB, nA,reconcileBoundary);
            }
            else if (!nAisMovable && nBisMovable)
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 1c" << std::endl;);
                return contractSecondAndVirtual(nA, nB,reconcileBoundary);
            }
            else
            {// nA and nB cannot be moved to each other. The calculation of a third point is necessary
                VerboseNodeContraction(1, "DislocationNodeContraction case 1d" << std::endl;);
                const auto planesA(nA->glidePlanes());
                const auto planesB(nB->glidePlanes());
                Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> N(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,planesA.size()+planesB.size()));
                Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P(Eigen::Matrix<double,dim,Eigen::Dynamic>::Zero(dim,planesA.size()+planesB.size()));
                
                int k=0;
                for(const auto& plane : planesA)
                {
                    N.col(k)=plane->unitNormal;
                    P.col(k)=plane->P;
                    ++k;
                }
                for(const auto& plane : planesB)
                {
                    N.col(k)=plane->unitNormal;
                    P.col(k)=plane->P;
                    ++k;
                }
                
                PlanesIntersection<dim> pInt(N,P,FLT_EPSILON);
                const std::pair<bool,VectorDim> snapped(pInt.snap(0.5*(nA->get_P()+nB->get_P())));
                if(snapped.first)
                {
                    const double maxRange = 4.0 * (nA->get_P() - nB->get_P()).norm();
                    return contractToPosition(nA,nB,snapped.second,maxRange,reconcileBoundary);
                }
                else
                {
                    return false;
                }
            }
        }
        else
        {
            if ((nA->get_P() - nB->get_P()).norm() <= FLT_EPSILON)
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 2a" << std::endl;);
                return contractYoungest(nA, nB,reconcileBoundary);
            }
            else
            {
                VerboseNodeContraction(1, "DislocationNodeContraction case 2b" << std::endl;);
                return false;
            }
        }
    }

template class DislocationNodeContraction<3>;
}
#endif

