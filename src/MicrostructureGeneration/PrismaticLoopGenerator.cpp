/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PrismaticLoopGenerator_cpp_
#define model_PrismaticLoopGenerator_cpp_

#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <numbers>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <limits>
#include <cmath>
#include <numbers> // std::numbers

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationDynamics/DislocationNodeIO.h>
#include <DislocationDynamics/DislocationLoopIO.h>
#include <DislocationDynamics/DislocationLoopLinkIO.h>
#include <DislocationDynamics/DislocationLoopNodeIO.h>
#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <DefectiveCrystal/DefectiveCrystalAuxIO.h>

#include <DislocationDynamics/DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationDynamics/DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <PrismaticLoopGenerator.h>
#include <PlanesIntersection.h>
#include <LatticePlane.h>

namespace model
{

//PrismaticLoopGenerator::PrismaticLoopGenerator(const std::string& fileName) :
///* init */ MicrostructureGeneratorBase(fileName)
//{
//    
//}

PrismaticLoopGenerator::PrismaticLoopGenerator(const PrismaticLoopDensitySpecification& spec,MicrostructureGenerator& mg)
{
    
    std::cout<<magentaBoldColor<<"Generating prismatic loop density"<<defaultColor<<std::endl;
    if(spec.targetDensity>0.0)
    {
        std::normal_distribution<double> radiusDistribution(spec.radiusDistributionMean/mg.ddBase.poly.b_SI,spec.radiusDistributionStd/mg.ddBase.poly.b_SI);
        std::mt19937 generator;
        
        std::uniform_int_distribution<> allowedSlipSystemDist(0,spec.allowedSlipSystemIDs.size()-1);

        std::set<int> allowedGrainIDsSet;
        for(const auto& gID : spec.allowedGrainIDs)
        {
            allowedGrainIDsSet.insert(gID);
        }
        const double allowAllGrains(allowedGrainIDsSet.size() ? (*allowedGrainIDsSet.begin())<0 : true);
        
        double density=0.0;
        while(density<spec.targetDensity)
        {
            const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
            const LatticeVector<3> L0=rp.first;
            const size_t grainID=rp.second;
                        
            if(allowAllGrains || allowedGrainIDsSet.find(grainID)!=allowedGrainIDsSet.end())
            {
                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID)->slipSystems().size()-1);
                
                const int allowedIndex(allowedSlipSystemDist(generator));
                const int allowedSlipID(spec.allowedSlipSystemIDs[allowedIndex]);
                const int rSS(allowedSlipID<0 ? ssDist(generator) : allowedSlipID); // a random SlipSystem
                const double radius(radiusDistribution(generator));
                try
                {
                    density+=generateSingle(mg,rSS,L0.cartesian(),radius,50.0)/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                    std::cout<<"prismatic loop density="<<density<<std::endl;
                }
                catch(const std::exception& e)
                {
                    std::cout<<e.what()<<std::endl;
                }
            }

        }
    }
}

PrismaticLoopGenerator::PrismaticLoopGenerator(const PrismaticLoopIndividualSpecification& spec,MicrostructureGenerator& mg)
{
    std::cout<<magentaBoldColor<<"Generating individual prismatic loops"<<defaultColor<<std::endl;
    if(spec.slipSystemIDs.size())
    {
        if(spec.slipSystemIDs.size()!=spec.loopRadii.size())
        {
            throw std::runtime_error("spec.slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO spec.loopRadii.size()="+std::to_string(spec.loopRadii.size()));
        }
        if(int(spec.slipSystemIDs.size())!=spec.loopCenters.rows())
        {
            throw std::runtime_error("spec.slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO spec.loopCenters.rows()="+std::to_string(spec.loopCenters.rows()));
        }
        if(spec.slipSystemIDs.size()!=spec.glideSteps.size())
        {
            throw std::runtime_error("spec.slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO spec.glideSteps.size()="+std::to_string(spec.glideSteps.size()));
        }
        
        for(size_t k=0;k<spec.slipSystemIDs.size();++k)
        {
            generateSingle(mg,spec.slipSystemIDs[k],spec.loopCenters.row(k),spec.loopRadii[k]/mg.ddBase.poly.b_SI,spec.glideSteps[k]);
        }
    }
}

double PrismaticLoopGenerator::generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& guessCenter,const double& radius,const double& L)
{
    std::pair<bool,const Simplex<3,3>*> found(mg.ddBase.mesh.search(guessCenter));
    if(found.first)
    {
        const int grainID(found.second->region->regionID);
        // Select the host grain from the center, also for polycrystalline meshes.
        const auto& grain(mg.ddBase.poly.grain(grainID));
        
        if(rSS>=0 && rSS<int(grain->slipSystems().size()))
        {
            // sort slip systems in the zone axis with incrasing angle from a reference
            const RationalLatticeDirection<3>  b(grain->slipSystems().at(rSS)->s);
            const VectorDimD unitAxis(b.cartesian().normalized());
            const VectorDimD refNormal(Plane<3>::getL2G(unitAxis).col(0));
            std::map<double,ReciprocalLatticeDirection<3>> ssMap;
            for(const auto& ss : grain->slipSystems())
            {
                if((ss.second->s-b).squaredNorm()<FLT_EPSILON)
                {
                    const double cosAngle(ss.second->unitNormal.dot(refNormal));
                    const double sinAngle(ss.second->unitNormal.dot(unitAxis.cross(refNormal)));
                    const double refAngle1(std::atan2(sinAngle,cosAngle));
                    const double refAngle2(std::atan2(-sinAngle,-cosAngle));
                    
                    if(std::fabs(refAngle1)<std::fabs(refAngle2))
                    {
                        ssMap.emplace(refAngle1,ss.second->n*(+1));
                    }
                    else
                    {
                        ssMap.emplace(refAngle2,ss.second->n*(-1));
                    }
                }
            }
            
            // Recenter and define sessile plane

            const VectorDimD center(grain->snapToLattice(guessCenter).cartesian());
            const ReciprocalLatticeDirection<3> axis(grain->reciprocalLatticeDirection(b.cartesian()));
            const int sessilePlaneHeight(axis.planeIndexOfPoint(center));
//            const LatticePlaneKey sessilePlaneKey(axis,sessilePlaneHeight,grain->grainID);
            const LatticePlaneKey sessilePlaneKey(sessilePlaneHeight,axis);
            const auto sessilePlane(mg.ddBase.dislocationDynamics().periodicGlidePlaneFactory.getFromKey(sessilePlaneKey));
            
            // Define prismatic planes.
            std::vector<LatticePlane> planes;
            
            for(const auto& pair : ssMap)
            {// first half of the prism
                const int planeHeight(pair.second.closestPlaneIndexOfPoint(center+pair.second.cartesian().normalized()*radius));
                planes.emplace_back(planeHeight,pair.second);
            }
            for(const auto& pair : ssMap)
            {// second half of the prism
                const int planeHeight(pair.second.closestPlaneIndexOfPoint(center-pair.second.cartesian().normalized()*radius));
                planes.emplace_back(planeHeight,pair.second);
            }
            
            // Intersect two consecutive prismatic planes and sessile plane to find polygon points
            std::vector<VectorDimD> sessileLoopNodePos;
            for(size_t k1=0;k1<planes.size();++k1)
            {
                size_t k2(k1==planes.size()-1? 0 : k1+1);
                const auto& plane1(planes[k1]);
                const auto& plane2(planes[k2]);
                
                Eigen::Matrix<double,3,3> N;
                N.col(0)=sessilePlane->referencePlane->unitNormal;
                N.col(1)=plane1.n.cartesian().normalized();
                N.col(2)=plane2.n.cartesian().normalized();
                
                Eigen::Matrix<double,3,3> P;
                P.col(0)=sessilePlane->referencePlane->P;
                P.col(1)=plane1.planeOrigin();
                P.col(2)=plane2.planeOrigin();
                
                PlanesIntersection<3> pi(N,P,FLT_EPSILON);
                const std::pair<bool,VectorDimD> intPair(pi.snap(VectorDimD::Zero()));
                if(intPair.first)
                {
                    sessileLoopNodePos.push_back(intPair.second);
                }
                else
                {
                    throw std::runtime_error("Cannot determine prismatic polygon.");
                }
                
            }
            
            // Map points along rays inside the convex prism. Crossing a physical
            // face rejects the prism; same-grain periodic faces translate the ray.
            const auto mappedShift=[&](const VectorDimD& target)
            {
                VectorDimD shift=VectorDimD::Zero();
                const VectorDimD direction=target-center;
                double parameter=0.0;
                for(size_t crossing=0;crossing<10000;++crossing)
                {
                    const PlanarMeshFace<3>* exitFace=nullptr;
                    double exitParameter=1.0;
                    for(const auto& entry : grain->region.faces())
                    {
                        const auto* face=entry.second.get();
                        const VectorDimD normal=face->outNormal(grainID);
                        const double speed=normal.dot(direction);
                        if(speed>FLT_EPSILON)
                        {
                            const double candidate=normal.dot(face->center()-center-shift)/speed;
                            if(candidate>=parameter-1.e-9 && candidate<exitParameter-1.e-9)
                            {
                                exitParameter=candidate;
                                exitFace=face;
                            }
                        }
                    }
                    if(!exitFace)
                    {
                        return shift;
                    }
                    const auto* partner=exitFace->periodicFacePair.second;
                    if(!partner || exitFace->regionIDs!=partner->regionIDs)
                    {
                        throw std::runtime_error("Prismatic loop crosses a physical grain boundary; entire prism rejected.");
                    }
                    shift+=partner->center()-exitFace->center();
                    parameter=exitParameter;
                }
                throw std::runtime_error("Prismatic loop exceeded periodic crossing limit; entire prism rejected.");
            };
            struct PreparedNode
            {
                VectorDimD position;
                VectorDimD shift;
                std::pair<short int,short int> edges;
            };
            struct PreparedLoop
            {
                std::shared_ptr<PeriodicGlidePlane<3>> plane;
                DislocationLoopIO<3>::DislocationLoopType type;
                std::vector<PreparedNode> nodes;
            };
            // Validate every face of the closed prism before modifying the IO.
            // The patch walk crosses only same-grain periodic faces. Any area
            // removed by its physical-boundary clipping rejects the whole prism.
            const auto prepare=[&](const std::vector<VectorDimD>& polygon,
                                   const std::shared_ptr<PeriodicGlidePlane<3>>& plane,
                                   const DislocationLoopIO<3>::DislocationLoopType type)
            {
                std::vector<PolyPoint> anchors(polygon.size());
                std::vector<std::pair<VectorDimD,const PolyPoint* const>> points;
                for(size_t k=0;k<polygon.size();++k)
                {
                    points.emplace_back(polygon[k],&anchors[k]);
                }
                const auto intersections=plane->polygonPatchIntersection(points,true);
                PreparedLoop result{plane,type,{}};
                for(const auto& point : intersections)
                {
                    result.nodes.push_back({plane->referencePlane->globalPosition(std::get<0>(point)),
                                            std::get<1>(point),std::get<2>(point)});
                }
                VectorDimD originalArea=VectorDimD::Zero();
                VectorDimD retainedArea=VectorDimD::Zero();
                for(size_t k=0;k<polygon.size();++k)
                {
                    originalArea+=0.5*(polygon[k]-polygon.front()).cross(polygon[(k+1)%polygon.size()]-polygon.front());
                }
                for(size_t k=0;k<result.nodes.size();++k)
                {
                    retainedArea+=0.5*(result.nodes[k].position-polygon.front()).cross(
                        result.nodes[(k+1)%result.nodes.size()].position-polygon.front());
                }
                if(result.nodes.size()<3 || originalArea.norm()<=FLT_EPSILON ||
                   (retainedArea-originalArea).norm()>1.e-8*originalArea.norm())
                {
                    throw std::runtime_error("Prismatic loop crosses a physical grain boundary; entire prism rejected.");
                }
                return result;
            };
            std::vector<PreparedLoop> prepared;
            prepared.push_back(prepare(sessileLoopNodePos,sessilePlane,DislocationLoopIO<3>::SESSILELOOP));
            const VectorDimD step(L*b.cartesian());
            std::vector<VectorDimD> top;
            for(const auto& point : sessileLoopNodePos)
            {
                mappedShift(point);
                mappedShift(point+step);
                top.push_back(point+step);
            }
            const VectorDimD topShift=mappedShift(center+step);
            for(auto& point : top)
            {
                point+=topShift;
            }
            const auto topPlane=mg.ddBase.dislocationDynamics().periodicGlidePlaneFactory.getFromKey(
                LatticePlaneKey(axis.planeIndexOfPoint(center+step+topShift),axis));
            // The cap is validation-only: inserting it would cancel the physical loop.
            prepare(top,topPlane,DislocationLoopIO<3>::SESSILELOOP);
            double loopLength=0.0;
            for(size_t k1=0;k1<planes.size();++k1)
            {
                const size_t k2=(k1+1)%planes.size();
                loopLength+=(sessileLoopNodePos[k2]-sessileLoopNodePos[k1]).norm();
                std::vector<VectorDimD> side{sessileLoopNodePos[k2],sessileLoopNodePos[k1],
                                                  sessileLoopNodePos[k1]+step,sessileLoopNodePos[k2]+step};
                const VectorDimD sideShift=mappedShift(0.5*(sessileLoopNodePos[k1]+sessileLoopNodePos[k2]+step));
                for(auto& point : side)
                {
                    point+=sideShift;
                }
                const LatticePlaneKey key(planes[k2].n.planeIndexOfPoint(side.front()),planes[k2].n);
                prepared.push_back(prepare(side,mg.ddBase.dislocationDynamics().periodicGlidePlaneFactory.getFromKey(key),
                                           DislocationLoopIO<3>::GLISSILELOOP));
            }
            // Commit all prepared loops together; also roll back unexpected
            // insertion failures, including the position-to-network-node index.
            const size_t oldLoops=mg.config().loops().size();
            const size_t oldNodes=mg.config().nodes().size();
            const size_t oldLoopNodes=mg.config().loopNodes().size();
            const size_t oldLinks=mg.config().loopLinks().size();
            auto oldNodeMap=mg.uniqueNetworkNodeMap;
            try
            {
                for(const auto& loop : prepared)
                {
                    const size_t id=mg.insertLoop(b.cartesian(),loop.plane->referencePlane->unitNormal,
                                                  loop.plane->referencePlane->P,grainID,loop.type);
                    std::vector<size_t> nodeIDs;
                    for(const auto& node : loop.nodes)
                    {
                        const VectorDimD position=node.position+node.shift;
                        auto foundNode=mg.uniqueNetworkNodeMap.find(position);
                        if(foundNode==mg.uniqueNetworkNodeMap.end())
                        {
                            const size_t networkID=mg.insertNetworkNode(position);
                            foundNode=mg.uniqueNetworkNodeMap.emplace(position,networkID).first;
                        }
                        nodeIDs.push_back(mg.insertLoopNode(id,node.position,foundNode->second,node.shift,node.edges));
                    }
                    mg.insertLoopLinks(id,nodeIDs);
                }
            }
            catch(...)
            {
                mg.config().loopLinks().resize(oldLinks);
                mg.config().loopNodes().resize(oldLoopNodes);
                mg.config().nodes().resize(oldNodes);
                mg.config().loops().resize(oldLoops);
                mg.uniqueNetworkNodeMap.swap(oldNodeMap);
                throw;
            }
            return loopLength;
        }
        else
        {
            if(rSS<0)
            {
                std::cout<<"Skipping slip system "<<rSS<<std::endl;
            }
            else
            {
                throw std::runtime_error("slipSystem "+std::to_string(rSS)+" not found, skipping.");
            }
        }
    }
    else
    {
        std::cout<<"Center outside mesh, skipping prismatic loop"<<std::endl;
    }
    return 0.0;
}

}
#endif
