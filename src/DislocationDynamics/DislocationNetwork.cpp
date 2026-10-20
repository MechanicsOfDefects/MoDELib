/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetwork_CPP_
#define model_DislocationNetwork_CPP_

#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <PeriodicContourTopology.h>
#include <numbers>

#include <DislocationNetwork.h>
#include <DislocationNucleation.h>
#include <algorithm>


namespace model
{


    template <int dim>
    DislocationNetwork<dim>::DislocationNetwork(MicrostructureContainerType& mc) :
    /* init */ MicrostructureBase<dim>("DislocationDynamics",mc)
    /* init */,glideStepsSinceLastClimb(0)
    /* init */,ddBase(this->microstructures.ddBase)
    /* init */,networkRemesher(*this)
    /* init */,junctionsMaker(*this)
    /* init */,crossSlipModel(DislocationCrossSlip<dim>::getModel(ddBase.poly,ddBase.simulationParameters.traitsIO))
    /* init */,crossSlipMaker(*this)
    /* init */,nodeContractor(*this)
    /* init */,timeStepper(*this)
    /* init */,stochasticForceGenerator(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("useStochasticForce",true)? new StochasticForceGenerator(ddBase.simulationParameters.traitsIO) : nullptr)
    /* init */,bulkNucleationModel(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("bulkNucleationModel",true))
    /* init */,surfaceNucleationModel(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("surfaceNucleationModel",true))
    /* init */,computeDDinteractions(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("computeDDinteractions",true))
    /* init */,outputQuadraturePoints(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("outputQuadraturePoints",true))
    /* init */,outputLinkingNumbers(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("outputLinkingNumbers",true))
    /* init */,outputLoopLength(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("outputLoopLength",true))
    /* init */,outputSegmentPairDistances(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("outputSegmentPairDistances",true))
    /* init */,outputPlasticDistortionPerSlipSystem(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("outputPlasticDistortionPerSlipSystem",true))
    /* init */,outputDislocationDensityPerSlipSystem(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("outputDislocationDensityPerSlipSystem",true))
    /* init */,computeElasticEnergyPerLength(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("computeElasticEnergyPerLength",true))
    /* init */,alphaLineTension(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("alphaLineTension",true))
    /* init */,use_velocityFilter(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("use_velocityFilter",true))
    /* init */,velocityReductionFactor(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("velocityReductionFactor",true))
    /* init */,nodalVelocityConstraints(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readMatrixCols<double>("nodalVelocityConstraints",dim,true))
    /* init */,subcyclingBins(getSubCyclingSet(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readArray<int>("subcyclingBins",true)))
    /* init */,verboseDislocationNode(TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("verboseDislocationNode",true))
    {
        assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
        LoopNetworkType::verboseLevel=TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("verboseLoopNetwork",true);
        verboseDislocationNetwork=TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<int>("verboseDislocationNetwork",true);
    }

template <int dim>
std::set<int> DislocationNetwork<dim>::getSubCyclingSet(const std::vector<int> &inpVector)
{
    std::set<int> temp;
    for (const auto &iv : inpVector)
    {
        temp.insert(iv);
    }
    return (temp.size()>1)? temp : std::set<int>();
}

    template <int dim>
    void DislocationNetwork<dim>::initializeConfiguration(const DefectiveCrystalConfigIO<dim>& configIO,const std::ofstream&,const std::ofstream&)
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        LoopType::initFromFile(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt");
        LoopNodeType::initFromFile(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt");
        LoopLinkType::initFromFile(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt");
        NetworkLinkType::initFromFile(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt");
        DislocationFieldBase<dim>::initFromFile(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt");
        
        glideSolver=DislocationGlideSolverFactory<dim>::getGlideSolver(*this,TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readString("glideSolverType",false));
        climbSolver=DislocationClimbSolverFactory<dim>::getClimbSolver(*this,TextFileParser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readString("climbSolverType",false));
        _inclusions=this->microstructures.template getUniqueTypedMicrostructure<InclusionMicrostructure<dim>>();
        
        setConfiguration(configIO);
    }

    template <int dim>
    void DislocationNetwork<dim>::addConfiguration(const DefectiveCrystalConfigIO<dim>& evl)
    {
        //    DislocationNode<3>::force_count(0);
        //    DislocationLoopNode<3>::force_count(0);
        //    DislocationLoop<3>::force_count(0);
        //    EshelbyInclusionBase<3>::force_count(0);
        //
        //    this->loopLinks().clear(); // erase base network to clear current config
        
        const size_t initialLoopCount(DislocationLoop<3>::get_count());
        const size_t initialNetworkNodeCount(NetworkNodeType::get_count());
        const size_t initialLoopNodeCount(DislocationLoopNode<3>::get_count());
        
        // Create Loops
        std::deque<std::shared_ptr<LoopType>> tempLoops; // keep loops alive during setConfiguration
        size_t loopNumber=1;
        for(const auto& loop : evl.loops())
        {
            const bool faulted(ddBase.poly.grain(loop.grainID)->rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false);
            VerboseDislocationNetwork(1,"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;);
            
            //std::cout<<"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
            const size_t loopIDtoUse(loop.sID+initialLoopCount);
            LoopType::set_count(loopIDtoUse);
            
            GlidePlaneKey<dim> loopPlaneKey(loop.P, ddBase.poly.grain(loop.grainID)->reciprocalLatticeDirection(loop.N));
            tempLoops.push_back(this->loops().create(loop.B, ddBase.dislocationDynamics().glidePlaneFactory.getFromKey(loopPlaneKey)));
            assert(this->loops().get(loopIDtoUse)->sID == loopIDtoUse);
            loopNumber++;
        }
        
        // Create NetworkNodes
        std::deque<std::shared_ptr<NetworkNodeType>> tempNetNodes; // keep loops alive during setConfiguration
        size_t netNodeNumber=1;
        for(const auto& node : evl.nodes())
        {
            VerboseDislocationNetwork(1,"Creating DislocationNode "<<node.sID<<" ("<<netNodeNumber<<" of "<<evl.nodes().size()<<")"<<std::endl;);
            const size_t nodeIDtoUse(node.sID+initialNetworkNodeCount);
            NetworkNodeType::set_count(nodeIDtoUse);
            tempNetNodes.push_back(this->networkNodes().create(node.P,node.V,node.climbVelocityScalar,node.velocityReduction));
            assert(this->networkNodes().get(nodeIDtoUse)->sID==nodeIDtoUse);
            netNodeNumber++;
        }
        
        // Create LoopNodes
        std::deque<std::shared_ptr<LoopNodeType>> tempLoopNodes; // keep loops alive during setConfiguration
        size_t loopNodeNumber=1;
        for(const auto& node : evl.loopNodes())
        {
            VerboseDislocationNetwork(1,"Creating DislocationLoopNode "<<node.sID<<" ("<<loopNodeNumber<<" of "<<evl.loopNodes().size()<<")"<<std::endl;);
            const size_t nodeIDtoUse(node.sID+initialLoopNodeCount);
            LoopNodeType::set_count(nodeIDtoUse);
            const auto loop(this->loops().get(node.loopID+initialLoopCount));
            const auto netNode(this->networkNodes().get(node.networkNodeID+initialNetworkNodeCount));
            const auto periodicPatch(loop->periodicGlidePlane? loop->periodicGlidePlane->patches().getFromKey(node.periodicShift) : nullptr);
            const auto periodicPatchEdge((periodicPatch && node.edgeIDs.first>=0)? (node.edgeIDs.second>=0 ? std::make_pair(periodicPatch->edges()[node.edgeIDs.first],
                                                                                                                            periodicPatch->edges()[node.edgeIDs.second]):
                                                                                    std::make_pair(periodicPatch->edges()[node.edgeIDs.first],nullptr)):std::make_pair(nullptr,nullptr));
            tempLoopNodes.push_back(this->loopNodes().create(loop,netNode,node.P,periodicPatch,periodicPatchEdge));
            assert(this->loopNodes().get(nodeIDtoUse)->sID==nodeIDtoUse);
            loopNodeNumber++;
        }
        
        // Insert Loops
        std::map<size_t,std::map<size_t,size_t>> loopMap;
        for(const auto& looplink : evl.loopLinks())
        {// Collect LoopLinks by loop IDs
            loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
        }
        assert(loopMap.size()==evl.loops().size());
        
        for(const auto& loop : evl.loops())
        {// for each loop in the DefectiveCrystalConfigIO<dim> object
            
            const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
            assert(loopFound!=loopMap.end());
            std::vector<std::shared_ptr<LoopNodeType>> loopNodes;
            loopNodes.push_back(this->loopNodes().get(loopFound->second.begin()->first+initialLoopNodeCount));
            for(size_t k=0;k<loopFound->second.size();++k)
            {
                const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
                if(k<loopFound->second.size()-1)
                {
                    loopNodes.push_back(this->loopNodes().get(nodeFound->second+initialLoopNodeCount));
                }
                else
                {
                    assert(nodeFound->second+initialLoopNodeCount==loopNodes[0]->sID);
                }
            }
            //        std::cout<<" Inserting loop "<<loop.sID<<std::endl;
            this->insertLoop(this->loops().get(loop.sID+initialLoopCount),loopNodes);
            // Recover contour bends from older files without changing the IO
            // format. Geometry is still committed at this point.
            retainPeriodicContourBends(loopNodes,1.e-8*(ddBase.mesh.xMax()-ddBase.mesh.xMin()).norm());
        }
        updateGeometry();
    }

    template <int dim>
    void DislocationNetwork<dim>::setConfiguration(const DefectiveCrystalConfigIO<dim>& evl)
    {
        DislocationNode<3>::force_count(0);
        DislocationLoopNode<3>::force_count(0);
        DislocationLoop<3>::force_count(0);
        //        EshelbyInclusionBase<3>::force_count(0);
        
        this->loopLinks().clear(); // erase base network to clear current config
        addConfiguration(evl);
        
        //        // Create Loops
        //        std::deque<std::shared_ptr<LoopType>> tempLoops; // keep loops alive during setConfiguration
        //        size_t loopNumber=1;
        //        for(const auto& loop : evl.loops())
        //        {
        //            const bool faulted(ddBase.poly.grain(loop.grainID)->rationalLatticeDirection(loop.B).rat.asDouble()!=1.0? true : false);
        //            VerboseDislocationNetwork(1,"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;);
        //
        //            //std::cout<<"Creating DislocationLoop "<<loop.sID<<" ("<<loopNumber<<" of "<<evl.loops().size()<<"), type="<<loop.loopType<<", faulted="<<faulted<<", |b|="<<loop.B.norm()<<std::endl;
        //            const size_t loopIDinFile(loop.sID);
        //            LoopType::set_count(loopIDinFile);
        //
        //            GlidePlaneKey<dim> loopPlaneKey(loop.P, ddBase.poly.grain(loop.grainID)->reciprocalLatticeDirection(loop.N));
        //            tempLoops.push_back(this->loops().create(loop.B, ddBase.dislocationDynamics().glidePlaneFactory.getFromKey(loopPlaneKey)));
        //            assert(this->loops().get(loopIDinFile)->sID == loopIDinFile);
        //            loopNumber++;
        //        }
        //
        //        // Create NetworkNodes
        //        std::deque<std::shared_ptr<NetworkNodeType>> tempNetNodes; // keep loops alive during setConfiguration
        //        size_t netNodeNumber=1;
        //        for(const auto& node : evl.nodes())
        //        {
        //            VerboseDislocationNetwork(1,"Creating DislocationNode "<<node.sID<<" ("<<netNodeNumber<<" of "<<evl.nodes().size()<<")"<<std::endl;);
        //            const size_t nodeIDinFile(node.sID);
        //            NetworkNodeType::set_count(nodeIDinFile);
        //            tempNetNodes.push_back(this->networkNodes().create(node.P,node.V,node.climbVelocityScalar,node.velocityReduction));
        //            assert(this->networkNodes().get(nodeIDinFile)->sID==nodeIDinFile);
        //            netNodeNumber++;
        //        }
        //
        //        // Create LoopNodes
        //        std::deque<std::shared_ptr<LoopNodeType>> tempLoopNodes; // keep loops alive during setConfiguration
        //        size_t loopNodeNumber=1;
        //        for(const auto& node : evl.loopNodes())
        //        {
        //            VerboseDislocationNetwork(1,"Creating DislocationLoopNode "<<node.sID<<" ("<<loopNodeNumber<<" of "<<evl.loopNodes().size()<<")"<<std::endl;);
        //            const size_t nodeIDinFile(node.sID);
        //            LoopNodeType::set_count(nodeIDinFile);
        //            const auto loop(this->loops().get(node.loopID));
        //            const auto netNode(this->networkNodes().get(node.networkNodeID));
        //            const auto periodicPatch(loop->periodicGlidePlane? loop->periodicGlidePlane->patches().getFromKey(node.periodicShift) : nullptr);
        //            const auto periodicPatchEdge((periodicPatch && node.edgeIDs.first>=0)? (node.edgeIDs.second>=0 ? std::make_pair(periodicPatch->edges()[node.edgeIDs.first],
        //                                                                                                                            periodicPatch->edges()[node.edgeIDs.second]):
        //                                                                                    std::make_pair(periodicPatch->edges()[node.edgeIDs.first],nullptr)):std::make_pair(nullptr,nullptr));
        //            tempLoopNodes.push_back(this->loopNodes().create(loop,netNode,node.P,periodicPatch,periodicPatchEdge));
        //            assert(this->loopNodes().get(nodeIDinFile)->sID==nodeIDinFile);
        //            loopNodeNumber++;
        //        }
        //
        //        // Insert Loops
        //        std::map<size_t,std::map<size_t,size_t>> loopMap;
        //        for(const auto& looplink : evl.loopLinks())
        //        {// Collect LoopLinks by loop IDs
        //            loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
        //        }
        //        assert(loopMap.size()==evl.loops().size());
        //
        //        for(const auto& loop : evl.loops())
        //        {// for each loop in the DefectiveCrystalConfigIO<dim> object
        //
        //            const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
        //            assert(loopFound!=loopMap.end());
        //            std::vector<std::shared_ptr<LoopNodeType>> loopNodes;
        //            loopNodes.push_back(this->loopNodes().get(loopFound->second.begin()->first));
        //            for(size_t k=0;k<loopFound->second.size();++k)
        //            {
        //                const auto nodeFound=loopFound->second.find(loopNodes.back()->sID);
        //                if(k<loopFound->second.size()-1)
        //                {
        //                    loopNodes.push_back(this->loopNodes().get(nodeFound->second));
        //                }
        //                else
        //                {
        //                    assert(nodeFound->second==loopNodes[0]->sID);
        //                }
        //            }
        //            //        std::cout<<" Inserting loop "<<loop.sID<<std::endl;
        //            this->insertLoop(this->loops().get(loop.sID),loopNodes);
        //        }
        //        updateGeometry();
    }

    template <int dim>
    void DislocationNetwork<dim>::updateGeometry()
    {
        VerboseDislocationNetwork(2,"DislocationNetwork::updateGeometry"<<std::endl;);
        for(auto& loop : this->loops())
        {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
            loop.second.lock()->updateGeometry();
        }
        // updatePlasticDistortionRateFromAreas();
        VerboseDislocationNetwork(3,"DislocationNetwork::updateGeometry DONE"<<std::endl;);
    }

    template <int dim>
    typename DislocationNetwork<dim>::MatrixDim DislocationNetwork<dim>::averagePlasticDistortion() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& loop : this->loops())
        {
            temp+= loop.second.lock()->averagePlasticDistortion();
        }
        return temp;
    }

    template <int dim>
    typename DislocationNetwork<dim>::MatrixDim DislocationNetwork<dim>::averagePlasticStrain() const
    {/*!\returns the plastic strain rate tensor generated during the last time step.
      */
        const MatrixDim apd(averagePlasticDistortion());
        return 0.5*(apd+apd.transpose());
    }

    template <int dim>
    std::map<std::pair<int,int>,double> DislocationNetwork<dim>::slipSystemAveragePlasticDistortion() const
    {
        std::map<std::pair<int,int>,double> temp; // <grainID,slipSystemID>
        for(const auto& grain : ddBase.poly.grains)
        {
            for(const auto& slipSystem : grain.second->slipSystems())
            {
                const std::pair<int,int> key(std::make_pair(grain.first,slipSystem.second->sID));
                temp[key]=0.0;
            }
        }
        
        for(const auto& weakloop : this->loops())
        {
            const auto loop(weakloop.second.lock());
            if(loop->slipSystem())
            {
                const double val = loop->slipSystem()->unitSlip.transpose()*loop->averagePlasticDistortion()*loop->slipSystem()->unitNormal;
                const std::pair<int,int> key(std::make_pair(loop->grain.grainID,loop->slipSystem()->sID));
                auto iter(temp.find(key));
                if(iter!=temp.end())
                {// iter found
                    iter->second += val ;
                }
                else
                {// iter not found
                    temp.emplace(key,val);
                }
            }
        }
        return temp;
    }

    template <int dim>
    typename DislocationNetwork<dim>::MatrixDim DislocationNetwork<dim>::averagePlasticDistortionRate() const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& loop : this->loops())
        {
            temp+= loop.second.lock()->averagePlasticDistortionRate();
        }
        return temp;
    }

    template <int dim>
    typename DislocationNetwork<dim>::MatrixDim DislocationNetwork<dim>::averagePlasticStrainRate() const
    {/*!\returns the plastic strain rate tensor generated during the last time step.
      */
        const MatrixDim apdr(averagePlasticDistortionRate());
        return 0.5*(apdr+apdr.transpose());
    }

    template <int dim>
    std::vector<std::tuple<double,double,double,double>> DislocationNetwork<dim>::networkLengthPerSlipSystem() const
    {
        std::vector<std::tuple<double,double,double,double>> temp;
        
        if(ddBase.poly.grains.size())
        {
            temp.resize(ddBase.poly.grains.begin()->second->slipSystems().size(),std::make_tuple(0.0,0.0,0.0,0.0));
            //            std::vector<std::tuple<double,double,double,double>> temp(ddBase.poly.grains.begin()->second->slipSystems().size(),std::make_tuple(0.0,0.0,0.0,0.0));
            for(auto& loop : this->loops())
            {
                if(loop.second.lock()->slipSystem())
                {
                    const size_t ssID(loop.second.lock()->slipSystem()->sID);
                    double& bulkGlissileLength(std::get<0>(temp[ssID]));
                    double& bulkSessileLength(std::get<1>(temp[ssID]));
                    double& boundaryLength(std::get<2>(temp[ssID]));
                    double& grainBoundaryLength(std::get<3>(temp[ssID]));
                    
                    for(const auto& loopLink : loop.second.lock()->loopLinks())
                    {
                        if(loopLink->networkLink())
                        {
                            if(!loopLink->networkLink()->hasZeroBurgers())
                            {
                                if(loopLink->networkLink()->isGrainBoundarySegment())
                                {
                                    grainBoundaryLength+=loopLink->networkLink()->chord().norm();
                                }
                                else if(loopLink->networkLink()->isBoundarySegment())
                                {
                                    boundaryLength+=loopLink->networkLink()->chord().norm();
                                }
                                else
                                {
                                    if(loopLink->networkLink()->isSessile())
                                    {
                                        bulkSessileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                                    }
                                    else
                                    {
                                        bulkGlissileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return temp;
    }

    template <int dim>
    std::tuple<double,double,double,double> DislocationNetwork<dim>::networkLength() const
    {/*!\returns the total line length of the DislocationNetwork. The return
      * value is a tuple, where the first value is the length of bulk glissile
      * dislocations, the second value is the length of bulk sessile
      * dislocations, and the third value is the length accumulated on
      * the mesh boundary.
      */
        double bulkGlissileLength(0.0);
        double bulkSessileLength(0.0);
        double boundaryLength(0.0);
        double grainBoundaryLength(0.0);
        
        for(auto& loop : this->loops())
        {
            for(const auto& loopLink : loop.second.lock()->loopLinks())
            {
                if(loopLink->networkLink())
                {
                    if(!loopLink->networkLink()->hasZeroBurgers())
                    {
                        if(loopLink->networkLink()->isGrainBoundarySegment())
                        {
                            grainBoundaryLength+=loopLink->networkLink()->chord().norm();
                        }
                        else if(loopLink->networkLink()->isBoundarySegment())
                        {
                            boundaryLength+=loopLink->networkLink()->chord().norm();
                        }
                        else
                        {
                            if(loopLink->networkLink()->isSessile())
                            {
                                bulkSessileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                            }
                            else
                            {
                                bulkGlissileLength+=loopLink->networkLink()->chord().norm()/loopLink->networkLink()->loopLinks().size();
                            }
                        }
                    }
                }
            }
        }
        return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
    }


    template <int dim>
    const std::shared_ptr<InclusionMicrostructure<dim>>& DislocationNetwork<dim>::inclusions() const
    {
        return _inclusions;
    }

    template <int dim>
    bool DislocationNetwork<dim>::contract(std::shared_ptr<NetworkNodeType> nA,
                                           std::shared_ptr<NetworkNodeType> nB)
    {
        return nodeContractor.contract(nA,nB);
    }

    template <int dim>
    typename DislocationNetwork<dim>::VectorDim DislocationNetwork<dim>::displacement(const VectorDim& x,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        VectorDim temp(VectorDim::Zero());
        
        for(const auto& loop : this->loops())
        {// sum solid angle of each loop
            for(const auto& shift : ddBase.periodicShifts)
            {
                temp-=loop.second.lock()->solidAngle(x+shift)/4.0/std::numbers::pi*loop.second.lock()->burgers();
            }
        }
        
        for(const auto& link : this->networkLinks())
        {// sum line-integral part of displacement field per segment
            if(   !link.second.lock()->hasZeroBurgers())
            {
                for(const auto& shift : ddBase.periodicShifts)
                {
                    temp+=link.second.lock()->straight.displacement(x+shift);
                }
            }
        }
        
        return temp;
    }

    template <int dim>
    typename DislocationNetwork<dim>::MatrixDim DislocationNetwork<dim>::averageStress() const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        return MatrixDim::Zero();
    }

    template <int dim>
    typename DislocationNetwork<dim>::VectorMSize DislocationNetwork<dim>::mobileConcentration(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        VectorMSize temp(VectorMSize::Zero());
        if(climbSolver)
        {
            const auto pointGrains(this->pointGrains(x,node,ele,guess));
            //            std::set<const Grain<dim>*> pointGrains;
            //            if(node)
            //            {
            //                for(const auto& nodeEle : *node)
            //                {
            //                    pointGrains.emplace(&ddBase.poly.grain(nodeEle->simplex.region->regionID));
            //                }
            //            }
            //            else
            //            {
            //                if(ele)
            //                {
            //                    pointGrains.emplace(&ddBase.poly.grain(ele->simplex.region->regionID));
            //                }
            //                else
            //                {
            //                    const std::pair<bool,const Simplex<dim,dim>*> found(ddBase.mesh.searchWithGuess(x,guess));
            //                    if(found.first)
            //                    {
            //                        pointGrains.emplace(&ddBase.poly.grain(found.second->region->regionID));
            //                    }
            //                }
            //            }
            
            if(pointGrains.size())
            {
                for(const auto& link : this->networkLinks())
                {// sum stress field per segment
                    if(   !link.second.lock()->hasZeroBurgers()
                       )
                    {
                        
                        const auto segment(link.second.lock());
                        const auto segGrains(segment->grains());
                        std::set<const Grain<dim>*> intersect;
                        std::set_intersection(pointGrains.begin(), pointGrains.end(), segGrains.begin(), segGrains.end(),
                                              std::inserter(intersect, intersect.begin()));
                        if(intersect.size()==1)
                        {
                            //                            const int grainID((*intersect.begin())->region.regionID);
                            //                            StressStraight<3> ss(ddBase.poly,segment->source->get_P(),segment->sink->get_P(),segment->burgers(),ddBase.dislocationDynamics().EwaldLength);
                            //                            for(const auto& shift : ddBase.periodicShifts)
                            //                            {
                            temp+=segment->clusterConcentration(x,climbSolver->CD->cdp);
                            
                            //                                temp+=ss.clusterConcentration(x+shift,grainID, segment->source->climbDirection(), segment->source->climbVelocityScalar , segment->sink->climbDirection(), segment->sink->climbVelocityScalar, climbSolver->CD->cdp);
                            //                            }
                        }
                    }
                }
            }
        }
        return temp;
    }


    template<int dim>
    typename DislocationNetwork<dim>::VectorISize DislocationNetwork<dim>::immobileClusters(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorISize::Zero();
    }
    

    template <int dim>
    typename DislocationNetwork<dim>::MatrixDim DislocationNetwork<dim>::stress(const VectorDim& x,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {/*!\param[in] P position vector
      * \returns The stress field generated by the DislocationNetwork at P
      *
      * Note:
      */
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& link : this->networkLinks())
        {// sum stress field per segment
            if(   !link.second.lock()->hasZeroBurgers()
               )
            {
                for(const auto& shift : ddBase.periodicShifts)
                {
                    temp+=link.second.lock()->straight.stress(x+shift);
                }
            }
        }
        return temp;
    }

    template <int dim>
    typename DislocationNetwork<dim>::VectorDim DislocationNetwork<dim>::inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const
    {
        return VectorDim::Zero();
    }

    template <int dim>
    void DislocationNetwork<dim>::updateConfiguration()
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        moveNodes(ddBase.simulationParameters.dt);
        executeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
        if(isClimbStep())
        {
            glideStepsSinceLastClimb=0;
        }
        else
        {
            glideStepsSinceLastClimb++;
        }
        updateGeometry();
    }

    template <int dim>
    double DislocationNetwork<dim>::getDt() const
    {
        if(isClimbStep())
        {
            return timeStepper.getDt(climbSolver->vClimbRef,1.0e-7);
            
        }
        else
        {
            //        climbSolver->CD->cdp
            return timeStepper.getDt(ddBase.poly.cs,1.0e-7);
        }
    }

    template <int dim>
    bool DislocationNetwork<dim>::isClimbStep() const
    {
        if(glideSolver)
        {
            return climbSolver?  (glideStepsSinceLastClimb && averagePlasticDistortionRate().norm()<climbSolver->glideEquilibriumRate): false;
        }
        else
        {
            if(climbSolver)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }

template <int dim>
void DislocationNetwork<dim>::reSolve()
{
    
    const bool isClimbingStep(isClimbStep());
    if(isClimbingStep)
    {
        std::cout<<" climbStep"<<std::flush;
    }
    else
    {
        std::cout<<" glideStep"<<std::flush;
    }
    double maxVelocity = 0.0;
    for (const auto &nodeIter : this->networkNodes())
    {
        const double vNorm(nodeIter.second.lock()->get_V().norm());
        if (vNorm > maxVelocity)
        {
            maxVelocity = vNorm;
        }
    }
    
#ifdef _OPENMP
    const size_t nThreads = omp_get_max_threads();
#else
    const size_t nThreads = 1;
#endif
    
    //! -1 Compute the interaction StressField between dislocation particles
    std::map<int, int> velocityBinMap{{1,0}};
    for (const auto &binVal : subcyclingBins)
    {
        velocityBinMap.emplace(binVal, 0);
    }
    
    
    std::cout <<" ,updating updateForcesAndVelocities (" << nThreads << " threads) " << std::flush;
#ifdef _OPENMP
    {
        // Topology is fixed during field evaluation. Preserve container order
        // and shared ownership, but make OpenMP's indexed access constant-time.
        // Release the snapshot before any subsequent topology changes.
        std::vector<std::shared_ptr<NetworkLinkType>> fieldLinks;
        fieldLinks.reserve(this->networkLinks().size());
        for(const auto& entry : this->networkLinks())
        {
            fieldLinks.push_back(entry.second.lock());
        }
#pragma omp parallel for
        for (size_t k = 0; k < fieldLinks.size(); ++k)
        {
            const auto& link=fieldLinks[k];
            const int velGroup((subcyclingBins.size() && !isClimbingStep) ? link->velocityGroup(maxVelocity, subcyclingBins) : 1);

            if ((ddBase.simulationParameters.runID % velGroup) == 0)
            {
                for(auto& qp : link->quadraturePoints())
                {
                    qp.updateForcesAndVelocities(*link,isClimbingStep);
                }
            }
            //            else
            //            {
            //                link->assembleGlide(false);
            //            }
        }
    }

#else
    for (auto &linkIter : this->networkLinks())
    {
        const auto link=linkIter.second.lock();
        const int velGroup(subcyclingBins.size() ? linkIter.second.lock()->velocityGroup(maxVelocity, subcyclingBins) : 1);
        
        if ((ddBase.simulationParameters.runID % velGroup) == 0)
        {
            for(auto& qp : link->quadraturePoints())
            {
                qp.updateForcesAndVelocities(*link,isClimbingStep);
            }
        }
        //            else
        //            {
        //                linkIter.second.lock()->assembleGlide(false);
        //            }
    }
#endif
    
    Eigen::VectorXd X(Eigen::VectorXd::Zero(0));
    if(isClimbingStep)
    {
        climbSolver->computeClimbScalarVelocities(false);
        size_t k=0;
        for (auto& networkNode : this->networkNodes())
        {
            networkNode.second.lock()->climbVelocityScalar=climbSolver->scalarVelocities()[k];
            ++k;
        }
        X=climbSolver->getNodeVelocities();
    }
    else
    {
        if(glideSolver)
        {
            size_t k=0;
            for (auto& networkNode : this->networkNodes())
            {
                networkNode.second.lock()->climbVelocityScalar.setZero();
                ++k;
            }
            X=glideSolver->getNodeVelocities();
        }
    }
    if(int(NdofXnode*this->networkNodes().size())==X.size())
    {
        size_t k=0;
        for (auto& networkNode : this->networkNodes())
        {
            networkNode.second.lock()->set_V(X.segment(NdofXnode*k,NdofXnode),isClimbingStep); // double cast to remove some numerical noise
            ++k;
        }
    }
    else
    {
        std::cout<<"NdofXnode*this->networkNodes().size()="<<NdofXnode*this->networkNodes().size()<<std::endl;
        std::cout<<"vSolver->getNodeVelocities().size()="<<X.size()<<std::endl;
        throw std::runtime_error("vSolver returned wrong velocity vector size.");
    }
    
    VerboseDislocationNetwork(2,"DislocationNetwork::updateRates"<<std::endl;);
    for(auto& loop : this->loops())
    {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
        loop.second.lock()->updateRates();
    }
    // updatePlasticDistortionRateFromAreas();
    VerboseDislocationNetwork(3,"DislocationNetwork::updateRates DONE"<<std::endl;);
    
//    storeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
    
}

    template <int dim>
    void DislocationNetwork<dim>::solve()
    {
        
        const bool isClimbingStep(isClimbStep());
        if(isClimbingStep)
        {
            std::cout<<" climbStep"<<std::flush;
        }
        else
        {
            std::cout<<" glideStep"<<std::flush;
        }
        double maxVelocity = 0.0;
        for (const auto &nodeIter : this->networkNodes())
        {
            const double vNorm(nodeIter.second.lock()->get_V().norm());
            if (vNorm > maxVelocity)
            {
                maxVelocity = vNorm;
            }
        }
        
    #ifdef _OPENMP
        const size_t nThreads = omp_get_max_threads();
    #else
        const size_t nThreads = 1;
    #endif
        
        //! -1 Compute the interaction StressField between dislocation particles
        std::map<int, int> velocityBinMap{{1,0}};
        for (const auto &binVal : subcyclingBins)
        {
            velocityBinMap.emplace(binVal, 0);
        }
        
        std::cout <<" creating qPoints "<< std::flush;
        for (const auto &links : this->networkLinks())
        {
            
            const int velGroup((subcyclingBins.size() && !isClimbingStep) ? links.second.lock()->velocityGroup(maxVelocity, subcyclingBins) : 1);
            auto velocityBinIter(velocityBinMap.find(velGroup));
            assert(velocityBinIter != velocityBinMap.end());
            velocityBinIter->second++;
            
            if ((ddBase.simulationParameters.runID % velGroup) == 0)
            {
                links.second.lock()->createQuadraturePoints(isClimbingStep);
            }
        }
        
        std::cout <<" ,updating qPoints (" << nThreads << " threads) " << std::flush;
    #ifdef _OPENMP
    {
        // Topology is fixed during field evaluation. Preserve container order
        // and shared ownership, but make OpenMP's indexed access constant-time.
        // Release the snapshot before any subsequent topology changes.
        std::vector<std::shared_ptr<NetworkLinkType>> fieldLinks;
        fieldLinks.reserve(this->networkLinks().size());
        for(const auto& entry : this->networkLinks())
        {
            fieldLinks.push_back(entry.second.lock());
        }
    #pragma omp parallel for
        for (size_t k = 0; k < fieldLinks.size(); ++k)
        {
            const auto& link=fieldLinks[k];
            const int velGroup((subcyclingBins.size() && !isClimbingStep) ? link->velocityGroup(maxVelocity, subcyclingBins) : 1);
            
            if ((ddBase.simulationParameters.runID % velGroup) == 0)
            {
                link->updateQuadraturePoints(isClimbingStep);
            }
            //            else
            //            {
            //                link->assembleGlide(false);
            //            }
        }
    }

    #else
        for (auto &linkIter : this->networkLinks())
        {
            const int velGroup(subcyclingBins.size() ? linkIter.second.lock()->velocityGroup(maxVelocity, subcyclingBins) : 1);
            
            if ((ddBase.simulationParameters.runID % velGroup) == 0)
            {
                linkIter.second.lock()->updateQuadraturePoints(isClimbingStep);
            }
            //            else
            //            {
            //                linkIter.second.lock()->assembleGlide(false);
            //            }
        }
    #endif
        
        Eigen::VectorXd X(Eigen::VectorXd::Zero(0));
        if(isClimbingStep)
        {
            climbSolver->computeClimbScalarVelocities(true);
            size_t k=0;
            for (auto& networkNode : this->networkNodes())
            {
                networkNode.second.lock()->climbVelocityScalar=climbSolver->scalarVelocities()[k];
                ++k;
            }
            X=climbSolver->getNodeVelocities();
        }
        else
        {
            if(glideSolver)
            {
                size_t k=0;
                for (auto& networkNode : this->networkNodes())
                {
                    networkNode.second.lock()->climbVelocityScalar.setZero();
                    ++k;
                }
                X=glideSolver->getNodeVelocities();
            }
        }
        if(int(NdofXnode*this->networkNodes().size())==X.size())
        {
            size_t k=0;
            for (auto& networkNode : this->networkNodes())
            {
                networkNode.second.lock()->set_V(X.segment(NdofXnode*k,NdofXnode),isClimbingStep); // double cast to remove some numerical noise
                ++k;
            }
        }
        else
        {
            std::cout<<"NdofXnode*this->networkNodes().size()="<<NdofXnode*this->networkNodes().size()<<std::endl;
            std::cout<<"vSolver->getNodeVelocities().size()="<<X.size()<<std::endl;
            throw std::runtime_error("vSolver returned wrong velocity vector size.");
        }
        
        VerboseDislocationNetwork(2,"DislocationNetwork::updateRates"<<std::endl;);
        for(auto& loop : this->loops())
        {// copmute slipped areas and right-handed normal // TODO: PARALLELIZE THIS LOOP
            loop.second.lock()->updateRates();
        }
        // updatePlasticDistortionRateFromAreas();
        VerboseDislocationNetwork(3,"DislocationNetwork::updateRates DONE"<<std::endl;);
        
        storeSingleGlideStepDiscreteEvents(ddBase.simulationParameters.runID);
        
    }

    template <int dim>
    void DislocationNetwork<dim>::moveNodes(const double & dt_in)
    {/*! Moves all nodes in the DislocationNetwork using the stored glide velocity and current dt
      */
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"Moving DislocationNodes (dt="<<dt_in<< ")... "<<std::flush;
        // Thermal hazards use this step's accepted time interval, not the
        // previous dt still present during the velocity solve.
        if(crossSlipModel && crossSlipModel->isThermal() && !isClimbStep())
        {
            crossSlipMaker.findCrossSlipSegments(dt_in);
            crossSlipMaker.execute();
        }
        std::vector<std::pair<std::shared_ptr<LoopNodeType>,bool>> savedFlags;
        for(const auto& entry : this->loopNodes())
        {
            const auto node=entry.second.lock();
            savedFlags.emplace_back(node,node->periodicContourVertex);
        }
        const auto savedDangling=danglingBoundaryLoopNodes;
        for(const auto& entry : this->loops())
        {
            std::vector<std::shared_ptr<LoopNodeType>> nodes;
            for(const auto* link : entry.second.lock()->linkSequence())
            {
                nodes.push_back(link->source);
            }
            retainPeriodicContourBends(nodes,1.e-8*(ddBase.mesh.xMax()-ddBase.mesh.xMin()).norm());
        }
        danglingBoundaryLoopNodes.clear();
        std::vector<std::pair<std::shared_ptr<LoopNodeType>,typename LoopNodeType::SplineNodeType>> saved;
        for(const auto& entry : this->loopNodes())
        {
            const auto node=entry.second.lock();
            saved.emplace_back(node,static_cast<const typename LoopNodeType::SplineNodeType&>(*node));
        }
        trialNodeMotion=true;
        try
        {
            // Network positions remain unchanged throughout this loop, so
            // shared nodes cannot be advanced twice through different loops.
            for(const auto& entry : this->networkNodes())
            {
                const auto node=entry.second.lock();
                node->trySet_P(node->get_P()+node->get_V()*dt_in);
            }
            updateBoundaryNodes();
            trialNodeMotion=false;
        }
        catch(...)
        {
            for(const auto& entry : saved)
            {
                static_cast<typename LoopNodeType::SplineNodeType&>(*entry.first)=entry.second;
            }
            for(const auto& entry : savedFlags)
            {
                entry.first->periodicContourVertex=entry.second;
            }
            danglingBoundaryLoopNodes=savedDangling;
            trialNodeMotion=false;
            throw;
        }
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
    }

    template <int dim>
    void DislocationNetwork<dim>::storeSingleGlideStepDiscreteEvents(const long int&)
    {
        if(!crossSlipModel || !crossSlipModel->isThermal())
        {
            crossSlipMaker.findCrossSlipSegments();
            crossSlipMaker.execute();
        }
    }

    template <int dim>
    void DislocationNetwork<dim>::executeSingleGlideStepDiscreteEvents(const long int& runID)
    {
        
        //    crossSlipMaker.execute();
        //        //! 13- Node redistribution
        networkRemesher.remesh(runID);
        //        //! 12- Form Junctions
        junctionsMaker.formJunctions(3.0*networkRemesher.Lmin);
        //Calling remesh again so that any other topological changes created by junctions whihc are otherwise removable can be removed
        //        //! 13- Node redistribution
        //        this->io().output(runID);
        
        networkRemesher.remesh(runID);
        //        updateVirtualBoundaryLoops();
        
        DislocationNucleation<dim> nucleator(*this);
        nucleator.bulkNucleate(bulkNucleationModel);
        nucleator.surfaceNucleate(surfaceNucleationModel);
        
    }

    template <int dim>
    void DislocationNetwork<dim>::output(DefectiveCrystalConfigIO<dim>& configIO,DefectiveCrystalAuxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        
        for(const auto& loop : this->loops())
        {
            configIO.loops().emplace_back(*loop.second.lock());
        }
        
        // Store LoopNodes
        for(const auto& node : this->loopNodes())
        {
            configIO.loopNodes().emplace_back(*node.second.lock());
        }
        
        // Store LoopLinks
        for(const auto& link : this->loopLinks())
        {
            configIO.loopLinks().emplace_back(link.second);
        }
        
        // Store NetworkNodes
        for(const auto& node : this->networkNodes())
        {
            configIO.nodes().emplace_back(*node.second.lock());
        }
        
        // AuxIO
        if (this->outputQuadraturePoints)
        {
            for (const auto& link : this->networkLinks())
            {
                for(const auto& qPoint : link.second.lock()->quadraturePoints())
                {
                    auxIO.quadraturePoints().push_back(qPoint);
                }
            }
        }
        
        const std::tuple<double,double,double,double> length(this->networkLength());
        const double densityFactor(1.0/this->ddBase.mesh.volume()/std::pow(this->ddBase.poly.b_SI,2));
        f_file<<std::get<0>(length)*densityFactor<<" "<<std::get<1>(length)*densityFactor<<" "<<std::get<2>(length)*densityFactor<<" "<<std::get<3>(length)*densityFactor<<" ";
        if(ddBase.simulationParameters.runID==0)
        {
            F_labels<<"glissile density [m^-2]\n";
            F_labels<<"sessile density [m^-2]\n";
            F_labels<<"boundary density [m^-2]\n";
            F_labels<<"grain boundary density [m^-2]\n";
        }
        
        if(outputDislocationDensityPerSlipSystem)
        {
            const auto ssDen(this->networkLengthPerSlipSystem());
            for(const auto& tup : ssDen)
            {
                f_file<<std::get<0>(tup)*densityFactor<<" "<<std::get<1>(tup)*densityFactor<<" "<<std::get<2>(tup)*densityFactor<<" "<<std::get<3>(tup)*densityFactor<<" ";
            }
            
            if(ddBase.simulationParameters.runID==0)
            {
                for(size_t kd=0;kd<ssDen.size();++kd)
                {
                    F_labels<<"SlipSystem_"<<kd<<" glissile density [m^-2]\n";
                    F_labels<<"SlipSystem_"<<kd<<"sessile density [m^-2]\n";
                    F_labels<<"SlipSystem_"<<kd<<"boundary density [m^-2]\n";
                    F_labels<<"SlipSystem_"<<kd<<"grain boundary density [m^-2]\n";
                }
            }
            
        }
        
        if(this->outputPlasticDistortionPerSlipSystem)
        {
            const auto ssapd(slipSystemAveragePlasticDistortion());
            for(const auto& pair : ssapd)
            {
                f_file<<pair.second<<" ";
                if(ddBase.simulationParameters.runID==0)
                {
                    F_labels<<"betaP_g"+std::to_string(pair.first.first)+"_ss"+std::to_string(pair.first.second)+" [-]\n";
                }
            }
        }
        
        if(this->computeElasticEnergyPerLength)
        {
            double eE(0.0);
            double eC(0.0);
            
            for(const auto& linkIter : this->networkLinks())
            {// Collect LoopLinks by loop IDs
                const auto link(linkIter.second.lock());
                for(const auto& qPoint : link->quadraturePoints())
                {
                    eE+=qPoint.elasticEnergyPerLength*qPoint.dL;
                    eC+=qPoint.coreEnergyPerLength*qPoint.dL;
                }
            }
            f_file<<eE<<" "<<eC<<" ";
            if(ddBase.simulationParameters.runID==0)
            {
                F_labels<<"dislocation elastic energy [mu b^3]\n";
                F_labels<<"dislocation core energy [mu b^3]\n";
            }
        }
    }

    template <int dim>
    void DislocationNetwork<dim>::updateBoundaryNodes(const std::function<void()>& checkpoint)
    {
        using NodePtr=std::shared_ptr<LoopNodeType>;
        using NetPtr=std::shared_ptr<NetworkNodeType>;
        using Edge=PeriodicPlaneEdge<dim>;
        using Record=PolygonPatchIntersectionRecord<dim,LoopNodeType>;
        const double tolerance=1.e-8*(ddBase.mesh.xMax()-ddBase.mesh.xMin()).norm();
        struct Arc
        {
            NodePtr source,sink;
            VectorDim chord;
            bool reversed;
            std::set<const NetworkLinkType*> links;
            std::vector<NodePtr> subdivisions;
            size_t group;
        };
        struct Plan
        {
            std::shared_ptr<LoopType> loop;
            std::vector<NodePtr> oldNodes,anchors,newNodes;
            std::vector<size_t> arcs;
            std::vector<Record> contour;
        };
        std::vector<Plan> plans;
        std::vector<Arc> arcs;
        std::set<NodePtr> promoted;
        // Snapshot complete loops first. These strong references also protect
        // borrowed provenance pointers until all loops have been reconciled.
        for(const auto& entry : this->loops())
        {
            const auto loop=entry.second.lock();
            if(!loop->periodicGlidePlane)
            {
                continue;
            }
            Plan plan;
            plan.loop=loop;
            std::vector<size_t> anchorPositions;
            for(const auto* link : loop->linkSequence())
            {
                plan.oldNodes.push_back(link->source);
                // Physical boundary vertices define the retained slipped area;
                // only artificial periodic subdivisions may be omitted.
                if(!link->source->isPeriodicSubdivision())
                {
                    anchorPositions.push_back(plan.oldNodes.size()-1);
                    plan.anchors.push_back(link->source);
                }
            }
            if(plan.anchors.size()<3)
            {
                if(retracesContourSegment(plan.oldNodes,anchorPositions,tolerance))
                {
                    // Disconnect only this loop's memberships. Contraction has
                    // transported derived points along their original segments;
                    // both resulting paths have now been verified to retrace.
                    plans.push_back(std::move(plan));
                    continue;
                }
                // No collapse proof: retain the actual contour, including seam
                // bends. Never infer zero slipped area from the anchor count.
                plan.anchors.clear();
                anchorPositions.clear();
                for(size_t i=0;i<plan.oldNodes.size();++i)
                {
                    if(plan.anchors.empty() ||
                       (plan.oldNodes[i]->get_P()-plan.anchors.back()->get_P()).norm()>tolerance)
                    {
                        promoted.insert(plan.oldNodes[i]);
                        plan.anchors.push_back(plan.oldNodes[i]);
                        anchorPositions.push_back(i);
                    }
                }
                if(plan.anchors.size()>1 &&
                   (plan.anchors.front()->get_P()-plan.anchors.back()->get_P()).norm()<=tolerance)
                {
                    plan.anchors.pop_back();
                    anchorPositions.pop_back();
                }
                if(plan.anchors.size()<3)
                {
                    throw std::runtime_error("Cannot establish a retraced periodic contour.");
                }
            }
            for(size_t i=0; i<plan.anchors.size(); ++i)
            {
                Arc arc;
                arc.source=plan.anchors[i];
                arc.sink=plan.anchors[(i+1)%plan.anchors.size()];
                arc.chord=arc.sink->get_P()-arc.source->get_P();
                arc.reversed=arc.source->networkNode->sID>arc.sink->networkNode->sID;
                if(arc.source->networkNode==arc.sink->networkNode)
                {
                    for(int axis=0; axis<dim; ++axis)
                    {
                        if(std::abs(arc.chord[axis])>tolerance)
                        {
                            arc.reversed=arc.chord[axis]<0;
                            break;
                        }
                    }
                }
                if(arc.reversed)
                {
                    arc.chord=-arc.chord;
                }
                size_t pos=anchorPositions[i];
                const size_t end=anchorPositions[(i+1)%anchorPositions.size()];
                do
                {
                    if(const auto link=plan.oldNodes[pos]->next.second->networkLink())
                    {
                        arc.links.insert(link.get());
                    }
                    pos=(pos+1)%plan.oldNodes.size();
                    if(pos!=end)
                    {
                        arc.subdivisions.push_back(plan.oldNodes[pos]);
                    }
                }
                while(pos!=end);
                arc.group=arcs.size();
                plan.arcs.push_back(arcs.size());
                arcs.push_back(std::move(arc));
            }
            std::vector<std::pair<VectorLowerDim,const LoopNodeType* const>> polygon;
            for(const auto& node : plan.anchors)
            {
                polygon.emplace_back(loop->periodicGlidePlane->referencePlane->localPosition(node->get_P()),node.get());
            }
            plan.contour=loop->periodicGlidePlane->polygonPatchIntersection(polygon,false);
            plans.push_back(std::move(plan));
        }
        // Establish junction lineage from the pre-edit network, not from
        // coincident positions or the endpoints left after trimming.
        const auto root=[&](size_t i)
        {
            while(arcs[i].group!=i)
            {
                i=arcs[i].group;
            }
            return i;
        };
        std::map<std::pair<size_t,size_t>,std::vector<size_t>> endpointArcs;
        for(size_t i=0; i<arcs.size(); ++i)
        {
            const auto key=std::minmax(arcs[i].source->networkNode->sID,arcs[i].sink->networkNode->sID);
            auto& peers=endpointArcs[key];
            for(size_t j : peers)
            {
                if((arcs[i].chord-arcs[j].chord).norm()<=tolerance)
                {
                    for(const auto* link : arcs[i].links)
                    {
                        if(arcs[j].links.count(link))
                        {
                            arcs[root(i)].group=root(j);
                            break;
                        }
                    }
                }
            }
            peers.push_back(i);
        }
        struct Crossing
        {
            double parameter;
            NetPtr node;
        };
        std::map<size_t,std::vector<Crossing>> crossings;
        for(size_t i=0; i<arcs.size(); ++i)
        {
            const auto& arc=arcs[i];
            const VectorDim d=arc.sink->get_P()-arc.source->get_P();
            for(const auto& node : arc.subdivisions)
            {
                if(!danglingBoundaryLoopNodes.count(node.get()) && d.squaredNorm()>0)
                {
                    double t=d.dot(node->get_P()-arc.source->get_P())/d.squaredNorm();
                    crossings[root(i)].push_back({arc.reversed ? 1-t : t,node->networkNode});
                }
            }
        }
        // Select a representative per original junction before editing any
        // positions. Prefer its existing RVE copy when several copies survive.
        std::map<size_t,VectorDim> committedPositions;
        for(const auto& plan : plans)
        {
            for(const auto& record : plan.contour)
            {
                if(const auto index=std::get<6>(record).originalVertex)
                {
                    const auto node=plan.anchors.at(*index)->networkNode;
                    const VectorDim mapped=plan.loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(record))+std::get<1>(record);
                    if(!committedPositions.count(node->sID) || (mapped-node->get_P()).norm()<=tolerance)
                    {
                        committedPositions[node->sID]=mapped;
                    }
                }
            }
        }
        std::map<size_t,std::vector<NetPtr>> remappedOriginals;
        std::map<std::pair<size_t,size_t>,std::vector<NetPtr>> sharedBoundaryCorners;
        std::map<std::pair<std::vector<size_t>,size_t>,std::vector<NetPtr>> enclosedCorners;
        for(auto& plan : plans)
        {
            // Only fully shared old contour lineages may share an enclosed
            // replacement. Equal coordinates alone never establish ownership.
            std::vector<size_t> enclosureOwners;
            const bool enclosed=std::any_of(plan.contour.begin(),plan.contour.end(),[](const Record& record)
            {
                return std::get<6>(record).enclosureLayer.has_value();
            });
            if(enclosed)
            {
                for(size_t arc : plan.arcs)
                {
                    enclosureOwners.push_back(root(arc));
                }
                // Canonicalize only cyclic start and traversal direction, not
                // the order of visits: different sheet connectivity must stay
                // distinct even if it uses the same set of old junction arcs.
                auto best=enclosureOwners;
                for(bool reverse : {false,true})
                {
                    auto candidate=enclosureOwners;
                    if(reverse)
                    {
                        std::reverse(candidate.begin(),candidate.end());
                    }
                    for(size_t start=0;start<candidate.size();++start)
                    {
                        best=std::min(best,candidate);
                        std::rotate(candidate.begin(),candidate.begin()+1,candidate.end());
                    }
                }
                enclosureOwners=std::move(best);
            }
            std::set<const LoopNodeType*> reused;
            for(const auto& record : plan.contour)
            {
                const auto& info=std::get<6>(record);
                const auto patch=plan.loop->periodicGlidePlane->getPatch(std::get<1>(record));
                const VectorDim position=plan.loop->periodicGlidePlane->referencePlane->globalPosition(std::get<0>(record));
                const VectorDim mapped=position+patch->shift;
                std::pair<std::shared_ptr<Edge>,std::shared_ptr<Edge>> edges;
                if(std::get<2>(record).first>=0)
                {
                    edges.first=patch->edges().at(std::get<2>(record).first);
                    if(std::get<2>(record).second>=0)
                    {
                        edges.second=patch->edges().at(std::get<2>(record).second);
                    }
                }
                NodePtr chosen;
                NetPtr networkNode;
                if(info.originalVertex)
                {
                    const auto old=plan.anchors.at(*info.originalVertex);
                    networkNode=old->networkNode;
                    if((committedPositions.at(networkNode->sID)-mapped).norm()>tolerance)
                    {
                        auto& alternatives=remappedOriginals[networkNode->sID];
                        networkNode.reset();
                        for(const auto& other : alternatives)
                        {
                            if((other->get_P()-mapped).norm()<=tolerance)
                            {
                                networkNode=other;
                                break;
                            }
                        }
                        if(!networkNode)
                        {
                            networkNode=this->networkNodes().create(mapped,old->networkNode->get_V(),old->networkNode->climbVelocityScalar,old->networkNode->velocityReduction());
                            alternatives.push_back(networkNode);
                        }
                    }
                    // Recompute confinement: sliding can change which edge
                    // contains an anchor, especially at mixed corners.
                    edges={nullptr,nullptr};
                    for(const auto& edge : patch->edges())
                    {
                        if((edge->meshIntersection->snap(mapped)-mapped).norm()<=tolerance)
                        {
                            if(!edges.first)
                            {
                                edges.first=edge;
                            }
                            else
                            {
                                edges.second=edge;
                            }
                        }
                    }
                    const auto& previous=old->periodicPlaneEdge;
                    if(old->networkNode==networkNode && old->periodicPlanePatch()==patch &&
                       ((previous.first==edges.first && previous.second==edges.second) ||
                        (previous.first==edges.second && previous.second==edges.first)))
                    {
                        chosen=old;
                    }
                }
                else if(info.polygonLocation)
                {
                    const auto& location=*info.polygonLocation;
                    const size_t arcID=plan.arcs.at(location.segment.sourceIndex);
                    const double t=arcs[arcID].reversed ? 1-location.parameter : location.parameter;
                    auto& candidates=crossings[root(arcID)];
                    for(const auto& candidate : candidates)
                    {
                        if(std::abs(candidate.parameter-t)<=1.e-8 && (candidate.node->get_P()-mapped).norm()<=tolerance)
                        {
                            networkNode=candidate.node;
                            break;
                        }
                    }
                    if(!networkNode)
                    {
                        networkNode=this->networkNodes().create(mapped);
                        candidates.push_back({t,networkNode});
                    }
                    for(const auto& old : arcs[arcID].subdivisions)
                    {
                        if(old->networkNode==networkNode && old->periodicPlanePatch()==patch &&
                           (old->get_P()-position).norm()<=tolerance && !danglingBoundaryLoopNodes.count(old.get()))
                        {
                            const auto& oldEdges=old->periodicPlaneEdge;
                            if((oldEdges.first==edges.first && oldEdges.second==edges.second) ||
                               (oldEdges.first==edges.second && oldEdges.second==edges.first))
                            {
                                chosen=old;
                                break;
                            }
                        }
                    }
                }
                if(!networkNode && info.enclosureLayer)
                {
                    auto& candidates=enclosedCorners[{enclosureOwners,*info.enclosureLayer}];
                    for(const auto& candidate : candidates)
                    {
                        if((candidate->get_P()-mapped).norm()<=tolerance)
                        {
                            networkNode=candidate;
                            break;
                        }
                    }
                    if(!networkNode)
                    {
                        networkNode=this->networkNodes().create(mapped);
                        candidates.push_back(networkNode);
                    }
                }
                // A boundary-following corner inherits a shared junction only
                // when BOTH ends of its boundary walk originate on shared arcs.
                // Matching coordinates alone must not join unrelated loops.
                if(!networkNode && info.onPhysicalBoundary && !plan.arcs.empty())
                {
                    const size_t index=&record-plan.contour.data();
                    std::optional<size_t> before,after;
                    for(size_t step=1;step<plan.contour.size();++step)
                    {
                        if(!before)
                        {
                            const auto& metadata=std::get<6>(plan.contour[(index+plan.contour.size()-step)%plan.contour.size()]);
                            if(metadata.polygonLocation)
                            {
                                before=root(plan.arcs.at(metadata.polygonLocation->segment.sourceIndex));
                            }
                        }
                        if(!after)
                        {
                            const auto& metadata=std::get<6>(plan.contour[(index+step)%plan.contour.size()]);
                            if(metadata.polygonLocation)
                            {
                                after=root(plan.arcs.at(metadata.polygonLocation->segment.sourceIndex));
                            }
                        }
                        if(before && after)
                        {
                            break;
                        }
                    }
                    if(before && after)
                    {
                        auto& candidates=sharedBoundaryCorners[std::minmax(*before,*after)];
                        for(const auto& candidate : candidates)
                        {
                            if((candidate->get_P()-mapped).norm()<=tolerance)
                            {
                                networkNode=candidate;
                                break;
                            }
                        }
                        if(!networkNode)
                        {
                            networkNode=this->networkNodes().create(mapped);
                            candidates.push_back(networkNode);
                        }
                    }
                }
                if(!networkNode)
                {
                    networkNode=this->networkNodes().create(mapped);
                }
                if(!chosen || reused.count(chosen.get()))
                {
                    chosen=this->loopNodes().create(plan.loop,networkNode,position,patch,edges);
                }
                if(info.originalVertex && (chosen->periodicPlaneEdge.first || plan.anchors.at(*info.originalVertex)->periodicContourVertex || promoted.count(plan.anchors.at(*info.originalVertex))))
                {
                    promoted.insert(chosen);
                }
                reused.insert(chosen.get());
                plan.newNodes.push_back(chosen);
            }
            if(!plan.newNodes.empty() && plan.newNodes.size()<3)
            {
                throw std::runtime_error("Boundary update produced a degenerate retained contour.");
            }
        }
        // Prepare the edit journal before mutating any committed objects.
        using LinkMap=typename std::remove_reference_t<decltype(this->loopLinks())>;
        using Key=typename LoopLinkType::KeyType;
        std::vector<Key> obsolete,created;
        std::vector<std::tuple<NodePtr,NodePtr,std::shared_ptr<LoopType>>> additions;
        std::vector<typename LinkMap::node_type> removed;
        struct PositionState
        {
            NetPtr node;
            typename NetworkNodeType::SplineNodeType spline;
            const Simplex<dim,dim>* simplex;
        };
        std::vector<PositionState> savedPositions;
        for(const auto& entry : committedPositions)
        {
            const auto node=this->networkNodes().get(entry.first);
            savedPositions.push_back({node,static_cast<const typename NetworkNodeType::SplineNodeType&>(*node),node->p_Simplex});
        }
        for(const auto& plan : plans)
        {
            std::set<Key> desired,replacements;
            for(size_t i=0;i<plan.newNodes.size();++i)
            {
                const auto a=plan.newNodes[i];
                const auto b=plan.newNodes[(i+1)%plan.newNodes.size()];
                const auto key=LoopLinkType::loopLinkKey(plan.loop,a,b);
                if(!desired.insert(key).second)
                {
                    throw std::runtime_error("Boundary contour repeats a loop-link identity.");
                }
                const auto found=this->loopLinks().find(key);
                if(found==this->loopLinks().end() || found->second.source!=a)
                {
                    additions.emplace_back(a,b,plan.loop);
                    replacements.insert(key);
                }
            }
            for(size_t i=0;i<plan.oldNodes.size();++i)
            {
                const auto a=plan.oldNodes[i];
                const auto b=plan.oldNodes[(i+1)%plan.oldNodes.size()];
                const auto key=LoopLinkType::loopLinkKey(plan.loop,a,b);
                if(!desired.count(key) || replacements.count(key))
                {
                    obsolete.push_back(key);
                }
            }
        }
        struct SegmentState
        {
            std::shared_ptr<NetworkLinkType> link;
            VectorDim burgers;
            double norm;
            std::shared_ptr<SlipSystem> slip;
        };
        std::map<const NetworkLinkType*,SegmentState> segmentStates;
        const auto saveSegment=[&](const std::shared_ptr<NetworkLinkType>& link)
        {
            if(link)
            {
                segmentStates.try_emplace(link.get(),SegmentState{link,link->Burgers,link->BurgersNorm,link->_slipSystem});
            }
        };
        for(const auto& key : obsolete)
        {
            saveSegment(this->loopLinks().at(key).networkLink());
        }
        for(const auto& item : additions)
        {
            saveSegment(this->networkLinks().get(NetworkLinkType::getKey(std::get<0>(item)->networkNode,std::get<1>(item)->networkNode)));
        }
        removed.reserve(obsolete.size());
        created.reserve(additions.size());
        try
        {
            for(const auto& entry : committedPositions)
            {
                this->networkNodes().get(entry.first)->set_P(entry.second);
                if(checkpoint)
                {
                    checkpoint();
                }
            }
            for(const auto& key : obsolete)
            {
                auto handle=this->loopLinks().extract(key);
                handle.mapped().detachTopology();
                removed.push_back(std::move(handle));
                if(checkpoint)
                {
                    checkpoint();
                }
            }
            for(const auto& item : additions)
            {
                const auto& [a,b,loop]=item;
                this->connect(a,b,loop);
                created.push_back(LoopLinkType::loopLinkKey(loop,a,b));
                if(checkpoint)
                {
                    checkpoint();
                }
            }
        }
        catch(...)
        {
            for(const auto& key : created)
            {
                this->loopLinks().erase(key);
            }
            for(auto& handle : removed)
            {
                handle.mapped().attachTopology();
                this->loopLinks().insert(std::move(handle));
            }
            for(const auto& entry : savedPositions)
            {
                static_cast<typename NetworkNodeType::SplineNodeType&>(*entry.node)=entry.spline;
                entry.node->p_Simplex=entry.simplex;
            }
            for(const auto& entry : savedPositions)
            {
                entry.node->updateGeometry();
            }
            for(const auto& entry : segmentStates)
            {
                const auto& state=entry.second;
                state.link->Burgers=state.burgers;
                state.link->BurgersNorm=state.norm;
                state.link->_slipSystem=state.slip;
                state.link->straight.updateGeometry();
            }
            throw;
        }
        for(const auto& node : promoted)
        {
            node->periodicContourVertex=true;
        }
        danglingBoundaryLoopNodes.clear();
    }

    template <int dim>
    int DislocationNetwork<dim>::verboseDislocationNetwork=0;

    template class DislocationNetwork<3>;
}
#endif
