/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamics_cpp_
#define model_ClusterDynamics_cpp_

#ifdef _OPENMP
#include <ClusterDynamics/ClusterDynamicsBase.h>
#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <omp.h>
#endif

#include <ClusterDynamics.h>
#include <exception>
//#include <ExternalAndInternalBoundary.h>
//#include <Fix.h>
#include <DislocationNetwork.h>
#include <MicrostructureGenerator.h>

namespace model
{

    //template <int dim>
    //struct BoundaryConcentration
    //{
    //    const IrradiationClimbParameters<dim>& icp;
    //    const Eigen::Matrix<double,dim,dim>& stress;
    //
    //    /**********************************************************************/
    //    BoundaryConcentration(const IrradiationClimbParameters<dim>& icp_in,const Eigen::Matrix<double,dim,dim>& stress_in) :
    //    /* */ icp(icp_in)
    //    /* */,stress(stress_in)
    //    {
    //
    //    }
    //
    //    /**********************************************************************/
    //    template <typename NodeType,int dofPerNode>
    //    Eigen::Matrix<double,dofPerNode,1> operator()(const NodeType& node,
    //                                                  Eigen::Matrix<double,dofPerNode,1>& val) const
    //    {
    //
    //        const auto v1(icp.equilibriumClusterConcentration(stress.trace()));
    //
    //        const auto outNormal(node.outNormal()); // used to compute traction
    //        const double tn = outNormal.dot(stress*outNormal);
    //
    //        val=(v1.array()*exp(-tn*icp.omega/icp.kB/icp.T*icp.speciesVector)).matrix().transpose();
    //
    //        return val;
    //
    //    }
    //
    //};

template <int dim>
typename ClusterDynamics<dim>::UniformControllerContainerType ClusterDynamics<dim>::getUniformControllers(const DefectiveCrystalBase<dim>& ddBase,const ClusterDynamicsParameters<dim>& cdp)
{
    typedef typename UniformControllerType::MatrixVoigt MatrixVoigt;
    typedef typename UniformControllerType::VectorVoigt VectorVoigt;
    
    UniformControllerContainerType temp;
    for(const auto& pair : cdp.D)
    {
        MatrixVoigt C(MatrixVoigt::Zero());
        const auto& grainID(pair.first);
        for(size_t k=0;k<pair.second.size();++k)
        {
            const auto& Dm(pair.second[k]);
            C.template block<dim,dim>(k*dim,k*dim)=Dm;
        }
        
        const auto f0(VectorVoigt::Zero()); // flux
        const auto f0Dot(VectorVoigt::Zero()); // flux rate
        const auto g0(VectorVoigt::Zero()); // Concentration gradient
        const auto g0Dot(VectorVoigt::Zero()); // Concentration gradient rate
        const auto stiffnessRatio(VectorVoigt::Ones()*1.0e20); // impose concentration gradient
        const double t0(0.0);

        temp.emplace(grainID,new UniformControllerType(t0,C,stiffnessRatio,g0,g0Dot,f0,f0Dot));
    }
    return temp;
}


    template<int dim>
    ClusterDynamics<dim>::ClusterDynamics(MicrostructureContainerType& mc) :
    /* init */ MicrostructureBase<dim>("ClusterDynamics",mc)
    /* init */,cdp(this->microstructures.ddBase.clusterDynamics())
//    /* init */,ClusterDynamicsBase<dim>(this->microstructures.ddBase)
    // /* init */,useClusterDynamicsFEM(this->microstructures.ddBase.clusterDynamics().useClusterDynamicsFEM)
    /* init */,useClusterDynamicsFEM(this->microstructures.ddBase.clusterDynamics().useClusterDynamicsFEM)
//    /* init */,useClusterDynamicsFEM(this->microstructures.ddBase.clusterDynamics().useClusterDynamicsFEM)
    /* init */,clusterDynamicsFEM(useClusterDynamicsFEM?  new ClusterDynamicsFEM<dim>(this->microstructures.ddBase,cdp) : nullptr)
    /* init */,uniformControllers(useClusterDynamicsFEM? UniformControllerContainerType() : getUniformControllers(this->microstructures.ddBase,this->cdp))
//    /* init */,mobileClustersIncrement(this->microstructures.ddBase.fe->template trial<'d',mSize>())
//    /* init */,dV(this->microstructures.ddBase.fe->template domain<EntireDomain,dVorder,GaussLegendre>())
//    /* init */,mBWF((test(this->mobileGrad),-this->microstructures.ddBase.poly.Omega*this->mobileFlux)*dV)
//    /* init */,dmBWF((test(grad(mobileClustersIncrement)),-this->microstructures.ddBase.poly.Omega*(FluxMatrix<dim>(this->cdp)*grad(mobileClustersIncrement)))*dV)
//    /* init */,mSolver(true,FLT_EPSILON)
////    /* init */,mSolver(mBWF,true,FLT_EPSILON)
//    /* init */,cascadeGlobalProduction(((test(clusterDynamicsFEM->mobileClusters),make_constant(this->cdp.G))*dV)->globalVector())
    ///* init */,cascadeGlobalProduction(((test(clusterDynamicsFEM->mobileClusters),make_constant(this->cdp.G.transpose().eval()))*dV)->globalVector())
    {
//        mobileClustersIncrement->setConstant(Eigen::Matrix<double,mSize,1>::Zero());
//        std::cout<<cascadeGlobalProduction<<std::endl;
    }

    

    template<int dim>
    void ClusterDynamics<dim>::initializeConfiguration(const DefectiveCrystalConfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        if(clusterDynamicsFEM)
        {
            clusterDynamicsFEM->initializeConfiguration(configIO,f_file,F_labels);
        }
        else
        {// nothing to initialize for uniformControllers
            
        }
    }

template<int dim>
void ClusterDynamics<dim>::applyBoundaryConditions()
{
 
    auto& fem=*clusterDynamicsFEM;
    Eigen::VectorXd offset=Eigen::VectorXd::Zero(fem.mobileClusters->gSize());
    std::map<size_t,double> prescribed;
    std::set<size_t> periodicNodes;
    for(const auto& pair : fem.mobileClusters->fe().periodicNodePairs())
    {
        periodicNodes.insert(pair.source);
        periodicNodes.insert(pair.sink);
    }
    const auto dn=this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>();
    const bool useInfluence=solvedQuantityBoundaryCache && dn && dn->climbSolver;
    Eigen::MatrixXd discreteConcentration;
    if(useInfluence)
    {
        if(solvedQuantityInfluence.empty())
        {
            std::map<size_t,size_t> columns;
            for(const auto& entry : dn->networkNodes())
            {
                columns.emplace(entry.first,solvedQuantitySourceIDs.size());
                solvedQuantitySourceIDs.push_back(entry.first);
            }
            for(const auto& node : fem.mobileClusters->fe().nodes())
            {
                if(cdp.computeReactions || fem.sinkNormals.count(node.gID) || periodicNodes.count(node.gID))
                {
                    solvedQuantityTargetRows.emplace(node.gID,solvedQuantityTargetRows.size());
                }
            }
            solvedQuantityInfluence.assign(mSize,Eigen::MatrixXd::Zero(solvedQuantityTargetRows.size(),columns.size()));
            std::vector<const NodeType*> targets(solvedQuantityTargetRows.size());
            for(const auto& node : fem.mobileClusters->fe().nodes())
            {
                const auto row=solvedQuantityTargetRows.find(node.gID);
                if(row!=solvedQuantityTargetRows.end())
                {
                    targets[row->second]=&node;
                }
            }
            std::vector<std::exception_ptr> failures(targets.size());
            // Each worker writes a distinct H row, preserving segment summation
            // order and the shared network-node columns within that row.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(targets.size());++row)
            {
                try
                {
                    const auto& node=*targets[row];
                    const auto pointGrains=dn->pointGrains(node.P0,&node,nullptr,nullptr);
                    for(const auto& entry : dn->networkLinks())
                    {
                        const auto segment=entry.second.lock();
                        if(segment->hasZeroBurgers())
                        {
                            continue;
                        }
                        const auto grains=segment->grains();
                        size_t common=0;
                        for(const auto* grain : grains)
                        {
                            common+=pointGrains.count(grain);
                        }
                        if(common==1)
                        {
                            const auto coefficients=segment->concentrationMatrices(node.P0,dn->climbSolver->CD->cdp);
                            for(int species=0;species<mSize;++species)
                            {
                                solvedQuantityInfluence[species](row,columns.at(segment->source->sID))+=coefficients(species,0);
                                solvedQuantityInfluence[species](row,columns.at(segment->sink->sID))+=coefficients(species,1);
                            }
                        }
                    }
                }
                catch(...)
                {
                    failures[row]=std::current_exception();
                }
            }
            for(const auto& failure : failures)
            {
                if(failure)
                {
                    solvedQuantityInfluence.clear();
                    solvedQuantitySourceIDs.clear();
                    solvedQuantityTargetRows.clear();
                    std::rethrow_exception(failure);
                }
            }
            ++solvedQuantityInfluenceBuilds;
        }
        if(dn->networkNodes().size()!=solvedQuantitySourceIDs.size())
        {
            throw std::runtime_error("DD topology changed during concentration coupling");
        }
        Eigen::MatrixXd velocities(solvedQuantitySourceIDs.size(),mSize);
        for(size_t column=0;column<solvedQuantitySourceIDs.size();++column)
        {
            const auto found=dn->networkNodes().find(solvedQuantitySourceIDs[column]);
            if(found==dn->networkNodes().end() || found->second.expired())
            {
                throw std::runtime_error("DD source changed during concentration coupling");
            }
            velocities.row(column)=found->second.lock()->climbVelocityScalar.matrix();
        }
        discreteConcentration.resize(solvedQuantityTargetRows.size(),mSize);
        for(int species=0;species<mSize;++species)
        {
            discreteConcentration.col(species).noalias()=solvedQuantityInfluence[species]*velocities.col(species);
        }
    }
    // Evaluate only uncached sink stresses in parallel. Commit to the map on
    // the calling thread after every evaluation succeeds.
    std::vector<const NodeType*> sinkTargets;
    for(const auto& node : fem.mobileClusters->fe().nodes())
    {
        if(fem.sinkNormals.count(node.gID) &&
           (!solvedQuantityBoundaryCache || !solvedQuantitySinkConcentrations.count(node.gID)))
        {
            sinkTargets.push_back(&node);
        }
    }
    std::vector<MatrixDim> sinkStresses(sinkTargets.size());
    std::vector<std::exception_ptr> sinkFailures(sinkTargets.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(!sinkTargets.empty())
#endif
    for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(sinkTargets.size());++row)
    {
        try
        {
            const auto& node=*sinkTargets[row];
            sinkStresses[row]=this->microstructures.stress(node.P0,&node,nullptr,nullptr);
        }
        catch(...)
        {
            sinkFailures[row]=std::current_exception();
        }
    }
    for(const auto& failure : sinkFailures)
    {
        if(failure)
        {
            std::rethrow_exception(failure);
        }
    }
    std::map<size_t,std::vector<size_t>> sinkClasses;
    for(size_t row=0;row<sinkTargets.size();++row)
    {
        sinkClasses[fem.mobileClusters->fe().periodicRepresentative(sinkTargets[row]->gID)].push_back(row);
    }
    for(const auto& entry : sinkClasses)
    {
        MatrixDim sigma=MatrixDim::Zero();
        for(const size_t row : entry.second)
        {
            sigma+=sinkStresses[row]/double(entry.second.size());
        }
        const auto& normals=fem.sinkNormals.at(sinkTargets[entry.second.front()]->gID);
        double normalTraction=0.0;
        for(const auto& normal : normals)
        {
            normalTraction+=normal.dot(sigma*normal)/double(normals.size());
        }
        const VectorMSize concentration=cdp.boundaryMobileConcentration(sigma.trace(),normalTraction);
        for(const size_t row : entry.second)
        {
            solvedQuantitySinkConcentrations[sinkTargets[row]->gID]=concentration;
        }
    }
    for(const auto& node : fem.mobileClusters->fe().nodes())
    {
        const auto sink=fem.sinkNormals.find(node.gID);
        // Nonlinear positivity checks require the known concentration throughout the domain.
        if(!cdp.computeReactions && sink==fem.sinkNormals.end() && !periodicNodes.count(node.gID))
        {
            continue;
        }
        VectorMSize otherConcentration=VectorMSize::Zero();
        for(const auto& microstructure : this->microstructures)
        {
            if(microstructure.get()!=static_cast<const MicrostructureBase<dim>*>(this))
            {
                if(useInfluence && microstructure.get()==dn.get())
                {
                    otherConcentration+=discreteConcentration.row(solvedQuantityTargetRows.at(node.gID)).transpose();
                }
                else
                {
                    otherConcentration+=microstructure->mobileConcentration(node.P0,&node,nullptr,nullptr);
                }
            }
        }
        offset.template segment<mSize>(mSize*node.gID)=-otherConcentration;
        if(sink!=fem.sinkNormals.end())
        {
            const auto& concentration=solvedQuantitySinkConcentrations.at(node.gID);
            for(int c=0;c<mSize;++c)
            {
                prescribed.emplace(mSize*node.gID+c,concentration(c)-otherConcentration(c));
            }
        }
    }
    fem.mobileClusters->setPeriodicNodalOffset(offset);
    fem.mobileClusters->dirichletConditions()=std::move(prescribed);
}

    template<int dim>
    void ClusterDynamics<dim>::registerSolvedQuantities(SolvedQuantityRegistry& registry)
    {
        SolvedQuantityContribution concentration;
        concentration.provider=this->tag;
        concentration.name="MobileConcentration";
        concentration.quantities={"MobileConcentration"};
        concentration.dependencies={{"Stress","",false,{}},{"MobileConcentration","",false,{}}};
        concentration.update=[this](const SolvedQuantityContext& context)
        {
            if(context.iteration==0)
            {
                solvedQuantitySinkConcentrations.clear();
                solvedQuantityInfluence.clear();
                solvedQuantitySourceIDs.clear();
                solvedQuantityTargetRows.clear();
            }
            solvedQuantityBoundaryCache=true;
            try
            {
                if(clusterDynamicsFEM)
                {
                    if(!clusterDynamicsFEM->solverInitialized)
                    {
                        clusterDynamicsFEM->initializeSolver();
                    }
                    applyBoundaryConditions();
                    const auto dn=this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>();
                    clusterDynamicsFEM->solveMobileClusters(dn && !dn->loops().empty());
                }
            }
            catch(...)
            {
                solvedQuantityBoundaryCache=false;
                throw;
            }
            solvedQuantityBoundaryCache=false;
        };
        concentration.state=[this]() -> Eigen::VectorXd
        {
            if(clusterDynamicsFEM)
            {
                Eigen::VectorXd result=clusterDynamicsFEM->mobileClusters->dofVector();
                const auto equilibrium=cdp.equilibriumMobileConcentration(0.0);
                for(Eigen::Index k=0;k<result.size();++k)
                {
                    result(k)-=equilibrium(k%mSize);
                }
                return result;
            }
            return Eigen::VectorXd();
        };
        concentration.restore=[this](const Eigen::VectorXd& state)
        {
            if(clusterDynamicsFEM)
            {
                Eigen::VectorXd values=state;
                const auto equilibrium=cdp.equilibriumMobileConcentration(0.0);
                for(Eigen::Index k=0;k<values.size();++k)
                {
                    values(k)+=equilibrium(k%mSize);
                }
                clusterDynamicsFEM->mobileClusters->setFullDofVector(values);
            }
        };
        TextFileParser input(this->microstructures.ddBase.simulationParameters.traitsIO.dcFile);
        const double fraction=input.readStringVector("solvedQuantityConcentrationAbsoluteFraction").empty()? 1.e-10 :
            input.readScalar<double>("solvedQuantityConcentrationAbsoluteFraction",true);
        if(!std::isfinite(fraction) || fraction<=0.0)
        {
            throw std::runtime_error("solvedQuantityConcentrationAbsoluteFraction must be finite and positive.");
        }
        concentration.absoluteTolerance=[this,fraction]() -> Eigen::VectorXd
        {
            if(!clusterDynamicsFEM)
            {
                return Eigen::VectorXd();
            }
            Eigen::VectorXd tolerance(clusterDynamicsFEM->mobileClusters->gSize());
            const auto equilibrium=cdp.equilibriumMobileConcentration(0.0);
            for(Eigen::Index k=0;k<tolerance.size();++k)
            {
                tolerance(k)=fraction*std::max(std::abs(equilibrium(k%mSize)),1.e-100);
            }
            return tolerance;
        };
        registry.add(std::move(concentration));
        SolvedQuantityContribution rates;
        rates.provider=this->tag;
        rates.name="ImmobileRates";
        rates.quantities={"ImmobileRates"};
        rates.dependencies={{"MobileConcentration","",false,{}}};
        rates.update=[this](const SolvedQuantityContext&)
        {
            if(clusterDynamicsFEM)
            {
                const auto dn=this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>();
                if(iSize>0 && (!dn || dn->loops().empty()))
                {
                    clusterDynamicsFEM->solveImmobileClusters();
                }
                else
                {
                    clusterDynamicsFEM->immobileClusterRate.setZero();
                }
            }
        };
        registry.add(std::move(rates));
    }

    template<int dim>
    void ClusterDynamics<dim>::solve()
    {
        if(clusterDynamicsFEM)
        {
            if(!clusterDynamicsFEM->solverInitialized)
            {
                clusterDynamicsFEM->initializeSolver();
            }
            std::cout<<", mobile BCs"<<std::flush;
            applyBoundaryConditions();

            auto DN(this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>());
            const bool hasDiscreteLoops(DN? (DN->loops().size()>0? true : false) : false);
            clusterDynamicsFEM->solve(hasDiscreteLoops); 
        }
        else
        {

        }
    }


    template<int dim>
    void ClusterDynamics<dim>::updateConfiguration()
    {
        if(clusterDynamicsFEM)
        {
            auto DN(this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>());
            const bool hasDiscreteLoops(DN? (DN->loops().size()>0 ? true : false) : false);
            const double dt(this->microstructures.getDt());

            if constexpr (iSize > 0)
            {
                if(!hasDiscreteLoops)
                {
                    clusterDynamicsFEM->updateImmobileClusters(dt);
                }

                // if(cdp.discretizationTime <= this->microstructures.ddBase.simulationParameters.totalTime && !hasDiscreteLoops)
                // {
                //     std::cout << " initializing discrete climb loops, " << std::endl;
                //     initializeDiscreteClimbLoops();
                // }
            }
        }
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
    }

    template<int dim>
    double ClusterDynamics<dim>::getDt() const
    {
        return this->microstructures.ddBase.simulationParameters.dtMax;
    }

    template<int dim>
    void ClusterDynamics<dim>::output(DefectiveCrystalConfigIO<dim>& configIO,DefectiveCrystalAuxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        if(clusterDynamicsFEM)
        {
            const size_t nNodes(clusterDynamicsFEM->mobileClusters->fe().nodes().size());
            configIO.cdMatrix().resize(nNodes,mSize+iSize);
            configIO.cdMatrix().block(0,0,nNodes,mSize)=clusterDynamicsFEM->mobileClusters->dofVector().reshaped(mSize,nNodes).transpose();
            configIO.cdMatrix().block(0,mSize,nNodes,iSize)=clusterDynamicsFEM->immobileClusters->dofVector().reshaped(iSize,nNodes).transpose();
        }
        else
        {
            
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorDim ClusterDynamics<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(clusterDynamicsFEM)
        {
            return clusterDynamicsFEM->inelasticDisplacementRate(x,node,ele,guess);
        }
        else
        {
            return VectorDim::Zero();
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averagePlasticDistortion() const
    {
        if(clusterDynamicsFEM)
        {
            return clusterDynamicsFEM->averageBetaP()+clusterDynamicsFEM->averageBetaV();
        }
        else
        {
            return MatrixDim::Zero();
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averagePlasticDistortionRate() const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorDim ClusterDynamics<dim>::displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::stress(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorMSize ClusterDynamics<dim>::mobileConcentration(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(uniformControllers.size())
        {
            const auto pointGrains(this->pointGrains(x,node,ele,guess));
            VectorMSize averageVal(VectorMSize::Zero());
            for(const auto& grain : pointGrains)
            {
                const auto grainID(grain->grainID);
                const auto& controller(uniformControllers.at(grainID));
                averageVal += controller->grad(this->microstructures.ddBase.simulationParameters.totalTime).reshaped(dim,mSize).transpose()*x;
            }
            return averageVal/pointGrains.size();
        }
        else
        {
            if(node)
            {
                return (*eval(clusterDynamicsFEM->mobileClusters))(*node);
            }
            else
            {
                if(ele)
                {
                    return (*eval(clusterDynamicsFEM->mobileClusters))(*ele,ele->simplex.pos2bary(x));
                }
                else
                {
                    return (*eval(clusterDynamicsFEM->mobileClusters))(x,guess);
                }
            }
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorISize ClusterDynamics<dim>::immobileClusters(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(node)
        {
            return (*eval(clusterDynamicsFEM->immobileClusters))(*node);
        }
        else
        {
            if(ele)
            {
                return (*eval(clusterDynamicsFEM->immobileClusters))(*ele,ele->simplex.pos2bary(x));
            }
            else
            {
                return (*eval(clusterDynamicsFEM->immobileClusters))(x,guess);
            }
        }
    }


    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averageStress() const
    {
        return MatrixDim::Zero();
    }


    // template<int dim>
    // std::set<const Simplex<dim,dim>*> ClusterDynamics<dim>::vertexSetNeighbors(const std::set<const Simplex<dim,dim>*>& inSet, const std::set<const Simplex<dim,dim>*>& usedEle) const
    // {
    //     std::set<const Simplex<dim,dim>*> outSet;
    //     for(const auto& tet : inSet)
    //     {
    //         const auto neighbors(tet->vertexNeighbors());
    //         for(const auto& neighbor : neighbors)
    //         {
    //             if(usedEle.find(neighbor)==usedEle.end())
    //             {
    //                 outSet.insert(neighbor);
    //             }
    //         }
    //     }
    //     return outSet;
    // }

    // template<int dim>
    // std::pair<double,double> ClusterDynamics<dim>::groupNi(const std::set<const Simplex<dim,dim>*>& inSet, const int &k) const
    // {
    //     double nIT(0.0);
    //     double cIT(0.0);
    //     for(const auto& tet : inSet)
    //     {
    //         const auto& ele(this->microstructures.ddBase.fe->elements().at(tet->xID));
    //         Eigen::Matrix<double,dim+1,1> bary(Eigen::Matrix<double,dim+1,1>::Constant(0.25));
    //         const Eigen::Matrix<double,iSize,1> sinkValue((*eval(clusterDynamicsFEM->immobileClusters))(ele,bary));
    //         // ADDED A CORRECTED SINK VALUE HERE -> SOME BARYCENTERS HAVE NEGATIVE VALUES, BUT NODES WERE NOT NEGATIVE
    //         Eigen::Matrix<double,iSize,1> correctedSinkValue(sinkValue);
    //         for(int dof=0;dof<iSize/2;dof++)
    //         { // Corrected sink strength is to adjust for negative values at the barycenter of the quadratic elements where we have a sharp gradient of sink values between two end nodes
    //             const double cmin = cdp.n_min(dof)*cdp.omega*correctedSinkValue(dof);
    //             if(correctedSinkValue(iSize/2+dof)<cmin)
    //             { // defect size cannot be smaller than a minimal value
    //                 correctedSinkValue(iSize/2+dof)=cmin;
    //             }
    //         }
    //         const auto cI(correctedSinkValue.template block<iSize/2,1>(iSize/2,0));
    //         const auto nI(correctedSinkValue.template block<iSize/2,1>(0,0));
    //         // const auto cI(sinkValue.template block<iSize/2,1>(iSize/2,0));
    //         // const auto nI(sinkValue.template block<iSize/2,1>(0,0));
    //         nIT+=nI(k)*tet->vol0;
    //         cIT+=cI(k)*tet->vol0/cdp.omega;
    //     }
    //     return std::make_pair(nIT,cIT);
    // }

    // template<int dim>
    // void ClusterDynamics<dim>::initializeDiscreteClimbLoops()
    // {
    //     if constexpr (iSize == 0)
    //     {
    //         return;
    //     }
    //     else 
    //     {
    //         std::default_random_engine generator;
    //         auto DN(this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>());

    //         for(int k=0; k < iSize/2; ++k)
    //         { //Iteration species-wise
    //             // for(const auto& ele : this->microstructures.ddBase.fe->elements())
    //             // {
    //             //     Eigen::Matrix<double,dim+1,1> bary(Eigen::Matrix<double,dim+1,1>::Constant(0.25));
    //             //     const Eigen::Matrix<double,iSize,1> sinkValue((*eval(clusterDynamicsFEM->immobileClusters))(ele.second,bary));
    //             //     const auto cI(sinkValue.template block<iSize/2,1>(iSize/2,0));
    //             //     if(cI(0)>0)
    //             //     {
    //             //         std::cout<< "ClusterDynamics: Negative concentration for species " << k << " in element " << ele.first << ". Concentration: " << cI(k) << std::endl;
    //             //     }
    //             // }
    //             double tetSize(0.0);
    //             std::multimap<double,const ElementType* const> CiMap; //Create Map
    //             for(const auto& ele : this->microstructures.ddBase.fe->elements())
    //             {
    //                 Eigen::Matrix<double,dim+1,1> bary(Eigen::Matrix<double,dim+1,1>::Constant(0.25));
    //                 const Eigen::Matrix<double,iSize,1> sinkValue((*eval(clusterDynamicsFEM->immobileClusters))(ele.second,bary));
    //                 const auto cI(sinkValue.template block<iSize/2,1>(iSize/2,0));
    //                 const auto cT(cI*ele.second.simplex.vol0/cdp.omega);
    //                 CiMap.emplace(cT(k),&ele.second);
    //             }
    //             std::set<const Simplex<dim,dim>*> usedEle;
    //             for(typename std::multimap<double,const ElementType* const>::reverse_iterator rIter=CiMap.rbegin(); rIter!=CiMap.rend(); rIter++)
    //             {
    //                 const ElementType* const ele(rIter->second);
    //                 if(usedEle.find(&ele->simplex)==usedEle.end())
    //                 {
    //                     std::set<const Simplex<dim,dim>*> currentTets;
    //                     currentTets.insert(&(ele->simplex));
    //                     auto NiCiK(groupNi(currentTets,k));
    //                     size_t currentTetsOldSize(currentTets.size());
    //                     while(NiCiK.first<cdp.discretizationFactor)
    //                     { //If not enough defects to create a loop, grab neighboring elements
    //                         currentTets=vertexSetNeighbors(currentTets,usedEle);
    //                         NiCiK=groupNi(currentTets,k);
    //                         if(currentTets.size()==currentTetsOldSize)
    //                         {
    //                             break;
    //                         }
    //                             else
    //                         {
    //                             currentTetsOldSize=currentTets.size();
    //                         }
    //                     }

    //                     for(const auto &tet : currentTets)
    //                     {
    //                         usedEle.insert(tet);
    //                     }
    //                     tetSize+=currentTets.size();

    //                     const double radius(sqrt(NiCiK.second*cdp.omega/(M_PI*cdp.b*cdp.immobileSpeciesBurgersMagnitude(k))));
    //                     Eigen::Matrix<double,dim+1,1> bary(Eigen::Matrix<double,dim+1,1>::Constant(0.25));
    //                     const auto &grain(this->microstructures.ddBase.poly.grains.at(ele->simplex.region->regionID));
    //                     std::uniform_int_distribution<int> tetDistribution(0,currentTets.size()-1);
    //                     auto tetIter(currentTets.begin());
    //                     std::advance(tetIter,tetDistribution(generator));
    //                     const VectorDim b(grain->latticeBasis*cdp.immobileSpeciesBurgers.col(k));
    //                     const ReciprocalLatticeDirection<dim> n(grain->reciprocalLatticeDirection(b));
    //                     VectorDim tetCenter((*tetIter)->bary2pos(bary).transpose());

    //                     bool insertedLoop(false);
    //                     while(!insertedLoop)
    //                     {
    //                         const long int planeIndex(n.closestPlaneIndexOfPoint(tetCenter));
    //                         GlidePlaneKey<dim> glidePlaneKey(planeIndex,n);
    //                         std::shared_ptr<PeriodicGlidePlane<dim>> periodicGlidePlane(this->microstructures.ddBase.dislocationDynamics().periodicGlidePlaneFactory.get(glidePlaneKey));
    //                         const VectorDim loopCenter(periodicGlidePlane->referencePlane->snapToPlane(tetCenter));
    //                         std::vector<VectorDim> pts;
    //                         if(radius > cdp.minimumLoopSize)
    //                         {
    //                             const int loopSides(12);
    //                             for(int j=0; j<loopSides; ++j)
    //                             {
    //                                 Eigen::Matrix<double,2,1> pLocal(std::cos(2.0*M_PI/loopSides*j),std::sin(2.0*M_PI/loopSides*j));
    //                                 Eigen::Matrix<double,3,1> pGlobal(periodicGlidePlane->referencePlane->globalPosition(radius*pLocal) - periodicGlidePlane->referencePlane->P + loopCenter);
    //                                 if(this->microstructures.ddBase.poly.mesh.searchRegionWithGuess(pGlobal,&ele->simplex).first)
    //                                 {
    //                                     pts.push_back(pGlobal);
    //                                 }
    //                                 else
    //                                 {
    //                                     break;
    //                                 }
    //                             }  
    //                             if(pts.size()==loopSides)
    //                             {
    //                                 insertedLoop=true;
    //                                 VectorDim nA(VectorDim::Zero());
    //                                 for(size_t k=0; k<pts.size(); ++k)
    //                                 {
    //                                     nA += 0.5 * (pts[k]-pts.front()).cross((pts[(k+1)%pts.size()]-pts[k]));
    //                                 }
    //                                 const double trBetaLoop(-b.dot(nA.normalized()));
    //                                 const int signedb(cdp.immobileSpeciesVector(k)/abs(cdp.immobileSpeciesVector(k))); //negative for vacancy, positive for interstitial loops
    //                                 const VectorDim loopBurgers((trBetaLoop * signedb > 0)? b : -b);

    //                                 auto newLoop(DN->loops().create(loopBurgers, this->microstructures.ddBase.dislocationDynamics().glidePlaneFactory.getFromKey(glidePlaneKey)));
    //                                 std::vector<std::shared_ptr<DislocationLoopNode<dim>>> loopNodes;
    //                                 for(const auto &pGlobal:pts)
    //                                 {
    //                                     auto newNetNode(DN->networkNodes().create(pGlobal));
    //                                     const auto periodicPatch(newLoop->periodicGlidePlane? newLoop->periodicGlidePlane->patches().getFromKey(VectorDim::Zero()) : nullptr);
    //                                     const auto periodicPlaneEdge(std::make_pair(nullptr,nullptr));
    //                                     loopNodes.push_back(DN->loopNodes().create(newLoop,newNetNode,pGlobal,periodicPatch,periodicPlaneEdge));
    //                                 }
    //                                 // const int bdn(loopBurgers.dot(nA.normalized()));
    //                                 // std::cout<<"Species: "<<k<<"b dot n: "<< bdn <<std::endl;
    //                                 DN->insertLoop(newLoop, loopNodes);
    //                             }                                                      
    //                             else
    //                             {
    //                                 // std::cout<< "Some Nodes outside of Mesh: "<<std::endl;
    //                                 const VectorDim regionCenter(grain->region.regionCenter());
    //                                 tetCenter+=(0.1)*(regionCenter-tetCenter);
    //                             }
    //                         }
    //                     }
    //                 }   
    //             }
    //         }
    //         // std::cout<<"Removing values for continuous immobile clusters"<<std::endl;
    //         for(size_t nodeID=0; nodeID<clusterDynamicsFEM->immobileClusters->nodeSize(); nodeID++)
    //         {
    //             for(int dof=0; dof<iSize; dof++)
    //             {
    //             //    clusterDynamicsFEM->immobileClusters->dofVector()(iSize*nodeID+dof)=1e-33;
    //             TrialBase<typename ClusterDynamicsFEM<dim>::ImmobileTrialType>::dofVector()(iSize*nodeID+dof)=1e-33;
    //             }
    //         }

    //     }
        
    // }



    template struct ClusterDynamics<3>;
}
#endif


//template<int dim>
//void ClusterDynamics<dim>::solveDiffusiveDisplacement()
//{
//    std::cout<<"Solving diffusiveDisplacementRate..."<<std::flush;
//    const auto t0= std::chrono::system_clock::now();
//    diffusiveDisplacementRate=Eigen::VectorXd::Zero(this->diffusiveDisplacement.gSize());
//    for(const auto& node: this->microstructures.ddBase.fe->nodes())
//    {
//        const Eigen::Matrix<double,dim*mSize,1> speciesFlux((*eval(this->mobileFlux))(node));
//        Eigen::Matrix<double,dim,1> netFlux(Eigen::Matrix<double,dim,1>::Zero());
//        for(int i=0; i<mSize; i++)
//        {
//            const int mSgn(this->cdp.msVector(i)/std::abs(this->cdp.msVector(i)));
//            netFlux+= speciesFlux.template block<dim,1>(i*dim,0)*mSgn;
//        }
//        diffusiveDisplacementRate.template segment<dim>(dim*node.gID)=netFlux;
//    }
//    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//
//}
