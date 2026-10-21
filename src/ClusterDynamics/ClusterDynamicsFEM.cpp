/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsFEM_cpp_
#define model_ClusterDynamicsFEM_cpp_

#include <cmath>
#include <ClusterDynamicsFEM.h>
#include <Eigen/SVD>
#include <Fix.h>

namespace model
{

    template<int dim>
    FluxMatrix<dim>::FluxMatrix(const ClusterDynamicsParameters<dim>& cdp_in) :
    /* init */cdp(cdp_in)
    {
        
    }

    template<int dim>
    const typename FluxMatrix<dim>::MatrixType FluxMatrix<dim>::operator() (const ElementType& elem, const BaryType&) const
    {

        const size_t grainID(elem.simplex.region->regionID);
        MatrixType D(MatrixType::Zero());
        for(size_t k=0;k<mSize;++k)
        {
            D.template block<dim,dim>(k*dim,k*dim)=-cdp.D.at(grainID)[k]/cdp.omega;
        }
        
        // for(size_t k=0;k<mSize;++k)
        // {
        //     std::cout
        //         << " species " << k
        //         << " D trace = "
        //         << cdp.D.at(grainID)[k].trace()
        //         << std::endl;
        // }

        return D;
    }

    template struct FluxMatrix<3>;

template<int dim>
InvDscaling<dim>::InvDscaling(const ClusterDynamicsParameters<dim>& cdp_in) :
/* init */cdp(cdp_in)
{
    
}

template<int dim>
const typename InvDscaling<dim>::MatrixType InvDscaling<dim>::operator() (const ElementType& elem, const BaryType&) const
{
    const size_t grainID(elem.simplex.region->regionID);
    MatrixType InvDscaling(MatrixType::Zero());
    for(size_t k=0;k<mSize;++k)
    {
        InvDscaling(k,k)=3.0/(cdp.D.at(grainID)[k].trace()/cdp.omega);
//        InvDscaling(k,k)=(cdp.D.at(grainID)[k].trace()/cdp.omega)/3.0;
//        InvDscaling(k,k)=1.0;

    }
    return InvDscaling;
}

template struct InvDscaling<3>;


    template<int dim>
    ClusterDynamicsFEM<dim>::ClusterDynamicsFEM(DefectiveCrystalBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in) :
    /* init */ ddBase(ddBase_in)
    /* init */,cdp(cdp_in)
    /* init */,iDs(std::make_shared<InvDscaling<dim>>(cdp))
    /* init */,mobileClusters(ddBase.fe->template trial<'m',mSize>())
    /* init */,mobileGrad(grad(mobileClusters))
    /* init */,mobileFlux(std::make_shared<FluxMatrix<dim>>(cdp)*mobileGrad)
    /* init */,immobileClusters(ddBase.fe->template trial<'i',iSize>())
    /* init */,sinkNormals(clusterBoundarySinks(*ddBase.fe))
    /* init */,mobileClustersIncrement(ddBase.fe->template trial<'d',mSize>())
    /* init */,dV(ddBase.fe->template domain<EntireDomain,dVorder,GaussLegendre>())
//    /* init */,mBWF((test(this->mobileGrad),-ddBase.poly.Omega*this->mobileFlux)*dV)
    /* init */,mBWF((test(grad(iDs*mobileClusters)),-ddBase.poly.Omega*this->mobileFlux)*dV)
    /* init */,dmBWF((test(grad(iDs*mobileClustersIncrement)),-ddBase.poly.Omega*(std::make_shared<FluxMatrix<dim>>(this->cdp)*grad(mobileClustersIncrement)))*dV)
    /* init */,immobileRateTrial(ddBase.fe->template trial<'r',iSize>())
    /* init */,solverInitialized(false)
    /* init */,cascadeGlobalProduction(((test(iDs*this->mobileClusters),make_constant(this->cdp.G))*dV)->globalVector())
    // /* init */,cascadeGlobalProduction(((test(this->mobileClusters),make_constant(this->cdp.G))*dV)->globalVector())

    {
        mobileClustersIncrement->setConstant(Eigen::Matrix<double,mSize,1>::Zero());
        mobileClusters->setConstant(cdp.equilibriumMobileConcentration(0.0).matrix().transpose());
        // immobileClusters->setConstant(Eigen::Matrix<double,iSize,1>::Zero());
        immobileClusters->setConstant(cdp.initLoopSinks);
        immobileClusterRate = Eigen::VectorXd::Zero(immobileClusters->gSize());
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::enforceMobileConstraints()
    {
        auto values=mobileClusters->dofVector();
        std::map<size_t,double> fixed;
        const auto& offset=mobileClusters->periodicNodalOffset();
        const auto lift=[&](size_t node,int component)
        {
            const size_t root=ddBase.fe->periodicRepresentative(node);
            return offset.size()? offset(node*mSize+component)-offset(root*mSize+component) : 0.0;
        };
        for(const auto& bc : mobileClusters->dirichletConditions())
        {
            const size_t node=bc.first/mSize;
            const int component=bc.first%mSize;
            const size_t root=ddBase.fe->periodicRepresentative(node)*mSize+component;
            const double value=bc.second-lift(node,component);
            const auto inserted=fixed.emplace(root,value);
            if(!inserted.second && std::abs(inserted.first->second-value)>1.e-8*std::max({1.e-50,std::abs(value),std::abs(inserted.first->second)}))
            {
                throw std::runtime_error("ClusterDynamics: incompatible concentrations on periodic sink nodes.");
            }
        }
        for(size_t node=0;node<ddBase.fe->nodeSize();++node)
        {
            for(int c=0;c<mSize;++c)
            {
                const size_t root=ddBase.fe->periodicRepresentative(node)*mSize+c;
                values(node*mSize+c)=(fixed.count(root)? fixed.at(root) : mobileClusters->dofVector()(root))+lift(node,c);
            }
        }
        mobileClusters->setFullDofVector(values);
    }

    template<int dim>
    Eigen::MatrixXd ClusterDynamicsFEM<dim>::conservedRows(const bool useImmobileClusters) const
    {
        if(!sinkNormals.empty())
        {
            return Eigen::MatrixXd(0,mobileClusters->gSize());
        }
        std::vector<Eigen::Matrix<double,mSize,1>> directions;
        const auto append=[&](const Eigen::Matrix<double,mSize,mSize>& matrix)
        {
            for(int c=0;c<mSize;++c)
            {
                if(matrix.col(c).norm()>0.0)
                {
                    directions.push_back(matrix.col(c).normalized());
                }
            }
        };
        if(cdp.computeReactions)
        {
            append(cdp.R1cd);
            for(int j=0;j<mSize;++j)
            {
                Eigen::Matrix<double,mSize,mSize> matrix;
                for(int row=0;row<mSize;++row)
                {
                    for(int k=0;k<mSize;++k)
                    {
                        matrix(row,k)=std::abs(cdp.msVector(row))*cdp.R2[row](j,k)/std::abs(cdp.msVector(j)*cdp.msVector(k));
                    }
                }
                append(matrix);
            }
            if(useImmobileClusters)
            {
                FirstOrderReaction<MobileTrialType,ImmobileTrialType> sink(immobileClusters,cdp,ddBase.poly);
                for(const auto& entry : ddBase.fe->elements())
                {
                    // Include every interpolation node, not only element centers.
                    for(int n=0;n<ElementType::nodesPerElement;++n)
                    {
                        append(sink(entry.second,ElementType::baryStarsAndBars.row(n).template cast<double>().transpose()/ElementType::order));
                    }
                }
            }
        }
        Eigen::MatrixXd invariants=Eigen::MatrixXd::Identity(mSize,mSize);
        if(!directions.empty())
        {
            Eigen::MatrixXd matrix(mSize,directions.size());
            for(size_t k=0;k<directions.size();++k)
            {
                matrix.col(k)=directions[k];
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix,Eigen::ComputeFullU);
            svd.setThreshold(1.e-10);
            invariants=svd.matrixU().rightCols(mSize-svd.rank());
        }
        const Eigen::VectorXd production=cdp.G;
        if((invariants.transpose()*production).norm()>1.e-10*production.norm())
        {
            throw std::runtime_error("ClusterDynamics: no steady state in the closed periodic domain; production changes a conserved species inventory.");
        }
        return invariants.transpose()*meanRows;
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solveMobileClusters(const bool hasDiscreteLoops)
    {
        std::cout<<", mobile solver, "<<std::flush;
        const Eigen::VectorXd saved=mobileClusters->dofVector();
        const Eigen::VectorXd savedIncrement=mobileClustersIncrement->dofVector();
        try
        {
            enforceMobileConstraints();
            const bool useImmobileClusters(!hasDiscreteLoops && iSize>0);
            const Eigen::MatrixXd conserved=conservedRows(useImmobileClusters);
            const Eigen::VectorXd inventory=conserved*mobileClusters->dofVector();
            if(!cdp.computeReactions)
            {
                if(mSolver->factorizationCount()==0)
                {
                    mSolver->setLinearConstraints(conserved);
                }
                mSolver->setLinearConstraintValues(inventory);
                mSolver->solveWithMasterSlaveConstraints();
                return;
            }
            // Do not solve singular diffusion alone in a closed periodic domain.
            const auto R1=cdp.R1cd.eval();
            auto sink=std::make_shared<FirstOrderReaction<MobileTrialType,ImmobileTrialType>>(immobileClusters,cdp,ddBase.poly);
            auto R2=std::make_shared<SecondOrderReaction<MobileTrialType>>(mobileClusters,cdp);
            auto j1=(test(iDs*mobileClustersIncrement),R1*(-1.0*mobileClustersIncrement))*dV;
            auto j2=(test(iDs*mobileClustersIncrement),R2*(-1.0*mobileClustersIncrement))*dV;
            auto js=(test(iDs*mobileClustersIncrement),sink*(-1.0*mobileClustersIncrement))*dV;
            auto f1=(test(iDs*mobileClustersIncrement),eval(R1*mobileClusters))*dV;
            auto f2=(test(iDs*mobileClustersIncrement),eval(R2*(0.5*mobileClusters)))*dV;
            auto fs=(test(iDs*mobileClustersIncrement),eval(sink*mobileClusters))*dV;
            const auto residual=[&]() -> Eigen::VectorXd
            {
                Eigen::VectorXd result=cascadeGlobalProduction-diffusionMatrix*mobileClusters->dofVector()+f1->globalVector()+f2->globalVector();
                if(useImmobileClusters)
                {
                    result+=fs->globalVector();
                }
                return result;
            };
            const auto reducedNorm=[&](const Eigen::VectorXd& value)
            {
                Eigen::VectorXd reduced=Eigen::VectorXd::Zero(value.size());
                std::set<size_t> fixed;
                for(const auto& bc : mobileClusters->dirichletConditions())
                {
                    fixed.insert(ddBase.fe->periodicRepresentative(bc.first/mSize)*mSize+bc.first%mSize);
                }
                for(size_t n=0;n<ddBase.fe->nodeSize();++n)
                {
                    for(int c=0;c<mSize;++c)
                    {
                        const size_t root=ddBase.fe->periodicRepresentative(n)*mSize+c;
                        if(!fixed.count(root))
                        {
                            reduced(root)+=value(n*mSize+c);
                        }
                    }
                }
                return reduced.norm();
            };
            rSolver->setLinearConstraints(conserved);
            double reference=0.0;
            for(int iteration=0;iteration<60;++iteration)
            {
                reactionLoad->value=residual();
                const double error=reducedNorm(reactionLoad->value);
                if(iteration==0)
                {
                    reference=std::max({error,reducedNorm(cascadeGlobalProduction),1.e-50});
                }
                if(error<=1.e-8*reference)
                {
                    return;
                }
                reactionMatrix->entries=diffusionTriplets;
                const auto append=[&](const auto& form)
                {
                    const auto entries=form->globalTriplets();
                    reactionMatrix->entries.insert(reactionMatrix->entries.end(),entries.begin(),entries.end());
                };
                append(j1);
                append(j2);
                if(useImmobileClusters)
                {
                    append(js);
                }
                rSolver->setLinearConstraintValues(inventory-conserved*mobileClusters->dofVector());
                rSolver->reassembleMatrix();
                rSolver->solveWithMasterSlaveConstraints();
                const Eigen::VectorXd old=mobileClusters->dofVector();
                bool accepted=false;
                for(int backtrack=0;backtrack<24;++backtrack)
                {
                    const double step=std::ldexp(1.0,-backtrack);
                    mobileClusters->setFullDofVector(old+step*mobileClustersIncrement->dofVector());
                    Eigen::VectorXd total=mobileClusters->dofVector();
                    if(mobileClusters->periodicNodalOffset().size())
                    {
                        total-=mobileClusters->periodicNodalOffset();
                    }
                    if(total.minCoeff()< -1.e-12*std::max(1.e-50,total.cwiseAbs().maxCoeff()))
                    {
                        continue;
                    }
                    const double trialError=reducedNorm(residual());
                    if(std::isfinite(trialError) && (trialError<error || trialError<=1.e-8*reference))
                    {
                        accepted=true;
                        break;
                    }
                }
                if(!accepted)
                {
                    throw std::runtime_error("ClusterDynamics: reaction iteration failed to reduce the residual.");
                }
            }
            throw std::runtime_error("ClusterDynamics: reaction iteration exceeded 60 steps.");
        }
        catch(...)
        {
            mobileClusters->setFullDofVector(saved);
            mobileClustersIncrement->setFullDofVector(savedIncrement);
            throw;
        }
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solveImmobileClusters()
    {
        if constexpr(iSize>0)
        {
            auto rate=std::make_shared<ImmobileSinkRate<MobileTrialType,ImmobileTrialType>>(mobileClusters,immobileClusters,cdp,ddBase.poly);
            auto rhs=(test(immobileRateTrial),rate)*dV;
            if(!iSolver)
            {
                auto mass=(test(immobileRateTrial),immobileRateTrial)*dV;
                WeakProblemOptions options;
                options.useDirectSolver=false;
                iSolver=weakProblem<Eigen::SparseLU<SparseMatrixType>,Eigen::BiCGSTAB<SparseMatrixType>>(mass,rhs,options);
            }
            else
            {
                iSolver->setRightHandSide(rhs);
            }
            iSolver->solveWithMasterSlaveConstraints();
            immobileClusterRate=immobileRateTrial->dofVector();
        }
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::updateImmobileClusters(const double dt)
    {
        immobileClusters->dofVector()+=immobileClusterRate * dt;
        for(size_t nodeID=0; nodeID<immobileClusters->nodeSize(); nodeID++)
        {
            for(int dof=0;dof<iSize/2;dof++)
            { // Corrected sink strength is to adjust for negative values at the barycenter of the quadratic elements where we have a sharp gradient of sink values between two end nodes
                // rmin = r_min * omega * N
                const double cmin = this->cdp.n_min(dof)*this->cdp.omega*immobileClusters->dofVector()(iSize*nodeID+dof);
                if(immobileClusters->dofVector()(iSize*nodeID+iSize/2+dof)<cmin) //  Nc, Na1, Na2, Na3, c, ca1, ca2, ca3
                { // Loop size cannot be smaller than a minimal value
                    immobileClusters->dofVector()(iSize*nodeID+iSize/2+dof)=cmin;
                }
            }
        }
        for(size_t node=0;node<ddBase.fe->nodeSize();++node)
        {
            const size_t root=ddBase.fe->periodicRepresentative(node);
            if(root!=node)
            {
                immobileClusters->dofVector().template segment<iSize>(iSize*node)=
                    immobileClusters->dofVector().template segment<iSize>(iSize*root);
            }
        }
    }


    template<int dim>
    void ClusterDynamicsFEM<dim>::solve(const bool hasDiscreteLoops)
    {
        const Eigen::VectorXd mobile=mobileClusters->dofVector();
        const Eigen::VectorXd increment=mobileClustersIncrement->dofVector();
        const Eigen::VectorXd rate=immobileClusterRate;
        const Eigen::VectorXd rateTrial=immobileRateTrial->dofVector();
        try
        {
            solveMobileClusters(hasDiscreteLoops);
            if(!hasDiscreteLoops && iSize>0)
            {
                solveImmobileClusters();
            }
            else
            {
                immobileClusterRate.setZero();
            }
        }
        catch(...)
        {
            mobileClusters->setFullDofVector(mobile);
            mobileClustersIncrement->setFullDofVector(increment);
            immobileClusterRate=rate;
            immobileRateTrial->setFullDofVector(rateTrial);
            throw;
        }
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::initializeSolver()
    {
        for(const auto& entry : sinkNormals)
        {
            for(int c=0;c<mSize;++c)
            {
                mobileClusters->dirichletConditions()[entry.first*mSize+c]=0.0;
                mobileClustersIncrement->dirichletConditions()[entry.first*mSize+c]=0.0;
            }
        }
        diffusionTriplets=dmBWF->globalTriplets();
        diffusionMatrix.resize(mobileClusters->gSize(),mobileClusters->gSize());
        diffusionMatrix.setFromTriplets(diffusionTriplets.begin(),diffusionTriplets.end());
        mobileLoad=std::make_shared<Load<MobileTrialType>>(Load<MobileTrialType>{mobileClusters,cascadeGlobalProduction});
        reactionLoad=std::make_shared<Load<MobileIncrementTrialType>>(Load<MobileIncrementTrialType>{mobileClustersIncrement,Eigen::VectorXd::Zero(mobileClusters->gSize())});
        reactionMatrix=std::make_shared<MatrixForm<MobileIncrementTrialType>>(MatrixForm<MobileIncrementTrialType>{mobileClustersIncrement,diffusionTriplets});
        mSolver=weakProblem<Eigen::SparseLU<SparseMatrixType>,Eigen::BiCGSTAB<SparseMatrixType>>(mBWF,mobileLoad);
        WeakProblemOptions reactionOptions;
        // Preserve iterative reaction solves on surface/GB-bounded domains.
        // Closed domains may need an indefinite conservation-constrained solve.
        reactionOptions.useDirectSolver=sinkNormals.empty();
        rSolver=weakProblem<Eigen::SparseLU<SparseMatrixType>,Eigen::BiCGSTAB<SparseMatrixType>>(reactionMatrix,reactionLoad,reactionOptions);
        meanRows=Eigen::MatrixXd::Zero(mSize,mobileClusters->gSize());
        for(int c=0;c<mSize;++c)
        {
            Eigen::Matrix<double,mSize,1> unit=Eigen::Matrix<double,mSize,1>::Zero();
            unit(c)=1.0;
            meanRows.row(c)=((test(mobileClusters),make_constant(unit))*dV)->globalVector().transpose()/ddBase.mesh.volume();
        }
        solverInitialized=true;
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::initializeConfiguration(const DefectiveCrystalConfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {    
        if(size_t(configIO.cdMatrix().size())==mobileClusters->gSize()+immobileClusters->gSize())
        {
            const size_t nNodes(mobileClusters->fe().nodes().size());
            const Eigen::VectorXd mobile=configIO.cdMatrix().block(0,0,nNodes,mSize).transpose().reshaped(mobileClusters->gSize(),1);
            const Eigen::VectorXd immobile=configIO.cdMatrix().block(0,mSize,nNodes,iSize).transpose().reshaped(immobileClusters->gSize(),1);
            for(const auto& pair : ddBase.fe->periodicNodePairs())
            {
                const auto a=immobile.template segment<iSize>(iSize*pair.source);
                const auto b=immobile.template segment<iSize>(iSize*pair.sink);
                if((a-b).norm()>1.e-8*std::max({1.e-50,a.norm(),b.norm()}))
                {
                    throw std::runtime_error("ClusterDynamics: restart immobile concentrations violate active periodic pairs.");
                }
            }
            mobileClusters->setFullDofVector(mobile);
            immobileClusters->setFullDofVector(immobile);

        }
        else
        {
            if(configIO.cdMatrix().size())
            {
                throw std::runtime_error("ClusterDynamics: TrialFunctions initializatoin size mismatch");
            }
        }        
    }


    template<int dim>
    typename ClusterDynamicsFEM<dim>::VectorDim ClusterDynamicsFEM<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        Eigen::Matrix<double,dim*mSize,1> speciesFlux(Eigen::Matrix<double,dim*mSize,1>::Zero());
        if(node)
        {
            speciesFlux=(*eval(this->mobileFlux))(*node);
        }
        else
        {
            if(ele)
            {
                speciesFlux=(*eval(this->mobileFlux))(*ele,ele->simplex.pos2bary(x));
            }
            else
            {
                speciesFlux=(*eval(this->mobileFlux))(x,guess);
            }
        }
        Eigen::Matrix<double,dim,1> netFlux(Eigen::Matrix<double,dim,1>::Zero());
        for(int i=0; i<mSize; i++)
        {
            const int mSgn(this->cdp.msVector(i)/std::abs(this->cdp.msVector(i)));
            netFlux += speciesFlux.template block<dim,1>(i*dim,0)*mSgn;
        }
        return netFlux;
    }


    template<int dim>
    std::vector<Eigen::Matrix<double,dim,dim>> ClusterDynamicsFEM<dim>::betaP(const ElementType& ele, const BaryType& bary) const
    {
        std::vector<Eigen::Matrix<double,dim,dim>> strainVector;
        const Eigen::Matrix<double,iSize,1> sinkValue((*eval(immobileClusters))(ele,bary));
        Eigen::Matrix<double,iSize,1> correctedSinkValue(sinkValue);
        for(int dof=0;dof<iSize/2;dof++)
        {
            const double cmin = cdp.n_min(dof)*cdp.omega*correctedSinkValue(dof);
            if(correctedSinkValue(iSize/2+dof)<cmin)
            { // Loop size cannot be smaller than a minimal value
                correctedSinkValue(iSize/2+dof)=cmin;
            }
        }
        const size_t gID(ele.simplex.region->regionID);
        const Eigen::Matrix<double,dim,dim> lat(ddBase.poly.grains.at(gID)->latticeBasis);
        const Eigen::Matrix<double,dim,iSize/2> localBurgers(lat*cdp.immobileSpeciesBurgers);
        const Eigen::Array<double,1,iSize/2> cr(cdp.clusterRadius(correctedSinkValue.template block<iSize/2,1>(iSize/2,0),correctedSinkValue.template block<iSize/2,1>(0,0)));
        const Eigen::Matrix<double,1,dim> localx(lat.col(0).normalized());
        const Eigen::Matrix<double,1,dim> localy((2.0*lat.col(1)-lat.col(0)).normalized());
        const Eigen::Matrix<double,1,dim> localz(lat.col(2).normalized());
        for(int k=0; k<iSize/2; k++)
        { // N A b ^b*^b
            Eigen::Matrix<double,dim,dim> tempStrain(Eigen::Matrix<double,dim,dim>::Zero());
            const double Area=M_PI*std::pow(cr(k),2);
            const Eigen::Matrix<double,dim,dim> highValue(cdp.immobileSpeciesVector(k)*correctedSinkValue(k)*Area*localBurgers.col(k)*localBurgers.col(k).normalized().transpose());
            const Eigen::Matrix<double,dim,dim> lowValue(cdp.immobileSpeciesVector(k)*correctedSinkValue(k)*cdp.delVPyramid*(localx.transpose()*localx+localy.transpose()*localy+localz.transpose()*localz)/3.0);
            tempStrain+=cdp.sigmoidalMatrixInterpolation(correctedSinkValue.template block<iSize/2,1>(iSize/2,0),correctedSinkValue.template block<iSize/2,1>(0,0),lowValue,highValue,k).matrix();
            strainVector.push_back(tempStrain);
        }
        return strainVector;
    }


    template<int dim>
    std::vector<Eigen::Matrix<double,dim,dim>> ClusterDynamicsFEM<dim>::betaV(const ElementType& ele, const BaryType& bary) const
    {
        std::vector<Eigen::Matrix<double,dim,dim>> strainVector;
        const Eigen::Matrix<double,iSize,1> sinkValue((*eval(immobileClusters))(ele,bary));
        Eigen::Matrix<double,iSize,1> correctedSinkValue(sinkValue);
        for(int dof=0;dof<iSize/2;dof++)
        {
        const double cmin = cdp.n_min(dof)*cdp.omega*correctedSinkValue(dof);
            if(correctedSinkValue(iSize/2+dof)<cmin)
            { // Loop size cannot be smaller than a minimal valu
                correctedSinkValue(iSize/2+dof)=cmin;
            }
        }
        const size_t gID(ele.simplex.region->regionID);
        const Eigen::Matrix<double,dim,dim> lat(ddBase.poly.grains.at(gID)->latticeBasis);
        const Eigen::Matrix<double,dim,iSize/2> localBurgers(lat*cdp.immobileSpeciesBurgers);
        const Eigen::Matrix<double,1,dim> localx(lat.col(0).normalized());
        const Eigen::Matrix<double,1,dim> localy((2.0*lat.col(1)-lat.col(0)).normalized());
        const Eigen::Matrix<double,1,dim> localz(lat.col(2).normalized());
        const Eigen::Array<double,1,iSize/2> cr(cdp.clusterRadius(correctedSinkValue.template block<iSize/2,1>(iSize/2,0),correctedSinkValue.template block<iSize/2,1>(0,0)));
        for(int k=0; k<iSize/2; k++)
        { // N A b ^b*^b
            Eigen::Matrix<double,dim,dim> tempStrain(Eigen::Matrix<double,dim,dim>::Zero());
            const double Area=M_PI*std::pow(cr(k),2);
            const Eigen::Matrix<double,dim,dim> highValue(cdp.immobileSpeciesRelRelaxVol(k)*cdp.immobileSpeciesVector(k)*correctedSinkValue(k)*Area*localBurgers.col(k)*localBurgers.col(k).normalized().transpose());
            const Eigen::Matrix<double,dim,dim> lowValue(cdp.immobileSpeciesRelRelaxVol(k)*cdp.immobileSpeciesVector(k)*correctedSinkValue(k)*cdp.delVPyramid*(localx.transpose()*localx+localy.transpose()*localy+localz.transpose()*localz)/3.0);
            tempStrain+=cdp.sigmoidalMatrixInterpolation(correctedSinkValue.template block<iSize/2,1>(iSize/2,0),correctedSinkValue.template block<iSize/2,1>(0,0),lowValue,highValue,k).matrix();
            strainVector.push_back(tempStrain);
        }
       // Add up the contribution from the mobile species
        Eigen::Matrix<double,dim,dim> tempStrain(Eigen::Matrix<double,dim,dim>::Zero());
        const Eigen::Matrix<double,mSize,1> inFactors((1.0/cdp.msVector.abs()).matrix());
        const Eigen::Matrix<double,mSize,1> localConcentration(inFactors.asDiagonal()*(*eval(mobileClusters))(ele,bary)); // Local concentration
        for(int k=0; k<mSize; k++)
        {
            tempStrain += 1.0/3.0*cdp.msVector(k)*localConcentration(k)*cdp.msRelRelaxVol(k)*Eigen::Matrix<double,dim,dim>::Identity();
        }
        strainVector.push_back(tempStrain);
        return strainVector;
    }
    

    /* plastic distortion */
    template<int dim>
    typename ClusterDynamicsFEM<dim>::MatrixDim ClusterDynamicsFEM<dim>::averageBetaPKernel(const Eigen::Matrix<double,dim,1>& a1, const ElementType& ele) const
    {
        const auto bary(BarycentricTraits<dim>::x2l(a1));

        const std::vector<Eigen::Matrix<double,dim,dim>> bps(betaP(ele,bary));
        Eigen::Matrix<double,dim,dim> temp(Eigen::Matrix<double,dim,dim>::Zero());  
        for(const Eigen::Matrix<double,dim,dim>& bp : bps)
        {
            temp += bp*ele.absJ(bary);
        }
        return temp; 
        // return betaP(ele,bary)[0]*ele.absJ(bary);
    }
    template<int dim>
    typename ClusterDynamicsFEM<dim>::MatrixDim ClusterDynamicsFEM<dim>::averageBetaP() const
    {
        Eigen::Matrix<double,dim,dim> temp(Eigen::Matrix<double,dim,dim>::Zero());
        // dV.integrate(this,temp,&ClusterDynamicsFEM<dim>::averageBetaPKernel);    // integrate through the volume
        for(size_t k=0; k<dV.size(); ++k)
        {
            VolumeIntegrationDomainType::QuadratureType::integrate(this,temp,&ClusterDynamicsFEM<dim>::averageBetaPKernel,dV.element(k));
        }
        return temp/ddBase.poly.mesh.volume();
    }

    /* volumetric distortion */   
    template<int dim>
    typename ClusterDynamicsFEM<dim>::MatrixDim ClusterDynamicsFEM<dim>::averageBetaVKernel(const Eigen::Matrix<double,dim,1>& a1, const ElementType& ele) const
    {
        const auto bary(BarycentricTraits<dim>::x2l(a1));

        const std::vector<Eigen::Matrix<double,dim,dim>> bvs(betaV(ele,bary));
        Eigen::Matrix<double,dim,dim> temp(Eigen::Matrix<double,dim,dim>::Zero());  
        for(const Eigen::Matrix<double,dim,dim>& bv : bvs)
        {
            temp += bv*ele.absJ(bary);
        }
        return temp; 
        // return betaP(ele,bary)[0]*ele.absJ(bary);
    }
    template<int dim>
    typename ClusterDynamicsFEM<dim>::MatrixDim ClusterDynamicsFEM<dim>::averageBetaV() const
    {
        Eigen::Matrix<double,dim,dim> temp(Eigen::Matrix<double,dim,dim>::Zero());
        // dV.integrate(this,temp,&ClusterDynamicsFEM<dim>::averageBetaVKernel);    // integrate through the volume
        for(size_t k=0; k<dV.size(); ++k)
        {
            VolumeIntegrationDomainType::QuadratureType::integrate(this,temp,&ClusterDynamicsFEM<dim>::averageBetaVKernel,dV.element(k));
        }
        return temp/ddBase.poly.mesh.volume();
    }


    template struct ClusterDynamicsFEM<3>;

}
#endif



//template<int dim>
//void ClusterDynamicsFEM<dim>::applyBoundaryConditions()
//{
//
//    const auto& nodesInternalExternal(mobileClusters->fe().nodeList(nodeListInternalExternal));
//
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//    for(size_t k=0;k<nodesInternalExternal.size();++k)
//    {
//        const auto& node(nodesInternalExternal[k]);
//        const auto outNormal(node->outNormal()); // used to compute traction
//        const MatrixDim sigma(microstructures.stress(node->P0,node,nullptr,nullptr));
//        const double normalTraction(outNormal.dot(sigma*outNormal));
//        const auto bndConcentration(this->cdp.boundaryMobileConcentration(sigma.trace(),normalTraction));
//
//        VectorMSize otherConcentration(VectorMSize::Zero());
//        for(const auto& microstructure : this->microstructures)
//        {
//            if(microstructure.get()!=static_cast<const MicrostructureBase<dim>* const>(this))
//            {// not the ClusterDynamics physics
//                otherConcentration += microstructure->mobileConcentration(node->P0,node,nullptr,nullptr);
//            }
//        }
//        for(int k=0;k<mSize;++k)
//        {
//            mobileClusters->dirichletConditions().at(mSize*node->gID+k) = bndConcentration(k) - otherConcentration(k);
//        }
//    }
//}
