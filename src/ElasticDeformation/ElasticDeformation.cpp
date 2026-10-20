/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ElasticDeformation_cpp_
#define model_ElasticDeformation_cpp_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ElasticDeformation.h>

namespace model
{

    bool kd(const int& i,const int& j)
    {
        return i==j;
    }

    template <int dim>
    std::unique_ptr<UniformController<SymmetricVoigtTraits<dim>::voigtSize>> getUniformEDcontroller(const DislocationDynamicsBase<dim>& ddBase)
    {
        
        typedef UniformController<SymmetricVoigtTraits<dim>::voigtSize> UniformControllerType;
        typedef typename UniformControllerType::MatrixVoigt MatrixVoigt;
        typedef typename UniformControllerType::MatrixVoigt MatrixVoigt;
        
        // Build stiffness matrix
        const double mu(ddBase.poly.mu);
        const double nu(ddBase.poly.nu);
        const double lam(2.0*nu*mu/(1.0-2.0*nu));
        MatrixVoigt C(MatrixVoigt::Zero()); // C_ijkl=lam*d_ij*d_kl+mu*(d_ik*d_jl+d_il*d_jk)
        for(int vI=0;vI<SymmetricVoigtTraits<dim>::voigtSize;++vI)
        {
            const auto& i(ddBase.voigtTraits.tensorIndex(vI,0));
            const auto& j(ddBase.voigtTraits.tensorIndex(vI,1));
            for(int vJ=0;vJ<SymmetricVoigtTraits<dim>::voigtSize;++vJ)
            {
                const auto& k(ddBase.voigtTraits.tensorIndex(vJ,0));
                const auto& l(ddBase.voigtTraits.tensorIndex(vJ,1));
                C(vI,vJ)= lam*kd(i,j)*kd(k,l)+mu*(kd(i,k)*kd(j,l)+kd(i,l)*kd(j,k));
            }
        }
        
        TextFileParser parser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/ElasticDeformation.txt");
        
        const auto f0(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStress0",true));
        const auto f0Dot(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStressRate",true));
        const auto g0(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStrain0",true));
        const auto g0Dot(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStrainRate",true));
        const auto stiffnessRatio(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("stiffnessRatio",true));
        const double t0(0.0);
        return std::unique_ptr<UniformControllerType>(new UniformControllerType(t0,C,stiffnessRatio,g0,g0Dot,f0,f0Dot));
    }

    template<int dim>
    ElasticDeformation<dim>::ElasticDeformation(MicrostructureContainerType& mc) :
    /* init */ MicrostructureBase<dim>("ElasticDeformation",mc)
    /* init */,useElasticDeformationFEM(this->microstructures.ddBase.fe? bool(TextFileParser(this->microstructures.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/ElasticDeformation.txt").readScalar<int>("useElasticDeformationFEM",true)) : false )
    /* init */,elasticDeformationFEM(useElasticDeformationFEM?  new ElasticDeformationFEM<dim>(this->microstructures.ddBase) : nullptr)
    /* init */,uniformLoadController(useElasticDeformationFEM? std::unique_ptr<UniformControllerType>(nullptr) : getUniformEDcontroller(this->microstructures.ddBase))
    /* init */,inertiaReliefPenaltyFactor(uniformLoadController? 0.0 : TextFileParser(this->microstructures.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/ElasticDeformation.txt").readScalar<double>("inertiaReliefPenaltyFactor",true))
    /* init */,solverInitialized(false)
    {
        if(useElasticDeformationFEM)
        {
            updatePeriodicGradient();
        }
    }

    template<int dim>
    void ElasticDeformation<dim>::updatePeriodicGradient()
    {
        // Free periodic jumps are solved jointly with the nodal field. Reset the
        // prescribed offset before solving; the solved gradient is stored afterwards.
        elasticDeformationFEM->u->setPeriodicGradient(MatrixDim::Zero());
    }

    template<int dim>
    void ElasticDeformation<dim>::initializeSolver()
    {
        auto settings=getUniformEDcontroller(this->microstructures.ddBase);
        femLoadController=std::make_shared<FEMControllerType>(elasticDeformationFEM->u,*settings);
        auto penalty=std::make_shared<GlobalRigidMotionPenalty<typename ElasticDeformationFEM<dim>::TrialFunctionType>>(elasticDeformationFEM->u,inertiaReliefPenaltyFactor,true);
        auto total=elasticDeformationFEM->bWF+penalty+femLoadController;
        auto rhs=femLoadController->rightHandSide(Eigen::VectorXd::Zero(elasticDeformationFEM->u->gSize()),0.0,FEMControllerType::Vector6::Zero());
        weakProblemSolver=weakProblem<DirectSolverType,IterativeSolverType>(total,rhs);
        weakProblemSolver->setAffineModes(femLoadController->affineModes());
        weakProblemSolver->setLinearConstraints(femLoadController->constraintRows());
        solverInitialized=true;
    }

    template<int dim>
    void ElasticDeformation<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        if(uniformLoadController)
        {
            const MatrixDim apd(this->microstructures.averagePlasticDistortion());
            const MatrixDim aps(0.5*(apd+apd.transpose()));
            uniformLoadController->gs=this->microstructures.ddBase.voigtTraits.m2v(aps,true);
        }
        else
        {
            updatePeriodicGradient();
            if(size_t(configIO.displacementMatrix().size())==elasticDeformationFEM->u->gSize()+elasticDeformationFEM->z->gSize())
            {
                std::cout<<" from configIO,"<<std::flush;
                const size_t nNodes(elasticDeformationFEM->u->fe().nodes().size());
                *elasticDeformationFEM->u=configIO.displacementMatrix().block(0,0,nNodes,dim).transpose().reshaped(elasticDeformationFEM->u->gSize(),1);
                *elasticDeformationFEM->z=configIO.displacementMatrix().block(0,dim,nNodes,dim).transpose().reshaped(elasticDeformationFEM->z->gSize(),1);
            }
            else
            {
                if(configIO.displacementMatrix().size())
                {
                    throw std::runtime_error("ElasticDeformation: TrialFunction u initializatoin size mismatch");
                }
            }
            zDot.setZero(elasticDeformationFEM->z->gSize());
        }
    }

template<int dim>
void ElasticDeformation<dim>::reSolve()
{
    if(uniformLoadController)
    {
//        const MatrixDim apd(this->microstructures.averagePlasticDistortion());
//        const MatrixDim aps(0.5*(apd+apd.transpose()));
//        uniformLoadController->gs=this->microstructures.ddBase.voigtTraits.m2v(aps,true);
    }
    else
    {
        std::cout<<", zDot "<<std::flush;
        for(const auto& node : elasticDeformationFEM->z->fe().nodes())
        {
            zDot.template segment<dim>(dim*node.gID)=this->microstructures.inelasticDisplacementRate(node.P0,&node,nullptr,nullptr);
        }
    }
}

    template<int dim>
    void ElasticDeformation<dim>::solve()
    {
        if(uniformLoadController)
        {
            const MatrixDim apd(this->microstructures.averagePlasticDistortion());
            const MatrixDim aps(0.5*(apd+apd.transpose()));
            uniformLoadController->gs=this->microstructures.ddBase.voigtTraits.m2v(aps,true);
        }
        else
        {
            updatePeriodicGradient();
            if(!solverInitialized)
            {
                initializeSolver();
                
            }
            const auto& base=this->microstructures.ddBase;
            Eigen::VectorXd source=Eigen::VectorXd::Zero(elasticDeformationFEM->u->gSize());
            typename FEMControllerType::Vector6 otherStress=FEMControllerType::Vector6::Zero();
            using Q=Quadrature<dim,4,GaussLegendre>;
            // Volume virtual work includes both physical-face tractions and the
            // work conjugate to free periodic jumps; face-only loading misses the latter.
            for(const auto& entry : base.fe->elements())
            {
                const auto& e=entry.second;
                for(size_t k=0;k<Q::quadratureOrder;++k)
                {
                    const auto bary=BarycentricTraits<dim>::x2l(Q::abscissa(k));
                    MatrixDim stress=MatrixDim::Zero();
                    for(const auto& physics : this->microstructures)
                    {
                        if(physics.get()!=static_cast<const MicrostructureBase<dim>*>(this))
                        {
                            stress+=physics->stress(e.position(bary),nullptr,&e,nullptr);
                        }
                    }
                    const auto voigt=base.voigtTraits.m2v(stress,false);
                    const double weight=e.absJ(bary)*Q::weight(k);
                    otherStress+=weight/base.mesh.volume()*voigt;
                    const auto local=(-weight*elasticDeformationFEM->u->sfmDef(e,bary).transpose()*voigt).eval();
                    for(int node=0;node<ElementType::nodesPerElement;++node)
                    {
                        source.template segment<dim>(dim*e.node(node).gID)+=local.template segment<dim>(dim*node);
                    }
                }
            }
            const MatrixDim plastic=this->microstructures.averagePlasticDistortion();
            const typename FEMControllerType::Vector6 strainOffset=base.voigtTraits.m2v(0.5*(plastic+plastic.transpose()),true)+elasticDeformationFEM->C.ldlt().solve(otherStress);
            const double time=base.simulationParameters.totalTime;
            weakProblemSolver->setRightHandSide(femLoadController->rightHandSide(source,time,strainOffset));
            weakProblemSolver->setLinearConstraintValues(femLoadController->constraintValues(time,strainOffset));
            weakProblemSolver->solveWithMasterSlaveConstraints();
            elasticDeformationFEM->u->setPeriodicGradient(femLoadController->gradient(weakProblemSolver->affineCoefficients()));

            // Compute zDot
            std::cout<<", zDot "<<std::flush;
//            zDot.setZero(elasticDeformationFEM->z->gSize());
            for(const auto& node : elasticDeformationFEM->z->fe().nodes())
            {
                zDot.template segment<dim>(dim*node.gID)=this->microstructures.inelasticDisplacementRate(node.P0,&node,nullptr,nullptr);
            }
        }
    }

    template<int dim>
    void ElasticDeformation<dim>::updateConfiguration()
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        if(uniformLoadController)
        {// already updated in solve()
        }
        else
        {
            // Update z
            *elasticDeformationFEM->z += zDot*this->timeSinceLastUpdate();
        }
    }

    template<int dim>
    double ElasticDeformation<dim>::getDt() const
    {
        return this->microstructures.ddBase.simulationParameters.dtMax;
    }

    template<int dim>
    void ElasticDeformation<dim>::output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        if(uniformLoadController)
        {
            const auto epsil(uniformLoadController->grad(this->microstructures.ddBase.simulationParameters.totalTime));
            const auto sigma(uniformLoadController->flux(this->microstructures.ddBase.simulationParameters.totalTime));
            
            f_file<<epsil.transpose()<<" "<<sigma.transpose()<<" ";
                        
            if(this->microstructures.ddBase.simulationParameters.runID==0)
            {
                for(int k=0;k<SymmetricVoigtTraits<dim>::voigtSize;++k)
                {
                    const auto& i(this->microstructures.ddBase.voigtTraits.tensorIndex(k,0)+1); // 1-based notation
                    const auto& j(this->microstructures.ddBase.voigtTraits.tensorIndex(k,1)+1); // 1-based notation
                    const std::string lab("e_"+std::to_string(i)+std::to_string(j));
                    F_labels<<lab<<"\n";
                }
                for(int k=0;k<SymmetricVoigtTraits<dim>::voigtSize;++k)
                {
                    const auto& i(this->microstructures.ddBase.voigtTraits.tensorIndex(k,0)+1); // 1-based notation
                    const auto& j(this->microstructures.ddBase.voigtTraits.tensorIndex(k,1)+1); // 1-based notation
                    const std::string lab("s_"+std::to_string(i)+std::to_string(j));
                    F_labels<<lab<<"\n";
                }
            }
        }
        else
        {
            const size_t nNodes(elasticDeformationFEM->u->fe().nodes().size());
            configIO.displacementMatrix().resize(nNodes,dim+dim);
            configIO.displacementMatrix().block(0,0,nNodes,dim)=elasticDeformationFEM->u->dofVector().reshaped(dim,nNodes).transpose();
            configIO.displacementMatrix().block(0,dim,nNodes,dim)=elasticDeformationFEM->z->dofVector().reshaped(dim,nNodes).transpose();
        }
    }

    template<int dim>
    typename ElasticDeformation<dim>::VectorDim ElasticDeformation<dim>::inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const
    {
        return VectorDim::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::MatrixDim ElasticDeformation<dim>::averagePlasticDistortion() const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::MatrixDim ElasticDeformation<dim>::averagePlasticDistortionRate() const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::VectorDim ElasticDeformation<dim>::displacement(const VectorDim& x,const NodeType* const  node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(uniformLoadController)
        {
            return this->microstructures.ddBase.voigtTraits.v2m(uniformLoadController->grad(this->microstructures.ddBase.simulationParameters.totalTime),true)*x;
        }
        else
        {
            if(node)
            {
                return (*eval(elasticDeformationFEM->u))(*node);
//                return (*eval(elasticDeformationFEM->u))(*node)+(*eval(elasticDeformationFEM->z))(*node);
            }
            else
            {
                if(ele)
                {
                    return (*eval(elasticDeformationFEM->u))(*ele,ele->simplex.pos2bary(x));
//                    return (*eval(elasticDeformationFEM->u))(*ele,ele->simplex.pos2bary(x))+(*eval(elasticDeformationFEM->z))(*ele,ele->simplex.pos2bary(x));

                }
                else
                {
                    return (*eval(elasticDeformationFEM->u))(x,guess);
//                    return (*eval(elasticDeformationFEM->u))(x,guess)+(*eval(elasticDeformationFEM->z))(x,guess);
                }
            }
        }
    }

    template<int dim>
    typename ElasticDeformation<dim>::MatrixDim ElasticDeformation<dim>::stress(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(uniformLoadController)
        {
            return this->microstructures.ddBase.voigtTraits.v2m(uniformLoadController->flux(this->microstructures.ddBase.simulationParameters.totalTime),false);
        }
        else
        {
            if(node)
            {
                return this->microstructures.ddBase.voigtTraits.v2m((*eval(elasticDeformationFEM->s))(*node),false);
            }
            else
            {
                if(ele)
                {
                    return this->microstructures.ddBase.voigtTraits.v2m((*eval(elasticDeformationFEM->s))(*ele,ele->simplex.pos2bary(x)),false);
                }
                else
                {
                    return this->microstructures.ddBase.voigtTraits.v2m((*eval(elasticDeformationFEM->s))(x,guess),false);
                }
            }
        }
    }

    template<int dim>
    typename ElasticDeformation<dim>::VectorMSize ElasticDeformation<dim>::mobileConcentration(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorMSize::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::VectorISize ElasticDeformation<dim>::immobileClusters(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorISize::Zero();
    }

    template<int dim>
    typename ElasticDeformation<dim>::MatrixDim ElasticDeformation<dim>::averageStress() const
    {
        if(uniformLoadController)
        {
            return this->microstructures.ddBase.voigtTraits.v2m(uniformLoadController->flux(this->microstructures.ddBase.simulationParameters.totalTime),false);
        }
        else
        {
            if(!femLoadController)
            {
                return MatrixDim::Zero();
            }
            return this->microstructures.ddBase.voigtTraits.v2m(elasticDeformationFEM->C*femLoadController->averageStrain(),false);
        }
    }

    template struct ElasticDeformation<3>;
}
#endif


//    template <int dim>
//    std::unique_ptr<ExternalLoadControllerBase<dim>> selectExternalLoadController(const DislocationDynamicsBase<dim>& ddBase)
//    {
//        if(ddBase.simulationParameters.externalLoadControllerName=="UniformExternalLoadController")
//        {
//            return std::unique_ptr<ExternalLoadControllerBase<dim>>(new UniformExternalLoadController<dim>(ddBase));
//        }
//        else
//        {
//            std::cout<<"Unknown externalLoadController name "<<ddBase.simulationParameters.externalLoadControllerName<<"No controller applied."<<std::endl;
//            return std::unique_ptr<ExternalLoadControllerBase<dim>>(nullptr);
//        }
//    }


//template <int dim>
//typename VoigtTraits<dim,dim>::VoigtStorageType getVoigtOrder()
//{
//    typename VoigtTraits<dim,dim>::VoigtStorageType voigtOrder;
//    for(int i=0;i<dim;++i)
//    {
//        for(int j=0;j<dim;++j)
//        {
//            voigtOrder(dim*i+j,0)=i;
//            voigtOrder(dim*i+j,1)=j;
//        }
//    }
//    return voigtOrder;
//}

//template <int dim>
//std::unique_ptr<UniformController<dim,dim>> getUniformController(const DislocationDynamicsBase<dim>& ddBase, const VoigtTraits<dim,dim>& vt)
//{
//    typedef typename UniformController<dim,dim>::MatrixGrad MatrixGrad;
//    typedef typename UniformController<dim,dim>::MatrixGrad MatrixGrad;
//
//    const double mu(ddBase.poly.mu);
//    const double nu(ddBase.poly.nu);
//    const double lam(2.0*nu*mu/(1.0-2.0*nu));
//
//    MatrixGrad C(MatrixGrad::Zero()); // C_ijkl=lam*d_ij*d_kl+mu*(d_ik*d_jl+d_il*d_jk)
//    for(int vI=0;vI<VoigtTraits<dim,dim>::voigtSize;++vI)
//    {
//        const auto& i(vt.tensorIndex(vI,0));
//        const auto& j(vt.tensorIndex(vI,1));
//        for(int vJ=0;vJ<VoigtTraits<dim,dim>::voigtSize;++vJ)
//        {
//            const auto& k(vt.tensorIndex(vJ,0));
//            const auto& l(vt.tensorIndex(vJ,1));
//            C(vI,vJ)= lam*kd(i,j)*kd(k,l)+mu*(kd(i,k)*kd(j,l)+kd(i,l)*kd(j,k));
//        }
//    }
//
//    TextFileParser parser(ddBase.simulationParameters.traitsIO.inputFilesFolder+"/uniformExternalLoadControllerNew.txt");
//
//    SymmetricVoigtTraits<dim> svt(parser.readMatrix<size_t,SymmetricVoigtTraits<dim>::voigtSize,2>("voigtOrder",true));
//    const auto f0(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStress0",true),false)));
//    const auto f0Dot(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStressRate",true),false)));
//    const auto g0(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStrain0",true),true)));
//    const auto g0Dot(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("ExternalStrainRate",true),true)));
//    const auto stiffnessRatio(vt.m2v(svt.v2m(parser.readMatrix<double,SymmetricVoigtTraits<dim>::voigtSize,1>("stiffnessRatio",true),false)));
//    const double t0(0.0);
//    std::cout<<f0<<std::endl<<std::endl;
//    std::cout<<f0Dot<<std::endl<<std::endl;
//    std::cout<<g0<<std::endl<<std::endl;
//    std::cout<<g0Dot<<std::endl<<std::endl;
//
//    return std::unique_ptr<UniformController<dim,dim>>(new UniformController<dim,dim>(t0,C,stiffnessRatio,g0,g0Dot,f0,f0Dot));
//}
