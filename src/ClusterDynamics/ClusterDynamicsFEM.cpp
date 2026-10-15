/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsFEM_cpp_
#define model_ClusterDynamicsFEM_cpp_

#include <cmath>
#include <chrono>
#include <algorithm>
#include <ClusterDynamicsFEM.h>
#include <ExternalAndInternalBoundary.h>
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
    ClusterDynamicsFEM<dim>::ClusterDynamicsFEM(DislocationDynamicsBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in) :
    /* init */ ddBase(ddBase_in)
    /* init */,cdp(cdp_in)
    /* init */,iDs(cdp)
    /* init */,mobileClusters(ddBase.fe->template trial<'m',mSize>())
    /* init */,mobileGrad(grad(mobileClusters))
    /* init */,mobileFlux(FluxMatrix<dim>(cdp)*mobileGrad)
    /* init */,immobileClusters(ddBase.fe->template trial<'i',iSize>())
    /* init */,immobileFloors(getImmobileFloors())
    /* init */,nodeListInternalExternal(ddBase.isPeriodicDomain ? -1 : ddBase.fe->template createNodeList<ExternalAndInternalBoundary>())
    /* init */,mobileClustersIncrement(ddBase.fe->template trial<'d',mSize>())
    /* init */,dV(ddBase.fe->template domain<EntireDomain,dVorder,GaussLegendre>())
    /* init */,R1((cdp.mobileDissociationRate+cdp.mobileOtherSinkRate).eval())
    /* init */,mBWF((test(grad(iDs*mobileClusters)),-cdp.omega*this->mobileFlux)*dV+(test(iDs*mobileClusters),R1*(-1.0*mobileClusters))*dV)
    /* init */,dmBWF((test(grad(iDs*mobileClustersIncrement)),-cdp.omega*(FluxMatrix<dim>(this->cdp)*grad(mobileClustersIncrement)))*dV)   
    /* init */,mSolver(true,FLT_EPSILON)
    /* init */,solverInitialized(false)
    /* init */,mobileFieldInitialized(false)
    /* init */,cascadeGlobalProduction(((test(iDs*this->mobileClusters),make_constant(this->cdp.G))*dV).globalVector())
    /* init */,rateEquations(cdp)
    /* init */,immobileSpeciesSolver(rateEquations,1.0e-8,1.0e-24)
    /* init */,nucleationFractions(1,rateEquations.nucleationWeights())
    {
        mobileClustersIncrement.setConstant(Eigen::Matrix<double,mSize,1>::Zero());
        mobileClusters.setConstant(cdp.equilibriumMobileConcentration(0.0).matrix().transpose());
        // immobileClusters.setConstant(Eigen::Matrix<double,iSize,1>::Zero());
        Eigen::Array<double,1,iSize> initialValues(cdp.initialImmobileClusters); // initialImmobileClusters_SI: number densities per unit volume, then contents per atom
        initialValues.head(iSize/2)*=cdp.omega; // clusters per atom
        immobileClusters.setConstant(initialValues.max(immobileFloors.transpose().array())); // a family can start empty. The floor keeps the ratio content/density defined
        if(iSize>0 && !ImmobileSpeciesSolver<dim>::available())
        {// the immobile species are integrated by CVODE
            throw std::runtime_error("ClusterDynamics: the immobile species need a build of MoDELib with SUNDIALS 7 (CVODE). Install SUNDIALS and configure again.");
        }
    }


    template<int dim>
    void ClusterDynamicsFEM<dim>::solveMobileReactions(const bool useImmobileClusters,const Eigen::VectorXd& production)
    {/*!@param[in] useImmobileClusters true if the absorption at the immobile clusters (LoopSinkReaction) is part of the equations
      * @param[in] production the global vector of the sources that do not depend on the mobile concentrations
      *
      * Newton iterations of the steady reaction-diffusion equations of the mobile species.
      */
        const double cTol(1e-5);
        double cError(1.0);

        const auto DR((this->cdp.mobileDissociationRate).eval());
        const auto OS((this->cdp.mobileOtherSinkRate).eval());
        auto bWF_DR((test(iDs*mobileClustersIncrement),DR*(-1.0*mobileClustersIncrement))*dV);
        auto bWF_OS((test(iDs*mobileClustersIncrement),OS*(-1.0*mobileClustersIncrement))*dV);
        
        while(cError>cTol)
        {
            SecondOrderReaction<MobileTrialType> R2(mobileClusters,this->cdp);
            auto bWF_R2((test(iDs*mobileClustersIncrement),R2*(-1.0*mobileClustersIncrement))*dV);
            auto lWF_R2((test(iDs*mobileClustersIncrement),eval(R2*(0.5*mobileClusters)))*dV);

            MobileReactionSolverType rSolver(false,FLT_EPSILON);

            if(useImmobileClusters)
            {// absorption at the immobile clusters, from the same rate equations that ImmobileSpeciesSolver integrates
                LoopSinkReaction<MobileTrialType,ImmobileTrialType> IS(immobileClusters,rateEquations);
                auto bWF_IS((test(iDs*mobileClustersIncrement),IS*(-1.0*mobileClustersIncrement))*dV);
                auto lWF_IS((test(iDs*mobileClustersIncrement),eval(IS*mobileClusters))*dV);
                rSolver.compute(dmBWF+bWF_DR+bWF_OS+bWF_R2+bWF_IS);
                mobileClustersIncrement=rSolver.solve(production-mSolver.getA()*mobileClusters.dofVector()+(lWF_R2+lWF_IS).globalVector());
            }
            else
            {// no immobile fields (no immobile species, or discrete loops): the same Newton step without their absorption
                rSolver.compute(dmBWF+bWF_DR+bWF_OS+bWF_R2);
                mobileClustersIncrement=rSolver.solve(production-mSolver.getA()*mobileClusters.dofVector()+lWF_R2.globalVector());
            }

            Eigen::MatrixXd cOld(mobileClusters.dofVector());
            cOld.resize(mSize,mobileClusters.gSize()/mSize);
            mobileClusters += mobileClustersIncrement.dofVector();

            Eigen::MatrixXd cNew(mobileClusters.dofVector());
            cNew.resize(mSize,mobileClusters.gSize()/mSize);

            const Eigen::VectorXd absErr((cNew-cOld).rowwise().norm());
            const Eigen::VectorXd cNewNorm((cNew.rowwise().norm().array()+1.e-50).matrix());
            const Eigen::VectorXd relErr((absErr.array()/cNewNorm.array()).matrix());

            cError=relErr.maxCoeff();
            std::cout<<"        convergenceError="<<cError<<std::endl;
        }
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solveMobileClusters(const bool hasDiscreteLoops)
    {/*! The mobile species are solved at steady state, for the immobile fields
      * currently held. The method is that of v2.0.0. The immobile clusters
      * enter through the sink term and through the vacancies that the vacancy
      * clusters release.
      */
        std::cout<<"    mobile solver"<<std::endl;
        if(!this->cdp.computeReactions || !mobileFieldInitialized)
        {// linear solve: diffusion and R1. It is the solution without reactions, or the first guess of the Newton iterations
            mobileClusters=mSolver.solve(cascadeGlobalProduction);
            mobileFieldInitialized=true;
        }
        else
        {// warm start: the converged field of the previous step, with the boundary values of this step
            Eigen::VectorXd& c(TrialBase<MobileTrialType>::dofVector());
            for(const auto& bc : TrialBase<MobileTrialType>::dirichletConditions())
            {
                c(bc.first)=bc.second;
            }
        }
        const bool useImmobileClusters(!hasDiscreteLoops && iSize > 0);

        if(this->cdp.computeReactions)
        {
            if(useImmobileClusters)
            {// the immobile clusters release defects (ICS) and absorb them (in solveMobileReactions)
                ImmobileClusterSource<MobileTrialType,ImmobileTrialType> ICS(immobileClusters,rateEquations);
                const Eigen::VectorXd production(cascadeGlobalProduction+((test(iDs*this->mobileClusters),ICS)*dV).globalVector());
                solveMobileReactions(useImmobileClusters,production);
            }
            else
            {// no immobile fields: the Newton iterations still solve diffusion, DR, OS and R2
                solveMobileReactions(useImmobileClusters,cascadeGlobalProduction);
            }
        }
    }



    template<int dim>
    void ClusterDynamicsFEM<dim>::updateImmobileClusters(const double dt)
    {
        // implicit integration of the rate equations at each node, with the mobile field of the last solve
        std::cout<<", immobile ODE solver (CVODE)"<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        Eigen::VectorXd& immobileDofs(TrialBase<ImmobileTrialType>::dofVector());
        immobileSpeciesSolver.integrate(mobileClusters.dofVector(),immobileDofs,nucleationFractions,dt);
        applyImmobileFloors(immobileDofs);
        std::cout<<", "<<immobileSpeciesSolver.lastInternalSteps()<<" internal steps, at most "<<immobileSpeciesSolver.lastMaxInternalSteps()<<" at a node"
        /*     */<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }


    template<int dim>
    void ClusterDynamicsFEM<dim>::scaleImmobileFamily(const int& family,const double& factor)
    {/*!@param[in] family the index of an immobile family
      * @param[in] factor the fraction of the family that remains in the fields
      *
      * Scales the number density and the defect content of a family by the
      * same factor, which leaves the size of its clusters unchanged.
      */
        for(size_t nodeID=0; nodeID<immobileClusters.nodeSize(); nodeID++)
        {
            double& N(TrialBase<ImmobileTrialType>::dofVector()(iSize*nodeID+family));
            double& c(TrialBase<ImmobileTrialType>::dofVector()(iSize*nodeID+iSize/2+family));
            N=std::max(N*factor,immobileFloors(family)); // clusters per atom
            c=std::max(c*factor,immobileFloors(iSize/2+family));
        }
    }

    template<int dim>
    auto ClusterDynamicsFEM<dim>::getImmobileFloors() const -> Eigen::Matrix<double,iSize,1>
    {/*!\returns the floor of each immobile field at a node: immobileFloor per unit volume for the
      * number of clusters (immobileFloor*Omega clusters per atom), immobileFloor for the contents per atom
      */
        Eigen::Matrix<double,iSize,1> temp;
        for(int k=0;k<iSize/2;++k)
        {
            temp(k)=immobileFloor*cdp.omega;
            temp(iSize/2+k)=immobileFloor;
        }
        return temp;
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::applyImmobileFloors(Eigen::VectorXd& immobileDofs) const
    {/*!@param[in,out] immobileDofs the nodal values of the immobile fields, raised to their floors
      */
        for(long int i=0;i<immobileDofs.size();++i)
        {
            immobileDofs(i)=std::max(immobileDofs(i),immobileFloors(i%iSize));
        }
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::solve(const bool hasDiscreteLoops)
    {
        solveMobileClusters(hasDiscreteLoops);
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::initializeSolver()
    {
        std::cout<<" decomposing:"<<std::endl;
        std::array<bool,mSize> allComps;
        for(int k=0;k<mSize;++k)
        {
            allComps[k]=true;
        }
        
        mobileClusters.addDirichletCondition(nodeListInternalExternal,Fix(),allComps); // apply to the four components of c
        mobileClustersIncrement.addDirichletCondition(nodeListInternalExternal,Fix(),allComps); // apply to the four components of c
        mSolver.compute(mBWF,true); // call this after assigning the BCs; prints the times of the stages
        solverInitialized=true;
    }

    template<int dim>
    void ClusterDynamicsFEM<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {    
        if(size_t(configIO.cdMatrix().size())==mobileClusters.gSize()+immobileClusters.gSize())
        {
            const size_t nNodes(mobileClusters.fe().nodes().size());
            mobileClusters=configIO.cdMatrix().block(0,0,nNodes,mSize).transpose().reshaped(mobileClusters.gSize(),1);
            immobileClusters=configIO.cdMatrix().block(0,mSize,nNodes,iSize).transpose().reshaped(immobileClusters.gSize(),1);
            for(size_t node=0;node<nNodes;++node)
            {// the evl file holds the number densities per unit volume: clusters per atom in the field
                for(int k=0;k<iSize/2;++k)
                {
                    TrialBase<ImmobileTrialType>::dofVector()(iSize*node+k)*=cdp.omega;
                }
            }

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
            speciesFlux=eval(this->mobileFlux)(*node);
        }
        else
        {
            if(ele)
            {
                speciesFlux=eval(this->mobileFlux)(*ele,ele->simplex.pos2bary(x));
            }
            else
            {
                speciesFlux=eval(this->mobileFlux)(x,guess);
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
    std::vector<Eigen::Matrix<double,dim,dim>> ClusterDynamicsFEM<dim>::immobileStrain(const Eigen::Matrix<double,iSize,1>& value,const size_t& gID,const bool& relaxation) const
    {/*!@param[in] value the immobile fields at one node: clusters per atom, then defects per atom
      * @param[in] gID the grain, which orients the Burgers vectors
      * @param[in] relaxation true: weighted by immobileSpeciesRelRelaxVol (betaV). false: betaP
      *\returns the strain of each immobile family (report Eqs. 34-35, summed over the clusters):
      *
      *   p_k*c_k*[(1-S_k)*I/3 + S_k*bhat_k^bhat_k]   (times immobileSpeciesRelRelaxVol for betaV)
      *
      * Each cluster changes the volume by -Omega per vacancy (+Omega per interstitial) in both shapes:
      * isotropically as a bi-pyramid, along its Burgers vector as a loop (Christiaen et al.,
      * Acta Mater. 179 (2019) 93, Eqs. 10-11 and Figs. 7-8). N*m*Omega = c, so the content c_k
      * (defects per atom) sets the volume change and the sigmoid only its anisotropy.
      */
        std::vector<Eigen::Matrix<double,dim,dim>> strainVector;
        const Eigen::Array<double,1,iSize/2> n(value.template block<iSize/2,1>(0,0)); // clusters per atom
        const Eigen::Array<double,1,iSize/2> c(value.template block<iSize/2,1>(iSize/2,0)); // defects per atom
//      const Eigen::Array<double,1,iSize/2> m(c/n); // defects per cluster
        const Eigen::Array<double,1,iSize/2> m(cdp.populated(n,c).select(c/n,1.0)); // defects per cluster (1 for an absent family, whose strain is zero)
//      const Eigen::Array<double,1,iSize/2> N(n/cdp.omega); // clusters per unit volume
        const Eigen::Matrix<double,dim,dim> lat(ddBase.poly.grains.at(gID)->latticeBasis);
        const Eigen::Matrix<double,dim,iSize/2> localBurgers(lat*cdp.immobileSpeciesBurgers);
//      const Eigen::Array<double,1,iSize/2> cr(cdp.effectiveClusterRadius(n,c));
        for(int k=0; k<iSize/2; k++)
        { // N A b ^b*^b
            const double scale(relaxation? cdp.immobileSpeciesRelRelaxVol(k) : 1.0);
//          const double Area=M_PI*std::pow(cr(k),2);
//          const Eigen::Matrix<double,dim,dim> highValue(scale*cdp.immobileSpeciesVector(k)*N(k)*Area*localBurgers.col(k)*localBurgers.col(k).normalized().transpose());
//          const Eigen::Matrix<double,dim,dim> lowValue(scale*cdp.immobileSpeciesVector(k)*N(k)*cdp.pyramidRelaxationVolume*Eigen::Matrix<double,dim,dim>::Identity()/3.0); // isotropic, in any crystal
            const Eigen::Matrix<double,dim,1> bHat(localBurgers.col(k).normalized());
            const Eigen::Matrix<double,dim,dim> highValue(scale*cdp.immobileSpeciesVector(k)*c(k)*bHat*bHat.transpose()); // loop: volume change along the Burgers vector
            const Eigen::Matrix<double,dim,dim> lowValue(scale*cdp.immobileSpeciesVector(k)*c(k)*Eigen::Matrix<double,dim,dim>::Identity()/3.0); // bi-pyramid: the same volume change, isotropic, in any crystal
            strainVector.push_back(cdp.shapeInterpolation(m,lowValue,highValue,k).matrix());
        }
        return strainVector;
    }

    template<int dim>
    std::vector<Eigen::Matrix<double,dim,dim>> ClusterDynamicsFEM<dim>::betaP(const ElementType& ele, const BaryType& bary) const
    {/*! Strain of the immobile families at bary. It is computed at the nodes of the element, from the nodal values
      * that ImmobileSpeciesSolver evolves, and interpolated with the shape functions. The immobile fields are never
      * evaluated between nodes, where the quadratic shape functions can make them negative.
      */
        const auto sf(ElementType::sf(bary));
        const size_t gID(ele.simplex.region->regionID);
        std::vector<Eigen::Matrix<double,dim,dim>> strainVector(iSize/2,Eigen::Matrix<double,dim,dim>::Zero());
        for(int n=0;n<ElementType::nodesPerElement;++n)
        {
            const std::vector<Eigen::Matrix<double,dim,dim>> nodalStrain(immobileStrain(immobileClusters.dofs(ele.node(n)),gID,false));
            for(int k=0;k<iSize/2;++k)
            {
                strainVector[k]+=sf(n)*nodalStrain[k];
            }
        }
        return strainVector;
    }



    template<int dim>
    std::vector<Eigen::Matrix<double,dim,dim>> ClusterDynamicsFEM<dim>::betaV(const ElementType& ele, const BaryType& bary) const
    {/*! Volumetric strain at bary. The immobile families are computed at the nodes and interpolated, as in betaP.
      * The mobile species are linear in the concentrations and are evaluated at bary.
      */
        const auto sf(ElementType::sf(bary));
        const size_t gID(ele.simplex.region->regionID);
        std::vector<Eigen::Matrix<double,dim,dim>> strainVector(iSize/2,Eigen::Matrix<double,dim,dim>::Zero());
        for(int n=0;n<ElementType::nodesPerElement;++n)
        {
            const std::vector<Eigen::Matrix<double,dim,dim>> nodalStrain(immobileStrain(immobileClusters.dofs(ele.node(n)),gID,true));
            for(int k=0;k<iSize/2;++k)
            {
                strainVector[k]+=sf(n)*nodalStrain[k];
            }
        }
       // Add up the contribution from the mobile species
        Eigen::Matrix<double,dim,dim> tempStrain(Eigen::Matrix<double,dim,dim>::Zero());
        const Eigen::Matrix<double,mSize,1> localConcentration(eval(mobileClusters)(ele,bary)); // content, defects per atom
        for(int k=0; k<mSize; k++)
        {
            tempStrain += 1.0/3.0*(cdp.msVector(k)>0.0? 1.0 : -1.0)*localConcentration(k)*cdp.msRelRelaxVol(k)*Eigen::Matrix<double,dim,dim>::Identity(); // content space
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
//    const auto& nodesInternalExternal(mobileClusters.fe().nodeList(nodeListInternalExternal));
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
//            mobileClusters.dirichletConditions().at(mSize*node->gID+k) = bndConcentration(k) - otherConcentration(k);
//        }
//    }
//}
