/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsFEM_H_
#define model_ClusterDynamicsFEM_H_

#ifdef MODELIB_CHOLMOD // SuiteSparse Cholmod module
#include <Eigen/CholmodSupport>
#endif

#ifdef MODELIB_UMFPACK // SuiteSparse UMFPACK module
#include <Eigen/UmfPackSupport>
#endif


#include <ClusterDynamicsParameters.h>
#include <DislocationDynamicsBase.h>
#include <EvalFunction.h>
#include <FixedDirichletSolver.h>
#include <SpatialODESolver.h>
#include <DDconfigIO.h>
#include <MicrostructureBase.h>
#include <SecondOrderReaction.h>
#include <FirstOrderReaction.h>
#include <ImmobileSinkRate.h>
#include <MicrostructureContainer.h>

namespace model
{

    template <int dim>
    struct FluxMatrix : public EvalFunction<FluxMatrix<dim>>
    {
        typedef typename DislocationDynamicsBase<dim>::ElementType ElementType;
        typedef Eigen::Matrix<double,dim+1,1> BaryType;
        constexpr static int mSize=ClusterDynamicsParameters<dim>::mSize;
        constexpr static int rows=dim*mSize;
        constexpr static int cols=rows;
        typedef Eigen::Matrix<double,rows,cols> MatrixType;

        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::map<size_t,std::vector<MatrixDim>> DiffusionContainerType;

        const ClusterDynamicsParameters<dim>& cdp;
        const DiffusionContainerType& D; // diffusion tensors of each species in each grain, at the current external stress
        
        FluxMatrix(const ClusterDynamicsParameters<dim>& cdp_in,const DiffusionContainerType& D_in);
        const MatrixType operator() (const ElementType& elem, const BaryType& bary) const;
    };

    template <int dim>
    struct InvDscaling : public EvalFunction<InvDscaling<dim>>
    {
        typedef typename DislocationDynamicsBase<dim>::ElementType ElementType;
        typedef Eigen::Matrix<double,dim+1,1> BaryType;
        constexpr static int mSize=ClusterDynamicsParameters<dim>::mSize;
        constexpr static int rows=mSize;
        constexpr static int cols=rows;
        typedef Eigen::Matrix<double,rows,cols> MatrixType;

        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::map<size_t,std::vector<MatrixDim>> DiffusionContainerType;

        const ClusterDynamicsParameters<dim>& cdp;
        const DiffusionContainerType& D; // diffusion tensors of each species in each grain, at the current external stress
        
        InvDscaling(const ClusterDynamicsParameters<dim>& cdp_in,const DiffusionContainerType& D_in);
        const MatrixType operator() (const ElementType& elem, const BaryType& bary) const;
    };

    template<int dim>
    struct ClusterDynamicsFEM
    {
        typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
        typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
        typedef typename MicrostructureBase<dim>::ElementType ElementType;
        typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
        typedef typename MicrostructureBase<dim>::NodeType NodeType;
        typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;
        typedef FiniteElement<ElementType> FiniteElementType;
        static constexpr int dVorder=4;
        typedef IntegrationDomain<FiniteElementType,0,dVorder,GaussLegendre> VolumeIntegrationDomainType;
        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
        static constexpr int iSize=ClusterDynamicsParameters<dim>::iSize;
        typedef TrialFunction<'m',mSize,FiniteElementType> MobileTrialType;
        typedef TrialFunction<'i',iSize,FiniteElementType> ImmobileTrialType;
        typedef TrialFunction<'z',dim,FiniteElementType> DiffusiveTrialType;
        typedef TrialGrad<MobileTrialType> MobileGradType;
        typedef TrialProd<FluxMatrix<dim>,MobileGradType> MobileFluxType;
        
        typedef TrialProd<InvDscaling<dim>,MobileTrialType> MobileTestType;
        typedef TrialGrad<MobileTestType> MobileTestGradType;
        typedef BilinearForm<MobileTestGradType,TrialProd<Constant<double,1,1>,MobileFluxType>> MobileBilinearFormType;

//        typedef BilinearForm<MobileGradType,TrialProd<Constant<double,1,1>,MobileFluxType>> MobileBilinearFormType;
        typedef BilinearWeakForm<MobileBilinearFormType,VolumeIntegrationDomainType> MobileBilinearWeakFormType;

        typedef TrialFunction<'d',mSize,FiniteElementType> MobileIncrementTrialType;
        // Drift-diffusion (activity) formulation: mobileClusters stores the activity W=c*exp(E^e/kB/T).
        // The physical concentration c is reconstructed node-by-node with mobileConcentrationFactors=exp(-E^e/kB/T).
        typedef NodalTrialScaling<MobileTrialType> MobileConcentrationType;
        typedef NodalTrialScaling<MobileIncrementTrialType> MobileConcentrationIncrementType;
        typedef TrialProd<InvDscaling<dim>,MobileIncrementTrialType> MobileIncrementTestType;
        typedef TrialGrad<MobileIncrementTestType> MobileIncrementTestGradType;

        typedef TrialGrad<MobileIncrementTrialType> MobileIncrementGradType;
        typedef TrialProd<FluxMatrix<dim>,MobileIncrementGradType> MobileIncrementFluxType;
        typedef BilinearForm<MobileIncrementTestGradType,TrialProd<Constant<double,1,1>,MobileIncrementFluxType>> MobileIncrementBilinearFormType;
        typedef BilinearWeakForm<MobileIncrementBilinearFormType,VolumeIntegrationDomainType> MobileIncrementBilinearWeakFormType;

        // typedef TrialGrad<MobileTrialType> MobileGradType;
        // typedef TrialProd<FluxMatrix<dim>,MobileGradType> MobileFluxType;

        // typedef TrialGrad<MobileTrialType> MobileTestGradType;
        // typedef BilinearForm<MobileTestGradType,TrialProd<Constant<double,1,1>,MobileFluxType>> MobileBilinearFormType;
        // typedef BilinearWeakForm<MobileBilinearFormType,VolumeIntegrationDomainType> MobileBilinearWeakFormType;

        // typedef TrialFunction<'d',mSize,FiniteElementType> MobileIncrementTrialType;
        // typedef TrialGrad<MobileIncrementTrialType> MobileIncrementTestGradType;

        // typedef TrialGrad<MobileIncrementTrialType> MobileIncrementGradType;
        // typedef TrialProd<FluxMatrix<dim>,MobileIncrementGradType> MobileIncrementFluxType;
        // typedef BilinearForm<MobileIncrementTestGradType,TrialProd<Constant<double,1,1>,MobileIncrementFluxType>> MobileIncrementBilinearFormType;
        // typedef BilinearWeakForm<MobileIncrementBilinearFormType,VolumeIntegrationDomainType> MobileIncrementBilinearWeakFormType;
        
        typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMatrixType;
#ifdef CHOLMOD_H // SuiteSparse Cholmod (LLT) module
    typedef Eigen::CholmodSupernodalLLT<SparseMatrixType> LltSolverType;
#else
    typedef Eigen::SimplicialLLT<SparseMatrixType> LltSolverType;
#endif
        
#ifdef UMFPACK_H // SuiteSparse Umfpack (LU) module
    typedef Eigen::UmfPackLU<SparseMatrixType> LuSolverType;
#else
    typedef Eigen::SparseLU<SparseMatrixType> LuSolverType;
#endif
                
        typedef FixedDirichletSolver<LltSolverType,Eigen::ConjugateGradient<SparseMatrixType>> MobileSolverType;
        typedef FixedDirichletSolver<LuSolverType,Eigen::BiCGSTAB<SparseMatrixType>> MobileReactionSolverType;
        typedef std::map<size_t,std::vector<MatrixDim>> DiffusionContainerType;
        typedef std::map<size_t,std::vector<MatrixDim>> InverseDiffusionContainerType;
        typedef std::map<size_t,Eigen::Array<double,mSize,1>> DiffusionDeterminantContainerType;
        typedef std::map<size_t,Eigen::Matrix<double,mSize,mSize>> FirstOrderReactionContainerType;
        typedef std::map<size_t,std::vector<Eigen::Matrix<double,mSize,mSize>>> SecondOrderReactionContainerType;

        const DislocationDynamicsBase<dim>& ddBase;
        const MicrostructureContainer<dim>& microstructures;
        const ClusterDynamicsParameters<dim>& cdp;
        MatrixDim externalStress;                          // macroscopic external stress entering the diffusion tensors
        DiffusionContainerType currentD;                   // diffusion tensors at externalStress, used by the FEM operators
        InverseDiffusionContainerType analyticalInvD;      // same, inverse, used by the analytical DD concentration kernel
        DiffusionDeterminantContainerType analyticalDetD;  // same, determinant, used by the analytical DD concentration kernel
        FirstOrderReactionContainerType reactionR1cd;      // first-order reaction coefficients of each grain at externalStress
        SecondOrderReactionContainerType reactionR2;       // second-order reaction coefficients of each grain at externalStress
        const InvDscaling<dim> iDs;
        
        const Eigen::Matrix<double,mSize,mSize> invTrD;
        
        MobileTrialType mobileClusters;
        Eigen::VectorXd mobileConcentrationFactors;
        MobileGradType mobileGrad;
        MobileFluxType mobileFlux;
        ImmobileTrialType immobileClusters;

        const int nodeListInternalExternal;
        MobileIncrementTrialType mobileClustersIncrement;
        VolumeIntegrationDomainType dV;
        MobileBilinearWeakFormType mBWF;
        MobileIncrementBilinearWeakFormType dmBWF;

        MobileSolverType mSolver;
        bool dirichletConditionsInitialized;
        bool solverInitialized;

        Eigen::VectorXd cascadeGlobalProduction;
        Eigen::VectorXd immobileClusterRate;

//         ClusterDynamicsFEM(const DislocationDynamicsBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in);
//         void solveMobileClusters();
//         void solveImmobileClusters();
//         void solve();
// //        void applyBoundaryConditions();
//         void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels);
//         void initializeSolver();
//         VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const;

        ClusterDynamicsFEM(DislocationDynamicsBase<dim>& ddBase_in,const MicrostructureContainer<dim>& microstructures_in,const ClusterDynamicsParameters<dim>& cdp_in);
        MatrixDim nodalStress(const NodeType& node) const;
        void updateDiffusionTensors();
        void updateReactionCoefficients();
        void updateMobileConcentrationFactors();
        void setExternalStress(const MatrixDim& externalStress_in);
        void solveMobileClusters(const bool hasDiscreteLoops);
        void solveImmobileClusters();
        void updateImmobileClusters(const double dt);
        void solve(const bool hasDiscreteLoops);
//        void applyBoundaryConditions();
        void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels);
        void initializeSolver();
        VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const;

        // Added all this below : 
        typedef Eigen::Matrix<double,dim+1,1> BaryType;
        std::vector<Eigen::Matrix<double,dim,dim>> betaP(const ElementType& ele, const BaryType& bary) const;
        std::vector<Eigen::Matrix<double,dim,dim>> betaV(const ElementType& ele, const BaryType& bary) const;
        MatrixDim averageBetaP() const;
        MatrixDim averageBetaPKernel(const Eigen::Matrix<double,dim,1>& a1, const ElementType& ele) const;
        MatrixDim averageBetaV() const;
        MatrixDim averageBetaVKernel(const Eigen::Matrix<double,dim,1>& a1, const ElementType& ele) const;

    };
    
}
#endif

