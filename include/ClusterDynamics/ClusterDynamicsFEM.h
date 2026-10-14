/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamicsFEM_H_
#define model_ClusterDynamicsFEM_H_

//#ifdef MODELIB_CHOLMOD // SuiteSparse Cholmod module
//#include <Eigen/CholmodSupport>
//#endif

#ifdef MODELIB_UMFPACK // SuiteSparse UMFPACK module
#include <Eigen/UmfPackSupport>
#endif


#include <ClusterDynamicsParameters.h>
#include <DislocationDynamicsBase.h>
#include <EvalFunction.h>
#include <FixedDirichletSolver.h>
// #include <SpatialODESolver.h> // explicit (Euler) integrator, replaced by ImmobileSpeciesSolver in 2.2.2 (SpatialODESolver_old.h)
#include <DDconfigIO.h>
#include <MicrostructureBase.h>
#include <SecondOrderReaction.h>
// #include <FirstOrderReaction.h> // explicit (Euler) integrator, replaced by ImmobileSpeciesSolver in 2.2.2 (FirstOrderReaction_old.h)
// #include <ImmobileSinkRate.h> // explicit (Euler) integrator, replaced by ImmobileSpeciesSolver in 2.2.2 (ImmobileSinkRate_old.h)
#include <ImmobileRateEquations.h>
#include <ImmobileSpeciesSolver.h>
#include <LoopSinkReaction.h>
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

        const ClusterDynamicsParameters<dim>& cdp;
        
        FluxMatrix(const ClusterDynamicsParameters<dim>& cdp_in);
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

        const ClusterDynamicsParameters<dim>& cdp;
        
        InvDscaling(const ClusterDynamicsParameters<dim>& cdp_in);
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
        static constexpr double immobileFloor=1.0e-33; // smallest value of the immobile fields: of N per unit volume, and of the contents per atom
//        static constexpr int dVprojectionOrder=14; // exact for the mass matrix of the quadratic elements
//        typedef IntegrationDomain<FiniteElementType,0,dVprojectionOrder,GaussLegendre> ProjectionIntegrationDomainType;
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
//#ifdef CHOLMOD_H // SuiteSparse Cholmod (LLT) module
//  typedef Eigen::CholmodSupernodalLLT<SparseMatrixType> LltSolverType;
//#else
//  typedef Eigen::SimplicialLLT<SparseMatrixType> LltSolverType;
//#endif
        
#ifdef UMFPACK_H // SuiteSparse Umfpack (LU) module
    typedef Eigen::UmfPackLU<SparseMatrixType> LuSolverType;
#else
    typedef Eigen::SparseLU<SparseMatrixType> LuSolverType;
#endif
                
//      typedef FixedDirichletSolver<LltSolverType,Eigen::ConjugateGradient<SparseMatrixType>> MobileSolverType;
        typedef FixedDirichletSolver<LuSolverType,Eigen::BiCGSTAB<SparseMatrixType>> MobileSolverType; // diffusion and R1 are not symmetric: LU
        typedef FixedDirichletSolver<LuSolverType,Eigen::BiCGSTAB<SparseMatrixType>> MobileReactionSolverType;
        typedef TrialProd<Constant<Eigen::Matrix<double,mSize,mSize>,mSize,mSize>,TrialProd<Constant<double,1,1>,MobileTrialType>> MobileLinearReactionType; // R1*(-1.0*c)
        typedef BilinearWeakForm<BilinearForm<MobileTestType,MobileLinearReactionType>,VolumeIntegrationDomainType> MobileLinearReactionBilinearWeakFormType;
        typedef BilinearWeakSum<MobileBilinearWeakFormType,MobileLinearReactionBilinearWeakFormType> MobileLinearBilinearWeakFormType; // diffusion and R1

        const DislocationDynamicsBase<dim>& ddBase;
        const ClusterDynamicsParameters<dim>& cdp;
        const InvDscaling<dim> iDs;
        
        const Eigen::Matrix<double,mSize,mSize> invTrD;
        
        MobileTrialType mobileClusters;
        MobileGradType mobileGrad;
        MobileFluxType mobileFlux;
        ImmobileTrialType immobileClusters; // per node: clusters per atom n_k, then defects per atom c_k
        const Eigen::Matrix<double,iSize,1> immobileFloors; // floor of each immobile field: immobileFloor*Omega for n_k, immobileFloor for c_k

        const int nodeListInternalExternal;
        MobileIncrementTrialType mobileClustersIncrement;
        VolumeIntegrationDomainType dV;
//        ProjectionIntegrationDomainType dVprojection; // used to project the rates of the immobile species
//      MobileBilinearWeakFormType mBWF;
        const Eigen::Matrix<double,mSize,mSize> R1; // DR+OS: the constant first-order terms of the mobile species, content space
        MobileLinearBilinearWeakFormType mBWF; // diffusion and R1
        MobileIncrementBilinearWeakFormType dmBWF;

        MobileSolverType mSolver;
        bool solverInitialized;
        bool mobileFieldInitialized; // true once the mobile field has been solved: the next solves start from it

        const Eigen::VectorXd cascadeGlobalProduction;
//        Eigen::VectorXd immobileClusterRate;

        // Immobile clusters, integrated by CVODE at each node
        typedef typename ImmobileRateEquations<dim>::ArrayChi ArrayChi;
        const ImmobileRateEquations<dim> rateEquations; // rate equations of the immobile clusters at a point
//        const bool useImmobileODESolver; // true: ImmobileSpeciesSolver. false: projection of the rates and explicit step
        ImmobileSpeciesSolver<dim> immobileSpeciesSolver;
        std::vector<ArrayChi> nucleationFractions; // one entry per node, or one for all nodes

//        static bool getUseImmobileODESolver(const DislocationDynamicsBase<dim>& ddBase);
        template<typename SinkType>
//      void solveMobileReactions(const SinkType& R1sink,const bool useImmobileClusters,const Eigen::VectorXd& production);
        void solveMobileReactions(const SinkType& IS,const bool useImmobileClusters,const Eigen::VectorXd& production);

//         ClusterDynamicsFEM(const DislocationDynamicsBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in);
//         void solveMobileClusters();
//         void solveImmobileClusters();
//         void solve();
// //        void applyBoundaryConditions();
//         void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels);
//         void initializeSolver();
//         VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const;

        ClusterDynamicsFEM(DislocationDynamicsBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in);
        void solveMobileClusters(const bool hasDiscreteLoops);
//        void solveImmobileClusters();
        void updateImmobileClusters(const double dt);
        void scaleImmobileFamily(const int& family,const double& factor);
        Eigen::Matrix<double,iSize,1> getImmobileFloors() const;
        void applyImmobileFloors(Eigen::VectorXd& immobileDofs) const;
        void solve(const bool hasDiscreteLoops);
//        void applyBoundaryConditions();
        void initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels);
        void initializeSolver();
        VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const;

        // Added all this below : 
        typedef Eigen::Matrix<double,dim+1,1> BaryType;
        std::vector<Eigen::Matrix<double,dim,dim>> betaP(const ElementType& ele, const BaryType& bary) const;
        std::vector<Eigen::Matrix<double,dim,dim>> betaV(const ElementType& ele, const BaryType& bary) const;
        std::vector<Eigen::Matrix<double,dim,dim>> immobileStrain(const Eigen::Matrix<double,iSize,1>& value,const size_t& gID,const bool& relaxation) const; // strain of each immobile family from the values at one node
        MatrixDim averageBetaP() const;
        MatrixDim averageBetaPKernel(const Eigen::Matrix<double,dim,1>& a1, const ElementType& ele) const;
        MatrixDim averageBetaV() const;
        MatrixDim averageBetaVKernel(const Eigen::Matrix<double,dim,1>& a1, const ElementType& ele) const;

    };
    
}
#endif

