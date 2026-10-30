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


#include <ClusterDynamics/ClusterDynamicsParameters.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <EvalFunction.h>
#include <WeakProblem.h>
#include <ClusterBoundarySinks.h>
#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <DefectiveCrystal/MicrostructureBase.h>
#include <SecondOrderReaction.h>
#include <FirstOrderReaction.h>
#include <ImmobileSinkRate.h>
#include <DefectiveCrystal/MicrostructureContainer.h>

namespace model
{

    template <int dim>
    struct FluxMatrix : public EvalFunction<FluxMatrix<dim>>
    {
        typedef typename DefectiveCrystalBase<dim>::ElementType ElementType;
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
        typedef typename DefectiveCrystalBase<dim>::ElementType ElementType;
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

        using SparseMatrixType=Eigen::SparseMatrix<double>;
        using Problem=WeakProblem<DefaultWeakProblemDirectSolver<SparseMatrixType>,Eigen::BiCGSTAB<SparseMatrixType>>;
        template<class Trial>
        struct Load
        {
            std::shared_ptr<Trial> trial;
            Eigen::VectorXd value;
            const auto& field() const
            {
                return trial->field();
            }
            Eigen::VectorXd globalVector() const
            {
                return value;
            }
        };
        template<class Trial>
        struct MatrixForm
        {
            std::shared_ptr<Trial> trial;
            std::vector<Eigen::Triplet<double>> entries;
            const auto& field() const
            {
                return trial->field();
            }
            auto globalTriplets() const
            {
                return entries;
            }
        };
        using ImmobileRateTrialType=TrialFunction<'r',iSize,FiniteElementType>;
        const DefectiveCrystalBase<dim>& ddBase;
        const ClusterDynamicsParameters<dim>& cdp;
        const std::shared_ptr<InvDscaling<dim>> iDs;
        
        const Eigen::Matrix<double,mSize,mSize> invTrD;
        
        std::shared_ptr<MobileTrialType> mobileClusters;
        std::shared_ptr<MobileGradType> mobileGrad;
        std::shared_ptr<MobileFluxType> mobileFlux;
        std::shared_ptr<ImmobileTrialType> immobileClusters;

        const std::map<size_t,std::vector<VectorDim>> sinkNormals;
        std::shared_ptr<MobileIncrementTrialType> mobileClustersIncrement;
        VolumeIntegrationDomainType dV;
        std::shared_ptr<MobileBilinearWeakFormType> mBWF;
        std::shared_ptr<MobileIncrementBilinearWeakFormType> dmBWF;

        std::shared_ptr<Problem> mSolver,rSolver,iSolver;
        std::shared_ptr<Load<MobileTrialType>> mobileLoad;
        std::shared_ptr<Load<MobileIncrementTrialType>> reactionLoad;
        std::shared_ptr<MatrixForm<MobileIncrementTrialType>> reactionMatrix;
        std::shared_ptr<ImmobileRateTrialType> immobileRateTrial;
        SparseMatrixType diffusionMatrix;
        std::vector<Eigen::Triplet<double>> diffusionTriplets;
        Eigen::MatrixXd meanRows;
        Eigen::MatrixXd conservedRows(const bool useImmobileClusters) const;
        void enforceMobileConstraints();
        bool solverInitialized;

        const Eigen::VectorXd cascadeGlobalProduction;
        Eigen::VectorXd immobileClusterRate;

//         ClusterDynamicsFEM(const DefectiveCrystalBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in);
//         void solveMobileClusters();
//         void solveImmobileClusters();
//         void solve();
// //        void applyBoundaryConditions();
//         void initializeConfiguration(const DefectiveCrystalConfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels);
//         void initializeSolver();
//         VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const;

        ClusterDynamicsFEM(DefectiveCrystalBase<dim>& ddBase_in,const ClusterDynamicsParameters<dim>& cdp_in);
        void solveMobileClusters(const bool hasDiscreteLoops);
        void solveImmobileClusters();
        void updateImmobileClusters(const double dt);
        void solve(const bool hasDiscreteLoops);
//        void applyBoundaryConditions();
        void initializeConfiguration(const DefectiveCrystalConfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels);
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

