/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ElasticDeformation_H_
#define model_ElasticDeformation_H_

#ifdef MODELIB_CHOLMOD // SuiteSparse Cholmod module
#include <Eigen/CholmodSupport>
#endif

#include <DefectiveCrystal/MicrostructureContainer.h>
#include <ElasticDeformation/ElasticDeformationFEM.h>
#include <IntegrationDomain.h>
#include <GaussLegendre.h>
#include <DefectiveCrystal/UniformController.h>
#include <WeakProblem.h>
#include <GlobalRigidMotionPenalty.h>
#include <FEMLoadController.h>

namespace model
{
    template<int dim>
    struct ElasticDeformation : public MicrostructureBase<dim>
{
    typedef typename MicrostructureBase<dim>::ElementType ElementType;
    typedef typename MicrostructureBase<dim>::FiniteElementType FiniteElementType;
    typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
    typedef typename MicrostructureBase<dim>::NodeType NodeType;
    typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;
    typedef typename MicrostructureBase<dim>::VectorISize VectorISize;
    typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
    typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
    typedef MicrostructureContainer<dim> MicrostructureContainerType;
    typedef UniformController<SymmetricVoigtTraits<dim>::voigtSize> UniformControllerType;
    
    using DirectSolverType=DefaultWeakProblemDirectSolver<typename ElasticDeformationFEM<dim>::SparseMatrixType>;
    using IterativeSolverType=Eigen::BiCGSTAB<typename ElasticDeformationFEM<dim>::SparseMatrixType>;
    using WeakProblemType=WeakProblem<DirectSolverType,IterativeSolverType>;

    const bool useElasticDeformationFEM;
    const std::unique_ptr<ElasticDeformationFEM<dim>> elasticDeformationFEM;
    const std::unique_ptr<UniformControllerType> uniformLoadController;
    const double inertiaReliefPenaltyFactor;
    using FEMControllerType=FEMLoadController<typename ElasticDeformationFEM<dim>::TrialFunctionType>;
    std::shared_ptr<FEMControllerType> femLoadController;
    std::shared_ptr<GlobalRigidMotionPenalty<typename ElasticDeformationFEM<dim>::TrialFunctionType>> rigidMotionPenalty;
    std::shared_ptr<WeakProblemType> weakProblemSolver;
    Eigen::VectorXd zDot;
    bool solverInitialized;
    
    ElasticDeformation(MicrostructureContainerType& mc);
    void initializeConfiguration(const DefectiveCrystalConfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels) override;
    void registerSolvedQuantities(SolvedQuantityRegistry&) override;
    void solveStress();
    void solve() override;
    void solveDisplacementRate();
    double getDt() const override;
    void output(DefectiveCrystalConfigIO<dim>& configIO,DefectiveCrystalAuxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const override;
    void updateConfiguration() override;
    MatrixDim averagePlasticDistortion() const override ;
    MatrixDim averagePlasticDistortionRate() const override;
    VectorDim displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override;
    MatrixDim stress(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override;
    MatrixDim averageStress() const override;
    VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
    VectorMSize mobileConcentration(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
    VectorISize immobileClusters(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
    void initializeSolver();
    void updatePeriodicGradient();
};
}
#endif

