/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// valgrind --leak-check=full --show-leak-kinds=all ./DDomp

#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_

#ifdef _MODEL_MPI_
#define _MODEL_DD_MPI_
#endif

// template header cpp
// https://www.codeproject.com/Articles/48575/How-to-Define-a-Template-Class-in-a-h-File-and-Imp

#ifdef _OPENMP
#include <optional>
#include <omp.h>
#endif

#include <vector>
#include <functional>
#include <chrono>
#include <map>
#include <memory>
#include <numbers>

#include <Eigen/Dense>

#include <TerminalColors.h>
#include <DislocationNetworkTraits.h>
#include <DislocationDynamicsModule.h>
#include <LoopNetwork.h>
//#include <DislocationNetworkComponent.h>
#include <DislocationLoop.h>
#include <DislocationLoopNode.h>
#include <DislocationLoopLink.h>
#include <DislocationNode.h>
#include <DislocationSegment.h>
#include <Hermite.h>

#include <DislocationNetworkRemesh.h>
#include <DislocationJunctionFormation.h>
#include <DislocationDynamics/DislocationFieldBase.h>
#include <DDtimeStepper.h>
#include <DefectiveCrystal/DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <GlidePlaneModule.h>

#include <DislocationNodeContraction.h>
#include <InclusionMicrostructure/EshelbyInclusionBase.h>
#include <InclusionMicrostructure/SphericalInclusion.h>
#include <InclusionMicrostructure/PolyhedronInclusion.h>

#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <DislocationGlideSolverFactory.h>
#include <DislocationClimbSolverFactory.h>
#include <CrossSlipModels.h>
#include <DefectiveCrystal/MicrostructureBase.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <InclusionMicrostructure.h>

#ifndef NDEBUG
#define VerboseDislocationNetwork(N,x) if(verboseDislocationNetwork>=N){std::cout<<x;}
#else
#define VerboseDislocationNetwork(N,x)
#endif

namespace model
{
    template <int dim>
    class DislocationNetwork : public MicrostructureBase<dim>
    /*                      */,public LoopNetwork<DislocationNetwork<dim> >
    {
        
    public:
        
        
        typedef TypeTraits<DislocationNetwork<dim>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::VectorLowerDim VectorLowerDim;
        typedef typename TraitsType::MatrixDim MatrixDim;
        typedef MicrostructureContainer<dim> MicrostructureContainerType;
        typedef typename MicrostructureBase<dim>::ElementType ElementType;
        typedef typename MicrostructureBase<dim>::NodeType NodeType;
        typedef typename MicrostructureBase<dim>::SimplexDim SimplexDim;
        typedef typename MicrostructureBase<dim>::VectorMSize VectorMSize;
        typedef typename MicrostructureBase<dim>::MobileClusterFluxType MobileClusterFluxType;
        typedef typename MicrostructureBase<dim>::VectorISize VectorISize;


        constexpr static int NdofXnode=NetworkNodeType::NdofXnode;
        static int verboseDislocationNetwork;
        size_t glideStepsSinceLastClimb;
        // During a time step, loop positions are provisional until reconciliation.
        bool trialNodeMotion=false;
        void moveNodes(const double & dt_in);
        void storeSingleGlideStepDiscreteEvents(const long int& runID);
        void executeSingleGlideStepDiscreteEvents(const long int& runID);
        // Optional checkpoint is for fault-injection regression tests. It runs
        // inside the rollback-protected commit after each topology operation.
        void updateBoundaryNodes(const std::function<void()>& checkpoint={});
        bool contract(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB);
        
        std::shared_ptr<DislocationGlideSolverBase<dim>> glideSolver;
        std::shared_ptr<DislocationClimbSolverBase<dim>> climbSolver;

    private:
        
        static std::set<int> getSubCyclingSet(const std::vector<int> &inpVector);
        std::shared_ptr<InclusionMicrostructure<dim>> _inclusions;

    public:

        DefectiveCrystalBase<dim>& ddBase;
        DislocationNetworkRemesh<dim> networkRemesher;
        DislocationJunctionFormation<dim> junctionsMaker;
        const std::shared_ptr<BaseCrossSlipModel<dim>> crossSlipModel;
        DislocationCrossSlip<dim> crossSlipMaker;
        DislocationNodeContraction<dim> nodeContractor;
        DDtimeStepper<dim> timeStepper;
        std::shared_ptr<StochasticForceGenerator> stochasticForceGenerator;
        int ddSolverType;
        int bulkNucleationModel;
        int surfaceNucleationModel;
        bool computeDDinteractions;
        bool outputQuadraturePoints;
        bool outputLinkingNumbers;
        bool outputLoopLength;
        bool outputSegmentPairDistances;
        const bool outputPlasticDistortionPerSlipSystem;
        const bool outputDislocationDensityPerSlipSystem;
        const bool computeElasticEnergyPerLength;
        double alphaLineTension;
        std::set<const LoopNodeType*> danglingBoundaryLoopNodes;
        const bool use_velocityFilter;
        const double velocityReductionFactor;
        const Eigen::Matrix<double,Eigen::Dynamic,dim> nodalVelocityConstraints;
        const std::set<int> subcyclingBins;
        
        const int verboseDislocationNode;
            
        DislocationNetwork(MicrostructureContainerType& mc);
                
        void initializeConfiguration(const DefectiveCrystalConfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels) override;
        void registerSolvedQuantities(SolvedQuantityRegistry&) override;
        void solve() override;
        std::optional<bool> solvedQuantityClimbStep;
        void solveCoupledClimb(bool first);
        Eigen::VectorXd solvedQuantityState() const;
        void restoreSolvedQuantityState(const Eigen::VectorXd&);
        double getDt() const override;
        void output(DefectiveCrystalConfigIO<dim>& configIO,DefectiveCrystalAuxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const override;
        void updateConfiguration() override;
        MatrixDim averagePlasticDistortion() const override ;
        MatrixDim averagePlasticDistortionRate() const override;
        VectorDim displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override;
        MatrixDim stress(const VectorDim&,const NodeType* const,const ElementType* const ele,const SimplexDim* const guess) const override;
        MatrixDim averageStress() const override;
        VectorDim inelasticDisplacementRate(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
        VectorMSize mobileConcentration(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
        MobileClusterFluxType mobileClusterFlux(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;
        VectorISize immobileClusters(const VectorDim&, const NodeType* const, const ElementType* const,const SimplexDim* const) const override;

        void setConfiguration(const DefectiveCrystalConfigIO<dim>&);
        void addConfiguration(const DefectiveCrystalConfigIO<dim>&);

        MatrixDim averagePlasticStrain() const;
        std::map<std::pair<int,int>,double> slipSystemAveragePlasticDistortion() const;
        MatrixDim averagePlasticStrainRate() const;
        void updateGeometry();//
        std::tuple<double,double,double,double> networkLength() const;
        std::vector<std::tuple<double,double,double,double>> networkLengthPerSlipSystem() const;
        bool isClimbStep() const;
        const std::shared_ptr<InclusionMicrostructure<dim>>&  inclusions() const;

    };
    
}
#endif
