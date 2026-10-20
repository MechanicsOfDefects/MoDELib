/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GalerkinGlideSolver_h_
#define model_GalerkinGlideSolver_h_

//#include <PlanarDislocationNode.h>


#ifndef NDEBUG
#define VerboseGalerkinGlideSolver(N,x) if(verboseGalerkinGlideSolver>=N){std::cout<<x;}
#else
#define VerboseGalerkinGlideSolver(N,x)
#endif

#include <DislocationDynamicsModule.h>
#include <DislocationGlideSolverFactory.h>
namespace model
{
    template<int dim> class DislocationNetwork;
    
    template <int dim>
//    template<int dim,int corder>
    class GalerkinGlideSolver : public DislocationGlideSolverBase<dim>
    {
//        typedef DislocationNetwork<dim,corder> DislocationNetwork<dim>;
        typedef typename DislocationNetwork<dim>::VectorDim   VectorDim;
        typedef typename DislocationNetwork<dim>::MatrixDim   MatrixDim;
        typedef typename DislocationNetwork<dim>::NetworkNodeType   NetworkNodeType;
        typedef typename DislocationNetwork<dim>::LoopNodeType   LoopNodeType;
        typedef typename DislocationNetwork<dim>::LoopType   LoopType;

        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        //        typedef std::vector<Eigen::Triplet<double> > TripletContainerType;
        typedef std::deque<Eigen::Triplet<double> > TripletContainerType;
//        enum {NdofXnode=NodeType::NdofXnode};
        constexpr static int NdofXnode=NetworkNodeType::NdofXnode;


//        void lumpedSolve();
//        void solve(const size_t& runID);
        Eigen::VectorXd lumpedSolve() const;
        size_t assembleNCtriplets(TripletContainerType& kqqT, Eigen::VectorXd& Fq) const;
//        void storeNodeSolution(const Eigen::VectorXd& X);
        size_t assemblePeriodicConstraints(TripletContainerType& zT) const;

        
        public:
        
        GalerkinGlideSolver(const DislocationNetwork<dim>& );
        Eigen::VectorXd getNodeVelocities() const override;

        
    };
    
}
#endif
