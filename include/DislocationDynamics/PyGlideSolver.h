/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PyGlideSolver_h_
#define model_PyGlideSolver_h_

//#include <PlanarDislocationNode.h>


#ifndef NDEBUG
#define VerbosePyGlideSolver(N,x) if(verbosePyGlideSolver>=N){std::cout<<x;}
#else
#define VerbosePyGlideSolver(N,x)
#endif

#include <iostream>

#include <DislocationDynamicsModule.h>
#include <DislocationGlideSolverFactory.h>

#ifdef _MODEL_PYBIND11_
#include <pybind11/embed.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#endif

namespace model
{
    template<int dim> class DislocationNetwork;
    
    template <int dim>
    class PyGlideSolver : public DislocationGlideSolverBase<dim>
    {
        typedef typename DislocationNetwork<dim>::VectorDim   VectorDim;
        typedef typename DislocationNetwork<dim>::MatrixDim   MatrixDim;
        typedef typename DislocationNetwork<dim>::NetworkNodeType   NetworkNodeType;
        typedef typename DislocationNetwork<dim>::LoopNodeType   LoopNodeType;
        typedef typename DislocationNetwork<dim>::LoopType   LoopType;
        constexpr static int NdofXnode=NetworkNodeType::NdofXnode;

        public:
        
        const std::string pyModuleName;
        
#ifdef _MODEL_PYBIND11_
        pybind11::module pyModule;
#endif

        Eigen::MatrixXd updateUnconnectedNodePythonVelocity(const Eigen::MatrixXd& nodePythonVelocity, const std::map<size_t, size_t>& removedNodeIDmap) const;
        Eigen::MatrixXd updateBoundaryNodePythonVelocity(const Eigen::MatrixXd& nodePythonVelocity, const std::map<size_t, size_t>& removedNodeIDmap) const;

        PyGlideSolver(const DislocationNetwork<dim>&);
        Eigen::VectorXd getNodeVelocities() const override;
        
    };
    
}
#endif
