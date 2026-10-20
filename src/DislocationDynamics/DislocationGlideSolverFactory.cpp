/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationGlideSolverFactory_cpp_
#define model_DislocationGlideSolverFactory_cpp_

//#include <PlanarDislocationNode.h>

#include <DislocationGlideSolverFactory.h>
#include <GalerkinGlideSolver.h>
#include <PyGlideSolver.h>

#include <DislocationNetwork.h>

namespace model
{
    
    template <int dim>
    DislocationGlideSolverBase<dim>::DislocationGlideSolverBase(const DislocationNetwork<dim>& DN) :
//    /* init */ DN(DN_in)
    /* init */ DislocationVelocitySolverBase<dim>(DN)
    {
        
    }

template <int dim>
std::shared_ptr<DislocationGlideSolverBase<dim>> DislocationGlideSolverFactory<dim>::getGlideSolver(const DislocationNetwork<dim>& DN,const std::string& solverType)
{
        if(solverType=="Galerkin" || solverType=="galerkin")
        {
            return std::shared_ptr<DislocationGlideSolverBase<dim>>(new GalerkinGlideSolver<dim>(DN));
        }
        else if(solverType=="Pybind11" || solverType=="pybind11")
        {
            return std::shared_ptr<DislocationGlideSolverBase<dim>>(new PyGlideSolver<dim>(DN));
        }
    else
    {
        std::cout<<redBoldColor<<"Unknown glide solver type "<<solverType<<". Glide disabled."<<defaultColor<<std::endl;
        return nullptr;
    }
    
}
    
    template struct DislocationGlideSolverBase<3>;
    template struct DislocationGlideSolverFactory<3>;

}
#endif
