/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationClimbSolverFactory_cpp_
#define model_DislocationClimbSolverFactory_cpp_

//#include <PlanarDislocationNode.h>

#include <DislocationClimbSolverFactory.h>
#include <GalerkinClimbSolver.h>
#include <TextFileParser.h>

#include <DislocationNetwork.h>

namespace model
{

    template <int dim>
    DislocationClimbSolverBase<dim>::DislocationClimbSolverBase(const DislocationNetwork<dim>& DN,const ClusterDynamics<dim>* const CD_in) :
    /* init */ DislocationVelocitySolverBase<dim>(DN)
    /* init */,glideEquilibriumRate(TextFileParser(DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("glideEquilibriumRate",true))
    /* init */,CD(CD_in)
    /* init */,vClimbRef(getVclimbRef())
    {
        if(!CD)
        {
            throw std::runtime_error("DislocationClimbSolverBase requires ClusterDynamics.");
        }
    }

    template <int dim>
    const typename DislocationClimbSolverBase<dim>::ScalarVelocitiesContainerType& DislocationClimbSolverBase<dim>::scalarVelocities() const
    {
        return *this;
    }

    template <int dim>
    typename DislocationClimbSolverBase<dim>::ScalarVelocitiesContainerType& DislocationClimbSolverBase<dim>::scalarVelocities()
    {
        return *this;
    }

    template <int dim>
    double DislocationClimbSolverBase<dim>::getVclimbRef() const
    {
        double vRef(0.0);
        const auto eqC(CD->cdp.equilibriumMobileConcentration(0.0));
        for(const auto& pair : CD->cdp.D)
        {
            const auto& dVec(pair.second);
            for(int m=0;m<CD->cdp.mSize;++m)
            {
                vRef=std::max(vRef,dVec[m].trace()*eqC[m]/3.0/this->DN.ddBase.poly.b);
            }
        }
        return  vRef;
    }

    template <int dim>
    std::shared_ptr<DislocationClimbSolverBase<dim>> DislocationClimbSolverFactory<dim>::getClimbSolver(const DislocationNetwork<dim>& DN,const std::string& solverType)
    {
        const auto CD(DN.microstructures.template getUniqueTypedMicrostructure<ClusterDynamics<dim>>());
        if(CD)
        {
            if(solverType=="Galerkin" || solverType=="galerkin")
            {
                
                return std::shared_ptr<DislocationClimbSolverBase<dim>>(new GalerkinClimbSolver<dim>(DN,CD.get()));
            }
            else
            {
                std::cout<<redBoldColor<<"Unknown climb solver type "<<solverType<<". Climb disabled."<<defaultColor<<std::endl;
                return nullptr;
            }
        }
        else
        {
            std::cout<<redBoldColor<<"ClusterDynamics not found. Climb disabled."<<defaultColor<<std::endl;
            return nullptr;
        }
    }

    template struct DislocationClimbSolverBase<3>;
    template struct DislocationClimbSolverFactory<3>;

}
#endif
