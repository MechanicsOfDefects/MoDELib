/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureBase_cpp_
#define model_MicrostructureBase_cpp_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fstream>
#include <map>
#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <DefectiveCrystal/DefectiveCrystalAuxIO.h>
#include <DefectiveCrystal/MicrostructureBase.h>
#include <DefectiveCrystal/MicrostructureContainer.h>

namespace model
{
    template<int dim>
    void MicrostructureBase<dim>::registerSolvedQuantities(SolvedQuantityRegistry& registry)
    {
        SolvedQuantityContribution contribution;
        contribution.provider=tag;
        contribution.name="State";
        contribution.quantities={"Stress"};
        contribution.update=[this](const SolvedQuantityContext&)
        {
            solve();
        };
        registry.add(std::move(contribution));
    }

    template <int dim>
    MicrostructureBase<dim>::MicrostructureBase(const std::string& tag_in,MicrostructureContainer<dim>& m_in):
    /* init */ tag(tag_in)
    /* init */,microstructures(m_in)
    {
        // The container's shared-base reference is not initialized during this constructor.
        if(tag!="MicrostructureContainer" && !microstructures.ddBase.hasPhysics(tag))
        {
            throw std::runtime_error(tag+" physics is not enabled.");
        }
    }

    template <int dim>
    double MicrostructureBase<dim>::timeSinceLastUpdate() const
    {
        return microstructures.ddBase.simulationParameters.totalTime-lastUpdateTime;
    }

    template <int dim>
    Eigen::Matrix<double,Eigen::Dynamic,dim> MicrostructureBase<dim>::displacement(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,dim>> points) const
    {
        Eigen::Matrix<double,Eigen::Dynamic,dim> temp(Eigen::Matrix<double,Eigen::Dynamic,3>::Zero(points.rows(),dim));
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(long int k=0;k<points.rows();++k)
        {
            temp.row(k)=displacement(points.row(k),nullptr,nullptr,nullptr);
        }
        return temp;
    }

template <int dim>
std::vector<Eigen::Matrix<double,dim,dim>> MicrostructureBase<dim>::stress(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,dim>> points) const
{
    std::vector<Eigen::Matrix<double,dim,dim>> temp(points.rows(),Eigen::Matrix<double,dim,dim>::Zero());
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(long int k=0;k<points.rows();++k)
    {
        temp[k]=stress(points.row(k),nullptr,nullptr,nullptr);
    }
    return temp;
}

template <int dim>
std::set<const Grain<dim>*> MicrostructureBase<dim>::pointGrains(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
{
    std::set<const Grain<dim>*> temp;
    if(node)
    {
        for(const auto& nodeEle : *node)
        {
            temp.emplace(microstructures.ddBase.poly.grain(nodeEle->simplex.region->regionID).get());
        }
    }
    else
    {
        if(ele)
        {
            temp.emplace(microstructures.ddBase.poly.grain(ele->simplex.region->regionID).get());
        }
        else
        {
            const std::pair<bool,const Simplex<dim,dim>*> found(microstructures.ddBase.mesh.searchWithGuess(x,guess));
            if(found.first)
            {
                temp.emplace(microstructures.ddBase.poly.grain(found.second->region->regionID).get());
            }
        }
    }
    return temp;
}

    template struct MicrostructureBase<3>;

}
#endif
