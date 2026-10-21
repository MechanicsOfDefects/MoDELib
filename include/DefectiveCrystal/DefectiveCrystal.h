/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_H_
#define model_DefectiveCrystal_H_

#include <fstream>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <DefectiveCrystal/DefectiveCrystalParameters.h>


namespace model
{
    template<int dim> class DislocationNetwork;
    template<int dim> struct ClusterDynamics;
    template<int dim> struct ElasticDeformation;
    template<int dim> struct InclusionMicrostructure;

    
    template <int _dim>
    class DefectiveCrystal : public MicrostructureContainer<_dim>
    {
        
        std::ofstream f_file;
        std::ofstream F_labels;

    public:
        static constexpr int dim=_dim; // make dim available outside class
        typedef typename MicrostructureBase<dim>::MatrixDim MatrixDim;
        typedef typename MicrostructureBase<dim>::VectorDim VectorDim;
        typedef MicrostructureContainer<_dim> MicrostructureContainerType;
        
        typedef InclusionMicrostructure<dim> InclusionMicrostructureType;
        typedef DislocationNetwork<dim> DislocationNetworkType;
        typedef ClusterDynamics<dim> ClusterDynamicsType;
        typedef ElasticDeformation<dim> ElasticDeformationType;
        
        DefectiveCrystal(DefectiveCrystalBase<_dim>& ddBase_in) ;
        using MicrostructureContainer<_dim>::initializeConfiguration;
        void initializeConfiguration(const DefectiveCrystalConfigIO<dim>& configIO);
        void runSteps();
        void runSingleStep();
        
        const DislocationNetwork<_dim>& dislocationNetwork() const;
        
    };
}
#endif
