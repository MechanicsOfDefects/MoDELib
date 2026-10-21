/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeContraction_H_
#define model_DislocationNodeContraction_H_

#include <memory>
#include <Eigen/Dense>
#include <TypeTraits.h>

namespace model
{
    template<int dim> class DislocationNetwork;
    template <int dim>
    class DislocationNodeContraction
    {
        typedef typename TypeTraits<DislocationNetwork<dim>>::NetworkNodeType NetworkNodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        std::vector<std::set<std::shared_ptr<NetworkNodeType>>> contractBoundary(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB);
        bool contractSecondAndVirtual(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB,const bool reconcileBoundary);
        bool contractYoungest(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB,const bool reconcileBoundary);
        bool contractToPosition(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB,const VectorDim& X, const double& maxRange,const bool reconcileBoundary);

        bool canReach(const std::shared_ptr<NetworkNodeType>& node,const VectorDim& target) const;
        bool commitContraction(std::shared_ptr<NetworkNodeType> a,std::shared_ptr<NetworkNodeType> b,const VectorDim& target,const bool reconcileBoundary);
        DislocationNetwork<dim>& DN;

    public:
        
        const int verboseNodeContraction;

        DislocationNodeContraction(DislocationNetwork<dim>& DN_in);
        bool contract(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB,const bool reconcileBoundary=true);
    };
}
#endif

