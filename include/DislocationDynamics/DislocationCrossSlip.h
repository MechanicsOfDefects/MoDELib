/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationCrossSlip_h_
#define model_DislocationCrossSlip_h_

#include <memory>

#include <DislocationDynamicsModule.h>
#include <CrossSlipModels.h>

#ifndef NDEBUG
#define VerboseCrossSlip(N,x) if(verboseCrossSlip>=N){std::cout<<x;}
#else
#define VerboseCrossSlip(N,x)
#endif

namespace model
{
    template<int dim> class DislocationNetwork;
        
    template <int dim>
    class DislocationCrossSlip
    {
        typedef TypeTraits<DislocationNetwork<dim>> TraitsType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        typedef std::deque<CrossSlipBranch<dim>> CrossSlipBranchContainerType;
        
        friend struct CrossSlipTestAccess;
        DislocationNetwork<dim>& DN; //! A reference to the DislocationNetwork
//        CrossSlipContainerType crossSlipDeq;
        CrossSlipBranchContainerType csNodes;
        
    public:
        
        const int verboseCrossSlip;

        DislocationCrossSlip(DislocationNetwork<dim>& DN_in) ;
        void findCrossSlipSegments(const double dt=-1.0);
        void execute();
        static std::shared_ptr<BaseCrossSlipModel<dim>> getModel(const PolycrystallineMaterialBase& material,const DefectiveCrystalTraitsIO& traitsIO);

        
    };
        
}
#endif
