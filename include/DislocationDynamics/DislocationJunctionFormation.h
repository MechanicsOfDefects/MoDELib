/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONJUNCTIONFORMATION_H_
#define model_DISLOCATIONJUNCTIONFORMATION_H_

#include <utility> // for std::pair
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <SegmentSegmentDistance.h>
#include <SweepPlane.h>
#include <TypeTraits.h>
#include <DislocationDynamics/StressStraight.h>

#include <DislocationNetworkRemesh.h>

namespace model
{
    template<int dim> class DislocationNetwork;
    
    template <int dim>
    class DislocationJunctionFormation
    {
        typedef typename TypeTraits<DislocationNetwork<dim>>::NetworkLinkType NetworkLinkType;
        typedef typename TypeTraits<DislocationNetwork<dim>>::NetworkNodeType NetworkNodeType;
        typedef typename TypeTraits<DislocationNetwork<dim>>::LoopNodeType LoopNodeType;
        typedef typename TypeTraits<DislocationNetwork<dim>>::LoopType LoopType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef typename NetworkLinkType::KeyType KeyType;
        typedef std::tuple<KeyType,KeyType,SegmentSegmentDistance<dim>,int> IntersectionType;
        typedef std::deque<IntersectionType> IntersectionTypeContainerType;
        
        void insertIntersection(std::deque<IntersectionTypeContainerType> &intersectionContainer,
                                const NetworkLinkType *const linkA,
                                const NetworkLinkType *const linkB,
                                const SegmentSegmentDistance<dim> &ssd,
                                const double &currentcCollisionTOL,
                                const bool &bndJunction, const bool &gbndJunction,const int event=0);
        friend struct JunctionFormationTestAccess;
        void findPairIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                                   const NetworkLinkType* linkA,const NetworkLinkType* linkB);
        void findIntersections(std::deque<IntersectionTypeContainerType>& intersectionContainer,
                               const size_t& nThreads);
        std::pair<std::shared_ptr<NetworkNodeType>, bool> junctionNode(const double &t,
                                                                       const VectorDim &x,
                                                                       const std::shared_ptr<NetworkLinkType> &L);
        size_t contractJunctions(const std::deque<IntersectionTypeContainerType>& intersectionContainer);
        void glissileJunctions(const double &dx);
        
        void breakJunctions();
        
        DislocationNetwork<dim>& DN;
        
    public:
                
        static double collisionTol;     //! The tolerance (in units of distance) used for collision detection
        const size_t maxJunctionIterations;
        const int verboseJunctions;
        const double infiniteLineLength;
        
        DislocationJunctionFormation(DislocationNetwork<dim>& DN_in);
        static void initFromFile(const std::string& fileName);
        void formJunctions(const double& dx);
    };
    
}
#endif
