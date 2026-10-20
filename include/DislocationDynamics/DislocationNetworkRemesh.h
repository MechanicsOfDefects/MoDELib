/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetworkRemesh_h_
#define model_DislocationNetworkRemesh_h_

#include <DislocationDynamicsModule.h>

namespace model
{
    template<int dim> class DislocationNetwork;
    
    /*! \brief Class template that handles the nodal remesh of the DislocationNetwork.
     */
    template <int dim>
    class DislocationNetworkRemesh
    {
        
        typedef TypeTraits<DislocationNetwork<dim>> TraitsType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::VectorLowerDim VectorLowerDim;

        
        DislocationNetwork<dim>& DN;
        friend struct RemeshTestAccess;
        bool zeroChordCandidate(const std::shared_ptr<typename TraitsType::NetworkLinkType>& link) const;
        size_t contractZeroChordCandidates(const std::vector<std::pair<size_t,size_t>>& candidates);
        
    public:
        
        const double Lmax;
        const double Lmin;
        const double absoluteAreaThreshold;
        const double relativeAreaThreshold;
        const short unsigned int remeshFrequency;
        
        DislocationNetworkRemesh(DislocationNetwork<dim>& DN_in);
        void remesh(const long int &runID);
        void remeshByRemoval();
        size_t simplifyPhysicalBoundaryChains();
        void remeshByExpansion();
        void contractBoundaryNodes(); //This function contracts boundary nodes if they are at the same position and share atleast one same neightbors
        void contract0chordSegments();
        void remove0AreaLoopAcrossBnd();
        void removeCollapsedLoops();
//        static double minMeshSize(const SimplicialMesh<dim> &mesh);
    };
    
} // namespace model
#endif

