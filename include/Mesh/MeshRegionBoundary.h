/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshRegionBoundary_H_
#define model_MeshRegionBoundary_H_

#include <set>
#include <deque>
#include <assert.h>

#include <Polyhedron.h>

namespace model
{

    template<int dim>
    struct MeshRegionBoundary : public std::set<const Simplex<dim,dim-1>*>
    /*                       */,public std::map<size_t,std::shared_ptr<Polyhedron<dim,dim-1>>>
    {
        typedef Simplex<dim,dim-1> SimplexType;
        typedef Polyhedron<dim,dim-1> PlanarMeshFaceType;
        typedef std::map<size_t,std::shared_ptr<PlanarMeshFaceType>> MeshFacesContainerType;
        
    private:
        
        PolyhedralMesh<dim>& topology;
        void buildSingleFace(const Simplex<dim,dim-1>* newStart,
                             std::shared_ptr<Polyhedron<dim,dim-1>> newFace,
                             std::set<const Simplex<dim,dim-1>*>& allSimplices);
        void buildFaces();
                
    public:
        
        const std::pair<int,int> regionBndID;

        MeshRegionBoundary(const std::pair<size_t,size_t>& rbndID, PolyhedralMesh<dim>& topologyIn);
        void update();
        const std::set<const SimplexType*>& simplices() const;
        std::set<const SimplexType*>& simplices();
        const MeshFacesContainerType& faces() const;
        MeshFacesContainerType& faces();
        
    };
    
}	// close namespace
#endif
