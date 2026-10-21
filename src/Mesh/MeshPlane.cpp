/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshPlane_cpp_
#define model_MeshPlane_cpp_

#include <numbers>
#include <cfloat>
#include <tuple>
#include <vector>
#include <Eigen/Dense>

#include <MeshModule.h>
#include <PolyhedralMesh.h>

namespace model
{
    
    /**********************************************************************/
    template <int dim>
    BoundingMeshSegments<dim> MeshPlane<dim>::getFaceBoundary(const Polyhedron<dim,dim-1>& face)
    {
        BoundingMeshSegments<dim> temp;
        for(const auto& incidence : face.perimeter)
        {
            temp.emplace_back(std::make_shared<MeshBoundarySegment<dim>>(
                incidence.source()->position,incidence.sink()->position,&face));
        }
        return temp;
    }
    
    /**********************************************************************/
    template <int dim>
    void MeshPlane<dim>::checkPlaneIntersections(const BoundingMeshSegments<dim>& temp)
    {
        
        if(temp.size()<3)
        {
            throw std::runtime_error("MeshPlane: Sorting MeshIntersections failed. temp.size()="+std::to_string(temp.size()));
        }
        
        for(const auto& seg : temp)
        {
            if(seg->hasZeroLength())
            {
                throw std::runtime_error("Plane-Face intersection has zero length");
            }
        }
    }
    
    /**********************************************************************/
    template <int dim>
    MeshPlane<dim>::MeshPlane(const SimplicialMesh<dim>& mesh,
                              const size_t& rID,
                              const VectorDim& p,
                              const VectorDim& n) :
    /* init */ Plane<dim>(p,n)
    /* init */,regionIDs(rID,rID)
    //        /* init */,meshIntersections(mesh,rID,*this)
    /* init */,meshIntersections(BoundingMeshSegments<dim>(mesh,rID,*this))
    {/*!\param[in] mesh
      * \param[in] rID the region ID where the plane is defined
      * \param[in] p position of the plane
      * \param[in] n normal to the plane
      * Constructor for plane internal to a mesh region
      */
        checkPlaneIntersections(meshIntersections);

    }
    
    /**********************************************************************/
    template <int dim>
    MeshPlane<dim>::MeshPlane(const Polyhedron<dim,dim-1>& face,
                              const size_t& rID1,
                              const size_t& rID2) :
    /* init */ Plane<dim>(face.asPlane())
    /* init */,regionIDs(rID1,rID2)
    //        /* init */ meshIntersections(getFaceBoundary(face)) // WARNING: CALLING meshIntersections with rID1
    /* init */,meshIntersections(getFaceBoundary(face))
    {
        checkPlaneIntersections(meshIntersections);
    }
    
    //        /**********************************************************************/
    //        template <class T>
    //        friend T& operator << (T& os, const MeshPlane<dim>& gp)
    //        {
    //            for (const auto& x : gp.meshIntersections)
    //            {
    //                os<<gp.sID<<" "<<*x;
    //            }
    //            return os;
    //        }
    
    template struct MeshPlane<3>;
    
}
#endif
