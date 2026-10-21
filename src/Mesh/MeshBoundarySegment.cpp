/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshBoundarySegment_cpp_
#define model_MeshBoundarySegment_cpp_

#include <cfloat>
#include <set>
#include <Eigen/Dense>
#include <FiniteLineSegment.h>

#include <MeshModule.h>
#include <PolyhedralMesh.h>
#include <limits>
#include <algorithm>

namespace model
{

    template <int dim>
    typename MeshBoundarySegment<dim>::VectorDim MeshBoundarySegment<dim>::periodicShift() const
    {
        VectorDim temp(VectorDim::Zero());
        for(const auto& face : faces)
        {
            if(face->periodicFacePair.second)
            {
                temp+=face->periodicFacePair.first;
            }
        }
        return temp;
    }


    template <int dim>
    typename MeshBoundarySegment<dim>::VectorDim MeshBoundarySegment<dim>::getBoundaryNormal(const typename MeshBoundarySegment<dim>::FaceContainerType& fcs)
    {
        if(fcs.size()==0)
        {
            throw std::runtime_error("EMPY FACE CONTAINER");
        }
        typename MeshBoundarySegment<dim>::VectorDim temp(MeshBoundarySegment<dim>::VectorDim::Zero());
        for(const auto& face : fcs)
        {
            temp+=face->outNormal();
        }
        const double tempNorm(temp.norm());
        return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
    }

template <int dim>
typename MeshBoundarySegment<dim>::VectorDim MeshBoundarySegment<dim>::regionBoundaryNormal(const int& k) const
{
 
    typename MeshBoundarySegment<dim>::VectorDim temp(MeshBoundarySegment<dim>::VectorDim::Zero());
    for(const auto& face : faces)
    {
        temp+=face->outNormal(k);
    }
    const double tempNorm(temp.norm());
    return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : VectorDim::Zero();
}



    template <int dim>
    MeshBoundarySegment<dim>::MeshBoundarySegment(const typename MeshBoundarySegment<dim>::VectorDim& p0,
                                                  const typename MeshBoundarySegment<dim>::VectorDim& p1,
                                                  const typename MeshBoundarySegment<dim>::PlanarMeshFaceType* const face) :
    /* init */ FiniteLineSegmentType(p0,p1)
    /* init */,faces(FaceContainerType({face}))
    /* init */,outNormal(getBoundaryNormal(faces))
    {
    }

    template <int dim>
    MeshBoundarySegment<dim>::MeshBoundarySegment(const typename MeshBoundarySegment<dim>::VectorDim& p0,
                                                  const typename MeshBoundarySegment<dim>::VectorDim& p1,
                                                  const typename MeshBoundarySegment<dim>::FaceContainerType& faces_in) :
    /* init */ FiniteLineSegmentType(p0,p1)
    /* init */,faces(faces_in)
    /* init */,outNormal(getBoundaryNormal(faces))
    {
    }

    namespace
    {
        using SectionKey=std::pair<int,std::vector<size_t>>;
        template<int dim>
        BoundingMeshSegments<dim> grainSection(const Polyhedron<dim,dim>& region,const Plane<dim>& plane)
        {
            using Vector=Eigen::Matrix<double,dim,1>;
            using Faces=typename MeshBoundarySegment<dim>::FaceContainerType;
            using Arc=std::pair<SectionKey,SectionKey>;
            std::map<std::vector<size_t>,double> distances;
            std::map<SectionKey,Vector> points;
            std::map<Arc,Faces> arcs;
            double scale=plane.P.norm();
            for(const auto& f : region.children)
            {
                for(const auto& incidence : f.second->perimeter)
                {
                    scale=std::max(scale,incidence.source()->position.norm());
                }
            }
            const double tolerance=std::max(double(FLT_EPSILON),64*std::numeric_limits<double>::epsilon()*scale);
            auto distance=[&](const std::shared_ptr<Polyhedron<dim,0>>& vertex)
            {
                auto found=distances.find(vertex->key);
                if(found==distances.end())
                {
                    double value=(vertex->position-plane.P).dot(plane.unitNormal);
                    if(std::abs(value)<=tolerance)
                    {
                        value=0;
                    }
                    found=distances.emplace(vertex->key,value).first;
                }
                return found->second;
            };
            std::set<size_t> coincidentFaces;
            for(const auto& f : region.children)
            {
                bool coincident=true;
                for(const auto& incidence : f.second->perimeter)
                {
                    coincident=coincident && distance(incidence.source())==0;
                }
                if(coincident)
                {
                    coincidentFaces.insert(f.first);
                }
            }
            auto addArc=[&](SectionKey a,SectionKey b,const Polyhedron<dim,dim-1>* face)
            {
                if(a==b)
                {
                    return;
                }
                if(b<a)
                {
                    std::swap(a,b);
                }
                arcs[{a,b}].insert(face);
            };
            for(const auto& f : region.children)
            {
                const auto& face=*f.second;
                std::set<SectionKey> roots;
                bool containsEdge=false;
                for(const auto& incidence : face.perimeter)
                {
                    const auto source=incidence.source();
                    const auto sink=incidence.sink();
                    const double d0=distance(source);
                    const double d1=distance(sink);
                    const SectionKey v0{0,source->key};
                    const SectionKey v1{0,sink->key};
                    if(d0==0)
                    {
                        points.emplace(v0,source->position);
                        roots.insert(v0);
                    }
                    if(d0==0 && d1==0)
                    {
                        points.emplace(v1,sink->position);
                        containsEdge=true;
                        const auto other=incidence.otherFace(region.key).face.lock();
                        // An edge internal to coplanar face patches is not part
                        // of the perimeter of the grain section.
                        if(!(coincidentFaces.count(f.first) && coincidentFaces.count(other->key)))
                        {
                            addArc(v0,v1,&face);
                        }
                    }
                    else if((d0<0 && d1>0) || (d0>0 && d1<0))
                    {
                        const SectionKey key{1,incidence.edge->key};
                        // Shared edge identity gives a single computed position.
                        if(!points.count(key))
                        {
                            points.emplace(key,source->position+(d0/(d0-d1))*(sink->position-source->position));
                        }
                        roots.insert(key);
                    }
                }
                if(!containsEdge)
                {
                    if(roots.size()==2)
                    {
                        addArc(*roots.begin(),*roots.rbegin(),&face);
                    }
                    else if(roots.size()>2)
                    {
                        throw std::runtime_error("Grain section requires convex faces: more than two crossings");
                    }
                }
            }
            BoundingMeshSegments<dim> result;
            if(arcs.empty())
            {
                return result;
            }
            std::map<SectionKey,std::vector<SectionKey>> adjacency;
            for(const auto& arc : arcs)
            {
                adjacency[arc.first.first].push_back(arc.first.second);
                adjacency[arc.first.second].push_back(arc.first.first);
            }
            // A tangent edge is a degenerate section, not an open polygon.
            if(arcs.size()==1)
            {
                return result;
            }
            for(const auto& entry : adjacency)
            {
                if(entry.second.size()!=2)
                {
                    // A chain of collinear on-plane edges can also be tangent.
                    const Vector base=points.at(arcs.begin()->first.first);
                    const Vector direction=points.at(arcs.begin()->first.second)-base;
                    bool collinear=true;
                    for(const auto& node : adjacency)
                    {
                        collinear=collinear && direction.cross(points.at(node.first)-base).norm()<=tolerance*direction.norm();
                    }
                    if(collinear)
                    {
                        return result;
                    }
                    throw std::runtime_error("Grain section has a non-manifold intersection vertex");
                }
            }
            const SectionKey start=adjacency.begin()->first;
            SectionKey current=start;
            SectionKey previous;
            std::set<Arc> visited;
            do
            {
                const auto& neighbors=adjacency.at(current);
                const SectionKey next=neighbors[0]==previous ? neighbors[1] : neighbors[0];
                const Arc arc=std::minmax(current,next);
                if(!visited.insert(arc).second)
                {
                    throw std::runtime_error("Grain section repeats an arc before closure");
                }
                result.emplace_back(std::make_shared<MeshBoundarySegment<dim>>(points.at(current),points.at(next),arcs.at(arc)));
                previous=current;
                current=next;
            } while(current!=start);
            if(visited.size()!=arcs.size())
            {
                throw std::runtime_error("Grain section has multiple contours; convex grain required");
            }
            Vector area=Vector::Zero();
            const Vector origin=result.front()->P0;
            for(const auto& segment : result)
            {
                area+=(segment->P0-origin).cross(segment->P1-origin);
            }
            if(area.dot(plane.unitNormal)<0)
            {
                std::reverse(result.begin(),result.end());
                for(auto& segment : result)
                {
                    segment=std::make_shared<MeshBoundarySegment<dim>>(segment->P1,segment->P0,segment->faces);
                }
            }
            // Remove only subdivision points on the same straight supported arc.
            // Connectivity above, never endpoint proximity, determines the walk.
            bool merged=true;
            while(merged && result.size()>2)
            {
                merged=false;
                for(size_t k=0;k<result.size();++k)
                {
                    const size_t next=(k+1)%result.size();
                    const auto a=result[k];
                    const auto b=result[next];
                    const Vector u=a->P1-a->P0;
                    const Vector v=b->P1-b->P0;
                    if(a->faces==b->faces && u.dot(v)>0 && u.cross(v).norm()<=tolerance*(u.norm()+v.norm()))
                    {
                        result[k]=std::make_shared<MeshBoundarySegment<dim>>(a->P0,b->P1,a->faces);
                        result.erase(result.begin()+next);
                        merged=true;
                        break;
                    }
                }
            }
            return result;
        }
    }

    template<int dim>
    BoundingMeshSegments<dim>::BoundingMeshSegments() = default;

    template<int dim>
    BoundingMeshSegments<dim>::BoundingMeshSegments(const SimplicialMesh<dim>& mesh,
                                                    const size_t& regionID,const Plane<dim>& plane)
    {
        const auto section=grainSection<dim>(mesh.polyhedralMesh().region(regionID),plane);
        this->assign(section.begin(),section.end());
    }

    /**********************************************************************/
    template <int dim>
    std::set<const MeshBoundarySegment<dim>*> BoundingMeshSegments<dim>::containingSegments(const VectorDim& P) const
    {
        std::set<const MeshBoundarySegment<dim>*> temp;
        
        for(const auto& seg : *this)
        {
            if(seg->contains(P))
            {
                temp.insert(seg.get());
            }
        }
        return temp;
    }

    /**********************************************************************/
    template <int dim>
    bool BoundingMeshSegments<dim>::contains(const VectorDim& P) const
    {
        return containingSegments(P).size();
    }

    /**********************************************************************/
    template <int dim>
    const BoundingMeshSegments<dim>& BoundingMeshSegments<dim>::boundingBoxSegments() const
    {
        return *this;
    }

    template <int dim>
    BoundingMeshSegments<dim>& BoundingMeshSegments<dim>::boundingBoxSegments()
    {
        return *this;
    }

    template class MeshBoundarySegment<3>;
    template struct BoundingMeshSegments<3>;

}
#endif
