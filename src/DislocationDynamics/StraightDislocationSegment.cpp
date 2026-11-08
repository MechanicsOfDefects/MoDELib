/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StraightDislocationSegment_cpp_
#define model_StraightDislocationSegment_cpp_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <limits>
#include <DislocationDynamics/StraightDislocationSegmentHelmholtz.h>
#include <DislocationDynamics/StraightDislocationSegmentFactory.h>
#include <DislocationDynamics/StraightDislocationSegment.h>

namespace model
{
    static_assert(_MODEL_NON_SINGULAR_DD_==1,"Analytical segments require the Cai core");
    
        template <int dim>
        StraightDislocationSegment<dim>::StraightDislocationSegment(const PolycrystallineMaterialBase& mat,
                                   const VectorDim& _P0,
                                   const VectorDim& _P1,
                                   const VectorDim& _b,
                                   const double& _length,
                                   const VectorDim& _t) :
        /* init  */ material(mat)
        /* init  */,P0(_P0)
        /* init  */,P1(_P1)
        /* init  */,b(_b)
        /* init  */,length(_length)
        /* init  */,t(_t)
        /* init  */,bCt(b.cross(t))
        {/*!\param[in] _P0 starting point of the segment
          * \param[in] _P0 ending point of the segment
          * \param[in] _b Burgers vector of the segment
          */
            setCore(CoreType::Cai,8,64,1);
        }

        template <int dim>
        StraightDislocationSegment<dim>::StraightDislocationSegment(const PolycrystallineMaterialBase& mat,
            const VectorDim& start,const VectorDim& end,const VectorDim& burgers)
            : ownedGeometry(std::make_shared<Geometry>(start,end,burgers)),material(mat),
              P0(ownedGeometry->p0),P1(ownedGeometry->p1),b(ownedGeometry->burgers),
              length(ownedGeometry->length),t(ownedGeometry->tangent),bCt(b.cross(t))
        {
            setCore(CoreType::Cai,8,64,1);
        }

        template <int dim>
        typename StraightDislocationSegment<dim>::MatrixDim StraightDislocationSegment<dim>::stress(
            const VectorDim& x,const helmholtzSegment::Evaluator& kernel) const
        {
            static_assert(dim==3,"Helmholtz segment stress is three-dimensional");
            return kernel.stress(x,{P0,P1,b});
        }

        template <int dim>
        void StraightDislocationSegment<dim>::setCore(CoreType type,int order,int resolution,double maxDistance,bool ewald,double ewaldWidth)
        {
            SegmentCoreParameters p;
            p.coreWidth=DislocationFieldBase<dim>::a;
            p.ewaldWidth=ewaldWidth>0?ewaldWidth:4*p.coreWidth;
            p.mu=material.mu;p.nu=material.nu;
            p.order=order;p.resolution=resolution;p.maxDistance=maxDistance;
            core_=makeStraightDislocationSegmentCore(type,ewald,p);
        }

        template <int dim>
        void StraightDislocationSegment<dim>::updateGeometry()
        {
            bCt=b.cross(t);
        }

        template <int dim>
        typename StraightDislocationSegment<dim>::MatrixDim StraightDislocationSegment<dim>::stress(const VectorDim& x) const
        {
        static_assert(dim==3,"Cai segment stress is three-dimensional");
        if(length<=std::numeric_limits<float>::epsilon())
        {
            return MatrixDim::Zero();
        }
        return core_->stress(x,{P0,P1,b});
        }
        
        template <int dim>
        typename StraightDislocationSegment<dim>::VectorDim StraightDislocationSegment<dim>::displacement(const VectorDim& x) const
        {/*!\returns the line-integral part of the displacement contribution of this straight segment.
          * Note: the return value  does NOT include the solid-angle contribution to the displacement field.
          */
            if(length<=std::numeric_limits<float>::epsilon())
            {
                return VectorDim::Zero();
            }
            return core_->displacement(x,{P0,P1,b});
        }
        
        template <int dim>
        double StraightDislocationSegment<dim>::elasticInteractionEnergy(const VectorDim& xA,const VectorDim& tA,const VectorDim& bA) const
        {
            return core_->elasticInteractionEnergy(xA,tA,bA,{P0,P1,b});
        }
        

template class StraightDislocationSegment<3>;
	
}
#endif

