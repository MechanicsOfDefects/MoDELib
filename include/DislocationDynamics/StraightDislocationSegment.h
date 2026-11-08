/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StraightDislocationSegment_H_
#define model_StraightDislocationSegment_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
#include <memory>
#include <DislocationDynamics/CoreType.h>
#include <DislocationDynamics/StraightDislocationSegmentBase.h>
#include <PolycrystallineMaterialBase.h>
#include <DislocationDynamics/DislocationFieldBase.h>

namespace model
{
    namespace helmholtzSegment { class Evaluator; }
    
    
    template <int dim>
    class StraightDislocationSegment
    {
        typedef double Scalar;
        typedef Eigen::Matrix<Scalar,dim,dim> MatrixDim;
        typedef Eigen::Matrix<Scalar,dim,1>   VectorDim;
        
        // Owned snapshots share immutable geometry so copies retain valid references.
        struct Geometry
        {
            VectorDim p0,p1,burgers,tangent;
            double length;
            Geometry(const VectorDim& start,const VectorDim& end,const VectorDim& burgersIn)
                : p0(start),p1(end),burgers(burgersIn),tangent(VectorDim::Zero()),length((end-start).norm())
            {
                if(length>0)
                {
                    tangent=(end-start)/length;
                }
            }
        };
        std::shared_ptr<const Geometry> ownedGeometry;

        
    public:
        
        const PolycrystallineMaterialBase& material;
        const VectorDim& P0;
        const VectorDim& P1;
        const VectorDim& b;
        const double& length;
        const VectorDim& t;
        
    private:

        VectorDim bCt;
        std::shared_ptr<const StraightDislocationSegmentBase> core_;
        
    public:

        // Standalone snapshot (geometry owned); material must outlive the segment.
        StraightDislocationSegment(const PolycrystallineMaterialBase& mat,
                                   const VectorDim& start,const VectorDim& end,const VectorDim& burgers);

        // Live network view; updateGeometry() refreshes cached products after motion.
        StraightDislocationSegment(const PolycrystallineMaterialBase& mat,
                                   const VectorDim& _P0,
                                   const VectorDim& _P1,
                                   const VectorDim& _b,
                                   const double& _length,
                                   const VectorDim& _t);

        void setCore(CoreType core,int order,int resolution,double maxDistance,bool ewald=false,double ewaldWidth=0);
        const StraightDislocationSegmentBase& core() const { return *core_; }
        void updateGeometry();
        MatrixDim stress(const VectorDim& x) const;
        // The supplied evaluator selects full/screened Helmholtz and shares cached tables.
        MatrixDim stress(const VectorDim& x,const helmholtzSegment::Evaluator& kernel) const;
        VectorDim displacement(const VectorDim& x) const;        
        double elasticInteractionEnergy(const VectorDim& xA,const VectorDim& tA,const VectorDim& bA) const;
        
	};	
	
}
#endif

