/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PrismaticLoopGenerator_H_
#define model_PrismaticLoopGenerator_H_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationDynamics/DislocationNodeIO.h>
#include <DislocationDynamics/DislocationLoopIO.h>
#include <DislocationDynamics/DislocationLoopLinkIO.h>
#include <DislocationDynamics/DislocationLoopNodeIO.h>
#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <DefectiveCrystal/DefectiveCrystalAuxIO.h>

#include <DislocationDynamics/DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationDynamics/DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
//#include <MicrostructureGeneratorBase.h>
#include <PrismaticLoopDensitySpecification.h>
#include <PrismaticLoopIndividualSpecification.h>

namespace model
{

    class PrismaticLoopGenerator //: public MicrostructureGeneratorBase
    {
        
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        
        static double generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& center,const double& radius,const double& step);
        
    public:
        
        PrismaticLoopGenerator(const PrismaticLoopDensitySpecification& spec,MicrostructureGenerator& mg);
        PrismaticLoopGenerator(const PrismaticLoopIndividualSpecification& spec,MicrostructureGenerator& mg);

//        PrismaticLoopGenerator(const std::string& fileName);
//
////        void generate(MicrostructureGenerator& mg) override;
//        void generateIndividual(MicrostructureGenerator& mg) override;
//        void generateDensity(MicrostructureGenerator& mg) override;

    };

}
#endif
