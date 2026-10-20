/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FrankLoopsGenerator_H_
#define model_FrankLoopsGenerator_H_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
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
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
//#include <MicrostructureGeneratorBase.h>
#include <FrankLoopsDensitySpecification.h>
#include <FrankLoopsIndividualSpecification.h>


namespace model
{
    class FrankLoopsGenerator //: public MicrostructureGeneratorBase
    {
        
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        
        static void generateSingle(MicrostructureGenerator& mg,const int& pID,const VectorDimD& center,const double& radius,const size_t& sides,const bool& isVacancyLoop);

    public:
        FrankLoopsGenerator(const FrankLoopsDensitySpecification&,MicrostructureGenerator&);
        FrankLoopsGenerator(const FrankLoopsIndividualSpecification&,MicrostructureGenerator&);

//        FrankLoopsGenerator(const std::string& fileName);
//        void generateIndividual(MicrostructureGenerator& mg) override;
//        void generateDensity(MicrostructureGenerator& mg) override;
    };
}
#endif
