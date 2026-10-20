/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDipoleGenerator_H_
#define model_PeriodicDipoleGenerator_H_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <limits>
#include <Eigen/LU>
#include <Eigen/Cholesky>

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
#include <PeriodicDipoleDensitySpecification.h>

namespace model
{
    class PeriodicDipoleGenerator
    {
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        static bool generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& dipolePoint,const int& exitFaceID,const int& dipoleHeight,const int& dipoleNodes,double glideStep);
    public:
        PeriodicDipoleGenerator(const PeriodicDipoleDensitySpecification& spec,MicrostructureGenerator& mg);
        PeriodicDipoleGenerator(const PeriodicDipoleIndividualSpecification& spec,MicrostructureGenerator& mg);
    };
}
#endif
