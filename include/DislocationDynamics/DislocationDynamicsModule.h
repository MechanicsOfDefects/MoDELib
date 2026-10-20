/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationDynamicsModule_h_
#define model_DislocationDynamicsModule_h_

//Generic includes
#include <iostream>
#include <iomanip>
#include <chrono>
#include <memory>
#include <utility> 
#include <assert.h>

//STL includes
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <algorithm>

//MPIInclude
//#include <MPIcout.h>

//Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

//MoDELib Utilities
#include <TypeTraits.h>
#include <TextFileParser.h>
#include <TerminalColors.h>

//MoDELib Mesh Module
#include <Simplex.h>
#include <PlanarMeshFace.h>

//MoDELib Lattice Module
#include <LatticeModule.h>






//DDD includes



//Topology layer
#include <LoopLink.h>
#include <Loop.h>
#include <LoopNode.h>
#include <NetworkNode.h>


//Geometry layer
#include <SplineNode.h>

//Math functions
#include <GramSchmidt.h>


//Geometry functions
#include <SweepPlane.h>
#include <SegmentSegmentDistance.h>
#include <FiniteLineSegment.h>
#include <LineLineIntersection.h>
#include <PlanePlaneIntersection.h>
#include <DislocationDynamics/DislocationLoopPatches.h>


namespace model
{

    template <int dim>
    class DislocationNetwork;

    template <int dim>
    struct DislocationVelocitySolverBase;

    template <int dim>
    struct DislocationGlideSolverBase;

    template <int dim>
    struct DislocationClimbSolverBase;

    template <int dim>
    class DislocationJunctionFormation;

    template <int dim>
    struct DDtimeStepper;

    template <int dim>
    class DefectiveCrystalConfigIO;

    struct MeshedDislocationLoop;

}


//Physics layer
#include <DislocationNetworkTraits.h>
#include <GlidePlaneModule.h>
#include <GlidePlane.h>
#include <Grain.h>
//#include <ConfinedDislocationObject.h>
#include <DislocationCrossSlip.h>
#include <DislocationJunctionFormation.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <DefectiveCrystal/DefectiveCrystalAuxIO.h>
#include <DislocationLoop.h>
#include <DislocationLoopNode.h>
#include <DislocationLoopLink.h>
#include <DislocationNode.h>
#include <DislocationSegment.h>
//#include <DislocationGlideSolver.h>
#include <DislocationNetworkRemesh.h>

//#include <DislocationNetworkRemesh.cpp>




#endif
