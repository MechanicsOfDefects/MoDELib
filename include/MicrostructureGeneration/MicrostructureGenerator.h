/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_H_
#define model_MicrostructureGenerator_H_

#ifdef _MODEL_PYBIND11_
    #undef slots
    #include <pybind11/embed.h>
    #include <pybind11/eigen.h>
    #include <pybind11/stl.h>
    #define slots Q_SLOTS
#endif

#include <filesystem>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>
#include <memory>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
#include <DefectiveCrystal/DefectiveCrystalTraitsIO.h>
#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <DefectiveCrystal/DefectiveCrystalAuxIO.h>

#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <DislocationDynamics/DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationDynamics/DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
//#include <MicrostructureGeneratorBase.h>

#include <ShearLoopDensitySpecification.h>
#include <ShearLoopIndividualSpecification.h>
#include <PeriodicDipoleDensitySpecification.h>
#include <PeriodicDipoleIndividualSpecification.h>
#include <PrismaticLoopDensitySpecification.h>
#include <PrismaticLoopIndividualSpecification.h>
#include <FrankLoopsDensitySpecification.h>
#include <FrankLoopsIndividualSpecification.h>
#include <StackingFaultTetrahedraDensitySpecification.h>
#include <StackingFaultTetrahedraIndividualSpecification.h>
#include <SphericalInclusionDensitySpecification.h>
#include <SphericalInclusionIndividualSpecification.h>
#include <EllipsoidalInclusionIndividualSpecification.h>
#include <EllipsoidalInclusionDensitySpecification.h>
#include <PolyhedronInclusionIndividualSpecification.h>
#include <PlanarLoopIndividualSpecification.h>
#include <CircularCrackIndividualSpecification.h>
#include <CrackMeshIndividualSpecification.h>

namespace model
{

struct PolyPoint
{
    
    std::shared_ptr<PeriodicPlanePatch<3>> periodicPlanePatch() const;



};

    class MicrostructureGenerator
    {
        constexpr static int dim=3;
//        typedef typename  MicrostructureGeneratorBase::VectorDimD VectorDimD;
//        typedef typename  MicrostructureGeneratorBase::DislocationLoopType DislocationLoopType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<int,dim,1> VectorDimI;
        typedef DislocationLoopIO<dim>::DislocationLoopType DislocationLoopType;

    public:
        
        DefectiveCrystalBase<3>& ddBase;
        DefectiveCrystalConfigIO<3> configIO;
        DefectiveCrystalAuxIO<3> auxIO;
        const bool outputBinary;
        const double minSize;
        const double maxSize;
        std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap;
        const DefectiveCrystalTraitsIO& traits() const;
        const DefectiveCrystalConfigIO<3>& config() const;
        const DefectiveCrystalAuxIO<3>& aux() const;
        DefectiveCrystalConfigIO<3>& config();
        DefectiveCrystalAuxIO<3>& aux();

        MicrostructureGenerator(DefectiveCrystalBase<3>& ddBase_in);
        void readMicrostructureFile();
        void addShearLoopDensity(const ShearLoopDensitySpecification& spec);
        void addShearLoopIndividual(const ShearLoopIndividualSpecification& spec);
        void addPeriodicDipoleDensity(const PeriodicDipoleDensitySpecification& spec);
        void addPeriodicDipoleIndividual(const PeriodicDipoleIndividualSpecification& spec);
        void addPrismaticLoopDensity(const PrismaticLoopDensitySpecification& spec);
        void addPrismaticLoopIndividual(const PrismaticLoopIndividualSpecification& spec);
        void addFrankLoopsDensity(const FrankLoopsDensitySpecification& spec);
        void addFrankLoopsIndividual(const FrankLoopsIndividualSpecification& spec);
        void addStackingFaultTetrahedraDensity(const StackingFaultTetrahedraDensitySpecification& spec);
        void addStackingFaultTetrahedraIndividual(const StackingFaultTetrahedraIndividualSpecification& spec);
        void addSphericalInclusionDensity(const SphericalInclusionDensitySpecification& spec);
        void addSphericalInclusionIndividual(const SphericalInclusionIndividualSpecification& spec);
        void addEllipsoidalInclusionDensity(const EllipsoidalInclusionDensitySpecification& spec);
        void addEllipsoidalInclusionIndividual(const EllipsoidalInclusionIndividualSpecification& spec);
        void addPolyhedronInclusionIndividual(const PolyhedronInclusionIndividualSpecification& spec);
        void addPlanarLoopIndividual(const PlanarLoopIndividualSpecification& spec);
        void addCircularCrackIndividual(const CircularCrackIndividualSpecification& spec);
        void addCrackMeshIndividual(const CrackMeshIndividualSpecification& spec);

        size_t insertLoop(const VectorDimD& b,const VectorDimD& unitNormal,const VectorDimD& P0,const size_t& grainID,const DislocationLoopType& loopType);
        size_t insertLoopNode(const size_t& loopID,const VectorDimD& loopNodePos,const size_t& networkNodeID,const VectorDimD& loopNodeShift,const std::pair<short int,short int>& periodicEdgeIDs);
        std::vector<size_t> insertLoopLinks(const size_t& loopID,const std::vector<size_t>& loopNodeIDs);
        size_t insertNetworkNode(const VectorDimD& networkNodePos);
        size_t insertInclusion(const VectorDimD& pos,const double& R, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type);
        size_t insertInclusion(const VectorDimD& pos,const Eigen::Matrix3d& axes,const Eigen::Matrix3d& eT,const double& vrc,const int& type);
        size_t insertInclusion(const std::map<size_t,Eigen::Vector3d>& nodes,const std::map<size_t,std::vector<size_t>>& faceMap, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type);
        void writeConfigFiles(const size_t& fileID);
        bool insertJunctionLoop(const std::vector<VectorDimD>& loopNodePos,
                                const std::shared_ptr<PeriodicGlidePlane<3>>& periodicPlane,
                                const VectorDimD& b,
                                const VectorDimD& unitNormal,
                                const VectorDimD& P0,
                                const size_t& grainID,
                                const DislocationLoopIO<dim>::DislocationLoopType& loopType);
        bool insertCrackMesh(const std::vector<VectorDimD>& nodes,const std::vector<VectorDimD>& burgers,const std::vector<VectorDimI>& triangles);
        bool allPointsInGrain(const std::vector<VectorDimD>& points,const int& grainID);
        bool allPointsInMesh(const std::vector<VectorDimD>& points);

    };
}
#endif
