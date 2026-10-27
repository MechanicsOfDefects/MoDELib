

add_subdirectory(ellipsoidalInclusion)

if(TARGET testSurfaceMesh)
    add_test(NAME surfaceMesh COMMAND testSurfaceMesh "${PROJECT_SOURCE_DIR}")
    set_property(TEST surfaceMesh APPEND PROPERTY LABELS "integration")
    add_test(NAME surfaceMeshWithoutDD COMMAND testSurfaceMesh "${PROJECT_SOURCE_DIR}" --without-dd)
    set_property(TEST surfaceMeshWithoutDD APPEND PROPERTY LABELS "integration")
    set_tests_properties(surfaceMeshWithoutDD PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 60)
    set_tests_properties(surfaceMesh PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 60)
endif()

if(TARGET testWrappingTabWidget)
    add_test(NAME wrappingTabWidget COMMAND testWrappingTabWidget)
    set_property(TEST wrappingTabWidget APPEND PROPERTY LABELS "integration")
    set_tests_properties(wrappingTabWidget PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 30)
endif()

if(TARGET testDDqtFrameWorkflow)
    add_test(NAME DDqtFrameWorkflow COMMAND testDDqtFrameWorkflow "${PROJECT_SOURCE_DIR}")
    set_property(TEST DDqtFrameWorkflow APPEND PROPERTY LABELS "integration")
    set_tests_properties(DDqtFrameWorkflow PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 60)
endif()

add_executable(testPeriodicMesh testPeriodicMesh.cpp)
target_link_libraries(testPeriodicMesh PRIVATE MoDELib)
add_test(NAME periodicMesh COMMAND testPeriodicMesh "${PROJECT_SOURCE_DIR}")
set_property(TEST periodicMesh APPEND PROPERTY LABELS "integration")
set_tests_properties(periodicMesh PROPERTIES TIMEOUT 120)

add_executable(testPeriodicFieldInitialization testPeriodicFieldInitialization.cpp)
target_link_libraries(testPeriodicFieldInitialization PRIVATE MoDELib)
add_test(NAME periodicFieldInitialization COMMAND testPeriodicFieldInitialization "${PROJECT_SOURCE_DIR}")
set_property(TEST periodicFieldInitialization APPEND PROPERTY LABELS "integration")
set_tests_properties(periodicFieldInitialization PROPERTIES TIMEOUT 120)

add_executable(testPeriodicPolygonWalk testPeriodicPolygonWalk.cpp)
target_link_libraries(testPeriodicPolygonWalk PRIVATE MoDELib)
add_test(NAME periodicPolygonWalk COMMAND testPeriodicPolygonWalk "${PROJECT_SOURCE_DIR}")
set_property(TEST periodicPolygonWalk APPEND PROPERTY LABELS "integration")
set_tests_properties(periodicPolygonWalk PROPERTIES TIMEOUT 120)

add_executable(testBoundaryNodeUpdate testBoundaryNodeUpdate.cpp)
target_link_libraries(testBoundaryNodeUpdate PRIVATE MoDELib)
add_test(NAME boundaryNodeUpdate COMMAND testBoundaryNodeUpdate "${PROJECT_SOURCE_DIR}")
set_property(TEST boundaryNodeUpdate APPEND PROPERTY LABELS "integration")
set_tests_properties(boundaryNodeUpdate PROPERTIES TIMEOUT 120)

add_executable(testNodeMotion testNodeMotion.cpp)
target_link_libraries(testNodeMotion PRIVATE MoDELib)
add_test(NAME nodeMotion COMMAND testNodeMotion "${PROJECT_SOURCE_DIR}")
set_property(TEST nodeMotion APPEND PROPERTY LABELS "integration")
set_tests_properties(nodeMotion PROPERTIES TIMEOUT 120)

add_executable(testJunctionSpline testJunctionSpline.cpp)
target_link_libraries(testJunctionSpline PRIVATE MoDELib)
add_test(NAME junctionSpline COMMAND testJunctionSpline "${PROJECT_SOURCE_DIR}")
set_property(TEST junctionSpline APPEND PROPERTY LABELS "numerical")
set_tests_properties(junctionSpline PROPERTIES TIMEOUT 120)

add_executable(testGlideVelocityDirection testGlideVelocityDirection.cpp)
target_link_libraries(testGlideVelocityDirection PRIVATE MoDELib)
add_test(NAME glideVelocityDirection COMMAND testGlideVelocityDirection "${PROJECT_SOURCE_DIR}")
set_property(TEST glideVelocityDirection APPEND PROPERTY LABELS "numerical")
set_tests_properties(glideVelocityDirection PROPERTIES TIMEOUT 120)

add_executable(testBicrystalPrismaticGlide testBicrystalPrismaticGlide.cpp)
target_link_libraries(testBicrystalPrismaticGlide PRIVATE MoDELib)
add_test(NAME bicrystalPrismaticGlide COMMAND testBicrystalPrismaticGlide "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalPrismaticGlide APPEND PROPERTY LABELS "integration")
set_tests_properties(bicrystalPrismaticGlide PROPERTIES TIMEOUT 120)

add_executable(testBicrystalBoundarySolver testBicrystalBoundarySolver.cpp)
target_link_libraries(testBicrystalBoundarySolver PRIVATE MoDELib)
add_test(NAME bicrystalBoundarySolver COMMAND testBicrystalBoundarySolver "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalBoundarySolver APPEND PROPERTY LABELS "integration")
set_tests_properties(bicrystalBoundarySolver PROPERTIES TIMEOUT 120)

add_executable(testBicrystalBoundaryTopology testBicrystalBoundaryTopology.cpp)
target_link_libraries(testBicrystalBoundaryTopology PRIVATE MoDELib)
add_test(NAME bicrystalBoundaryTopology COMMAND testBicrystalBoundaryTopology "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalBoundaryTopology APPEND PROPERTY LABELS "integration")
set_tests_properties(bicrystalBoundaryTopology PROPERTIES TIMEOUT 120)

if(TARGET testLoopAreaTriangulation)
    add_test(NAME loopAreaTriangulation COMMAND testLoopAreaTriangulation "${PROJECT_SOURCE_DIR}")
    set_property(TEST loopAreaTriangulation APPEND PROPERTY LABELS "numerical")
    set_tests_properties(loopAreaTriangulation PROPERTIES TIMEOUT 120)
endif()

if(TARGET testBicrystalNodeClassification)
    add_test(NAME bicrystalNodeClassification COMMAND testBicrystalNodeClassification "${PROJECT_SOURCE_DIR}")
    set_property(TEST bicrystalNodeClassification APPEND PROPERTY LABELS "integration")
    set_tests_properties(bicrystalNodeClassification PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 120)
endif()

add_executable(testBicrystalBoundaryRemesh testBicrystalBoundaryRemesh.cpp)
target_link_libraries(testBicrystalBoundaryRemesh PRIVATE MoDELib)
add_test(NAME bicrystalBoundaryRemesh COMMAND testBicrystalBoundaryRemesh "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalBoundaryRemesh APPEND PROPERTY LABELS "integration")
set_tests_properties(bicrystalBoundaryRemesh PROPERTIES TIMEOUT 120)

add_executable(testBicrystalOneSidedMotion testBicrystalOneSidedMotion.cpp)
target_link_libraries(testBicrystalOneSidedMotion PRIVATE MoDELib)
add_test(NAME bicrystalOneSidedMotion COMMAND testBicrystalOneSidedMotion "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalOneSidedMotion APPEND PROPERTY LABELS "integration")
set_tests_properties(bicrystalOneSidedMotion PROPERTIES TIMEOUT 120)

if(TARGET testMeshGrainBoundaries)
    add_test(NAME meshGrainBoundaries COMMAND testMeshGrainBoundaries "${PROJECT_SOURCE_DIR}")
    set_property(TEST meshGrainBoundaries APPEND PROPERTY LABELS "integration")
    set_tests_properties(meshGrainBoundaries PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 120)
endif()

if(TARGET testSurfaceInternalFaces)
    add_test(NAME surfaceInternalFaces COMMAND testSurfaceInternalFaces "${PROJECT_SOURCE_DIR}")
    set_property(TEST surfaceInternalFaces APPEND PROPERTY LABELS "integration")
    set_tests_properties(surfaceInternalFaces PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 120)
endif()

add_executable(testBicrystalJunctionCollapse testBicrystalJunctionCollapse.cpp)
target_link_libraries(testBicrystalJunctionCollapse PRIVATE MoDELib)
add_test(NAME bicrystalJunctionCollapse COMMAND testBicrystalJunctionCollapse "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalJunctionCollapse APPEND PROPERTY LABELS "integration")
add_test(NAME bicrystalJunction283 COMMAND testBicrystalJunctionCollapse "${PROJECT_SOURCE_DIR}" 283)
set_property(TEST bicrystalJunction283 APPEND PROPERTY LABELS "integration")
add_executable(testBicrystalCorner1537 testBicrystalCorner1537.cpp)
target_link_libraries(testBicrystalCorner1537 PRIVATE MoDELib)
add_test(NAME bicrystalCorner1537 COMMAND testBicrystalCorner1537 "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalCorner1537 APPEND PROPERTY LABELS "integration")
set_tests_properties(bicrystalJunctionCollapse PROPERTIES TIMEOUT 120)

add_executable(testBicrystalOpenBoundary testBicrystalOpenBoundary.cpp)
target_link_libraries(testBicrystalOpenBoundary PRIVATE MoDELib)
add_test(NAME bicrystalOpenBoundary COMMAND testBicrystalOpenBoundary "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalOpenBoundary APPEND PROPERTY LABELS "integration")
set_tests_properties(bicrystalOpenBoundary PROPERTIES TIMEOUT 120)

add_executable(testNonconvexMeshSearch testNonconvexMeshSearch.cpp)
target_link_libraries(testNonconvexMeshSearch PRIVATE MoDELib)
add_test(NAME nonconvexMeshSearch COMMAND testNonconvexMeshSearch "${PROJECT_SOURCE_DIR}")
set_property(TEST nonconvexMeshSearch APPEND PROPERTY LABELS "integration")
set_tests_properties(nonconvexMeshSearch PROPERTIES TIMEOUT 120)

add_executable(testPolycrystalEnclosedBoundary testPolycrystalEnclosedBoundary.cpp)
target_link_libraries(testPolycrystalEnclosedBoundary PRIVATE MoDELib)
add_test(NAME polycrystalEnclosedBoundary COMMAND testPolycrystalEnclosedBoundary "${PROJECT_SOURCE_DIR}")
set_property(TEST polycrystalEnclosedBoundary APPEND PROPERTY LABELS "integration")
set_tests_properties(polycrystalEnclosedBoundary PROPERTIES TIMEOUT 120)

add_executable(testPrismaticLoopAtomic testPrismaticLoopAtomic.cpp)
target_link_libraries(testPrismaticLoopAtomic PRIVATE MoDELib)
add_test(NAME prismaticLoopAtomic COMMAND testPrismaticLoopAtomic "${PROJECT_SOURCE_DIR}")
set_property(TEST prismaticLoopAtomic APPEND PROPERTY LABELS "integration")
set_tests_properties(prismaticLoopAtomic PROPERTIES TIMEOUT 120)

add_executable(testPolycrystalClosedBoundary testPolycrystalClosedBoundary.cpp)
target_link_libraries(testPolycrystalClosedBoundary PRIVATE MoDELib)
add_test(NAME polycrystalClosedBoundary COMMAND testPolycrystalClosedBoundary "${PROJECT_SOURCE_DIR}")
set_property(TEST polycrystalClosedBoundary APPEND PROPERTY LABELS "integration")
set_tests_properties(polycrystalClosedBoundary PROPERTIES TIMEOUT 120)

add_executable(testPolygonClipping testPolygonClipping.cpp)
target_link_libraries(testPolygonClipping PRIVATE MoDELib)
target_include_directories(testPolygonClipping PRIVATE "${PROJECT_SOURCE_DIR}/tools/DDqt/include")
add_test(NAME polygonClipping COMMAND testPolygonClipping "${PROJECT_SOURCE_DIR}")
set_property(TEST polygonClipping APPEND PROPERTY LABELS "numerical")
set_tests_properties(polygonClipping PROPERTIES TIMEOUT 120)

add_test(NAME partialPeriodicGlideSolver COMMAND testBicrystalBoundarySolver "${PROJECT_SOURCE_DIR}" partial)
set_property(TEST partialPeriodicGlideSolver APPEND PROPERTY LABELS "integration")
add_test(NAME nonPeriodicGlideSolver COMMAND testBicrystalBoundarySolver "${PROJECT_SOURCE_DIR}" physical)
set_property(TEST nonPeriodicGlideSolver APPEND PROPERTY LABELS "integration")
set_tests_properties(partialPeriodicGlideSolver nonPeriodicGlideSolver PROPERTIES TIMEOUT 120)

add_executable(testPartialPeriodicContraction testPartialPeriodicContraction.cpp)
target_link_libraries(testPartialPeriodicContraction PRIVATE MoDELib)
add_test(NAME partialPeriodicContraction COMMAND testPartialPeriodicContraction "${PROJECT_SOURCE_DIR}")
set_property(TEST partialPeriodicContraction APPEND PROPERTY LABELS "integration")
set_tests_properties(partialPeriodicContraction PROPERTIES TIMEOUT 120)

add_executable(testContractionSafety testContractionSafety.cpp)
target_link_libraries(testContractionSafety PRIVATE MoDELib)
add_test(NAME contractionSafety COMMAND testContractionSafety "${PROJECT_SOURCE_DIR}")
set_property(TEST contractionSafety APPEND PROPERTY LABELS "integration")
set_tests_properties(contractionSafety PROPERTIES TIMEOUT 120)

add_executable(testBoundaryJunctionFormation testBoundaryJunctionFormation.cpp)
target_link_libraries(testBoundaryJunctionFormation PRIVATE MoDELib)
add_test(NAME boundaryJunctionFormation COMMAND testBoundaryJunctionFormation "${PROJECT_SOURCE_DIR}")
set_property(TEST boundaryJunctionFormation APPEND PROPERTY LABELS "integration")
set_tests_properties(boundaryJunctionFormation PROPERTIES TIMEOUT 120)

add_executable(testJunctionFormationSafety testJunctionFormationSafety.cpp)
target_link_libraries(testJunctionFormationSafety PRIVATE MoDELib)
add_test(NAME junctionFormationSafety COMMAND testJunctionFormationSafety "${PROJECT_SOURCE_DIR}")
set_property(TEST junctionFormationSafety APPEND PROPERTY LABELS "integration")
set_tests_properties(junctionFormationSafety PROPERTIES TIMEOUT 120)

add_executable(testPartialPeriodicRemesh testPartialPeriodicRemesh.cpp)
target_link_libraries(testPartialPeriodicRemesh PRIVATE MoDELib)
add_test(NAME partialPeriodicRemesh COMMAND testPartialPeriodicRemesh "${PROJECT_SOURCE_DIR}")
set_property(TEST partialPeriodicRemesh APPEND PROPERTY LABELS "integration")
set_tests_properties(partialPeriodicRemesh PROPERTIES TIMEOUT 120)

add_executable(testRemeshSafety testRemeshSafety.cpp)
target_link_libraries(testRemeshSafety PRIVATE MoDELib)
add_test(NAME remeshSafety COMMAND testRemeshSafety "${PROJECT_SOURCE_DIR}")
set_property(TEST remeshSafety APPEND PROPERTY LABELS "integration")
set_tests_properties(remeshSafety PROPERTIES TIMEOUT 120)

add_executable(testJunctionLoopInsertion testJunctionLoopInsertion.cpp)
target_link_libraries(testJunctionLoopInsertion PRIVATE MoDELib)
add_test(NAME junctionLoopInsertion COMMAND testJunctionLoopInsertion "${PROJECT_SOURCE_DIR}")
set_property(TEST junctionLoopInsertion APPEND PROPERTY LABELS "integration")
set_tests_properties(junctionLoopInsertion PROPERTIES TIMEOUT 120)

add_executable(testCrossSlipTopology testCrossSlipTopology.cpp)
target_link_libraries(testCrossSlipTopology PRIVATE MoDELib)
add_test(NAME crossSlipTopology COMMAND testCrossSlipTopology "${PROJECT_SOURCE_DIR}")
set_property(TEST crossSlipTopology APPEND PROPERTY LABELS "integration")
set_tests_properties(crossSlipTopology PROPERTIES TIMEOUT 120)

add_executable(testThermalCrossSlipSampling testThermalCrossSlipSampling.cpp)
target_link_libraries(testThermalCrossSlipSampling PRIVATE MoDELib)
add_test(NAME thermalCrossSlipSampling COMMAND testThermalCrossSlipSampling)
set_property(TEST thermalCrossSlipSampling APPEND PROPERTY LABELS "numerical")
set_tests_properties(thermalCrossSlipSampling PROPERTIES TIMEOUT 120)

add_executable(testThermalCrossSlipTopology testThermalCrossSlipTopology.cpp)
target_link_libraries(testThermalCrossSlipTopology PRIVATE MoDELib)
add_test(NAME thermalCrossSlipTopology COMMAND testThermalCrossSlipTopology "${PROJECT_SOURCE_DIR}")
set_property(TEST thermalCrossSlipTopology APPEND PROPERTY LABELS "integration")
set_tests_properties(thermalCrossSlipTopology PROPERTIES TIMEOUT 120)

add_executable(testThermalCrossSlipBranches testThermalCrossSlipBranches.cpp)
target_link_libraries(testThermalCrossSlipBranches PRIVATE MoDELib)
add_test(NAME thermalCrossSlipBranches COMMAND testThermalCrossSlipBranches "${PROJECT_SOURCE_DIR}")
set_property(TEST thermalCrossSlipBranches APPEND PROPERTY LABELS "integration")
set_tests_properties(thermalCrossSlipBranches PROPERTIES TIMEOUT 120)

add_executable(testPeriodicGlissileJunction testPeriodicGlissileJunction.cpp)
target_link_libraries(testPeriodicGlissileJunction PRIVATE MoDELib)
if(MODELIB_ENABLE_ASSERTIONS)
    add_test(NAME periodicGlissileJunction COMMAND testPeriodicGlissileJunction "${PROJECT_SOURCE_DIR}" --require-assertions)
    set_property(TEST periodicGlissileJunction APPEND PROPERTY LABELS "integration")
else()
    add_test(NAME periodicGlissileJunction COMMAND testPeriodicGlissileJunction "${PROJECT_SOURCE_DIR}")
    set_property(TEST periodicGlissileJunction APPEND PROPERTY LABELS "integration")
endif()
set_tests_properties(periodicGlissileJunction PROPERTIES TIMEOUT 120)

add_executable(testBicrystalCorner1915 testBicrystalCorner1915.cpp)
target_link_libraries(testBicrystalCorner1915 PRIVATE MoDELib)
add_test(NAME bicrystalCorner1915 COMMAND testBicrystalCorner1915 "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalCorner1915 APPEND PROPERTY LABELS "integration")

add_executable(testPeriodicVelocityConstraints testPeriodicVelocityConstraints.cpp)
target_link_libraries(testPeriodicVelocityConstraints PRIVATE MoDELib)
add_test(NAME periodicVelocityConstraints COMMAND testPeriodicVelocityConstraints "${PROJECT_SOURCE_DIR}")
set_property(TEST periodicVelocityConstraints APPEND PROPERTY LABELS "numerical")

add_executable(testGlideIndexCaches testGlideIndexCaches.cpp)
target_link_libraries(testGlideIndexCaches PRIVATE MoDELib)
add_test(NAME glideIndexCaches COMMAND testGlideIndexCaches "${PROJECT_SOURCE_DIR}")
set_property(TEST glideIndexCaches APPEND PROPERTY LABELS "integration")
set_tests_properties(glideIndexCaches PROPERTIES TIMEOUT 120)

if(TARGET testActorIndexCaches)
    add_test(NAME actorIndexCaches COMMAND testActorIndexCaches "${PROJECT_SOURCE_DIR}")
    set_property(TEST actorIndexCaches APPEND PROPERTY LABELS "integration")
    set_tests_properties(actorIndexCaches PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 120)
endif()

add_executable(testCrossSlipOptimizations testCrossSlipOptimizations.cpp)
target_link_libraries(testCrossSlipOptimizations PRIVATE MoDELib)
add_test(NAME crossSlipOptimizations COMMAND testCrossSlipOptimizations "${PROJECT_SOURCE_DIR}")
set_property(TEST crossSlipOptimizations APPEND PROPERTY LABELS "integration")
set_tests_properties(crossSlipOptimizations PROPERTIES TIMEOUT 120)

add_executable(testBicrystalCollapse4244 testBicrystalCollapse4244.cpp)
target_link_libraries(testBicrystalCollapse4244 PRIVATE MoDELib)
add_test(NAME bicrystalCollapse4244 COMMAND testBicrystalCollapse4244 "${PROJECT_SOURCE_DIR}")
set_property(TEST bicrystalCollapse4244 APPEND PROPERTY LABELS "integration")
set_tests_properties(bicrystalCollapse4244 PROPERTIES TIMEOUT 120)

add_executable(testEnclosedJunctionOwnership testEnclosedJunctionOwnership.cpp)
target_link_libraries(testEnclosedJunctionOwnership PRIVATE MoDELib)
add_test(NAME enclosedJunctionOwnership COMMAND testEnclosedJunctionOwnership "${PROJECT_SOURCE_DIR}")
set_property(TEST enclosedJunctionOwnership APPEND PROPERTY LABELS "integration")

add_executable(testPeriodicFEMNodes testPeriodicFEMNodes.cpp)
target_link_libraries(testPeriodicFEMNodes PRIVATE MoDELib)
add_test(NAME periodicFEMNodes COMMAND testPeriodicFEMNodes "${PROJECT_SOURCE_DIR}")
set_property(TEST periodicFEMNodes APPEND PROPERTY LABELS "integration")
set_tests_properties(periodicFEMNodes PROPERTIES TIMEOUT 300)

add_executable(testFEMOwnership testFEMOwnership.cpp)
target_link_libraries(testFEMOwnership PRIVATE MoDELib)
add_test(NAME FEMOwnership COMMAND testFEMOwnership "${PROJECT_SOURCE_DIR}")
set_property(TEST FEMOwnership APPEND PROPERTY LABELS "integration")
set_tests_properties(FEMOwnership PROPERTIES TIMEOUT 120)

add_executable(testWeakProblem testWeakProblem.cpp)
target_link_libraries(testWeakProblem PRIVATE MoDELib)
add_test(NAME weakProblem COMMAND testWeakProblem "${PROJECT_SOURCE_DIR}")
set_property(TEST weakProblem APPEND PROPERTY LABELS "numerical")
set_tests_properties(weakProblem PROPERTIES TIMEOUT 120)

add_executable(testElasticGlobalPenalty testElasticGlobalPenalty.cpp)
target_link_libraries(testElasticGlobalPenalty PRIVATE MoDELib)
add_test(NAME elasticGlobalPenalty COMMAND testElasticGlobalPenalty "${PROJECT_SOURCE_DIR}")
set_property(TEST elasticGlobalPenalty APPEND PROPERTY LABELS "numerical")
set_tests_properties(elasticGlobalPenalty PROPERTIES TIMEOUT 120)

add_executable(testAffineDirichlet testAffineDirichlet.cpp)
target_link_libraries(testAffineDirichlet PRIVATE MoDELib)
add_test(NAME affineDirichlet COMMAND testAffineDirichlet "${PROJECT_SOURCE_DIR}")
set_property(TEST affineDirichlet APPEND PROPERTY LABELS "numerical")
set_tests_properties(affineDirichlet PROPERTIES TIMEOUT 120)

add_executable(testPlaneSectionMesh testPlaneSectionMesh.cpp)
target_link_libraries(testPlaneSectionMesh PRIVATE MoDELib)
add_test(NAME planeSectionMesh COMMAND testPlaneSectionMesh "${PROJECT_SOURCE_DIR}")
set_property(TEST planeSectionMesh APPEND PROPERTY LABELS "numerical")
set_tests_properties(planeSectionMesh PROPERTIES TIMEOUT 120)

if(MODELIB_CD_MSIZE STREQUAL "4" AND MODELIB_CD_ISIZE STREQUAL "8")
    add_executable(testClusterPeriodic testClusterPeriodic.cpp)
    target_link_libraries(testClusterPeriodic PRIVATE MoDELib)
    add_test(NAME clusterPeriodic COMMAND testClusterPeriodic "${PROJECT_SOURCE_DIR}")
    set_property(TEST clusterPeriodic APPEND PROPERTY LABELS "integration")
    set_tests_properties(clusterPeriodic PROPERTIES TIMEOUT 180)
endif()

add_subdirectory(spatialCDtest)

add_executable(testPhysicsBases testPhysicsBases.cpp)
target_link_libraries(testPhysicsBases PRIVATE MoDELib)
add_test(NAME physicsBases COMMAND testPhysicsBases "${PROJECT_SOURCE_DIR}")
set_property(TEST physicsBases APPEND PROPERTY LABELS "integration")
set_tests_properties(physicsBases PROPERTIES TIMEOUT 120)

if(TARGET pyMoDELib)
    add_test(NAME pythonMesh
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/testPythonMesh.py"
                "$<TARGET_FILE_DIR:pyMoDELib>")
    set_tests_properties(pythonMesh PROPERTIES LABELS "integration" TIMEOUT 120)
    add_test(NAME pythonPhysicsBases
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/testPhysicsBases.py"
                "$<TARGET_FILE_DIR:pyMoDELib>" "${PROJECT_SOURCE_DIR}")
    set_property(TEST pythonPhysicsBases APPEND PROPERTY LABELS "integration")
    set_tests_properties(pythonPhysicsBases PROPERTIES TIMEOUT 120)
endif()

# Separate translation units ensure each shared header is self-contained.
add_executable(testSharedHeaders
    sharedHeaders/main.cpp
    sharedHeaders/DefectiveCrystal.cpp
    sharedHeaders/DefectiveCrystalBase.cpp
    sharedHeaders/DefectiveCrystalParameters.cpp
    sharedHeaders/MicrostructureBase.cpp
    sharedHeaders/MicrostructureContainer.cpp
)
add_test(NAME sharedHeaders COMMAND testSharedHeaders)
set_property(TEST sharedHeaders APPEND PROPERTY LABELS "integration")

add_executable(testPhysicsHeaders
    physicsHeaders/main.cpp
    physicsHeaders/ClusterDynamicsBase.cpp
    physicsHeaders/ClusterDynamicsParameters.cpp
    physicsHeaders/ElasticDeformationBase.cpp
    physicsHeaders/ElasticDeformationFEM.cpp
    physicsHeaders/EllipsoidalInclusion.cpp
    physicsHeaders/EshelbyInclusionBase.cpp
    physicsHeaders/InclusionMicrostructureBase.cpp
    physicsHeaders/PolyhedronInclusion.cpp
    physicsHeaders/SphericalInclusion.cpp
    physicsHeaders/UniformController.cpp
    physicsHeaders/VoigtTraits.cpp
)
add_test(NAME physicsHeaders COMMAND testPhysicsHeaders)
set_property(TEST physicsHeaders APPEND PROPERTY LABELS "integration")

add_executable(testDDConcreteAPI testDDConcreteAPI.cpp)
target_link_libraries(testDDConcreteAPI PRIVATE MoDELib)
add_test(NAME DDConcreteAPI COMMAND testDDConcreteAPI)
set_property(TEST DDConcreteAPI APPEND PROPERTY LABELS "integration")

add_executable(testSourceSplits testSourceSplits.cpp)
target_link_libraries(testSourceSplits PRIVATE MoDELib)
add_test(NAME sourceSplits COMMAND testSourceSplits)
set_property(TEST sourceSplits APPEND PROPERTY LABELS "numerical;integration")

add_executable(testPolyhedralMesh testPolyhedralMesh.cpp)
target_link_libraries(testPolyhedralMesh PRIVATE MoDELib)
add_test(NAME polyhedralMesh COMMAND testPolyhedralMesh "${PROJECT_SOURCE_DIR}")
set_tests_properties(polyhedralMesh PROPERTIES LABELS "numerical;integration" TIMEOUT 120)

add_executable(testPolyhedralMeshDimensions testPolyhedralMeshDimensions.cpp)
target_link_libraries(testPolyhedralMeshDimensions PRIVATE MoDELib)
add_test(NAME polyhedralMeshDimensions COMMAND testPolyhedralMeshDimensions)
set_tests_properties(polyhedralMeshDimensions PROPERTIES LABELS "numerical;integration" TIMEOUT 120)

# Compile-only: the legacy mesh library currently links explicit 3D instantiations.
add_library(testPolyhedralMeshAdapters OBJECT polyhedralMeshAdapters.cpp)
target_link_libraries(testPolyhedralMeshAdapters PRIVATE MoDELib)

add_executable(testSimplexFactories testSimplexFactories.cpp)
target_link_libraries(testSimplexFactories PRIVATE MoDELib)
add_test(NAME simplexFactories COMMAND testSimplexFactories "${PROJECT_SOURCE_DIR}")
set_tests_properties(simplexFactories PROPERTIES LABELS "numerical;integration" TIMEOUT 120)

if(MODELIB_CD_MSIZE STREQUAL "4" AND MODELIB_CD_ISIZE STREQUAL "8")
    add_executable(testPeriodicClimbSolver testPeriodicClimbSolver.cpp)
    target_link_libraries(testPeriodicClimbSolver PRIVATE MoDELib)
    add_test(NAME periodicClimbSolver COMMAND testPeriodicClimbSolver "${PROJECT_SOURCE_DIR}")
    set_tests_properties(periodicClimbSolver PROPERTIES LABELS "integration" TIMEOUT 120)
endif()

add_executable(testBoxTree testBoxTree.cpp)
target_link_libraries(testBoxTree PRIVATE MoDELib)
add_test(NAME boxTree COMMAND testBoxTree)
set_tests_properties(boxTree PROPERTIES LABELS "numerical" TIMEOUT 60)

if(TARGET testBoxTreeActor)
    add_test(NAME boxTreeActor COMMAND testBoxTreeActor "${PROJECT_SOURCE_DIR}")
    set_tests_properties(boxTreeActor PROPERTIES LABELS "integration"
        ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 60)
endif()

add_executable(testPlaneBoxIntersection testPlaneBoxIntersection.cpp)
target_link_libraries(testPlaneBoxIntersection PRIVATE MoDELib)
add_test(NAME planeBoxIntersection COMMAND testPlaneBoxIntersection)
set_tests_properties(planeBoxIntersection PROPERTIES LABELS "numerical" TIMEOUT 60)
