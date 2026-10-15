

add_subdirectory(ellipsoidalInclusion)

if(TARGET testSurfaceMesh)
    add_test(NAME surfaceMesh COMMAND testSurfaceMesh "${PROJECT_SOURCE_DIR}")
    set_tests_properties(surfaceMesh PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 60)
endif()

if(TARGET testWrappingTabWidget)
    add_test(NAME wrappingTabWidget COMMAND testWrappingTabWidget)
    set_tests_properties(wrappingTabWidget PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 30)
endif()

if(TARGET testDDqtFrameWorkflow)
    add_test(NAME DDqtFrameWorkflow COMMAND testDDqtFrameWorkflow "${PROJECT_SOURCE_DIR}")
    set_tests_properties(DDqtFrameWorkflow PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen" TIMEOUT 60)
endif()

add_executable(testPeriodicMesh testPeriodicMesh.cpp)
target_link_libraries(testPeriodicMesh PRIVATE MoDELib)
add_test(NAME periodicMesh COMMAND testPeriodicMesh "${PROJECT_SOURCE_DIR}")
set_tests_properties(periodicMesh PROPERTIES TIMEOUT 120)
