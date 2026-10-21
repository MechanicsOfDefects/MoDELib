cmake_minimum_required(VERSION 3.16)

# Configure this directory directly for a focused build, or include it from
# the main build with BUILD_TESTS=ON.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(EllipsoidalInclusionTests LANGUAGES CXX)
    enable_testing()
endif()
find_package(Threads REQUIRED)
get_filename_component(MODELIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)

if(TARGET MoDELib)
    set(ELLIPSOIDAL_TEST_LIBRARY MoDELib)
else()
    # This target compiles repository sources, never stand-in base classes.
    find_package(Eigen3 QUIET)
    if(NOT TARGET Eigen3::Eigen)
        find_path(EIGEN3_INCLUDE_DIR Eigen/Dense PATH_SUFFIXES eigen3)
        if(NOT EIGEN3_INCLUDE_DIR)
            message(FATAL_ERROR "Set EIGEN3_INCLUDE_DIR to the directory containing Eigen/Dense")
        endif()
    endif()
    find_path(ELLIPSOIDAL_BOOST_INCLUDE_DIR boost/math/special_functions/bessel.hpp)
    if(NOT ELLIPSOIDAL_BOOST_INCLUDE_DIR)
        message(FATAL_ERROR "Set ELLIPSOIDAL_BOOST_INCLUDE_DIR to your Boost headers directory")
    endif()
    set(MODELIB_EXTRA_INCLUDE_DIRS "" CACHE STRING "Additional MoDELib dependency include directories")
    add_library(ellipsoidalInclusionTestCore STATIC
        "${MODELIB_ROOT}/src/InclusionMicrostructure/EshelbyInclusionBase.cpp"
        "${MODELIB_ROOT}/src/InclusionMicrostructure/SphericalInclusion.cpp"
        "${MODELIB_ROOT}/src/InclusionMicrostructure/EllipsoidalInclusion.cpp"
        "${MODELIB_ROOT}/src/InclusionMicrostructure/EllipsoidalInclusionIO.cpp"
        "${MODELIB_ROOT}/src/DislocationDynamics/DislocationNodeIO.cpp"
        "${MODELIB_ROOT}/src/DislocationDynamics/DislocationLoopIO.cpp"
        "${MODELIB_ROOT}/src/DislocationDynamics/DislocationLoopNodeIO.cpp"
        "${MODELIB_ROOT}/src/DislocationDynamics/DislocationLoopLinkIO.cpp"
        "${MODELIB_ROOT}/src/DislocationDynamics/DislocationQuadraturePointIO.cpp"
        "${MODELIB_ROOT}/src/DislocationDynamics/DislocationSegmentIO.cpp"
        "${MODELIB_ROOT}/src/DefectiveCrystal/DefectiveCrystalConfigIO.cpp"
        "${MODELIB_ROOT}/src/DefectiveCrystal/DefectiveCrystalBaseIO.cpp"
        "${MODELIB_ROOT}/src/MicrostructureGeneration/MicrostructureSpecificationBase.cpp"
        "${MODELIB_ROOT}/src/MicrostructureGeneration/EllipsoidalInclusionIndividualSpecification.cpp")
    target_sources(ellipsoidalInclusionTestCore PRIVATE
        "${MODELIB_ROOT}/src/InclusionMicrostructure/SphericalInclusionIO.cpp"
        "${MODELIB_ROOT}/src/InclusionMicrostructure/PolyhedronInclusionIO.cpp"
        "${MODELIB_ROOT}/src/InclusionMicrostructure/PolyhedronInclusionNodeIO.cpp"
        "${MODELIB_ROOT}/src/InclusionMicrostructure/PolyhedronInclusionEdgeIO.cpp"
        "${MODELIB_ROOT}/src/IO/StrUtilities.cpp")
    target_compile_features(ellipsoidalInclusionTestCore PUBLIC cxx_std_20)
    target_precompile_headers(ellipsoidalInclusionTestCore PRIVATE <cfloat>)
    target_include_directories(ellipsoidalInclusionTestCore PUBLIC
        "${MODELIB_ROOT}/include" "${ELLIPSOIDAL_BOOST_INCLUDE_DIR}" ${MODELIB_EXTRA_INCLUDE_DIRS})
    file(GLOB model_include_dirs LIST_DIRECTORIES true "${MODELIB_ROOT}/include/*")
    foreach(dir IN LISTS model_include_dirs)
        if(IS_DIRECTORY "${dir}" AND NOT dir MATCHES "_NEW$")
            target_include_directories(ellipsoidalInclusionTestCore PUBLIC "${dir}")
        endif()
    endforeach()
    if(TARGET Eigen3::Eigen)
        target_link_libraries(ellipsoidalInclusionTestCore PUBLIC Eigen3::Eigen)
    else()
        target_include_directories(ellipsoidalInclusionTestCore PUBLIC "${EIGEN3_INCLUDE_DIR}")
    endif()
    target_link_libraries(ellipsoidalInclusionTestCore PUBLIC Threads::Threads)
    set(ELLIPSOIDAL_TEST_LIBRARY ellipsoidalInclusionTestCore)
endif()

add_executable(testEllipsoidalInclusion testEllipsoidalInclusion.cpp)
add_executable(testEllipsoidalInclusionIO testEllipsoidalInclusionIO.cpp)
add_executable(testDDconfigEllipsoidalIO testDDconfigEllipsoidalIO.cpp)
add_executable(testEllipsoidalSpecification testEllipsoidalSpecification.cpp)
add_executable(testEllipsoidalPlacement testEllipsoidalPlacement.cpp)
add_executable(benchmarkEllipsoidalInclusion benchmarkStress.cpp)
add_executable(studyEllipsoidalTolerance toleranceStudy.cpp)
foreach(target IN ITEMS testEllipsoidalInclusion testEllipsoidalInclusionIO testDDconfigEllipsoidalIO testEllipsoidalSpecification testEllipsoidalPlacement benchmarkEllipsoidalInclusion studyEllipsoidalTolerance)
    target_compile_features(${target} PRIVATE cxx_std_20)
    target_precompile_headers(${target} PRIVATE <cfloat>)
    target_link_libraries(${target} PRIVATE ${ELLIPSOIDAL_TEST_LIBRARY} Threads::Threads)
endforeach()
target_compile_definitions(testEllipsoidalInclusion PRIVATE COMPARE_SUPPLIED_SPHERE)
add_test(NAME ellipsoidalInclusion COMMAND testEllipsoidalInclusion)
set_property(TEST ellipsoidalInclusion APPEND PROPERTY LABELS "numerical")
add_test(NAME ellipsoidalInclusionIO COMMAND testEllipsoidalInclusionIO)
set_property(TEST ellipsoidalInclusionIO APPEND PROPERTY LABELS "integration")
add_test(NAME DDconfigEllipsoidalIO COMMAND testDDconfigEllipsoidalIO)
set_property(TEST DDconfigEllipsoidalIO APPEND PROPERTY LABELS "integration")
add_test(NAME ellipsoidalSpecification COMMAND testEllipsoidalSpecification
    "${MODELIB_ROOT}/Library/Microstructures/ellipsoidalInclusionIndividual.txt")
set_property(TEST ellipsoidalSpecification APPEND PROPERTY LABELS "integration")
set_tests_properties(ellipsoidalInclusion PROPERTIES TIMEOUT 60)

add_test(NAME ellipsoidalPlacement COMMAND testEllipsoidalPlacement)
set_property(TEST ellipsoidalPlacement APPEND PROPERTY LABELS "numerical")

# This integration test uses the full library and a real mesh/material fixture.
if(TARGET MoDELib)
    add_executable(testEllipsoidalGenerator testEllipsoidalGenerator.cpp)
    target_compile_features(testEllipsoidalGenerator PRIVATE cxx_std_20)
    target_precompile_headers(testEllipsoidalGenerator PRIVATE <cfloat>)
    target_link_libraries(testEllipsoidalGenerator PRIVATE MoDELib Threads::Threads)
    add_test(NAME ellipsoidalGenerator COMMAND testEllipsoidalGenerator "${MODELIB_ROOT}")
    set_property(TEST ellipsoidalGenerator APPEND PROPERTY LABELS "integration")
    set_tests_properties(ellipsoidalGenerator PROPERTIES TIMEOUT 60)
endif()

# Register after the root has enabled testing (tools are configured earlier).
if(TARGET testEllipsoidalInclusionPolyData)
    add_test(NAME ellipsoidalInclusionPolyData COMMAND testEllipsoidalInclusionPolyData)
    set_property(TEST ellipsoidalInclusionPolyData APPEND PROPERTY LABELS "integration")
endif()

# Timings are observations, not portable pass/fail performance thresholds.
add_test(NAME ellipsoidalStressBenchmark COMMAND benchmarkEllipsoidalInclusion)
set_tests_properties(ellipsoidalStressBenchmark PROPERTIES
    LABELS "benchmark" TIMEOUT 120 RUN_SERIAL TRUE)
