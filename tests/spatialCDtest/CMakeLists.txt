# Original surface-bounded 120-step tutorial regression.
if(MODELIB_CD_MSIZE STREQUAL "4" AND MODELIB_CD_ISIZE STREQUAL "8")
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_executable(testSpatialCD spatialCDtest.cpp)
    target_link_libraries(testSpatialCD PRIVATE MoDELib)
    add_test(NAME spatialCDtest COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/run.py"
        --executable "$<TARGET_FILE:testSpatialCD>"
        --reports "${CMAKE_CURRENT_SOURCE_DIR}/reports")
    set_tests_properties(spatialCDtest PROPERTIES
        LABELS "integration;slow;ClusterDynamics" TIMEOUT 1800 PROCESSORS 12 RUN_SERIAL TRUE)
endif()
