# The legacy MPI helpers are optional; normal shared-memory builds need no MPI.
option(MODELIB_ENABLE_MPI "Build MPI support" OFF)
if(MODELIB_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ModelMPIbase.cpp)
    target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
    target_compile_definitions(${PROJECT_NAME} PUBLIC _MODEL_MPI_)
    target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include/MPI)
endif()
