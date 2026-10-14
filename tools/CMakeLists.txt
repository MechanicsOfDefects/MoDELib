
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGenerator MicrostructureGenerator)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/DDomp DDomp)
if(BUILD_DDQT)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/DDqt DDqt)
endif()
if(pybind11_FOUND)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pyMoDELib pyMoDELib)
endif()

