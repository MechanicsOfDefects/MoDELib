
# MoDELib build file.
#
# Typical configuration:
#   cmake -S . -B build
#   cmake --build build -j
#
# Options (pass as -D<name>=<value>):
#   USE_PYBIND11        ON/OFF  build pyMoDELib and the Python solver hooks (needs Python 3 and pybind11)
#   BUILD_TOOLS         ON/OFF  build DDomp, microstructureGenerator, DDqt and pyMoDELib
#   BUILD_DDQT          ON/OFF  build the DDqt viewer (needs Qt 6 and VTK)
#   MODELIB_MARCH       value of -march, "native" by default; use e.g. x86-64-v2 for portable binaries
#   MODELIB_CD_MSIZE    number of mobile cluster-dynamics species (4 by default)
#   MODELIB_CD_ISIZE    number of immobile cluster-dynamics fields (8 by default, 0 for none)
#   EIGEN3_INCLUDE_DIRS path to the Eigen 3 headers, when they are not found automatically
#   CMAKE_PREFIX_PATH   where to find Qt 6, e.g. ~/Qt/6.7.3/macos/lib/cmake
#   Python_EXECUTABLE   Python interpreter to use for pybind11
#
# Example on macOS with MacPorts (clang 18, Python 3.11):
#   export LDFLAGS=-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#   export LDFLAGS="$LDFLAGS -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/system"
#   cmake .. -DCMAKE_CXX_COMPILER=clang++-mp-18 -DCMAKE_PREFIX_PATH=$HOME/Qt/6.7.3/macos/lib/cmake \
#            -DPython_EXECUTABLE=/opt/local/bin/python3.11 \
#            -Dpybind11_DIR=$(/opt/local/bin/python3.11 -m pybind11 --cmakedir)
# See also https://mcilloni.ovh/2021/02/09/cxx-cross-clang/

cmake_minimum_required(VERSION 3.20)
project(MoDELib VERSION 2.2.0 LANGUAGES C CXX)

option(USE_PYBIND11 "enables pybind11" ON)
option(BUILD_TOOLS "builds tools" ON)
option(BUILD_DDQT "builds the DDqt viewer (requires Qt 6 and VTK)" ON)
set(MODELIB_MARCH "native" CACHE STRING "value passed to -march")
set(MODELIB_CD_MSIZE "4" CACHE STRING "number of mobile cluster-dynamics species, e.g. Cv Ci C2i C3i")
set(MODELIB_CD_ISIZE "8" CACHE STRING "number of immobile cluster-dynamics fields, e.g. Nc Na1 Na2 Na3 cc ca1 ca2 ca3; 0 disables immobile species")
add_definitions(-DMODELIB_CD_MSIZE=${MODELIB_CD_MSIZE} -DMODELIB_CD_ISIZE=${MODELIB_CD_ISIZE})

### Eigen ###
# Order: EIGEN3_INCLUDE_DIRS given by the user, then the Eigen3 CMake package, then common install locations.
if(NOT EIGEN3_INCLUDE_DIRS)
	find_package(Eigen3 QUIET NO_MODULE)
	if(Eigen3_FOUND)
		get_target_property(EIGEN3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
	else()
		find_path(EIGEN3_INCLUDE_DIRS NAMES Eigen/Dense signature_of_eigen3_matrix_library
			PATHS /opt/local/include/eigen3 /opt/homebrew/include/eigen3 /usr/local/include/eigen3 /usr/include/eigen3)
	endif()
endif()
if(NOT EIGEN3_INCLUDE_DIRS)
	message(FATAL_ERROR "Eigen 3 not found. Pass -DEIGEN3_INCLUDE_DIRS=/path/to/eigen3")
endif()
message(STATUS "EIGEN3_INCLUDE_DIRS:" ${EIGEN3_INCLUDE_DIRS})

### Run-time path of the MacPorts gcc libraries, when present ###
foreach(RPATH /opt/local/lib/gcc14 /opt/local/lib/libgcc)
	if(EXISTS ${RPATH})
		list(APPEND CMAKE_BUILD_RPATH ${RPATH})
	endif()
endforeach()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=${MODELIB_MARCH} -fopenmp -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-parameter") # -Wno-class-memaccess
# set(CMAKE_C_FLAGS "-Wall -Wextra -march=${MODELIB_MARCH} -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-parameter")

#set(CMAKE_AUTOMOC ON)
#set(CMAKE_AUTORCC ON)
#set(CMAKE_AUTOUIC ON)

### Boost ### 
find_package(Boost)
if(Boost_FOUND)
	message(STATUS "Boost_VERSION:" ${Boost_VERSION})
	include_directories(${Boost_INCLUDE_DIRS})
else()
	message(WARNING "Boost not found")
endif()

### SuiteSparse ### 
find_package(CHOLMOD)
if(CHOLMOD_FOUND)
	message(STATUS "CHOLMOD_VERSION:" ${CHOLMOD_VERSION})
	message(STATUS "CHOLMOD_INCLUDE_DIR:" ${CHOLMOD_INCLUDE_DIR})
	message(STATUS "CHOLMOD_LIBRARIES:" ${CHOLMOD_LIBRARIES})
	include_directories(${CHOLMOD_INCLUDE_DIR})
	add_definitions(-D MODELIB_CHOLMOD)
else()
	message(WARNING "CHOLMOD not found")
endif()
find_package(UMFPACK)
if(UMFPACK_FOUND)
	message(STATUS "UMFPACK_VERSION:" ${UMFPACK_VERSION})
	message(STATUS "UMFPACK_INCLUDE_DIR:" ${UMFPACK_INCLUDE_DIR})
	message(STATUS "UMFPACK_LIBRARIES:" ${UMFPACK_LIBRARIES})
	include_directories(${UMFPACK_INCLUDE_DIR})
	add_definitions(-D MODELIB_UMFPACK)
else()
	message(WARNING "UMFPACK not found")
endif()

### SUNDIALS (optional) ###
# CVODE integrates the rate equations of the immobile cluster-dynamics species.
# Without it they are integrated by the explicit scheme, which has no dissolution, emission and coalescence.
option(USE_SUNDIALS "uses SUNDIALS (CVODE) for the immobile cluster-dynamics species, when found" ON)
if(USE_SUNDIALS)
	find_package(SUNDIALS 7.0 QUIET) # the interface of version 7 is used
	if(SUNDIALS_FOUND)
		message(STATUS "SUNDIALS_VERSION:" ${SUNDIALS_VERSION})
		add_definitions(-D_MODEL_SUNDIALS_)
	else()
		message(WARNING "SUNDIALS 7 not found: the CVODE integrator of cluster dynamics is disabled. Install it (apt-get install libsundials-dev, brew install sundials) or pass -DCMAKE_PREFIX_PATH=<SUNDIALS install>")
	endif()
endif()

### fftw3 ### find_package is currently broken. Workaround using find_library
find_library(FFTW3_LIBRARIES fftw3)  
if(FFTW3_LIBRARIES)
	cmake_path(GET FFTW3_LIBRARIES PARENT_PATH FFTW3_LIBRARY_DIRS)
	set(FFTW3_INCLUDE_DIRS ${FFTW3_LIBRARY_DIRS}/../include)
	message(STATUS "FFTW3_LIBRARIES:" ${FFTW3_LIBRARIES})
	message(STATUS "FFTW3_LIBRARY_DIRS:" ${FFTW3_LIBRARY_DIRS})
	message(STATUS "FFTW3_INCLUDE_DIRS:" ${FFTW3_INCLUDE_DIRS})
	include_directories(${BOOST_INCLUDE_DIRS})
	include_directories(${FFTW3_INCLUDE_DIRS})
	add_definitions(-D _MODEL_GLIDE_PLANE_NOISE_GENERATOR_)
else()
  message(WARNING "library ${FFTW3_LIBRARIES} not found, noise generator is disabled")
endif()

### pybind11 ###
if(USE_PYBIND11)
	cmake_policy(SET CMP0057 NEW)
	find_package(Python 3 COMPONENTS Interpreter Development)
	if(Python_FOUND)
		message(STATUS "Python_VERSION:" ${Python_VERSION})
		find_package(pybind11 CONFIG)
		if(NOT pybind11_FOUND)
			message(WARNING "pybind11 not found. Pass -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir), or -DUSE_PYBIND11=OFF")
		endif()
	else()
		message(WARNING "Python not found, pybind11 is disabled")
	endif()
endif()

if(pybind11_FOUND)
	message(STATUS "pybind11_VERSION:" ${pybind11_VERSION})
	set(pybind11_LIBRARIES pybind11::embed pybind11::headers pybind11::module)
	message(STATUS "pybind11_LIBRARIES: ${pybind11_LIBRARIES}")
	add_definitions(-D_MODEL_PYBIND11_)
	Python_add_library(MoDELib SHARED "")
else()
	message(WARNING "pybind11 is disabled")
	add_library(MoDELib "")
endif()

set(MoDELib_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories ( ${EIGEN3_INCLUDE_DIRS} )
include_directories(${MoDELib_INCLUDE_DIR}/Utilities/)
include_directories(${MoDELib_INCLUDE_DIR}/IO/)
include_directories(${MoDELib_INCLUDE_DIR}/Math/)
include_directories(${MoDELib_INCLUDE_DIR}/Geometry/)
include_directories(${MoDELib_INCLUDE_DIR}/Mesh/)
include_directories(${MoDELib_INCLUDE_DIR}/Lattices/)
include_directories(${MoDELib_INCLUDE_DIR}/FEM/)
include_directories(${MoDELib_INCLUDE_DIR}/GlidePlanes/)
include_directories(${MoDELib_INCLUDE_DIR}/Quadrature/)
include_directories(${MoDELib_INCLUDE_DIR}/PolycrystallineMaterials/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationDynamicsIO/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationMobilities/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationDynamicsBase/)
include_directories(${MoDELib_INCLUDE_DIR}/MicrostructureGeneration/)
include_directories(${MoDELib_INCLUDE_DIR}/LoopNetwork/)
include_directories(${MoDELib_INCLUDE_DIR}/ParticleInteraction/)
include_directories(${MoDELib_INCLUDE_DIR}/SpaceDecomposition/)
include_directories(${MoDELib_INCLUDE_DIR}/DiscreteCrackMechanics/)
include_directories(${MoDELib_INCLUDE_DIR}/ClusterDynamics/)
include_directories(${MoDELib_INCLUDE_DIR}/ElasticDeformation/)
include_directories(${MoDELib_INCLUDE_DIR}/InclusionMicrostructure/)
include_directories(${MoDELib_INCLUDE_DIR}/DislocationDynamics/)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/ src)
target_link_libraries(MoDELib PUBLIC ${FFTW3_LIBRARIES} ${pybind11_LIBRARIES} ${CHOLMOD_LIBRARIES} ${UMFPACK_LIBRARIES})
if(USE_SUNDIALS AND SUNDIALS_FOUND)
	target_link_libraries(MoDELib PUBLIC SUNDIALS::cvode SUNDIALS::nvecserial SUNDIALS::sunmatrixdense SUNDIALS::sunlinsoldense)
endif()

if(BUILD_TOOLS)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tools tools)
else()
endif()
