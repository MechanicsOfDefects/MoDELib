target_sources(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/InclusionMicrostructure.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/EllipsoidalInclusion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/EshelbyInclusionBase.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InclusionMicrostructureBase.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SphericalInclusion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/EllipsoidalInclusionIO.cpp
)

target_sources(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/SphericalInclusionIO.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusionIO.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusionNodeIO.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusionEdgeIO.cpp
)
