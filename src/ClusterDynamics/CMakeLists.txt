target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/ClusterDynamics.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ClusterDynamicsFEM.cpp
			)


target_sources(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ClusterDynamicsBase.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ClusterDynamicsParameters.cpp
)
