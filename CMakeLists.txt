cmake_minimum_required( VERSION 2.8.3 )
project(MedRob)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options( -std=c++11 )

find_package( Eigen3 3.3 REQUIRED NO_MODULE )
find_package( yaml-cpp REQUIRED )
find_package(Boost COMPONENTS system filesystem REQUIRED)
set(CMAKE_BUILD_TYPE=Release)

###########
## Build ##
###########

include_directories(
  include
  ${CMAKE_SOURCE_DIR}/vrep/remoteApi
  ${CMAKE_SOURCE_DIR}/vrep/include
  ${EIGEN3_INCLUDE_DIR}
  ${YAML_CPP_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
  ${VREP_DIR}
)
link_directories(${Boost_LIBRARY_DIR})

set( SOURCES ${CMAKE_SOURCE_DIR}
              src/kukaDynRCM.cpp
              src/externalForcesEstimation.cpp
              src/others.cpp
    )

add_executable( main ${CMAKE_SOURCE_DIR}/src/main.cpp ${SOURCES} )

target_link_libraries( main Eigen3::Eigen yaml-cpp ${Boost_LIBRARIES})
