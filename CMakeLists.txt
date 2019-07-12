cmake_minimum_required(VERSION 2.8.3)
project(MedRob)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++11)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package( yaml-cpp REQUIRED)

set(CMAKE_BUILD_TYPE=Release)


###########
## Build ##
###########

include_directories (${VREP_DIR})

include_directories(
  include
  ${CMAKE_SOURCE_DIR}/remoteApi
  ${EIGEN3_INCLUDE_DIR}
  ${YAML_CPP_INCLUDE_DIR}
)

set( SOURCES ${CMAKE_SOURCE_DIR}
              src/kukaDynRCM.cpp
              src/externalForcesEstimation.cpp
              src/others.cpp
#              src/readMeasurements.cpp
#              src/RoboteqCommand.cpp
#              src/RoboteqParser.cpp
#              src/joyMotCmd.cpp
#              src/g2o_ba.cpp
#              src/visualOdom.cpp
    )

add_executable( main ${CMAKE_SOURCE_DIR}/src/main.cpp ${SOURCES} )

target_link_libraries( main Eigen3::Eigen yaml-cpp )
