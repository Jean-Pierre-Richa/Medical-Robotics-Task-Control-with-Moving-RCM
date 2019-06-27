cmake_minimum_required(VERSION 2.8.3)
project(MedRob)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++11)

find_package( Eigen3 3.3 REQUIRED NO_MODULE )

###########
## Build ##
###########

include_directories(
  include
)

set( SOURCES ${CMAKE_SOURCE_DIR}
              src/kukaDynRCM.cpp
#              src/streamImages.cpp
#              src/readMeasurements.cpp
#              src/RoboteqCommand.cpp
#              src/RoboteqParser.cpp
#              src/joyMotCmd.cpp
#              src/g2o_ba.cpp
#              src/visualOdom.cpp
    )

add_executable( main ${CMAKE_SOURCE_DIR}/src/main.cpp ${SOURCES} )

target_link_libraries( main Eigen3::Eigen )
