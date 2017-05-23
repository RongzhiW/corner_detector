project(corner_detector)
cmake_minimum_required(VERSION 2.8)

find_package(catkin REQUIRED COMPONENTS
  image_transport
  message_generation
  roscpp
  rospy
  sensor_msgs
  std_msgs
  cv_bridge
)
set(OpenCV_DIR /usr/local/share/OpenCV)
find_package(OpenCV REQUIRED)

include_directories(
./
${EIGEN3_INCLUDE_DIR}
${GeographicLib_INCLUDE_DIRS}
)

#aux_source_directory(. SRC_LIST)

#add_executable(${PROJECT_NAME} ${SRC_LIST})

add_executable(corner_detector
                src/main.cpp
                src/detector.cpp
                src/detector.h

)

target_link_libraries(corner_detector
                ${OpenCV_LIBRARIES}
                ${catkin_LIBRARIES}

                )

#message(STATUS "OpenCV_LIBRARIES: " ${OpenCV_LIBRARIES})
set_target_properties( corner_detector
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_SOURCE_DIR}/bin
)



