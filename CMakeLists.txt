cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(jasp)
add_definitions(-D_USE_MATH_DEFINES)

AUX_SOURCE_DIRECTORY("src" DIR_SRCS)
add_executable(${PROJECT_NAME} ${DIR_SRCS})

find_package(MKL)
if(MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
endif()
