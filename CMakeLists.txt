cmake_minimum_required(VERSION 2.8)

project(jasp-1-0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")


AUX_SOURCE_DIRECTORY("src" DIR_SRCS)
add_executable(${PROJECT_NAME} ${DIR_SRCS})

