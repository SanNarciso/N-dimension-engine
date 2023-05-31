cmake_minimum_required(VERSION 3.15)
project(Engine)

set(CMAKE_CXX_STANDARD 20)
add_executable(${CMAKE_PROJECT_NAME}_run main.cpp lib/Engine/RayCast.h)

# set this flag when running coverage tests in Clion
#set(CMAKE_CXX_FLAGS "--coverage")

# generate debug information
set(CMAKE_CXX_FLAGS "-g")

include_directories(Engine)
add_subdirectory(Engine)

add_subdirectory(test/Math_tests)

target_link_libraries(${CMAKE_PROJECT_NAME}_run ${CMAKE_PROJECT_NAME}_lib)