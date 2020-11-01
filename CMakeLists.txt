cmake_minimum_required(VERSION 2.8)

set( PROJECT_NAME "tests" )

project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)

if (MSVC)
    add_compile_options(/Wall /std:c++17)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wnon-virtual-dtor -Woverloaded-virtual -Wno-unused-parameter)
endif()

add_executable(test1d_u demos/test1d_uniform.cpp)
add_executable(test1d_n demos/test1d_nonuniform.cpp)
add_executable(test2d_u demos/test2d_uniform.cpp)
add_executable(test2d_n demos/test2d_nonuniform.cpp)
add_executable(test3d_u demos/test3d_uniform.cpp)
add_executable(test3d_n demos/test3d_nonuniform.cpp)
add_executable(testNd_u demos/testNd_uniform.cpp)
add_executable(testNd_n demos/testNd_nonuniform.cpp)

# for using angle brackets for headers
set_property(TARGET test1d_u test1d_n test2d_u test2d_n test3d_u test3d_n testNd_u testNd_n PROPERTY INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR})

# moving executables to bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set_target_properties(test1d_u test1d_n test2d_u test2d_n test3d_u test3d_n testNd_u testNd_n PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)


