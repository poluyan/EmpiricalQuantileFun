cmake_minimum_required(VERSION 2.8)

## Use the variable PROJECT_NAME for changing the target name
set( PROJECT_NAME "tests" )

## Set our project name
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)

if (MSVC)
    add_compile_options(/Wall /std:c++17)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wnon-virtual-dtor -Woverloaded-virtual -Wno-unused-parameter -Wconversion)
endif()

add_executable(test1d test1d.cpp)
add_executable(test2d test2d.cpp)

set(foo_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/mveqf)

set_property(TARGET test1d test2d PROPERTY INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR})


