cmake_minimum_required(VERSION 3.0.0)

# set the project name and version
project(simulation VERSION 0.0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall")

add_subdirectory(components)

add_executable(simulation main.cpp)

# add the library to the program
target_link_libraries(simulation PUBLIC components)

# add the binary tree to the search path for include files
target_include_directories(simulation PUBLIC "${PROJECT_SOURCE_DIR}/components")
