cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(gmshdemo)
add_executable(top top.cpp)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
	    IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so"
	        INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
	)
	set(GMSH_LIBRARIES gmsh)

target_link_libraries(top ${GMSH_LIBRARIES})
