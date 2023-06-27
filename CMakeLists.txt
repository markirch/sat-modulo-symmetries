cmake_minimum_required(VERSION 3.12)

project(SMS VERSION "1.0" LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if (CMAKE_BUILD_TYPE)
  message("Build type set to ${CMAKE_BUILD_TYPE}.")
else()
  set(CMAKE_BUILD_TYPE Release)
endif ()


# find_package(Z3 4.12.0 REQUIRED)

if (DEFINED GLASGOW)
	message("Linking with the Glasgow Subgraph Solver")
	set(CMAKE_CXX_STANDARD 20)
	add_compile_definitions(GLASGOW)
	include_directories(glasgow-subgraph-solver/src)
	link_libraries("${PROJECT_SOURCE_DIR}/glasgow-subgraph-solver/libcommon.a")
	link_libraries("pthread")
	link_libraries("stdc++fs")
	link_libraries("boost_container")
	link_libraries("boost_thread")
	link_libraries("boost_iostreams")
else ()
	message("Skipping Glasgow")
endif()

add_subdirectory(src)
# add_subdirectory(src_directed)
# add_subdirectory(src_hypher_graph)
# add_subdirectory(src_hypher_graph_incidence)

# -DZ3_INCLUDE_DIR=[path][/src/api
# -DZ3_LIBRARIES=[path]/debug/libz3.so

# cmake . -B./build -DZ3_INCLUDE_DIR=./z3-master/src/api -DZ3_LIBRARIES=./z3-master/build/libz3.so
