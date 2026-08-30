cmake_minimum_required(VERSION 3.12)

project(SMS VERSION "1.0" LANGUAGES C CXX)

# set(CMAKE_EXE_LINKER_FLAGS "-static")
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
	set(CMAKE_CXX_STANDARD 20)
	add_compile_definitions(GLASGOW)
	include_directories(glasgow-subgraph-solver)
	#link_libraries("${PROJECT_SOURCE_DIR}/glasgow-subgraph-solver/libcommon.a")
	link_libraries("${PROJECT_SOURCE_DIR}/glasgow-subgraph-solver/build/gss/libglasgow_subgraphs.a")
	link_libraries("pthread")
	if(NOT APPLE)
		link_libraries("stdc++fs")
	endif()
	# Glasgow uses GMP for arbitrary-precision integers
	find_path(GMP_INCLUDE_DIR NAMES gmp.h)
	find_library(GMP_LIBRARY NAMES gmp)
	find_library(GMPXX_LIBRARY NAMES gmpxx)
	if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
		message(FATAL_ERROR
			"GMP and GMPXX are required for GLASGOW builds. "
			"Install GMP or set CMAKE_PREFIX_PATH to its installation prefix.")
	endif()
	include_directories("${GMP_INCLUDE_DIR}")
	link_libraries("${GMPXX_LIBRARY}" "${GMP_LIBRARY}")
	find_package(Boost REQUIRED COMPONENTS container thread iostreams)
	include_directories(${Boost_INCLUDE_DIRS})
	link_libraries(${Boost_LIBRARIES})
	# link_libraries("boost_container")
	# link_libraries("boost_thread")
	# link_libraries("boost_iostreams")
	message("-- Linking with the Glasgow Subgraph Solver")
	unset(GLASGOW) # to suppress the warning about an unused variable in src/CMakeLists.txt
else ()
	message("-- Building without the Glasgow Subgraph Solver")
endif()

add_subdirectory(src)
# add_subdirectory(src_directed)
# add_subdirectory(src_hypher_graph)
# add_subdirectory(src_hypher_graph_incidence)

# -DZ3_INCLUDE_DIR=[path][/src/api
# -DZ3_LIBRARIES=[path]/debug/libz3.so

# cmake . -B./build -DZ3_INCLUDE_DIR=./z3-master/src/api -DZ3_LIBRARIES=./z3-master/build/libz3.so
