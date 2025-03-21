set(Boost_USE_STATIC_LIBS OFF)

# set(Boost_USE_MULTITHREADED ON)
# set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.74.0 REQUIRED COMPONENTS program_options graph)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

set(sources "minimalityCheckCommon.cpp" "sms.cpp" "graph.cpp" "useful.cpp" "options.cpp" "minimalityCheck.cpp" "minimalityCheck_dir.cpp" "cubing.cpp" "qbf/universal.cpp" "qbf/universal2.cpp" "qbf/qcirParser.cpp" )
if(GLASGOW)
	list(APPEND sources "other/forbiddenSubgraph.cpp")
endif()

# append graph propagators
file(GLOB prop_resources "graphPropagators/*.cpp")
message("Graph propagators: ${prop_resources}")

list(APPEND sources ${prop_resources})


set(glasgow_wrapper "glasgow_wrapper.cpp")

find_library(CADICAL_LIBRARY
NAMES cadical
HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../cadical_sms/build/"
)
if (NOT CADICAL_LIBRARY)
	message(FATAL_ERROR "CaDiCaL not found (")
else()
	message("CaDiCaL found at ${CADICAL_LIBRARY}")
endif()
add_library(cadical STATIC IMPORTED GLOBAL)
set_target_properties(cadical PROPERTIES IMPORTED_LOCATION ${CADICAL_LIBRARY})
link_libraries(cadical)

# cross-platform way to elevate the warning level
# -O3 not needed, because it is implied by CMAKE_BUILD_TYPE=Release (which is the default)
if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
    # additional warnings
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

if (DEFINED CLINGO)
	find_package(Clingo)
	if (Clingo_FOUND) # only build with Clingo if the library was explicitly requested and is available
		list(APPEND sources "${clingo_wrapper}")
		add_compile_definitions(INCLUDE_CLINGO)
		link_libraries(libclingo)
		message("Linking with Clingo")
	else()
		message("Could not find Clingo, continuing without it")
	endif()
endif()


#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "time")

add_library(sms SHARED ${sources})
add_library(sms_static STATIC ${sources})

set_target_properties(sms_static PROPERTIES OUTPUT_NAME sms)

# targets must be added *after* all calls to link_libraries()
add_executable(smsg "main.cpp")
target_link_libraries(smsg PRIVATE sms_static)

install(TARGETS
	smsg
	sms
	sms_static
)
install(FILES ${CADICAL_LIBRARY} DESTINATION lib RENAME libsmscadical.a)
install(FILES "sms.hpp" "cadical.hpp" "useful.h" "graphChecker.hpp" "minimalityCheck.hpp" DESTINATION include/sms)


# commands for adding Z3
# target_include_directories(graphgen PRIVATE ${Z3_CXX_INCLUDE_DIRS})
# target_link_libraries(graphgen PRIVATE ${Z3_LIBRARIES})
# target_compile_options(graphgen PRIVATE ${Z3_COMPONENT_CXX_FLAGS})

# target_include_directories(graphgendirected PRIVATE ${Z3_CXX_INCLUDE_DIRS})
# target_link_libraries(graphgendirected PRIVATE ${Z3_LIBRARIES})
# target_compile_options(graphgendirected PRIVATE ${Z3_COMPONENT_CXX_FLAGS})

# add_executable(graphgenBoth "main.cpp" "minimalityCheck.cpp"  "solveCadical.cpp" "useful.cpp" "arguments.cpp" "path.cpp" "connectedComponents.cpp" "coloring.cpp" "solveGeneral.cpp")
# target_link_libraries(graphgenBoth PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../cadical-internal-theory-propagation/build/libcadical.a")
# target_compile_options(graphgenBoth PRIVATE "-Wall" "-Wextra" "-pedantic" "-O3")
