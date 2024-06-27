set(Boost_USE_STATIC_LIBS ON)

# set(Boost_USE_MULTITHREADED ON)
# set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.74.0 REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

set(sources "minimalityCheckCommon.cpp" "sms.cpp" "cadicalSMS.cpp" "useful.cpp" "connectedComponents.cpp" "coloring.cpp" "connectedChecker.cpp" "planarity.cpp" "coloringCheck.cpp"  "universal.cpp" "universal2.cpp" "qcirParser.cpp" "subgraphIsomorphism" "efx.cpp" "domination.cpp")
if(GLASGOW)
	list(APPEND sources "forbiddenSubgraph.cpp")
endif()

set(clingo_wrapper "clingoSMS.cpp")
set(glasgow_wrapper "glasgow_wrapper.cpp")
# set(cadical_static_lib "${CMAKE_CURRENT_SOURCE_DIR}/../cadical-rel-1.7.4/build/libcadical.a")
set(cadical_static_lib "${CMAKE_CURRENT_SOURCE_DIR}/../cadical/build/libcadical.a")
link_libraries("${cadical_static_lib}")

# -O3 not needed, because it is implied by CMAKE_BUILD_TYPE=Release (which is the default)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -Wextra -pedantic")

find_package(Clingo)
if(Clingo_FOUND) # only build with clingo if the library is available
	list(APPEND sources "${clingo_wrapper}")
	add_compile_definitions(INCLUDE_CLINGO)
	link_libraries(libclingo)
endif()

#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "time")

add_library(sms SHARED "minimalityCheck.cpp" ${sources})
add_library(smsdir SHARED "minimalityCheck_dir.cpp" ${sources})
add_library(sms_static STATIC "minimalityCheck.cpp" ${sources})
add_library(smsdir_static STATIC "minimalityCheck_dir.cpp" ${sources})

# targets must be added *after* all calls to link_libraries()
add_executable(smsg "main.cpp")
add_executable(smsd "main.cpp")

target_compile_definitions(smsd PRIVATE DIRECTED)
target_compile_definitions(smsdir PRIVATE DIRECTED)
target_compile_definitions(smsdir_static PRIVATE DIRECTED)

target_link_libraries(smsg PRIVATE sms_static)
target_link_libraries(smsd PRIVATE smsdir_static)

install(TARGETS smsg smsd sms smsdir sms_static smsdir_static)
install(FILES "${cadical_static_lib}" DESTINATION lib)
install(FILES "sms.hpp" "cadical.hpp" "cadicalSMS.hpp" "useful.h" "graphChecker.hpp" "minimalityCheck.hpp" DESTINATION include/sms)
if (Clinfo_FOUND)
	install(FILES "clingoSMS.hpp" DESTINATION include/sms)
endif()


# commands for adding Z3
# target_include_directories(graphgen PRIVATE ${Z3_CXX_INCLUDE_DIRS})
# target_link_libraries(graphgen PRIVATE ${Z3_LIBRARIES})
# target_compile_options(graphgen PRIVATE ${Z3_COMPONENT_CXX_FLAGS})

# target_include_directories(graphgendirected PRIVATE ${Z3_CXX_INCLUDE_DIRS})
# target_link_libraries(graphgendirected PRIVATE ${Z3_LIBRARIES})
# target_compile_options(graphgendirected PRIVATE ${Z3_COMPONENT_CXX_FLAGS})

# add_executable(graphgenBoth "main.cpp" "minimalityCheck.cpp"  "solveCadical.cpp" "useful.cpp" "arguments.cpp" "path.cpp" "connectedComponents.cpp" "coloring.cpp" "solveGeneral.cpp")
# target_link_libraries(graphgenBoth PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../cadical-internal-theory-propagation/build/libcadical.a")
# target_link_libraries(graphgenBoth PRIVATE libclingo)
# target_compile_options(graphgenBoth PRIVATE "-Wall" "-Wextra" "-pedantic" "-O3")
