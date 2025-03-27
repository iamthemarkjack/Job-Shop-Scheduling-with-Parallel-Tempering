cmake_minimum_required(VERSION 3.10)
project(JSSP_Parallel_Tempering CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost REQUIRED COMPONENTS graph)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(SOURCES
    src/main.cpp
    src/parallel_tempering.cpp
    src/jssp_problem.cpp
    src/solution.cpp
)

add_executable(jssp_solver ${SOURCES})

target_include_directories(jssp_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(jssp_solver PRIVATE Boost::graph)

target_link_libraries(tsp_solver PRIVATE OpenMP::OpenMP_CXX)

install(TARGETS jssp_solver DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})