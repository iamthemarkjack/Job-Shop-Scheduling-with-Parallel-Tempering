cmake_minimum_required(VERSION 3.10)
project(TSP_Parallel_Tempering CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMP package
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add optimization flags for Release build
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Set source files
set(SOURCES
    src/main.cpp
    src/parallel_tempering.cpp
    src/tsp_problem.cpp
    src/solution.cpp
)

# Add executable
add_executable(tsp_solver ${SOURCES})

# Include directories
target_include_directories(tsp_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Link OpenMP
target_link_libraries(tsp_solver PRIVATE OpenMP::OpenMP_CXX)

# Create a data directory in the build folder
add_custom_command(
    TARGET tsp_solver POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
    ${CMAKE_BINARY_DIR}/data
)

# Install
install(TARGETS tsp_solver DESTINATION bin)