cmake_minimum_required(VERSION 3.16)
project(Ex1 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")

# Optional OpenMP support on macOS
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(include)

# Source files (main will link them)
file(GLOB SOURCES "src/*.cpp")

add_executable(galaxy_grid ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(galaxy_grid PRIVATE OpenMP::OpenMP_CXX)
endif()