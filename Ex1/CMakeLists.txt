cmake_minimum_required(VERSION 3.16)
project(CompAstro_Exercise01 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")

# Point CMake to Homebrew's libomp on macOS
if(APPLE)
    if(EXISTS "/opt/homebrew/opt/libomp")
        set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
    elseif(EXISTS "/usr/local/opt/libomp")
        set(OpenMP_ROOT "/usr/local/opt/libomp")
    endif()
endif()

find_package(OpenMP)

include_directories(include)
file(GLOB SOURCES "src/*.cpp")
add_executable(galaxy_grid ${SOURCES})

if(OpenMP_CXX_FOUND)
    # OpenMP::OpenMP_CXX automatically attaches compile flags, includes, and linker flags
    target_link_libraries(galaxy_grid PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found! Compiling with multithreading support.")
else()
    message(STATUS "OpenMP NOT found. Compiling single-threaded.")
endif()