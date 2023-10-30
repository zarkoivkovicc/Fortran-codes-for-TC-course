cmake_minimum_required(VERSION 3.5)

project(homework1)
enable_language(Fortran)
# Compiler selection
include(CMakeForceCompiler) 
file(GLOB_RECURSE sources src/linear_algebra.f90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/Part1")
add_executable(part1 "src/kinds.f95" "src/linear_algebra.f95" "src/inout.f95" "src/part1.f95")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/Part2")
add_executable(part2 "src/kinds.f95" "src/math.f95" "src/inout.f95" "src/part2.f95")
