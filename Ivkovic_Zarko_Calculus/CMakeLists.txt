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
add_executable(part11 "src/kinds.f95" "src/integrals.f95" "src/part11.f95")
add_executable(part12 "src/kinds.f95" "src/integrals.f95" "src/part12.f95")
add_executable(part13 "src/kinds.f95" "src/integrals.f95" "src/part13.f95")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/Part2")
add_executable(part21 "src/kinds.f95" "src/optimization.f95" "src/part21.f95")
add_executable(part22 "src/kinds.f95" "src/optimization.f95" "src/part22.f95")
