###
# Compilation/Installation Instructions
#
# - mkdir build
# - cmake ./ -B ./build/ -G "Unix Makefiles"
# - cd build
# - make OR mingw32-make
# ./MCScattering
###

# Minimum CMAKE Version to use
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
# Project Name and Type
project(MonteCarloScattering Fortran)
enable_language(Fortran)
#add_compile_options(-fopenmp)
set (CMAKE_Fortran_FLAGS_RELEASE "-ffree-line-length-none")

# Set Fortran .mod Directory (where they should be outputted)
# set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/src)
# Include Source Code Directory
include_directories(src)
# Add certain libraries
add_library(mathConstants src/hwlib/mathConstants.f90)
add_library(speeds src/hwlib/speeds.f90)
add_library(directions src/hwlib/directions.f90)
add_library(mod_tests src/hwlib/testingMods.f90)
add_library(m_config src/hwlib/m_config.f90)

# Create Executable
add_executable(MCScattering src/MCScattering.f90)

# Link libraries, order is important so that relevant .mod files are created before compilation of specific modules
target_link_libraries (MCScattering mathConstants)
target_link_libraries (mod_tests mathConstants)
target_link_libraries (mod_tests directions)
target_link_libraries (speeds mathConstants)
target_link_libraries (directions speeds)
target_link_libraries (MCScattering speeds)
target_link_libraries (MCScattering directions)
target_link_libraries (MCScattering mod_tests)
target_link_libraries (MCScattering m_config)