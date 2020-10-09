###
# Compilation/Installation Instructions
#
# - mkdir build
# - cmake ./ -B ./build/ && cd build
# - make
# ./MCScattering
###

# Minimum CMAKE Version to use
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
# Project Name and Type
project(MonteCarloScattering Fortran)
enable_language(Fortran)

# Set Fortran .mod Directory (where they should be outputted)
# set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/Modules)
# Include Modules Directory
include_directories(Modules)
# Add certain libraries
add_library(mathConstants Modules/mathConstants.f90)
add_library(mod_tests Modules/testingMods.f90)

# Create Executable
add_executable(MCScattering MCScattering.f90)

# Link libraries, order is important so that relevant .mod files are created before compilation of specific modules
target_link_libraries (MCScattering mathConstants)
target_link_libraries (MCScattering mod_tests)