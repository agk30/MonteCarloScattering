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

# Set Fortran .mod Directory (where they should be outputted)
# set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/src)
# Include Source Code Directory
include_directories(src)
# Add certain libraries
add_library(mathConstants src/mathConstants.f90)
add_library(mod_tests src/testingMods.f90)

# Create Executable
add_executable(MCScattering src/MCScattering.f90)

# Link libraries, order is important so that relevant .mod files are created before compilation of specific modules
target_link_libraries (MCScattering mathConstants)
target_link_libraries (MCScattering mod_tests)

# Create image directories
add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory $"../Images"
    COMMAND ${CMAKE_COMMAND} -E make_directory $"../Images2"
    COMMAND ${CMAKE_COMMAND} -E make_directory $"../Images3")