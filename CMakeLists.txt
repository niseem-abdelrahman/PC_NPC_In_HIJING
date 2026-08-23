cmake_minimum_required(VERSION 3.0)
project(HijingFortranCXX)

# Compiler options
set(CMAKE_CXX_STANDARD 17)  # Set C++ standard to 14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_OSX_ARCHITECTURES "x86_64")

# Find Fortran compiler
enable_language(Fortran)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# Add executable
add_executable(HIJrun  hipyset1.36.f hijing1.420.f wrapper.cpp wrapper.h GeneralQCorrelator.h WrapperQCorrEvent.h)


# Set Fortran standard
set_target_properties(HIJrun PROPERTIES LINKER_LANGUAGE Fortran)

# Specify linker flags
target_link_libraries(HIJrun  PRIVATE gfortran ${ROOT_LIBRARIES})

# Set include directories
target_include_directories(HIJrun PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} PRIVATE ${ROOT_INCLUDE_DIRS})

# Workaround for ROOT with C++14
target_compile_definitions(HIJrun PRIVATE ROOT_FORCE_STL=1)
