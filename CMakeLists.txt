cmake_minimum_required(VERSION 3.4...3.18)
project(GenomeRearrangement)

add_subdirectory(pybind11)
ADD_LIBRARY(LibsModule 
    GenomeRearrangementSRC/errorMsg.cpp
    GenomeRearrangementSRC/FastZip.cpp
    GenomeRearrangementSRC/genome.cpp
    GenomeRearrangementSRC/GenomeRearrangements.cpp
    GenomeRearrangementSRC/genomes.cpp
    GenomeRearrangementSRC/logFile.cpp
    GenomeRearrangementSRC/RandomGenerators.cpp
    GenomeRearrangementSRC/readTree.cpp
    GenomeRearrangementSRC/simulator.cpp
    GenomeRearrangementSRC/tree.cpp
)

pybind11_add_module(GenomeRearrangement src/main.cpp)
target_link_libraries(GenomeRearrangement PRIVATE LibsModule)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(GenomeRearrangement
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall --std=c++11 -O3 -fPIC" )
set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -O3 -fPIC" )
