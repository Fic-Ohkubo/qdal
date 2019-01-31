

# exit if cmake version is below 2.8
# we need at least 2.8 for external projects support
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(QIC_MAX_L_VALUE "4" CACHE STRING "Max L value for QICgen")

# this will generate sources
# we call here another cmake instance
exec_program(
    ${CMAKE_COMMAND}
    ${CMAKE_BINARY_DIR}/interpreted
    ARGS -DMAIN_CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}
         -DMAIN_CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
         -DQIC_MAX_L_VALUE=${QIC_MAX_L_VALUE}
         ${CMAKE_SOURCE_DIR}/cmake/interpreted/)
exec_program(
    make
    ${CMAKE_BINARY_DIR}/interpreted
    ARGS generate_qic
    )

# set project name and languanges that are involved
project(K4MIRROR CXX)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# these are paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_BINARY_DIR}
    )

include(ConfigParentSettings)
include(ConfigSafeGuards)

include_directories(${CMAKE_SOURCE_DIR})

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-march=native -fopenmp")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -pg -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-w -Ofast")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS "-openmp")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -pg -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHOST")
endif()

add_library(
    newcontractions
    ${CMAKE_SOURCE_DIR}/integrals/newcontractions.cpp
    )

add_library(
    quimera
    ${CMAKE_SOURCE_DIR}/low/MPIwrap.cpp
    ${CMAKE_SOURCE_DIR}/low/chrono.cpp
    ${CMAKE_SOURCE_DIR}/low/plotters.cpp
    ${CMAKE_SOURCE_DIR}/math/gamma.cpp
    ${CMAKE_SOURCE_DIR}/math/angular.cpp
    ${CMAKE_SOURCE_DIR}/libquimera/libquimera.cpp
    ${CMAKE_SOURCE_DIR}/2eints/IIC64K4.cpp
    ${CMAKE_SOURCE_DIR}/2eints/IIC64MIRROR.cpp
    ${CMAKE_SOURCE_DIR}/2eints/IICinit.cpp
    ${CMAKE_SOURCE_DIR}/2eints/IIC32.cpp
    ${CMAKE_SOURCE_DIR}/2eints/quimera.cpp
    ${CMAKE_SOURCE_DIR}/2eints/VILIC.cpp
    ) 

add_library(
    echidna
    ${CMAKE_SOURCE_DIR}/basis/SPprototype.cpp
    ${CMAKE_SOURCE_DIR}/basis/shellpair.cpp
    ${CMAKE_SOURCE_DIR}/integrals/atomprod.cpp
    ${CMAKE_SOURCE_DIR}/integrals/ERIbatch.cpp
    ${CMAKE_SOURCE_DIR}/integrals/ERIblocks.cpp
    ${CMAKE_SOURCE_DIR}/integrals/fock2e.cpp
    ${CMAKE_SOURCE_DIR}/integrals/rotations.cpp
    ${CMAKE_SOURCE_DIR}/libechidna/libechidna.cpp
    ${CMAKE_SOURCE_DIR}/low/cache64.cpp
    ${CMAKE_SOURCE_DIR}/low/cache32.cpp
    ${CMAKE_SOURCE_DIR}/math/eigen.cpp
    ${CMAKE_SOURCE_DIR}/math/newsparse.cpp
    ${CMAKE_SOURCE_DIR}/math/tensors.cpp
    )

file(GLOB K2C_SOURCES "${CMAKE_BINARY_DIR}/k2c/*.cpp")
add_library(
    k2c_library
    ${K2C_SOURCES}
    )

install(TARGETS newcontractions quimera echidna ARCHIVE DESTINATION lib)
