
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Ofast -march=native -fopenmp -w")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS "-O3 -xHOST -openmp")
endif()

add_executable(
    QICgen
    ${MAIN_CMAKE_SOURCE_DIR}/ICgen/IICgen.cpp
    ${MAIN_CMAKE_SOURCE_DIR}/ICgen/IICgen_RRs.cpp
    ${MAIN_CMAKE_SOURCE_DIR}/ICgen/IICwrite.cpp
    ${MAIN_CMAKE_SOURCE_DIR}/ICgen/SKS.cpp
    ${MAIN_CMAKE_SOURCE_DIR}/ICgen/integrals.cpp
    ${MAIN_CMAKE_SOURCE_DIR}/math/angular.cpp
    ${MAIN_CMAKE_SOURCE_DIR}/ICgen/QICgen.cpp
    )

add_custom_target(
    generate_qic
    COMMAND mkdir -p ${MAIN_CMAKE_BINARY_DIR}/qic\;
            mkdir -p ${MAIN_CMAKE_BINARY_DIR}/k2c\;
            ./QICgen ${QIC_MAX_L_VALUE} ${MAIN_CMAKE_BINARY_DIR}/qic ${MAIN_CMAKE_BINARY_DIR}/k2c
    DEPENDS QICgen
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
