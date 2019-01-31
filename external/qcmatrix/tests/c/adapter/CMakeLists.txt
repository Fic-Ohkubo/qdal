# Required version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

# Sets the build mode: if the user specifies -DCMAKE_BUILD_TYPE on the command line,
# take their definition and dump it in the cache along with proper documentation,
# otherwise set CMAKE_BUILD_TYPE prior to calling PROJECT(), from
# http://www.cmake.org/pipermail/cmake/2008-September/023808.html
IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ELSE()
    SET(CMAKE_BUILD_TYPE None CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ENDIF()

# Project name, enable C, CXX and Fortran
PROJECT(QcMatrix C Fortran)

# Host program could include this file by setting LIB_QCMATRIX_PATH
IF("${LIB_QCMATRIX_PATH}" STREQUAL "")
    SET(LIB_QCMATRIX_PATH ${PROJECT_SOURCE_DIR})
ENDIF()

# Adds cmake directory to CMake module
SET(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Some options to control the properties of matrix library
OPTION(QCMATRIX_64BIT_INTEGER "Use 64 bit integer." OFF)
IF(QCMATRIX_64BIT_INTEGER)
    ADD_DEFINITIONS(-DQCMATRIX_64BIT_INTEGER)
ENDIF()
#
OPTION(QCMATRIX_BLAS_64BIT "Use 64 bit BLAS and LAPACK libraries." OFF)
IF(QCMATRIX_BLAS_64BIT)
    ADD_DEFINITIONS(-DQCMATRIX_BLAS_64BIT)
ENDIF()
#
OPTION(QCMATRIX_ZERO_BASED "Zero-based numbering." ON)
IF(QCMATRIX_ZERO_BASED)
    ADD_DEFINITIONS(-DQCMATRIX_ZERO_BASED)
ENDIF()
#
OPTION(QCMATRIX_ROW_MAJOR "Row major order for matrix elements" OFF)
IF(QCMATRIX_ROW_MAJOR)
    ADD_DEFINITIONS(-DQCMATRIX_ROW_MAJOR)
ENDIF()
#
OPTION(QCMATRIX_SINGLE_PRECISION "Use single precision for real numbers." OFF)
IF(QCMATRIX_SINGLE_PRECISION)
    ADD_DEFINITIONS(-DQCMATRIX_SINGLE_PRECISION)
ENDIF()
#
OPTION(QCMATRIX_STORAGE_MODE "Enable different matrix storage modes." OFF)
IF(QCMATRIX_STORAGE_MODE)
    ADD_DEFINITIONS(-DQCMATRIX_STORAGE_MODE)
ENDIF()
#
OPTION(QCMATRIX_ENABLE_VIEW "Enable matrix I/O." OFF)
IF(QCMATRIX_ENABLE_VIEW)
    ADD_DEFINITIONS(-DQCMATRIX_ENABLE_VIEW)
ENDIF()
#
OPTION(QCMATRIX_ENABLE_HDF5 "Enable the use of HDF5 library for matrix I/O." OFF)
IF(QCMATRIX_ENABLE_HDF5)
    ADD_DEFINITIONS(-DQCMATRIX_ENABLE_HDF5)
ENDIF()
#
#OPTION(QCMATRIX_ENABLE_MXML "Enable the use of Mini-XML library for matrix I/O." OFF)
#IF(QCMATRIX_ENABLE_MXML)
#    ADD_DEFINITIONS(-DQCMATRIX_ENABLE_MXML)
#ENDIF()
SET(QCMATRIX_ENABLE_MXML OFF)
#
# standard C functions for matrix I/O if any external I/O library disabled
IF(NOT QCMATRIX_ENABLE_HDF5 OR NOT QCMATRIX_ENABLE_MXML)
    ADD_DEFINITIONS(-DQCMATRIX_STANDARD_IO)
ENDIF()
#
OPTION(QCMATRIX_3M_METHOD "Enable 3M method for complex matrix-matrix multiplication." ON)
IF(QCMATRIX_3M_METHOD)
    ADD_DEFINITIONS(-DQCMATRIX_3M_METHOD)
ENDIF()
#
OPTION(QCMATRIX_STRASSEN_METHOD
       "Strassen's method for the square block complex matrix-matrix multiplication."
       ON)
IF(QCMATRIX_STRASSEN_METHOD)
    ADD_DEFINITIONS(-DQCMATRIX_STRASSEN_METHOD)
ENDIF()
#
OPTION(QCMATRIX_AUTO_ERROR_EXIT "Enable automatic exit on error." OFF)
IF(QCMATRIX_AUTO_ERROR_EXIT)
    ADD_DEFINITIONS(-DQCMATRIX_AUTO_ERROR_EXIT)
ENDIF()

# Sets the type of matrix library
OPTION(BLOCK_MATRIX "Square block matrix." OFF)
OPTION(COMPLEX_MATRIX "Complex matrix." OFF)

# Initializes the source codes
SET(QCMATRIX_SRCS)

# Sets the type of matrix implemented in the external library
SET(ADAPTER_BLOCK_CMPLX OFF)
SET(ADAPTER_BLOCK_REAL OFF)
SET(ADAPTER_CMPLX_MAT OFF)
SET(ADAPTER_REAL_MAT OFF)
IF(BLOCK_MATRIX)
    IF(COMPLEX_MATRIX)
        MESSAGE(FATAL_ERROR "Square block complex matrix not supported!")
    ELSE()
        SET(ADAPTER_MATRIX_TYPE "ADAPTER_BLOCK_REAL")
        SET(ADAPTER_BLOCK_REAL ON)
    ENDIF()
ELSE()
    IF(COMPLEX_MATRIX)
        SET(ADAPTER_MATRIX_TYPE "ADAPTER_CMPLX_MAT")
        SET(ADAPTER_CMPLX_MAT ON)
    ELSE()
        SET(ADAPTER_MATRIX_TYPE "ADAPTER_REAL_MAT")
        SET(ADAPTER_REAL_MAT ON)
    ENDIF()
ENDIF()
ADD_DEFINITIONS(-D${ADAPTER_MATRIX_TYPE})

# QcMatrix uses HDF5 library to save the structures of matrices
IF(QCMATRIX_ENABLE_VIEW AND QCMATRIX_ENABLE_HDF5)
    # In case that the host program has found the HDF5 library
    IF(NOT DEFINED LIB_HDF5_RELEASE)
        # On UNIX systems, the variable HDF5_USE_STATIC_LIBRARIES is used to determine whether
        # or not to prefer a static link to a dynamic link for HDF5 and all of it's dependencies
        IF(UNIX AND NOT DEFINED HDF5_USE_STATIC_LIBRARIES)
            OPTION(HDF5_USE_STATIC_LIBRARIES "Enable the static link for HDF5." ON)
        ENDIF()
        # To provide the module with a hint about where to find the HDF5 installation, one can
        # set the environment variable HDF5_ROOT
        IF(NOT DEFINED HDF5_ROOT)
            SET(HDF5_ROOT None CACHE STRING "Provide a hind about where to find the HDF5 installation.")
        ENDIF()
        # Tries to find the HDF5 library
        INCLUDE(FindHDF5)
        IF(HDF5_FOUND)
            # Location of the HDF5 includes
            INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
            # Required compiler definitions for HDF5
            ADD_DEFINITIONS(${HDF5_DEFINITIONS})
            # Required libraries for all requested bindings are defined in ${HDF5_LIBRARIES_RELEASE}
            FOREACH(_LIB ${HDF5_LIBRARIES_RELEASE})
                STRING(FIND ${_LIB} "libpthread.a" _LIB_PTHREAD)
                IF(${_LIB_PTHREAD} EQUAL -1)
                    SET(LIB_HDF5_RELEASE ${LIB_HDF5_RELEASE} ${_LIB})
                ELSE()
                    # To prevent, for instance undefined reference to `__syscall_error'
                    SET(LIB_HDF5_RELEASE ${LIB_HDF5_RELEASE} "-lpthread")
                ENDIF()
            ENDFOREACH()
            #ADD_LIBRARY(LIB_HDF5_RELEASE UNKNOWN IMPORTED)
            #SET_PROPERTY(TARGET LIB_HDF5_RELEASE PROPERTY IMPORTED_LOCATION "${HDF5_LIBRARIES_RELEASE}")
        ELSE()
            MESSAGE(FATAL_ERROR "HDF5 library is needed for QcMatrix!")
        ENDIF()
    ENDIF()
    # The HDF5-1.8 APIs are used, so it is necessary to define H5_NO_DEPRECATED_SYMBOLS on some systems
    ADD_DEFINITIONS(-DH5_NO_DEPRECATED_SYMBOLS)
    MESSAGE("HDF5 library: ${LIB_HDF5_RELEASE}")
ENDIF()

# Header files of QcMatrix
INCLUDE_DIRECTORIES(${LIB_QCMATRIX_PATH}/include)

# Adds the project binary directory to include directories
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# Source codes
IF(ADAPTER_BLOCK_REAL)
    INCLUDE(RealMat)
    INCLUDE(BlockReal)
ELSEIF(ADAPTER_CMPLX_MAT)
    INCLUDE(RealMat)
    INCLUDE(CmplxMat)
ELSEIF(ADAPTER_REAL_MAT)
    INCLUDE(RealMat)
ELSE()
    MESSAGE(FATAL_ERROR "Unknown matrix type of the adapter!")
ENDIF()

# Generates a configure file
INCLUDE(QcMatrixConfig.h)
#CONFIGURE_FILE(${LIB_QCMATRIX_PATH}/cmake/QcMatrixConfig.h.cmake
#               ${CMAKE_BINARY_DIR}/qcmatrix_config.h) 

# Creates the QcMatrix library
SET(QCMATRIX_SRCS
    ${QCMATRIX_SRCS}
    src/error_exit.F90
    src/lapack/qcmatrix_blas.F90
    src/lapack/qcmatrix_c_blas.F90)
ADD_LIBRARY(adapter_matrix ${QCMATRIX_SRCS})
IF(QCMATRIX_ENABLE_HDF5)
    TARGET_LINK_LIBRARIES(adapter_matrix
                          ${LIB_HDF5_RELEASE})
ENDIF()

#
INSTALL(TARGETS adapter_matrix ARCHIVE DESTINATION lib)
