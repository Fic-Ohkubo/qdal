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
PROJECT(QcMatrix C CXX Fortran)

# Host program could include this file by setting LIB_QCMATRIX_PATH
IF("${LIB_QCMATRIX_PATH}" STREQUAL "")
    SET(LIB_QCMATRIX_PATH ${PROJECT_SOURCE_DIR})
ENDIF()

# Adds cmake directory to CMake module
SET(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ConfigParentSettings.cmake configures directories where to place module files
# and where to find header files, this is important if there are module and header
# file dependencies between libraries and where the parent code needs to enforce
# a directory structure to compile correctly
INCLUDE(ConfigParentSettings)

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
#
SET(LIB_QCMATRIX_NAME "qcmatrix" CACHE STRING "Sets the name of the QcMatrix library.")

# Builds the test suite as excutables, otherwise the test suite is included into the QcMatrix library
OPTION(QCMATRIX_TEST_EXECUTABLE "Build the test suite excutables." ON)

OPTION(QCMATRIX_TEST_3M_METHOD "Build the test of efficiency of 3M method." OFF)

# A header file will be created for Fortran name mangling if the Fortran 90
# adapter and/or APIs are built
SET(QCMATRIX_FC_MANGLING OFF)

# Sets up external or internal matrix libraries
OPTION(QCMATRIX_BUILD_ADAPTER "Build the adapter for external matrix library." OFF)
IF(QCMATRIX_BUILD_ADAPTER)
    OPTION(EXTERNAL_BLOCK_MATRIX "Square block matrix implemented in the external library." OFF)
    OPTION(EXTERNAL_COMPLEX_MATRIX "Complex matrix implemented in the external library." OFF)
    SET(QCMATRIX_EXTERNAL_PATH None CACHE STRING "Sets the path of external libraries.")
    SET(QCMATRIX_EXTERNAL_LIBRARIES None CACHE STRING "Sets the external libraries (like -lxxxx).")
    IF(NOT ${QCMATRIX_EXTERNAL_PATH} STREQUAL "None")
        INCLUDE_DIRECTORIES(${QCMATRIX_EXTERNAL_PATH})
        LINK_DIRECTORIES(${QCMATRIX_EXTERNAL_PATH})
    ENDIF()
    # Makes a drop down list for the external library, inspired by
    # http://blog.bethcodes.com/cmake-tips-tricks-drop-down-list
    # and
    # https://clauandr.sdf-eu.org/blog/2012/04/06/cmake-select-among-exclusive-multiple-options
    SET(QCMATRIX_ADAPTER_LIST "C" "CXX" "F90" "F03")
    SET(QCMATRIX_ADAPTER_TYPE "C" CACHE STRING
        "Choose the type of the external library, valid entries are ${QCMATRIX_ADAPTER_LIST}.")
    SET_PROPERTY(CACHE QCMATRIX_ADAPTER_TYPE PROPERTY STRINGS ${QCMATRIX_ADAPTER_LIST})
    LIST(FIND QCMATRIX_ADAPTER_LIST ${QCMATRIX_ADAPTER_TYPE} QCMATRIX_ADAPTER_CHOSEN)
    IF(${QCMATRIX_ADAPTER_CHOSEN} EQUAL -1)
        MESSAGE(FATAL_ERROR
                "External '${QCMATRIX_ADAPTER_TYPE}' library not supported, valid entries are ${QCMATRIX_ADAPTER_LIST}.")
    ENDIF()
ELSE()
    # UNSET's are needed to clear and hide variables from the user
    UNSET(EXTERNAL_BLOCK_MATRIX)
    UNSET(EXTERNAL_COMPLEX_MATRIX)
    UNSET(QCMATRIX_ADAPTER_TYPE)
ENDIF()

# Let the user set different application programming interfaces (API)
OPTION(QCMATRIX_CXX_API "Build C++ API." OFF)
SET(QCMATRIX_Fortran_API None CACHE STRING "Build Fortran API, options are: None;F90;F03.")

# Initializes the source codes
SET(QCMATRIX_SRCS)

# Source code of error handling
UNSET(QCMATRIX_ERROR_SRC)

# Sets the adapter for the external library
IF(QCMATRIX_BUILD_ADAPTER)
    # Sets the type of matrix implemented in the external library
    SET(ADAPTER_BLOCK_CMPLX OFF)
    SET(ADAPTER_BLOCK_REAL OFF)
    SET(ADAPTER_CMPLX_MAT OFF)
    SET(ADAPTER_REAL_MAT OFF)
    SET(ADAPTER_C_LANG OFF)
    SET(ADAPTER_F90_LANG OFF)
    SET(ADAPTER_F03_LANG OFF)
    IF(EXTERNAL_BLOCK_MATRIX)
        IF(EXTERNAL_COMPLEX_MATRIX)
            SET(ADAPTER_MATRIX_TYPE "ADAPTER_BLOCK_CMPLX")
            SET(ADAPTER_BLOCK_CMPLX ON)
        ELSE()
            SET(ADAPTER_MATRIX_TYPE "ADAPTER_BLOCK_REAL")
            SET(ADAPTER_BLOCK_REAL ON)
        ENDIF()
    ELSE()
        IF(EXTERNAL_COMPLEX_MATRIX)
            SET(ADAPTER_MATRIX_TYPE "ADAPTER_CMPLX_MAT")
            SET(ADAPTER_CMPLX_MAT ON)
        ELSE()
            SET(ADAPTER_MATRIX_TYPE "ADAPTER_REAL_MAT")
            SET(ADAPTER_REAL_MAT ON)
        ENDIF()
    ENDIF()
    ADD_DEFINITIONS(-D${ADAPTER_MATRIX_TYPE})
    # C
    IF(${QCMATRIX_ADAPTER_CHOSEN} EQUAL 0)
        IF(ADAPTER_BLOCK_CMPLX)
            MESSAGE(FATAL_ERROR
                    "QcMatrix does not support external C square block complex matrix.")
        ENDIF()
        SET(ADAPTER_C_LANG ON)
        ADD_DEFINITIONS(-DADAPTER_C_LANG)
        INCLUDE(c_adapter)
        # Source codes for setting/getting external C matrix
        SET(QCMATRIX_SRCS
            ${QCMATRIX_SRCS}
            ${LIB_QCMATRIX_PATH}/src/qcmat/c/QcMatSetExternalMat.c
            ${LIB_QCMATRIX_PATH}/src/qcmat/c/QcMatGetExternalMat.c)
    # C++
    ELSEIF(${QCMATRIX_ADAPTER_CHOSEN} EQUAL 1)
      
    # Fortran 90
    ELSEIF(${QCMATRIX_ADAPTER_CHOSEN} EQUAL 2)
        SET(ADAPTER_F90_LANG ON)
        ADD_DEFINITIONS(-DADAPTER_F90_LANG)
        INCLUDE(f90_adapter)
        SET(QCMATRIX_FC_MANGLING ON)
        # APIs should be the same as the adapter
        SET(QCMATRIX_Fortran_API "F90")
    # Fortran 2003
    ELSEIF(${QCMATRIX_ADAPTER_CHOSEN} EQUAL 3)
        SET(ADAPTER_F03_LANG ON)
        ADD_DEFINITIONS(-DADAPTER_F03_LANG)
        INCLUDE(f03_adapter)
        # APIs should be the same as the adapter
        SET(QCMATRIX_Fortran_API "F03")
        SET(QCMATRIX_ERROR_SRC ${LIB_QCMATRIX_PATH}/src/qcmat/qcmatrix_error.F90)
    ENDIF()
ENDIF()

# Sets different APIs
IF(QCMATRIX_CXX_API)

ENDIF()
IF(${QCMATRIX_Fortran_API} STREQUAL "F90")
    INCLUDE(f90_api)
    SET(QCMATRIX_FC_MANGLING ON)
    SET(QCMATRIX_ERROR_SRC ${LIB_QCMATRIX_PATH}/src/qcmat/qcmatrix_error.F90)
ELSEIF(${QCMATRIX_Fortran_API} STREQUAL "F03")
    INCLUDE(f03_api)
    SET(QCMATRIX_ERROR_SRC ${LIB_QCMATRIX_PATH}/src/qcmat/qcmatrix_error.F90)
ENDIF()

# Adds source code of error handling
IF(DEFINED QCMATRIX_ERROR_SRC)
    SET(QCMATRIX_SRCS ${QCMATRIX_SRCS} ${QCMATRIX_ERROR_SRC})
ENDIF()

# Creates a header file f90_mangling.h with correct mangling for Fortran subroutines
IF(QCMATRIX_FC_MANGLING)
    INCLUDE(FortranCInterface)
    FortranCInterface_VERIFY()
    #FortranCInterface_VERIFY(CXX)
    FortranCInterface_HEADER(f90_mangling.h
                             MACRO_NAMESPACE "FC_"
                             SYMBOLS ${FC_MANGLING_SUB})
ENDIF()

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
INCLUDE(QcMatrixFunctions)
IF(QCMATRIX_BUILD_ADAPTER)
    IF(ADAPTER_BLOCK_REAL)
        INCLUDE(AdapterBlockReal)
        INCLUDE(CmplxMat)
    ELSEIF(ADAPTER_CMPLX_MAT)
        INCLUDE(AdapterCmplxMat)
        INCLUDE(BlockCmplx)
    ELSEIF(ADAPTER_REAL_MAT)
        INCLUDE(AdapterRealMat)
        INCLUDE(CmplxMat)
        INCLUDE(BlockCmplx)
    ELSEIF(NOT ADAPTER_BLOCK_CMPLX)
        MESSAGE(FATAL_ERROR "Unknown matrix type of the adapter!")
    ENDIF()
ELSE()
    INCLUDE(RealMat)
    INCLUDE(CmplxMat)
    INCLUDE(BlockCmplx)
ENDIF()
# Adds the source codes of the adapter
IF(DEFINED ADAPTER_SRCS)
    SET(QCMATRIX_SRCS ${QCMATRIX_SRCS} ${ADAPTER_SRCS})
ENDIF()
# Adds the source codes of APIs
IF(DEFINED API_SRCS)
    SET(QCMATRIX_SRCS ${QCMATRIX_SRCS} ${API_SRCS})
ENDIF()

# Source codes of C test suite
SET(QCMATRIX_C_TEST_SRCS
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatValues.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatDuplicate.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatZeroEntries.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatGetTrace.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatScale.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatAXPY.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatTranspose.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatGEMM.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatGetMatProdTrace.c
    ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatrix.c)
    IF(QCMATRIX_ENABLE_VIEW)
        SET(QCMATRIX_C_TEST_SRCS
            ${QCMATRIX_C_TEST_SRCS}
            ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatView.c)
    ENDIF()
    #IF(QCMATRIX_BUILD_ADAPTER AND ADAPTER_C_LANG)
    #    SET(QCMATRIX_C_TEST_SRCS
    #        ${QCMATRIX_C_TEST_SRCS}
    #        ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatSetExternalMat.c
    #        ${LIB_QCMATRIX_PATH}/tests/c/api/test_c_QcMatGetExternalMat.c)
    #ENDIF()
# Source codes of Fortran 90 test suite
IF(${QCMATRIX_Fortran_API} STREQUAL "F90")
    SET(QCMATRIX_F_TEST_SRCS
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatValues.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatDuplicate.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatZeroEntries.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatGetTrace.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatScale.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatAXPY.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatTranspose.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatGEMM.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatGetMatProdTrace.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatrix.F90)
    IF(QCMATRIX_ENABLE_VIEW)
        SET(QCMATRIX_F_TEST_SRCS
            ${QCMATRIX_F_TEST_SRCS}
            ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatView.F90)
    ENDIF()
    IF(QCMATRIX_BUILD_ADAPTER AND (ADAPTER_F90_LANG OR ADAPTER_F03_LANG))
        SET(QCMATRIX_F_TEST_SRCS
            ${QCMATRIX_F_TEST_SRCS}
            ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatSetExternalMat.F90
            ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatGetExternalMat.F90)
    ENDIF()
# Source codes of Fortran 2003 test suite
ELSEIF(${QCMATRIX_Fortran_API} STREQUAL "F03")
    SET(QCMATRIX_F_TEST_SRCS
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatValues.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatDuplicate.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatZeroEntries.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatGetTrace.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatScale.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatAXPY.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatTranspose.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatGEMM.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatGetMatProdTrace.F90
        ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatrix.F90)
    IF(QCMATRIX_ENABLE_VIEW)
        SET(QCMATRIX_F_TEST_SRCS
            ${QCMATRIX_F_TEST_SRCS}
            ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatView.F90)
    ENDIF()
    IF(QCMATRIX_BUILD_ADAPTER AND (ADAPTER_F90_LANG OR ADAPTER_F03_LANG))
        SET(QCMATRIX_F_TEST_SRCS
            ${QCMATRIX_F_TEST_SRCS}
            ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatSetExternalMat.F90
            ${LIB_QCMATRIX_PATH}/tests/f90/api/test_f_QcMatGetExternalMat.F90)
    ENDIF()
ENDIF()

# The test suite needs BLAS and LAPACK libraries
IF(QCMATRIX_TEST_EXECUTABLE OR QCMATRIX_TEST_3M_METHOD)
    IF(NOT DEFINED BLAS_LIBRARIES)
        INCLUDE(FindBLAS)
        IF(NOT BLAS_FOUND)
            MESSAGE(FATAL_ERROR "BLAS library is needed for the test suite of QcMatrix!")
        ENDIF()
    ENDIF()
    IF(NOT DEFINED LAPACK_LIBRARIES)
        INCLUDE(FindLAPACK)
        IF(NOT LAPACK_FOUND)
            MESSAGE(FATAL_ERROR "LAPACK library is needed for the test suite of QcMatrix!")
        ENDIF()
    ENDIF()
    MESSAGE("BLAS library: ${BLAS_LIBRARIES}")
    MESSAGE("LAPACK library: ${LAPACK_LIBRARIES}")
ENDIF()

# interface of BLAS and LAPACK libraries
SET(QCMATRIX_BLAS_SRCS
    ${LIB_QCMATRIX_PATH}/src/lapack/qcmatrix_blas.F90
    ${LIB_QCMATRIX_PATH}/src/lapack/qcmatrix_c_blas.F90)

# Generates a configure file
INCLUDE(QcMatrixConfig.h)
#CONFIGURE_FILE(${LIB_QCMATRIX_PATH}/cmake/QcMatrixConfig.h.cmake
#               ${CMAKE_BINARY_DIR}/qcmatrix_config.h) 

# Builds executables for the test suite
IF(QCMATRIX_TEST_EXECUTABLE)
    ADD_DEFINITIONS(-DQCMATRIX_TEST_EXECUTABLE)
    # Creates the QcMatrix library
    ADD_LIBRARY(${LIB_QCMATRIX_NAME}
                ${QCMATRIX_BLAS_SRCS}
                ${QCMATRIX_SRCS})
    IF(QCMATRIX_ENABLE_HDF5)
        TARGET_LINK_LIBRARIES(${LIB_QCMATRIX_NAME}
                              ${LIB_HDF5_RELEASE})
    ENDIF()
    # Creates C test suite
    ADD_EXECUTABLE(${LIB_QCMATRIX_NAME}_c_test ${QCMATRIX_C_TEST_SRCS})
    # Linked libraries for the C test suite
    TARGET_LINK_LIBRARIES(${LIB_QCMATRIX_NAME}_c_test
                          ${LIB_QCMATRIX_NAME})
    IF(QCMATRIX_BUILD_ADAPTER)
        TARGET_LINK_LIBRARIES(${LIB_QCMATRIX_NAME}_c_test
                              ${QCMATRIX_EXTERNAL_LIBRARIES})
    ENDIF()
    TARGET_LINK_LIBRARIES(${LIB_QCMATRIX_NAME}_c_test
                          ${LAPACK_LIBRARIES}
                          ${BLAS_LIBRARIES})
    SET_TARGET_PROPERTIES(${LIB_QCMATRIX_NAME}_c_test
                          PROPERTIES
                          APPEND_STRING PROPERTY LINK_FLAGS
                          "${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
    # We need the flag -nofor-main for ifort compiler to get rid of the error
    # multiple definition of `main'
    GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    IF(Fortran_COMPILER_NAME MATCHES "ifort.*")
        SET_TARGET_PROPERTIES(${LIB_QCMATRIX_NAME}_c_test
                              PROPERTIES
                              APPEND_STRING PROPERTY LINK_FLAGS
                              "-nofor-main")
    ENDIF()
    # Creates Fortran 90/2003 test suite
    IF(DEFINED QCMATRIX_F_TEST_SRCS)
        ADD_EXECUTABLE(${LIB_QCMATRIX_NAME}_f_test ${QCMATRIX_F_TEST_SRCS})
        # Linked libraries for the Fortran 90/2003 test suite
        TARGET_LINK_LIBRARIES(${LIB_QCMATRIX_NAME}_f_test
                              ${LIB_QCMATRIX_NAME})
        IF(QCMATRIX_BUILD_ADAPTER)
            TARGET_LINK_LIBRARIES(${LIB_QCMATRIX_NAME}_f_test
                                  ${QCMATRIX_EXTERNAL_LIBRARIES})
        ENDIF()
        TARGET_LINK_LIBRARIES(${LIB_QCMATRIX_NAME}_f_test
                              ${LAPACK_LIBRARIES}
                              ${BLAS_LIBRARIES})
        SET_TARGET_PROPERTIES(${LIB_QCMATRIX_NAME}_f_test
                              PROPERTIES
                              APPEND_STRING PROPERTY LINK_FLAGS
                              "${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
    ENDIF()
# Creates the QcMatrix library with the test suite
ELSE()
    ADD_LIBRARY(${LIB_QCMATRIX_NAME}
                ${QCMATRIX_BLAS_SRCS}
                ${QCMATRIX_SRCS}
                ${QCMATRIX_C_TEST_SRCS}
                ${QCMATRIX_F_TEST_SRCS})
    IF(QCMATRIX_ENABLE_HDF5)
        TARGET_LINK_LIBRARIES(${LIB_QCMATRIX_NAME}
                              ${LIB_HDF5_RELEASE})
    ENDIF()
ENDIF()

# Creates an executable to test the efficiency of 3M method by changing the option QCMATRIX_3M_METHOD
IF(QCMATRIX_TEST_3M_METHOD)
    ADD_EXECUTABLE(test_3M_method
                   ${LIB_QCMATRIX_PATH}/src/qcmat/qcmatrix_error.F90
                   ${LIB_QCMATRIX_PATH}/tests/f90/timer.F90
                   ${LIB_QCMATRIX_PATH}/tests/f90/test_3M_method.F90)
    TARGET_LINK_LIBRARIES(test_3M_method ${LIB_QCMATRIX_NAME})
    IF(QCMATRIX_BUILD_ADAPTER)
        TARGET_LINK_LIBRARIES(test_3M_method ${QCMATRIX_EXTERNAL_LIBRARIES})
    ENDIF()
    TARGET_LINK_LIBRARIES(test_3M_method
                          ${LAPACK_LIBRARIES}
                          ${BLAS_LIBRARIES})
    SET_TARGET_PROPERTIES(test_3M_method
                          PROPERTIES
                          APPEND_STRING PROPERTY LINK_FLAGS
                          "${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
ENDIF()

#
INSTALL(TARGETS ${LIB_QCMATRIX_NAME} ARCHIVE DESTINATION lib)
