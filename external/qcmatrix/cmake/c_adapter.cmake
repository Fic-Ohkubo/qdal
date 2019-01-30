# CMake file for C adapter
#
# Sets C compiler flags
SET(LANG_C_HEADER "lib_matrix.h" CACHE STRING "Name of header file of the external C library.")
SET(LANG_C_MATRIX "matrix_t" CACHE STRING "Name of external C matrix struct.")
ADD_DEFINITIONS(-DLANG_C_HEADER=${LANG_C_HEADER})
ADD_DEFINITIONS(-DLANG_C_MATRIX=${LANG_C_MATRIX})
# Source codes of the adapter for the external C library
SET(ADAPTER_SRCS ${LIB_QCMATRIX_PATH}/src/adapter/c_adapter.c)
