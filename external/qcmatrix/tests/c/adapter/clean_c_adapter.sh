#!/bin/bash

# removes files and directories
rm -rv ${PATH_LIB_QCMATRIX}/tests/c/adapter/build
rm -rv ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/api
rm -rv ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/lapack
rm -rv ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
rm -rv ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/utilities
rm -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/cmake/RealMat.cmake
rm -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/cmake/CmplxMat.cmake
rm -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/cmake/QcMatrixConfig.h.cmake
rm -rv ${PATH_LIB_QCMATRIX}/tests/c/adapter/src
