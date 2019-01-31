#!/bin/bash

# creates directories
mkdir -pv ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/api
mkdir -v ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/lapack
mkdir -v ${PATH_LIB_QCMATRIX}/tests/f90/adapter/cmake

# copies files
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_f_basic.h90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_f_mat_data.h90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_f_mat_duplicate.h90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_f_mat_operations.h90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_f_mat_storage.h90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_f_mat_symmetry.h90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_f_mat_view.h90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/lapack/qcmatrix_blas.h90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/include/lapack
cp -v ${PATH_LIB_QCMATRIX}/cmake/QcMatrixConfig.h.cmake ${PATH_LIB_QCMATRIX}/tests/f90/adapter/cmake
cp -v ${PATH_LIB_QCMATRIX}/src/lapack/qcmatrix_blas.F90 ${PATH_LIB_QCMATRIX}/tests/f90/adapter/src

# renames the names of subroutines in BLAS/LAPACK interface
sed -i "s/subroutine BLAS_/subroutine Real_BLAS_/g" ${PATH_LIB_QCMATRIX}/tests/f90/adapter/src/qcmatrix_blas.F90
sed -i "s/call QErrorExit/call error_exit/g" ${PATH_LIB_QCMATRIX}/tests/f90/adapter/src/qcmatrix_blas.F90
