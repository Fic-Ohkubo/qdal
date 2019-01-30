#!/bin/bash

# creates directories
mkdir -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/api
mkdir -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/lapack
mkdir -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
mkdir -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/utilities
mkdir -pv ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/lapack
mkdir -pv ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/real_mat
mkdir -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/cmplx_mat
mkdir -v ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/qcmat

# copies files
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_f_basic.h90 ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/api/qcmatrix_c_type.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/api
cp -v ${PATH_LIB_QCMATRIX}/include/lapack/qcmatrix_blas.h90 ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/lapack
cp -v ${PATH_LIB_QCMATRIX}/include/lapack/qcmatrix_c_blas.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/lapack
cp -v ${PATH_LIB_QCMATRIX}/include/types/mat_data.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
cp -v ${PATH_LIB_QCMATRIX}/include/types/mat_duplicate.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
cp -v ${PATH_LIB_QCMATRIX}/include/types/mat_operations.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
cp -v ${PATH_LIB_QCMATRIX}/include/types/mat_storage.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
cp -v ${PATH_LIB_QCMATRIX}/include/types/mat_symmetry.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
cp -v ${PATH_LIB_QCMATRIX}/include/types/mat_view.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
cp -v ${PATH_LIB_QCMATRIX}/include/types/qcmatrix_basic_types.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
cp -v ${PATH_LIB_QCMATRIX}/include/types/qcmatrix_mat_types.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/types
cp -v ${PATH_LIB_QCMATRIX}/include/utilities/qcmatrix_error.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/utilities
cp -v ${PATH_LIB_QCMATRIX}/include/utilities/qcmatrix_algebra.h ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/utilities
#
cp -v ${PATH_LIB_QCMATRIX}/cmake/RealMat.cmake ${PATH_LIB_QCMATRIX}/tests/c/adapter/cmake
cp -v ${PATH_LIB_QCMATRIX}/cmake/CmplxMat.cmake ${PATH_LIB_QCMATRIX}/tests/c/adapter/cmake
cp -v ${PATH_LIB_QCMATRIX}/cmake/QcMatrixConfig.h.cmake ${PATH_LIB_QCMATRIX}/tests/c/adapter/cmake
#
cp -v ${PATH_LIB_QCMATRIX}/tests/f90/adapter/src/error_exit.F90 ${PATH_LIB_QCMATRIX}/tests/c/adapter/src
cp -v ${PATH_LIB_QCMATRIX}/src/lapack/qcmatrix_blas.F90 ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/lapack
cp -v ${PATH_LIB_QCMATRIX}/src/lapack/qcmatrix_c_blas.F90 ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/lapack
cp -v ${PATH_LIB_QCMATRIX}/src/real_mat/RealMat*.c ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/real_mat
cp -v ${PATH_LIB_QCMATRIX}/src/cmplx_mat/CmplxMat*.c ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/cmplx_mat
cp -v ${PATH_LIB_QCMATRIX}/src/qcmat/QcMat*.c ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/qcmat

# renames the names of subroutines in BLAS/LAPACK interface
sed -i "s/C_BLAS_/Real_C_BLAS_/g" ${PATH_LIB_QCMATRIX}/tests/c/adapter/include/lapack/qcmatrix_c_blas.h
sed -i "s/subroutine BLAS_/subroutine Real_BLAS_/g" ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/lapack/qcmatrix_blas.F90
sed -i "s/call QErrorExit/call error_exit/g" ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/lapack/qcmatrix_blas.F90
sed -i "s/C_BLAS_/Real_C_BLAS_/g" ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/lapack/qcmatrix_c_blas.F90
sed -i "s/call BLAS_/call Real_BLAS_/g" ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/lapack/qcmatrix_c_blas.F90
# renames the names of struct and functions
cd ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/real_mat
for file in `ls RealMat*.c`
do
    sed -i "s/RealMat/Real_Mat_/g" ${file}
    sed -i "s/C_BLAS_/Real_C_BLAS_/g" ${file}
done
cd ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/cmplx_mat
for file in `ls CmplxMat*.c`
do
    sed -i "s/RealMat/Real_Mat_/g" ${file}
    sed -i "s/CmplxMat/Cmplx_Mat_/g" ${file}
done
# processes source codes in qcmat for square block real matrix
cd ${PATH_LIB_QCMATRIX}/tests/c/adapter/src/qcmat
for file in `ls QcMat*.c`
do
    sed -i "s/CmplxMat/Real_Mat_/g" ${file}
    sed -i "s/QcMat/Block_Real_/g" ${file}
done
