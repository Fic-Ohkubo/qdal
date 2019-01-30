!!  QcMatrix: an abstract matrix library
!!  Copyright 2012-2015 Bin Gao
!!
!!  QcMatrix is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  QcMatrix is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with QcMatrix. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file is the test suite of Fortran 90 APIs.
!!
!!  2014-03-13, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/test_3M_method.F90"

    program test_3M_method
        use qcmatrix_f, only: QINT,              &
                              QREAL,             &
                              QcMat,             &
                              QSYMMAT,           &
                              QcMatCreate_f,     &
                              QcMatSetRandMat_f, &
                              QcMatDestroy_f,    &
                              QCMPLXMAT,         &
                              QcMatGEMM_f,       &
                              MAT_NO_OPERATION
        implicit none
        integer(kind=4), parameter :: STDOUT = 6
        type(QcMat) A, B, C                                                  !matrices for tests
        integer(kind=QINT) dim_block                                         !dimension of blocks
        integer(kind=QINT) dim_mat                                           !dimension of each block
        real(kind=QREAL), parameter :: alpha(2) = (/0.01_QREAL,0.01_QREAL/)  !the scalar number
        real(kind=QREAL), parameter :: beta(2) = (/0.01_QREAL,0.01_QREAL/)   !the scalar number
        real(kind=QREAL) curr_time                                           !current time
        integer(kind=4) ierr                                                 !error information
        ! tests different numbers of blocks
        do dim_block = 1, MAX_DIM_BLOCK
            ! tests different dimension of each block
            do dim_mat = 100, 500, 100
                ierr = QcMatCreate_f(A)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatCreate_f(A) passed ..."
                ! generates a random matrix A according to its symmetry and data types
                ierr = QcMatSetRandMat_f(A,         &
                                         QSYMMAT,   &
                                         QCMPLXMAT, &
                                         dim_block, &
                                         dim_mat)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatSetRandMat_f(A) passed ..."
                ierr = QcMatCreate_f(B)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatCreate_f(B) passed ..."
                ! generates a random matrix B according to its symmetry and data types
                ierr = QcMatSetRandMat_f(B,         &
                                         QSYMMAT,   &
                                         QCMPLXMAT, &
                                         dim_block, &
                                         dim_mat)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatSetRandMat_f(B) passed ..."
                ierr = QcMatCreate_f(C)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatCreate_f(C) passed ..."
                ! generates a random matrix C according to its symmetry and data types
                ierr = QcMatSetRandMat_f(C,         &
                                         QSYMMAT,   &
                                         QCMPLXMAT, &
                                         dim_block, &
                                         dim_mat)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatSetRandMat_f(C) passed ..."
                write(STDOUT,100) "dimensions", dim_block, dim_mat
                call TimerGet(curr_time)
                ierr = QcMatGEMM_f(MAT_NO_OPERATION, &
                                   MAT_NO_OPERATION, &
                                   alpha,            &
                                   A,                &
                                   B,                &
                                   beta,             &
                                   C)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatGEMM_f(C) passed ..."
                call TimerView(curr_time, "QcMatGEMM_f", STDOUT)
                ! cleans
                ierr = QcMatDestroy_f(C)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatDestroy_f(C) passed ..."
                ierr = QcMatDestroy_f(B)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatDestroy_f(B) passed ..."
                ierr = QcMatDestroy_f(A)
                call QErrorCheckCode(STDOUT, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(STDOUT,100) "QcMatDestroy_f(A) passed ..."
            end do
        end do
100     format("test_3M_method>> ",A,2I8)
    end program test_3M_method

#undef QCMATRIX_F_TEST_SRC
