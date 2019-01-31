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

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatrix.F90"

#if defined(QCMATRIX_TEST_EXECUTABLE)
    program test_f_QcMatrix
#else
    subroutine test_f_QcMatrix(io_log)
#endif
        use qcmatrix_f, only: QINT,              &
                              QcMat,             &
                              QSYMMAT,           &
                              QANTISYMMAT,       &
                              QNONSYMMAT,        &
                              QcMatCreate_f,     &
                              QcMatSetRandMat_f, &
                              QcMatDestroy_f,    &
                              QREALMAT,          &
                              QIMAGMAT,          &
                              QCMPLXMAT
        implicit none
#if defined(QCMATRIX_TEST_EXECUTABLE)
        integer(kind=4), parameter :: io_log = 6
#else
        integer(kind=4), intent(in) :: io_log
#endif
        type(QcMat) A, B                              !matrices for tests
        integer(kind=QINT) dim_block                  !dimension of blocks
        integer(kind=QINT), parameter :: dim_mat = 6  !dimension of each block
        integer(kind=QINT), parameter :: sym_type(3) = (/QSYMMAT,QANTISYMMAT,QNONSYMMAT/)  !all symmetry types
        integer(kind=QINT), parameter :: data_type(3) = (/QREALMAT,QIMAGMAT,QCMPLXMAT/)    !all data types
        integer(kind=QINT) isym, jsym                 !incremental recorders for symmetry types
        integer(kind=QINT) idat, jdat                 !incremental recorders for data types
        integer(kind=4) ierr                          !error information
        ! tests different numbers of blocks
        do dim_block = 1, MAX_DIM_BLOCK
            ! tests QcMatSetValues() and QcMatGetValues()
            call test_f_QcMatValues(dim_block, io_log)
            ! tests different symmetry types (symmetric, anti-symmetric, non-symmetric) for matrix A
            do isym = 1, 3
                ! tests different data types (real, imaginary, complex) for matrix A
                do idat = 1, 3
                    ierr = QcMatCreate_f(A)
                    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                    write(io_log,100) "QcMatCreate_f(A) passed ..."
                    ! generates a random matrix A according to its symmetry and data types
                    ierr = QcMatSetRandMat_f(A,               &
                                             sym_type(isym),  &
                                             data_type(idat), &
                                             dim_block,       &
                                             dim_mat,         &
                                             dim_mat)
                    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                    write(io_log,100) "QcMatSetRandMat_f(A) passed ..."
                    ! tests QcMatDuplicate_f()
                    call test_f_QcMatDuplicate(A, io_log)
                    ! tests QcMatZeroEntries_f()
                    call test_f_QcMatZeroEntries(A, io_log)
                    ! tests QcMatGetTrace_f()
                    call test_f_QcMatGetTrace(A, io_log)
                    ! tests QcMatScale_f()
                    call test_f_QcMatScale(A, io_log)
                    ! tests QcMatAXPY_f()
                    ierr = QcMatCreate_f(B)
                    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                    call test_f_QcMatAXPY(A, B, io_log)
                    ierr = QcMatDestroy_f(B)
                    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                    ! tests QcMatTranspose_f() in-place and out-of-place
                    call test_f_QcMatTranspose(A, io_log)
#if defined(QCMATRIX_ENABLE_VIEW)
                    ! tests QcMatWrite_f() and QcMatRead_f()
                    call test_f_QcMatView(A,                        &
                                         "B"//char(dim_block+48)// &
                                         "S"//char(isym+48)//      &
                                         "D"//char(idat+48),       &
                                         io_log)
#endif
#if defined(ADAPTER_F90_LANG) || defined(ADAPTER_F03_LANG)
                    ! tests QcMatGetExternalMat_f()
                    call test_f_QcMatGetExternalMat(A, io_log)
                    ! tests QcMatSetExternalMat_f()
                    call test_f_QcMatSetExternalMat(A, io_log)
#endif
                    ! tests different symmetry types (symmetric, anti-symmetric, non-symmetric) for matrix B
                    do jsym = 1, 3
                        ! tests different data types (real, imaginary, complex) for matrix B
                        do jdat = 1, 3
                            ierr = QcMatCreate_f(B)
                            call QErrorCheckCode(io_log,   &
                                                 ierr,     &
                                                 __LINE__, &
                                                 QCMATRIX_F_TEST_SRC)
                            write(io_log,100) "QcMatCreate_f(B) passed ..."
                            ! generates a random matrix B according to its symmetry and data types
                            ierr = QcMatSetRandMat_f(B,               &
                                                     sym_type(jsym),  &
                                                     data_type(jdat), &
                                                     dim_block,       &
                                                     dim_mat,         &
                                                     dim_mat)
                            call QErrorCheckCode(io_log,   &
                                                 ierr,     &
                                                 __LINE__, &
                                                 QCMATRIX_F_TEST_SRC)
                            write(io_log,100) "QcMatSetRandMat_f(B) passed ..."
                            ! tests QcMatAXPY_f()
                            call test_f_QcMatAXPY(A, B, io_log)
                            ! tests QcMatGEMM_f()
                            call test_f_QcMatGEMM(A, B, io_log)
                            ! tests QcMatGetMatProdTrace_f()
                            call test_f_QcMatGetMatProdTrace(A, B, io_log)
                            ! cleans
                            ierr = QcMatDestroy_f(B)
                            call QErrorCheckCode(io_log,   &
                                                 ierr,     &
                                                 __LINE__, &
                                                 QCMATRIX_F_TEST_SRC)
                            write(io_log,100) "QcMatDestroy_f(B) passed ..."
                        end do
                    end do
                    ! cleans
                    ierr = QcMatDestroy_f(A)
                    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                    write(io_log,100) "QcMatDestroy_f(A) passed ..."
                end do
            end do
        end do
#if defined(QCMATRIX_TEST_EXECUTABLE)
100     format("test_f_QcMatrix>> ",A)
    end program test_f_QcMatrix
#else
        return
100     format("test_f_QcMatrix>> ",A)
    end subroutine test_f_QcMatrix
#endif

#undef QCMATRIX_F_TEST_SRC
