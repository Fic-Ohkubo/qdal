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
!!  This file tests the function QcMatTranspose_f().
!!
!!  2014-03-28, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatTranspose.F90"

    !% \brief tests the function QcMatTranspose_f()
    !  \author Bin Gao
    !  \date 2014-03-28
    !  \param[QcMat:type]{in} A the matrix
    !% \param[integer]{in} io_log IO of the logfile
    subroutine test_f_QcMatTranspose(A, io_log)
        use qcmatrix_f, only: QINT,                  &
                              QREAL,                 &
                              QcMat,                 &
                              QcMatIsAssembled_f,    &
                              QcMatCreate_f,         &
                              QcMatGetDimBlock_f,    &
                              QcMatGetDimMat_f,      &
                              QcMatGetAllValues_f,   &
                              QcMatTranspose_f,      &
                              MAT_NO_OPERATION,      &
                              MAT_TRANSPOSE,         &
                              MAT_HERM_TRANSPOSE,    &
                              MAT_COMPLEX_CONJUGATE, &
                              QcMatIsEqual_f,        &
                              QcMatCfArray_f,        &
                              QcMatDestroy_f,        &
#if defined(QCMATRIX_ENABLE_VIEW)
                              QcMatWrite_f,          &
                              ASCII_VIEW,            &
#endif
                              QSUCCESS
        implicit none
        type(QcMat), intent(in) :: A
        integer(kind=4), intent(in) :: io_log
        type(QcMat) B                                    !duplication of the matrix A
        logical(kind=4) assembled                        !indicates if the matrix is assembled or not
        integer(kind=QINT) dim_block                     !dimension of blocks
        integer(kind=QINT) dim_mat                       !dimension of each block
        integer(kind=QINT) size_values                   !number of elements in the matrix
        real(kind=QREAL), allocatable :: values_real(:)  !values of the real part
        real(kind=QREAL), allocatable :: values_imag(:)  !values of the imaginary part
        real(kind=QREAL) value_tmp                       !temporary value
        integer(kind=QINT), parameter :: all_mat_operations(3) = (/MAT_TRANSPOSE,      &
                                                                   MAT_HERM_TRANSPOSE, &
                                                                   MAT_COMPLEX_CONJUGATE/)
        logical(kind=4) :: cf_values = .true.            !comparing values
        logical(kind=4) :: row_major = .false.           !values in column major order
        logical(kind=4) is_equal                         !indicates if the matrix and array have the same values
        integer(kind=QINT) iop, icol, irow, ival, jval   !incremental recorders
        integer(kind=4) ierr                             !error information
        ! checks if the matrix is assembled
        ierr = QcMatIsAssembled_f(A, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix A is not assembled ..."
                write(io_log,100) "QcMatTranspose_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! duplicates the matrix A, and uses the duplication for the test
        ierr = QcMatCreate_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ierr = QcMatTranspose_f(MAT_NO_OPERATION, A, B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! checks if B = A
        ierr = QcMatIsEqual_f(A, B, cf_values, is_equal)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        if (is_equal) then
            write(io_log,100) "QcMatTranspose_f(MAT_NO_OPERATION, A) passed ..."
        else
            ! dumps results to check
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite_f(A, "QcMatTranspose_A", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ierr = QcMatWrite_f(B, "QcMatTranspose_B", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! gets the dimension of blocks
        ierr = QcMatGetDimBlock_f(A, dim_block)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimBlock_f(A) passed ..."
        ! gets the dimension of each block
        ierr = QcMatGetDimMat_f(A, dim_mat, dim_mat)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimMat_f(A) passed ..."
        size_values = dim_block*dim_block*dim_mat*dim_mat
        ! allocates memory for the elements of the matrix
        allocate(values_real(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(values_imag(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        do iop = 1, 3
            ! tests QcMatTranspose_f() out-of-place with option all_mat_operations(iop)
            ierr = QcMatTranspose_f(all_mat_operations(iop), A, B)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ! gets all the values of the matrix
            ierr = QcMatGetAllValues_f(A, row_major, size_values, values_real, values_imag)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            write(io_log,100) "QcMatGetAllValues_f(A) passed ..."
            ! manually transposes
            select case (all_mat_operations(iop))
            case (MAT_TRANSPOSE)
                ival = -dim_block*dim_mat
                do icol = 1, dim_block*dim_mat
                    ival = ival+dim_block*dim_mat
                    do irow = 1, icol-1
                        jval = (irow-1)*dim_block*dim_mat+icol
                        ! real part
                        value_tmp = values_real(ival+irow)
                        values_real(ival+irow) = values_real(jval)
                        values_real(jval) = value_tmp
                        ! imaginary part
                        value_tmp = values_imag(ival+irow)
                        values_imag(ival+irow) = values_imag(jval)
                        values_imag(jval) = value_tmp
                    end do
                end do
            case (MAT_HERM_TRANSPOSE)
                ival = -dim_block*dim_mat
                do icol = 1, dim_block*dim_mat
                    ival = ival+dim_block*dim_mat
                    do irow = 1, icol
                        jval = (irow-1)*dim_block*dim_mat+icol
                        ! real part
                        value_tmp = values_real(ival+irow)
                        values_real(ival+irow) = values_real(jval)
                        values_real(jval) = value_tmp
                        ! imaginary part
                        value_tmp = values_imag(ival+irow)
                        values_imag(ival+irow) = -values_imag(jval)
                        values_imag(jval) = -value_tmp
                    end do
                end do
            case (MAT_COMPLEX_CONJUGATE)
                values_imag = -values_imag
            case default
                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end select
            ierr = QcMatCfArray_f(B,           &
                                  row_major,   &
                                  size_values, &
                                  values_real, &
                                  values_imag, &
                                  is_equal)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            if (is_equal) then
                write(io_log,100) "QcMatTranspose_f(A)",   &
                                  all_mat_operations(iop), &
                                  " out-of-place passed ..."
            else
                ! dumps results to check
#if defined(QCMATRIX_ENABLE_VIEW)
                ierr = QcMatWrite_f(A, "QcMatTranspose_A", ASCII_VIEW)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                ierr = QcMatWrite_f(B, "QcMatTranspose_B", ASCII_VIEW)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
                write(io_log,100) "real part of B from hand coding"
                write(io_log,"(4Es24.12)") values_real
                write(io_log,100) "imaginary part of B from hand coding"
                write(io_log,"(4Es24.12)") values_imag
                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end if
            ! tests QcMatTranspose_f() in-place with option all_mat_operations(iop)
            ierr = QcMatTranspose_f(all_mat_operations(iop), B, B)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ! manually transposes
            select case (all_mat_operations(iop))
            case (MAT_TRANSPOSE)
                ival = -dim_block*dim_mat
                do icol = 1, dim_block*dim_mat
                    ival = ival+dim_block*dim_mat
                    do irow = 1, icol-1
                        jval = (irow-1)*dim_block*dim_mat+icol
                        ! real part
                        value_tmp = values_real(ival+irow)
                        values_real(ival+irow) = values_real(jval)
                        values_real(jval) = value_tmp
                        ! imaginary part
                        value_tmp = values_imag(ival+irow)
                        values_imag(ival+irow) = values_imag(jval)
                        values_imag(jval) = value_tmp
                    end do
                end do
            case (MAT_HERM_TRANSPOSE)
                ival = -dim_block*dim_mat
                do icol = 1, dim_block*dim_mat
                    ival = ival+dim_block*dim_mat
                    do irow = 1, icol
                        jval = (irow-1)*dim_block*dim_mat+icol
                        ! real part
                        value_tmp = values_real(ival+irow)
                        values_real(ival+irow) = values_real(jval)
                        values_real(jval) = value_tmp
                        ! imaginary part
                        value_tmp = values_imag(ival+irow)
                        values_imag(ival+irow) = -values_imag(jval)
                        values_imag(jval) = -value_tmp
                    end do
                end do
            case (MAT_COMPLEX_CONJUGATE)
                values_imag = -values_imag
            case default
                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end select
            ierr = QcMatCfArray_f(B,           &
                                  row_major,   &
                                  size_values, &
                                  values_real, &
                                  values_imag, &
                                  is_equal)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            if (is_equal) then
                write(io_log,100) "QcMatTranspose_f(B)",   &
                                  all_mat_operations(iop), &
                                  " in-place passed ..."
            else
                ! dumps results to check
#if defined(QCMATRIX_ENABLE_VIEW)
                ierr = QcMatWrite_f(B, "QcMatTranspose_B", ASCII_VIEW)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
                write(io_log,100) "real part of B from hand coding"
                write(io_log,"(4Es24.12)") values_real
                write(io_log,100) "imaginary part of B from hand coding"
                write(io_log,"(4Es24.12)") values_imag
                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end if
        end do
        ! cleans
        deallocate(values_real)
        deallocate(values_imag)
        ierr = QcMatDestroy_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        return
100     format("test_f_QcMatTranspose>> ",A,I4,A)
    end subroutine test_f_QcMatTranspose

#undef QCMATRIX_F_TEST_SRC
