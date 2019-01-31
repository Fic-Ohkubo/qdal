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
!!  This file tests the function QcMatZeroEntries_f().
!!
!!  2014-03-28, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatZeroEntries.F90"

    !% \brief tests the function QcMatZeroEntries_f()
    !  \author Bin Gao
    !  \date 2014-03-28
    !  \param[QcMat:type]{in} A the matrix
    !% \param[integer]{in} io_log IO of the logfile
    subroutine test_f_QcMatZeroEntries(A, io_log)
        use qcmatrix_f, only: QINT,                   &
                              QREAL,                  &
                              QcMat,                  &
                              QcMatIsAssembled_f,     &
                              QcMatCreate_f,          &
                              QcMatDuplicate_f,       &
                              QcMatGetDimBlock_f,     &
                              QcMatGetDimMat_f,       &
                              QcMatZeroEntries_f,     &
                              QcMatCfArray_f,         &
                              QcMatDestroy_f,         &
                              COPY_PATTERN_AND_VALUE, &
#if defined(QCMATRIX_ENABLE_VIEW)
                              QcMatWrite_f,           &
                              ASCII_VIEW,             &
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
        logical(kind=4) :: row_major = .false.           !values in column major order
        logical(kind=4) is_equal                         !indicates if the matrix and array have the same values
        integer(kind=4) ierr                             !error information
        ! checks if the matrix is assembled
        ierr = QcMatIsAssembled_f(A, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix A is not assembled ..."
                write(io_log,100) "QcMatZeroEntries_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! duplicates the matrix, and uses duplication for the test
        ierr = QcMatCreate_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ierr = QcMatDuplicate_f(A, COPY_PATTERN_AND_VALUE, B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! gets the dimension of blocks
        ierr = QcMatGetDimBlock_f(B, dim_block)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimBlock_f(B) passed ..."
        ! gets the dimension of each block
        ierr = QcMatGetDimMat_f(B, dim_mat, dim_mat)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimMat_f(B) passed ..."
        size_values = dim_block*dim_block*dim_mat*dim_mat
        ! allocates memory for the elements of the matrix
        allocate(values_real(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        values_real = 0.0_QREAL
        ! zeros all entries of the matrix by QcMatZeroEntries_f()
        ierr = QcMatZeroEntries_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ierr = QcMatCfArray_f(B,           &
                              row_major,   &
                              size_values, &
                              values_real, &
                              values_real, &
                              is_equal)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        if (is_equal) then
            write(io_log,100) "QcMatZeroEntries_f(B) passed ..."
        else
            ! dumps results to check
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite_f(B, "QcMatZeroEntries_B", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! cleans
        deallocate(values_real)
        ierr = QcMatDestroy_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        return
100     format("test_f_QcMatZeroEntries>> ",A)
    end subroutine test_f_QcMatZeroEntries

#undef QCMATRIX_F_TEST_SRC
