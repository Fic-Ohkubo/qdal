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
!!  This file tests the functions QcMatSetValues_f(), QcMatAddValues_f()
!!  and QcMatGetValues_f().
!!
!!  2014-08-22, Bin Gao:
!!  * adds the test of QcMatAddValues_f()
!!
!!  2014-03-25, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatValues.F90"

    !% \brief tests the functions QcMatSetValues_f(), QcMatAddValues_f()
    !         and QcMatGetValues_f()
    !  \author Bin Gao
    !  \date 2014-03-25
    !  \param[integer]{in} dim_block the dimension of blocks
    !% \param[integer]{in} io_log IO of the logfile
    subroutine test_f_QcMatValues(dim_block, io_log)
        use qcmatrix_f, only: QINT,               &
                              QREAL,              &
                              QcMat,              &
                              QcMatCreate_f,      &
                              QcMatBlockCreate_f, &
                              QcMatSetDataType_f, &
                              QcMatSetDimMat_f,   &
                              QcMatAssemble_f,    &
                              QcMatSetValues_f,   &
                              QcMatAddValues_f,   &
                              QcMatCfArray_f,     &
                              QcMatDestroy_f,     &
                              QCMPLXMAT,          &
#if defined(QCMATRIX_ENABLE_VIEW)
                              QcMatWrite_f,       &
                              ASCII_VIEW,         &
#endif
                              QSUCCESS
        implicit none
        integer(kind=QINT), intent(in) :: dim_block
        integer(kind=4), intent(in) :: io_log
        type(QcMat) A                                                !matrix for test
        integer(kind=QINT) num_blocks                                !number of blocks, as \var{dim_block}*\var{dim_block}
        integer(kind=QINT), allocatable :: idx_block_row(:)          !row indices of the blocks
        integer(kind=QINT), allocatable :: idx_block_col(:)          !column indices of the blocks
        integer(kind=QINT), allocatable :: data_type(:)              !data types of the blocks
        integer(kind=QINT), parameter :: dim_mat = 6                 !dimension of each block
        integer(kind=QINT), parameter :: size_mat = dim_mat*dim_mat  !number of elements in each block
#if defined(QCMATRIX_ZERO_BASED)
        integer(kind=QINT), parameter :: idx_first_row = 0           !index of the first row from which the values are set/added
        integer(kind=QINT), parameter :: num_row_set = dim_mat       !number of rows to set/add
        integer(kind=QINT), parameter :: idx_first_col = 0           !index of the first column from which the values are set/added
        integer(kind=QINT), parameter :: num_col_set = dim_mat       !number of columns to set/add
#else
        integer(kind=QINT), parameter :: idx_first_row = 1           !index of the first row from which the values are set/added
        integer(kind=QINT), parameter :: num_row_set = dim_mat       !number of rows to set/add
        integer(kind=QINT), parameter :: idx_first_col = 1           !index of the first column from which the values are set/added
        integer(kind=QINT), parameter :: num_col_set = dim_mat       !number of columns to set/add
#endif
        real(kind=QREAL), allocatable :: values_real(:)              !values of the real part
        logical(kind=4) :: row_major = .false.                       !values in column major order
        logical(kind=4) is_equal                                     !indicates if the matrix and array have the same values
        integer(kind=QINT) iblk, jblk, kblk                          !incremental recorders for blocks
        integer(kind=4) ierr                                         !error information
        ierr = QcMatCreate_f(A)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatCreate_f(A) passed ..."
        ierr = QcMatBlockCreate_f(A, dim_block)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatBlockCreate_f(A) passed ..."
        ! sets the data types of the blocks
        num_blocks = dim_block*dim_block
        allocate(idx_block_row(num_blocks), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(idx_block_col(num_blocks), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(data_type(num_blocks), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        kblk = 0
#if defined(QCMATRIX_ZERO_BASED)
        do iblk = 0, dim_block-1
            do jblk = 0, dim_block-1
#else
        do iblk = 1, dim_block
            do jblk = 1, dim_block
#endif
                kblk = kblk+1
                ! QcMatrix uses row major order for the blocks
                idx_block_row(kblk) = iblk
                idx_block_col(kblk) = jblk
                data_type(kblk) = QCMPLXMAT
            end do
        end do
        ierr = QcMatSetDataType_f(A, num_blocks, idx_block_row, idx_block_col, data_type)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        deallocate(data_type)
        ! sets the dimension of each block
        ierr = QcMatSetDimMat_f(A, dim_mat, dim_mat)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! assembles the matrix
        ierr = QcMatAssemble_f(A)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! allocates memory for setting the elements of each block
        allocate(values_real(num_blocks*size_mat), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! we sets all the elements as 1+i
        values_real = 1.0_QREAL
        ! sets the elements block by block
        do iblk = 1, num_blocks
            ierr = QcMatSetValues_f(A=A,                                 &
                                    idx_block_row=idx_block_row(iblk),   &
                                    idx_block_col=idx_block_col(iblk),   &
                                    idx_first_row=idx_first_row,         &
                                    num_row_set=num_row_set,             &
                                    idx_first_col=idx_first_col,         &
                                    num_col_set=num_col_set,             &
                                    values_real=values_real(1:size_mat), &
                                    values_imag=values_real(1:size_mat))
            if (ierr/=QSUCCESS) then
                write(io_log,100) "block (", idx_block_row(iblk), ",", &
                                  idx_block_col(iblk), ")"
                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end if
        end do
        ! checks the elements set by QcMatSetValues()
        ierr = QcMatCfArray_f(A,                   &
                              row_major,           &
                              num_blocks*size_mat, &
                              values_real,         &
                              values_real,         &
                              is_equal)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        if (is_equal) then
            write(io_log,100) "QcMatSetValues_f(A) passed ..."
        else
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite_f(A, "QcMatValues_A", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! adds the elements block by block
        do iblk = 1, num_blocks
            ierr = QcMatAddValues_f(A=A,                                 &
                                    idx_block_row=idx_block_row(iblk),   &
                                    idx_block_col=idx_block_col(iblk),   &
                                    idx_first_row=idx_first_row,         &
                                    num_row_add=num_row_set,             &
                                    idx_first_col=idx_first_col,         &
                                    num_col_add=num_col_set,             &
                                    values_real=values_real(1:size_mat), &
                                    values_imag=values_real(1:size_mat))
            if (ierr/=QSUCCESS) then
                write(io_log,100) "block (", idx_block_row(iblk), ",", &
                                  idx_block_col(iblk), ")"
                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end if
        end do
        ! all the elements become 2+2*i
        values_real = 2.0_QREAL
        ! checks the elements set by QcMatAddValues()
        ierr = QcMatCfArray_f(A,                   &
                              row_major,           &
                              num_blocks*size_mat, &
                              values_real,         &
                              values_real,         &
                              is_equal)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        if (is_equal) then
            write(io_log,100) "QcMatAddValues_f(A) passed ..."
        else
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite_f(A, "QcMatValues_A", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! cleans up
        deallocate(idx_block_row)
        deallocate(idx_block_col)
        deallocate(values_real)
        ierr = QcMatDestroy_f(A)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        return
100     format("test_f_QcMatValues>> ",A,I4,A,I4,A)
    end subroutine test_f_QcMatValues

#undef QCMATRIX_F_TEST_SRC
