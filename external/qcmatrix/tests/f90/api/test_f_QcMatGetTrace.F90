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
!!  This file tests the function QcMatGetTrace_f().
!!
!!  2014-03-28, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatGetTrace.F90"

    !% \brief tests the function QcMatGetTrace_f()
    !  \author Bin Gao
    !  \date 2014-03-28
    !  \param[QcMat:type]{in} A the matrix
    subroutine test_f_QcMatGetTrace(A, io_log)
        use qcmatrix_f, only: QINT,               &
                              QREAL,              &
                              QcMat,              &
                              QcMatGetDimBlock_f, &
                              QcMatGetDataType_f, &
                              QcMatGetDimMat_f,   &
                              QcMatGetValues_f,   &
                              QcMatGetTrace_f,    &
                              QSUCCESS,           &
                              QREALMAT,           &
                              QIMAGMAT,           &
                              QCMPLXMAT,          &
                              QZEROTHRSH
        implicit none
        type(QcMat), intent(in) :: A
        integer(kind=4), intent(in) :: io_log
        integer(kind=QINT) dim_block                         !dimension of blocks
        integer(kind=QINT), allocatable :: idx_block_row(:)  !row indices of the blocks
        integer(kind=QINT), allocatable :: idx_block_col(:)  !column indices of the blocks
        integer(kind=QINT), allocatable :: data_type(:)      !data types of the blocks
        integer(kind=QINT) dim_mat                           !dimension of each block
        integer(kind=QINT) size_mat                          !number of elements in each block
        integer(kind=QINT) idx_first_row                     !index of the first row to get values
        integer(kind=QINT) num_row_get                       !number of rows to get
        integer(kind=QINT) idx_first_col                     !index of the first column to get values
        integer(kind=QINT) num_col_get                       !number of columns to get
        real(kind=QREAL), allocatable :: values_real(:)      !values of the real part
        real(kind=QREAL), allocatable :: values_imag(:)      !values of the imaginary part
        real(kind=QREAL), allocatable :: cf_trace(:)         !trace to compare with
        real(kind=QREAL), allocatable :: trace(:)            !trace from QcMatGetTrace()
        integer(kind=QINT) iblk, jblk, kblk                  !incremental recorders for blocks
        integer(kind=QINT) ival, jval                        !incremental recorders for elements of each block
        integer(kind=4) ierr                                 !error information
        ! gets the dimension of blocks
        ierr = QcMatGetDimBlock_f(A, dim_block)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimBlock_f(A) passed ..."
        ! gets the data types of the diagonal blocks
        allocate(idx_block_row(dim_block), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(idx_block_col(dim_block), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(data_type(0:dim_block-1), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        do iblk = 1, dim_block
#if defined(QCMATRIX_ZERO_BASED)
            idx_block_row(iblk) = iblk-1
#else
            idx_block_row(iblk) = iblk
#endif
            idx_block_col(iblk) = idx_block_row(iblk)
        end do
        ierr = QcMatGetDataType_f(A, dim_block, idx_block_row, idx_block_col, data_type)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        deallocate(idx_block_row)
        deallocate(idx_block_col)
        ! gets the dimension of each block
        ierr = QcMatGetDimMat_f(A, dim_mat, dim_mat)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimMat_f(A) passed ..."
        size_mat = dim_mat*dim_mat
        ! allocates memory for getting the elements of each diagonal block
        allocate(values_real(size_mat), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(values_imag(size_mat), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! calculates the trace block by block
        allocate(cf_trace(2*dim_block), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(QCMATRIX_ZERO_BASED)
        idx_first_row = 0
        idx_first_col = 0
#else
        idx_first_row = 1
        idx_first_col = 1
#endif
        num_row_get = dim_mat
        num_col_get = dim_mat
        do iblk = 0, dim_block-1
            jblk = 2*iblk+1
            kblk = 2*iblk+2
            cf_trace(jblk) = 0  !real part
            cf_trace(kblk) = 0  !imaginary part
            select case (data_type(iblk))
            ! real block
            case (QREALMAT)
                ierr = QcMatGetValues_f(A=A,                         &
#if defined(QCMATRIX_ZERO_BASED)
                                        idx_block_row=iblk,          &
                                        idx_block_col=iblk,          &
#else
                                        idx_block_row=iblk+1,        &
                                        idx_block_col=iblk+1,        &
#endif
                                        idx_first_row=idx_first_row, &
                                        num_row_get=num_row_get,     &
                                        idx_first_col=idx_first_col, &
                                        num_col_get=num_col_get,     &
                                        values_real=values_real)
                if (ierr==QSUCCESS) then
                    jval = -dim_mat
                    do ival = 1, dim_mat
                        jval = jval+dim_mat+1
                        cf_trace(jblk) = cf_trace(jblk)+values_real(jval)
                    end do
                else
                    write(io_log,100) "block (", iblk, ",", iblk, ")"
                    call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
                end if
            ! imaginary block
            case (QIMAGMAT)
                ierr = QcMatGetValues_f(A=A,                         &
#if defined(QCMATRIX_ZERO_BASED)
                                        idx_block_row=iblk,          &
                                        idx_block_col=iblk,          &
#else
                                        idx_block_row=iblk+1,        &
                                        idx_block_col=iblk+1,        &
#endif
                                        idx_first_row=idx_first_row, &
                                        num_row_get=num_row_get,     &
                                        idx_first_col=idx_first_col, &
                                        num_col_get=num_col_get,     &
                                        values_imag=values_imag)
                if (ierr==QSUCCESS) then
                    jval = -dim_mat
                    do ival = 1, dim_mat
                        jval = jval+dim_mat+1
                        cf_trace(kblk) = cf_trace(kblk)+values_imag(jval)
                    end do
                else
                    write(io_log,100) "block (", iblk, ",", iblk, ")"
                    call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
                end if
            ! complex block
            case (QCMPLXMAT)
                ierr = QcMatGetValues_f(A=A,                         &
#if defined(QCMATRIX_ZERO_BASED)
                                        idx_block_row=iblk,          &
                                        idx_block_col=iblk,          &
#else
                                        idx_block_row=iblk+1,        &
                                        idx_block_col=iblk+1,        &
#endif
                                        idx_first_row=idx_first_row, &
                                        num_row_get=num_row_get,     &
                                        idx_first_col=idx_first_col, &
                                        num_col_get=num_col_get,     &
                                        values_real=values_real,     &
                                        values_imag=values_imag)
                if (ierr==QSUCCESS) then
                    jval = -dim_mat
                    do ival = 1, dim_mat
                        jval = jval+dim_mat+1
                        cf_trace(jblk) = cf_trace(jblk)+values_real(jval)
                        cf_trace(kblk) = cf_trace(kblk)+values_imag(jval)
                    end do
                else
                    write(io_log,100) "block (", iblk, ",", iblk, ")"
                    call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
                end if
            end select
        end do
        deallocate(data_type)
        deallocate(values_real)
        deallocate(values_imag)
        ! gets the trace by calling QcMatGetTrace()
        allocate(trace(2*dim_block), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ierr = QcMatGetTrace_f(A, dim_block, trace)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        do iblk = 1, 2*dim_block
            if (abs(trace(iblk)-cf_trace(iblk))>CF_THRESHOLD) then
                jblk = (iblk+1)/2
                if (mod(iblk,2)==1) then
                    write(io_log,100) "real part of block (", &
                                      jblk, ",", jblk, ")",   &
                                      trace(iblk), cf_trace(iblk), CF_THRESHOLD
                else
                    write(io_log,100) "imaginary part of block (", &
                                      jblk, ",", jblk, ")",        &
                                      trace(iblk), cf_trace(iblk), &
                                      CF_THRESHOLD
                end if
                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end if
        end do
        write(io_log,100) "QcMatGetTrace_f(A) passed ..."
        deallocate(cf_trace)
        deallocate(trace)
        return
100     format("test_f_QcMatGetTrace>> ",A,I4,A,I4,A,3Es16.8)
    end subroutine test_f_QcMatGetTrace

#undef QCMATRIX_F_TEST_SRC
