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
!!  This file tests the function QcMatGetExternalMat_f().
!!
!!  2014-03-13, Bin Gao:
!!  * first version

! checks if the external Fortran module and matrix type are given
#include "adapter/check_f_mod.h"

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatGetExternalMat.F90"

    !% \brief tests the function QcMatGetExternalMat_f()
    !  \author Bin Gao
    !  \date 2014-03-13
    !  \param[QcMat:type]{inout} A the matrix
    !% \param[integer]{in} io_log IO of the logfile
    subroutine test_f_QcMatGetExternalMat(A, io_log)
        use LANG_F_MODULE, only: LANG_F_MATRIX,      &
#if defined(ADAPTER_BLOCK_REAL)
                                 Matrix_IsAssembled, &
#endif
                                 Matrix_Scale
        use qcmatrix_f, only: QINT,                   &
                              QREAL,                  &
                              QcMat,                  &
                              QcMatIsAssembled_f,     &
                              QcMatCreate_f,          &
                              QcMatDuplicate_f,       &
#if defined(ADAPTER_CMPLX_MAT) || defined(ADAPTER_REAL_MAT)
                              QcMatGetDimBlock_f,     &
                              QcMatGetDataType_f,     &
#endif
                              QcMatGetExternalMat_f,  &
                              QNULLMAT,               &
                              QREALMAT,               &
                              QIMAGMAT,               &
                              QCMPLXMAT,              &
                              QcMatScale_f,           &
                              QcMatIsEqual_f,         &
                              QcMatDestroy_f,         &
                              COPY_PATTERN_AND_VALUE, &
#if defined(QCMATRIX_ENABLE_VIEW)
                              QcMatWrite_f,           &
                              ASCII_VIEW,             &
#endif
                              QSUCCESS
        implicit none
        type(QcMat), intent(inout) :: A
        integer(kind=4), intent(in) :: io_log
        type(QcMat) B                                           !duplication of the matrix A
        logical(kind=4) assembled                               !indicates if the matrix is assembled or not
#if defined(ADAPTER_CMPLX_MAT) || defined(ADAPTER_REAL_MAT)
        integer(kind=QINT) dim_block                            !dimension of the blocks
        integer(kind=QINT) num_blocks                           !number of blocks, as \var{dim_block}*\var{dim_block}
        integer(kind=QINT), allocatable :: idx_block_row(:)     !row indices of the blocks
        integer(kind=QINT), allocatable :: idx_block_col(:)     !column indices of the blocks
        integer(kind=QINT), allocatable :: data_type(:)         !data types of the blocks
        integer(kind=QINT) iblk, jblk, kblk
#endif
        type(LANG_F_MATRIX), pointer :: A_ext                   !external matrix
        real(kind=QREAL), parameter :: scal_number = 0.1_QREAL  !scaling number
        logical(kind=4) :: cf_values = .true.                   !comparing values
        logical(kind=4) is_equal                                !indicates if two matrices are equal (pattern and values)
        integer(kind=4) ierr                                    !error information
        ! checks if the matrix is assembled
        ierr = QcMatIsAssembled_f(A, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix A is not assembled ..."
                write(io_log,100) "QcMatGetExternalMat_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! duplicates the matrix A
        ierr = QcMatCreate_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatCreate_f(B) passed ..."
        ierr = QcMatDuplicate_f(A, COPY_PATTERN_AND_VALUE, B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatDuplicate_f(A,COPY_PATTERN_AND_VALUE) passed ..."
        ! scales the matrix B
        ierr = QcMatScale_f((/scal_number,0.0_QREAL/), B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatScale_f(B) passed ..."
#if defined(ADAPTER_CMPLX_MAT) || defined(ADAPTER_REAL_MAT)
        ! gets the dimension of blocks
        ierr = QcMatGetDimBlock_f(A, dim_block)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimBlock_f(A) passed ..."
        num_blocks = dim_block*dim_block
        ! gets the data type of each block
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
                idx_block_row(kblk) = iblk
                idx_block_col(kblk) = jblk
            end do
        end do
        ierr = QcMatGetDataType_f(A, num_blocks, idx_block_row, idx_block_col, data_type)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
#if defined(ADAPTER_BLOCK_CMPLX)
        ! tests if we could get the external matrix
        ierr = QcMatGetExternalMat_f(A=A, A_ext=A_ext)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetExternalMat_f(A) passed ..."
        ! tests Matrix_Scale() in the external matrix library
        call Matrix_Scale((/scal_number,0.0_QREAL/), A_ext)
        ! A_ext could be nullified after its use
        nullify(A_ext)
#elif defined(ADAPTER_BLOCK_REAL)
        ! gets the external matrix (real part)
        ierr = QcMatGetExternalMat_f(A=A,                &
                                     data_type=QREALMAT, &
                                     A_ext=A_ext)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetExternalMat_f(A) passed ..."
        ! tests Matrix_Scale() in the external matrix library
        call Matrix_IsAssembled(A_ext, assembled)
        if (assembled) then
            call Matrix_Scale(scal_number, A_ext)
        end if
        ! A_ext could be nullified after its use
        nullify(A_ext)
        ! gets the external matrix (imaginar part)
        ierr = QcMatGetExternalMat_f(A=A,                &
                                     data_type=QIMAGMAT, &
                                     A_ext=A_ext)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetExternalMat_f(A) passed ..."
        ! tests Matrix_Scale() in the external matrix library
        call Matrix_IsAssembled(A_ext, assembled)
        if (assembled) then
            call Matrix_Scale(scal_number, A_ext)
        end if
        ! A_ext could be nullified after its use
        nullify(A_ext)
#elif defined(ADAPTER_CMPLX_MAT)
        do iblk = 1, num_blocks
            if (data_type(iblk)/=QNULLMAT) then
                ! gets the external complex matrix for this block
                ierr = QcMatGetExternalMat_f(A=A,                               &
                                             idx_block_row=idx_block_row(iblk), &
                                             idx_block_col=idx_block_col(iblk), &
                                             A_ext=A_ext)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(io_log,100) "QcMatGetExternalMat_f(A) passed ..."
                ! tests Matrix_Scale() in the external matrix library
                call Matrix_Scale((/scal_number,0.0_QREAL/), A_ext)
                ! A_ext could be nullified after its use
                nullify(A_ext)
            end if
        end do
#elif defined(ADAPTER_REAL_MAT)
        do iblk = 1, num_blocks
            ! real part of the block
            if (data_type(iblk)==QREALMAT .or. data_type(iblk)==QCMPLXMAT) then
                ierr = QcMatGetExternalMat_f(A=A,                               &
                                             idx_block_row=idx_block_row(iblk), &
                                             idx_block_col=idx_block_col(iblk), &
                                             data_type=QREALMAT,                &
                                             A_ext=A_ext)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(io_log,100) "QcMatGetExternalMat_f(A) passed ..."
                ! tests Matrix_Scale() in the external matrix library
                call Matrix_Scale(scal_number, A_ext)
                nullify(A_ext)
            end if
            ! imaginary part of the block
            if (data_type(iblk)==QIMAGMAT .or. data_type(iblk)==QCMPLXMAT) then
                ierr = QcMatGetExternalMat_f(A=A,                               &
                                             idx_block_row=idx_block_row(iblk), &
                                             idx_block_col=idx_block_col(iblk), &
                                             data_type=QIMAGMAT,                &
                                             A_ext=A_ext)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                write(io_log,100) "QcMatGetExternalMat_f(A) passed ..."
                call Matrix_Scale(scal_number, A_ext)
                ! A_ext could be nullified after its use
                nullify(A_ext)
            end if
        end do
#endif
#if defined(ADAPTER_CMPLX_MAT) || defined(ADAPTER_REAL_MAT)
        ! cleans
        deallocate(idx_block_row)
        deallocate(idx_block_col)
        deallocate(data_type)
#endif
        ! tests if we really scaled the matrix A by using the subroutine in the external matrix library
        ierr = QcMatIsEqual_f(A, B, cf_values, is_equal)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        if (is_equal) then
            write(io_log,100) "QcMatIsEqual_f(A, B) passed ..."
        else
            ! dumps results to check
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite_f(A, "QcMatGetExternalMat_A", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ierr = QcMatWrite_f(B, "QcMatGetExternalMat_B", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! frees the space taken by the matrix B
        ierr = QcMatDestroy_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        return
100     format("test_f_QcMatGetExternalMat>> ",A)
    end subroutine test_f_QcMatGetExternalMat

#undef QCMATRIX_F_TEST_SRC
