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
!!  This file tests the function QcMatGetAdapterMat().
!!
!!  2014-03-13, Bin Gao:
!!  * first version

! checks if the external Fortran module and matrix type are given
#include "adapter/check_f_mod.h"

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

    !% \brief tests the function QcMatGetAdapterMat()
    !  \author Bin Gao
    !  \date 2014-03-13
    !% \param[QcMat:type]{inout} A the matrix
    subroutine test_f_QcMatGetAdapterMat(A)
        use LANG_F_MODULE, only: LANG_F_MATRIX,      &
#if defined(ADAPTER_BLOCK_REAL)
                                 Matrix_IsAssembled, &
#endif
                                 Matrix_Scale
        use qcmatrix, only: REALK,                  &
                           QcMat,                   &
                           QcMatIsAssembled,        &
                           QcMatCreate,             &
                           QcMatDuplicate,          &
#if defined(ADAPTER_CMPLX_MAT) || defined(ADAPTER_REAL_MAT)
                           QcMatGetDimBlock,        &
                           QcMatGetDataType,        &
#endif
                           QcMatGetAdapterMat,      &
                           QNULLMAT,               &
                           QREALMAT,               &
                           QIMAGMAT,               &
                           QCMPLXMAT,              &
                           QcMatScale,              &
                           QcMatIsEqual,            &
                           QcMatDestroy,            &
                           COPY_PATTERN_AND_VALUE, &
                           QSUCCESS,               &
                           QTRUE,                  &
                           QcMatWrite,              &
                           ASCII_VIEW
        implicit none
        type(QcMat), intent(inout) :: A
        type(QcMat) B                                       !duplication of the matrix A
        integer assembled                                  !indicates if the matrix is assembled or not
#if defined(ADAPTER_CMPLX_MAT) || defined(ADAPTER_REAL_MAT)
        integer dim_block                                  !dimension of the blocks
        integer num_blocks                                 !number of blocks, as \var{dim_block}*\var{dim_block}
        integer, allocatable :: idx_block_row(:)           !row indices of the blocks
        integer, allocatable :: idx_block_col(:)           !column indices of the blocks
        integer, allocatable :: data_type(:)               !data types of the blocks
        integer iblk, jblk, kblk
#endif
        type(LANG_F_MATRIX), pointer :: A_adapter          !adapter matrix for the external matrix library
        real(REALK), parameter :: scal_number = 0.1_REALK  !scaling number
        integer is_equal                                   !indicates if two matrices are equal (pattern and values)
        integer ierr                                       !error information
        ! checks if the matrix is assembled
        ierr = QcMatIsAssembled(A, assembled)
        if (ierr==QSUCCESS) then
            if (assembled/=QTRUE) then
                write(STDOUT,100) "matrix A is not assembled ..."
                write(STDOUT,100) "QcMatGetAdapterMat() will not be tested ..."
                return
            end if
        else
            call error_stop("test_f_QcMatScale",                  &
                            "failed to call QcMatIsAssembled(A)", &
                            ierr)
        end if
        ! duplicates the matrix A
        ierr = QcMatCreate(B)
        if (ierr==QSUCCESS) then
            write(STDOUT,100) "QcMatCreate(B) passed ..."
        else
            call error_stop("test_f_QcMatGetAdapterMat", "QcMatCreate(B) failed", ierr)
        end if
        ierr = QcMatDuplicate(A, COPY_PATTERN_AND_VALUE, B)
        if (ierr==QSUCCESS) then
            write(STDOUT,100) "QcMatDuplicate(A,COPY_PATTERN_AND_VALUE) passed ..."
        else
            call error_stop("test_f_QcMatGetAdapterMat",                                &
                            "failed to call QcMatDuplicate(A, COPY_PATTERN_AND_VALUE)", &
                            ierr)
        end if
        ! scales the matrix B
        ierr = QcMatScale((/scal_number,0.0_REALK/), B)
        if (ierr==QSUCCESS) then
            write(STDOUT,100) "QcMatScale(B) passed ..."
        else
            call error_stop("test_f_QcMatGetAdapterMat",    &
                            "failed to call QcMatScale(B)", &
                            ierr)
        end if
#if defined(ADAPTER_CMPLX_MAT) || defined(ADAPTER_REAL_MAT)
        ! gets the dimension of blocks
        ierr = QcMatGetDimBlock(A, dim_block)
        if (ierr==QSUCCESS) then
            write(STDOUT,100) "QcMatGetDimBlock(A) passed ..."
        else
            call error_stop("test_f_QcMatGetAdapterMat",          &
                            "failed to call QcMatGetDimBlock(A)", &
                            ierr)
        end if
        num_blocks = dim_block*dim_block
        ! gets the data type of each block
        allocate(idx_block_row(num_blocks), stat=ierr)
        if (ierr/=QSUCCESS) then
            call error_stop("test_f_QcMatGetAdapterMat",         &
                            "failed to allocate idx_block_row", &
                            ierr)
        end if
        allocate(idx_block_col(num_blocks), stat=ierr)
        if (ierr/=QSUCCESS) then
            call error_stop("test_f_QcMatGetAdapterMat",         &
                            "failed to allocate idx_block_col", &
                            ierr)
        end if
        allocate(data_type(num_blocks), stat=ierr)
        if (ierr/=QSUCCESS) then
            call error_stop("test_f_QcMatGetAdapterMat",     &
                            "failed to allocate data_type", &
                            ierr)
        end if
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
        ierr = QcMatGetDataType(A, num_blocks, idx_block_row, idx_block_col, data_type)
        if (ierr/=QSUCCESS) then
            call error_stop("test_f_QcMatGetAdapterMat",          &
                            "failed to call QcMatGetDataType(A)", &
                            ierr)
        end if
#endif
#if defined(ADAPTER_BLOCK_CMPLX)
        ! tests if we could get the adapter matrix
        ierr = QcMatGetAdapterMat(A=A, f_A=A_adapter)
        if (ierr/=QSUCCESS) then
            call error_stop("test_f_QcMatGetAdapterMat",            &
                            "failed to call QcMatGetAdapterMat(A)", &
                            ierr)
        else
            write(STDOUT,100) "QcMatGetAdapterMat(A) passed ..."
        end if
        ! tests Matrix_Scale() in the external matrix library
        call Matrix_Scale((/scal_number,0.0_REALK/), A_adapter)
        ! A_adapter could be nullified after its use
        nullify(A_adapter)
#elif defined(ADAPTER_BLOCK_REAL)
        ! gets the adapter matrix (real part)
        ierr = QcMatGetAdapterMat(A=A,                &
                                 data_type=QREALMAT, &
                                 f_A=A_adapter)
        if (ierr/=QSUCCESS) then
            call error_stop("test_f_QcMatGetAdapterMat",                     &
                            "failed to call QcMatGetAdapterMat(A,QREALMAT)", &
                            ierr)
        else
            write(STDOUT,100) "QcMatGetAdapterMat(A) passed ..."
        end if
        ! tests Matrix_Scale() in the external matrix library
        call Matrix_IsAssembled(A_adapter, assembled)
        if (assembled==QTRUE) then
            call Matrix_Scale(scal_number, A_adapter)
        end if
        ! A_adapter could be nullified after its use
        nullify(A_adapter)
        ! gets the adapter matrix (imaginar part)
        ierr = QcMatGetAdapterMat(A=A,                &
                                 data_type=QIMAGMAT, &
                                 f_A=A_adapter)
        if (ierr/=QSUCCESS) then
            call error_stop("test_f_QcMatGetAdapterMat",                     &
                            "failed to call QcMatGetAdapterMat(A,QIMAGMAT)", &
                            ierr)
        else
            write(STDOUT,100) "QcMatGetAdapterMat(A) passed ..."
        end if
        ! tests Matrix_Scale() in the external matrix library
        call Matrix_IsAssembled(A_adapter, assembled)
        if (assembled==QTRUE) then
            call Matrix_Scale(scal_number, A_adapter)
        end if
        ! A_adapter could be nullified after its use
        nullify(A_adapter)
#elif defined(ADAPTER_CMPLX_MAT)
        do iblk = 1, num_blocks
            if (data_type(iblk)/=QNULLMAT) then
                ! gets the adapter complex matrix for this block
                ierr = QcMatGetAdapterMat(A=A,                               &
                                         idx_block_row=idx_block_row(iblk), &
                                         idx_block_col=idx_block_col(iblk), &
                                         f_A=A_adapter)
                if (ierr/=QSUCCESS) then
                    call error_stop("test_f_QcMatGetAdapterMat",            &
                                    "failed to call QcMatGetAdapterMat(A)", &
                                    ierr)
                else
                    write(STDOUT,100) "QcMatGetAdapterMat(A) passed ..."
                end if
                ! tests Matrix_Scale() in the external matrix library
                call Matrix_Scale((/scal_number,0.0_REALK/), A_adapter)
                ! A_adapter could be nullified after its use
                nullify(A_adapter)
            end if
        end do
#elif defined(ADAPTER_REAL_MAT)
        do iblk = 1, num_blocks
            ! real part of the block
            if (data_type(iblk)==QREALMAT .or. data_type(iblk)==QCMPLXMAT) then
                ierr = QcMatGetAdapterMat(A=A,                               &
                                         idx_block_row=idx_block_row(iblk), &
                                         idx_block_col=idx_block_col(iblk), &
                                         data_type=QREALMAT,                &
                                         f_A=A_adapter)
                if (ierr/=QSUCCESS) then
                    call error_stop("test_f_QcMatGetAdapterMat",                     &
                                    "failed to call QcMatGetAdapterMat(A,QREALMAT)", &
                                    ierr)
                else
                    write(STDOUT,100) "QcMatGetAdapterMat(A) passed ..."
                end if
                ! tests Matrix_Scale() in the external matrix library
                call Matrix_Scale(scal_number, A_adapter)
                nullify(A_adapter)
            end if
            ! imaginary part of the block
            if (data_type(iblk)==QIMAGMAT .or. data_type(iblk)==QCMPLXMAT) then
                ierr = QcMatGetAdapterMat(A=A,                               &
                                         idx_block_row=idx_block_row(iblk), &
                                         idx_block_col=idx_block_col(iblk), &
                                         data_type=QIMAGMAT,                &
                                         f_A=A_adapter)
                if (ierr/=QSUCCESS) then
                    call error_stop("test_f_QcMatGetAdapterMat",                     &
                                    "failed to call QcMatGetAdapterMat(A,QIMAGMAT)", &
                                    ierr)
                else
                    write(STDOUT,100) "QcMatGetAdapterMat(A) passed ..."
                end if
                call Matrix_Scale(scal_number, A_adapter)
                ! A_adapter could be nullified after its use
                nullify(A_adapter)
            end if
        end do
#endif
#if defined(ADAPTER_CMPLX_MAT) || defined(ADAPTER_REAL_MAT)
        ! cleans
        deallocate(idx_block_row)
        deallocate(idx_block_col)
        deallocate(data_type)
#endif
        ! tests if we really scales the matrix A by using the subroutine in the external matrix library
        ierr = QcMatIsEqual(A, B, QTRUE, is_equal)
        if (ierr==QSUCCESS) then
            if (is_equal==QTRUE) then
                write(STDOUT,100) "QcMatIsEqual(A, B) passed ..."
            else
                ! dumps results to check
                ierr = QcMatWrite(A, "QcMatGetAdapterMat_A", ASCII_VIEW)
                if (ierr/=QSUCCESS) then
                    call error_stop("test_f_QcMatGetAdapterMat",    &
                                    "failed to call QcMatWrite(A)", &
                                    ierr)
                end if
                ierr = QcMatWrite(B, "QcMatGetAdapterMat_B", ASCII_VIEW)
                if (ierr/=QSUCCESS) then
                    call error_stop("test_f_QcMatGetAdapterMat",    &
                                    "failed to call QcMatWrite(B)", &
                                    ierr)
                end if
                call error_stop("test_f_QcMatGetAdapterMat", &
                                "QcMatIsEqual(A, B) failed", &
                                is_equal)
            end if
        else
            call error_stop("test_f_QcMatGetAdapterMat",         &
                            "failed to call QcMatIsEqual(A, B)", &
                            ierr)
        end if
        ! frees the space taken by the matrix B
        ierr = QcMatDestroy(B)
        if (ierr/=QSUCCESS) then
            call error_stop("test_f_QcMatGetAdapterMat", "QcMatDestroy(B) failed", ierr)
        else
            return
        end if
100     format("test_f_QcMatGetAdapterMat>> ",A)
    end subroutine test_f_QcMatGetAdapterMat
