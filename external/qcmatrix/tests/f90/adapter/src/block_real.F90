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
!!  This file is the module of the square block real matrix and its subroutines.
!!
!!  2014-06-15, Bin Gao:
!!  * first version

module real_matrix
    implicit none

    integer, parameter, public :: SYMMAT = 1       !symmetric (Hermitian) matrix
    integer, parameter, public :: ANTISYMMAT = -1  !anti-symmetric (anti-Hermitian) matrix
    integer, parameter, public :: NONSYMMAT = 0    !non-symmetric (non-Hermitian) matrix

    ! data type for real numbers
#if defined(QCMATRIX_SINGLE_PRECISION)
    integer, public, parameter :: QREAL = kind(1.0)
#else
    integer, public, parameter :: QREAL = kind(1.0D0)
#endif

    ! real matrix
    type, public :: real_mat_t
        private
        integer :: sym_type = NONSYMMAT        !symmetry of the matrix
        integer :: num_row = 0                 !number of rows
        integer :: num_col = 0                 !number of columns
        real(QREAL), allocatable :: values(:)  !numerical values
    end type real_mat_t

#define QCMATRIX_ENABLE_VIEW
#define QCMATRIX_STORAGE_MODE

    public :: Matrix_Create
    public :: Matrix_Duplicate
#if defined(QCMATRIX_ENABLE_VIEW)
    public :: Matrix_Read
#endif
    public :: Matrix_Destroy
    public :: Matrix_SetSymType
    public :: Matrix_SetDataType
#if defined(QCMATRIX_STORAGE_MODE)
    public :: Matrix_SetStorageMode
#endif
    public :: Matrix_SetDimBlock
    public :: Matrix_SetBlocks
    public :: Matrix_SetDimMat
    public :: Matrix_Assemble
    public :: Matrix_GetSymType
    public :: Matrix_GetDataType
#if defined(QCMATRIX_STORAGE_MODE)
    public :: Matrix_GetStorageMode
#endif
    public :: Matrix_GetDimBlock
    public :: Matrix_GetNumBlocks
    public :: Matrix_GetIdxBlocks
    public :: Matrix_GetDimMat
    public :: Matrix_IsAssembled
    public :: Matrix_ZeroEntries
    public :: Matrix_GetTrace
#if defined(QCMATRIX_ENABLE_VIEW)
    public :: Matrix_Write
#endif
    public :: Matrix_Scale
    public :: Matrix_AXPY
    public :: Matrix_HermTranspose
    public :: Matrix_MatMult
    public :: Matrix_MatHermTransMult

    contains

    subroutine Matrix_Create(A)
        type(real_mat_t), intent(out) :: A
    end subroutine Matrix_Create

    subroutine Matrix_Duplicate(A, duplicate_option, B)
        type(real_mat_t), intent(in) :: A
        integer, intent(in) :: duplicate_option
        type(real_mat_t), intent(out) :: B
    end subroutine Matrix_Duplicate

#if defined(QCMATRIX_ENABLE_VIEW)
    subroutine Matrix_Read(A, mat_label, view_option)
        type(real_mat_t), intent(out) :: A
        character*(*), intent(in) :: mat_label
        integer, intent(in) :: view_option
    end subroutine Matrix_Read
#endif

    subroutine Matrix_Destroy(A)
        type(real_mat_t), intent(inout) :: A
    end subroutine Matrix_Destroy

    subroutine Matrix_SetSymType(A, sym_type)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: sym_type
    end subroutine Matrix_SetSymType

    subroutine Matrix_SetDataType(A, data_type)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: data_type
    end subroutine Matrix_SetDataType

#if defined(QCMATRIX_STORAGE_MODE)
    subroutine Matrix_SetStorageMode(A, storage_mode)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: storage_mode
    end subroutine Matrix_SetStorageMode
#endif

    subroutine Matrix_SetDimBlock(A, dim_block)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: dim_block
    end subroutine Matrix_SetDimBlock

    subroutine Matrix_SetBlocks(A, num_blocks, row_idx, col_idx)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: num_blocks
        integer, intent(in) :: row_idx(num_blocks)
        integer, intent(in) :: col_idx(num_blocks)
    end subroutine Matrix_SetBlocks

    subroutine Matrix_SetDimMat(A, dim_mat)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: dim_mat
    end subroutine Matrix_SetDimMat

    subroutine Matrix_Assemble(A)
        type(real_mat_t), intent(inout) :: A
    end subroutine Matrix_Assemble

    subroutine Matrix_GetSymType(A, sym_type)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: sym_type
    end subroutine Matrix_GetSymType

    subroutine Matrix_GetDataType(A, data_type)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: data_type
    end subroutine Matrix_GetDataType

#if defined(QCMATRIX_STORAGE_MODE)
    subroutine Matrix_GetStorageMode(A, storage_mode)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: storage_mode
    end subroutine Matrix_GetStorageMode
#endif

    subroutine Matrix_GetDimBlock(A, dim_block)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: dim_block
    end subroutine Matrix_GetDimBlock

    subroutine Matrix_GetNumBlocks(A, num_blocks)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: num_blocks
    end subroutine Matrix_GetNumBlocks

    subroutine Matrix_GetIdxBlocks(A, num_blocks, row_idx, col_idx)
        type(real_mat_t), intent(in) :: A
        integer, intent(in) :: num_blocks
        integer, intent(out) :: row_idx(num_blocks)
        integer, intent(out) :: col_idx(num_blocks)
    end subroutine Matrix_GetIdxBlocks

    subroutine Matrix_GetDimMat(A, dim_mat)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: dim_mat
    end subroutine Matrix_GetDimMat

    subroutine Matrix_IsAssembled(A, assembled)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: assembled
    end subroutine Matrix_IsAssembled

    subroutine Matrix_ZeroEntries(A)
        type(real_mat_t), intent(inout) :: A
    end subroutine Matrix_ZeroEntries

    subroutine Matrix_GetTrace(A, num_blocks, trace)
        type(real_mat_t), intent(in) :: A
        integer, intent(in) :: num_blocks
        real(QREAL), intent(out) :: trace(2*num_blocks)
    end subroutine Matrix_GetTrace

#if defined(QCMATRIX_ENABLE_VIEW)
    subroutine Matrix_Write(A, mat_label, view_option)
        type(real_mat_t), intent(in) :: A
        character*(*), intent(in) :: mat_label
        integer, intent(in) :: view_option
    end subroutine Matrix_Write
#endif

    subroutine Matrix_Scale(scal_number, A)
        real(QREAL), intent(in) :: scal_number(2)
        type(real_mat_t), intent(inout) :: A
    end subroutine Matrix_Scale

    subroutine Matrix_AXPY(multiplier, X, Y)
        real(QREAL), intent(in) :: multiplier(2)
        type(real_mat_t), intent(in) :: X
        type(real_mat_t), intent(inout) :: Y
    end subroutine Matrix_AXPY

    subroutine Matrix_HermTranspose(A, B)
        type(real_mat_t), intent(in) :: A
        type(real_mat_t), intent(inout) :: B
    end subroutine Matrix_HermTranspose

    subroutine Matrix_MatMult(A, B, C)
        type(real_mat_t), intent(in) :: A
        type(real_mat_t), intent(in) :: B
        type(real_mat_t), intent(inout) :: C
    end subroutine Matrix_MatMult

    subroutine Matrix_MatHermTransMult(A, B, C)
        type(real_mat_t), intent(in) :: A
        type(real_mat_t), intent(in) :: B
        type(real_mat_t), intent(inout) :: C
    end subroutine Matrix_MatHermTransMult

end module real_matrix
