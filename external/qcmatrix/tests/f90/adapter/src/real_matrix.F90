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
!!  This file is the module of the real matrix and its subroutines.
!!
!!  2011-12-13, Bin Gao:
!!  * first version

module real_matrix

    implicit none

! parameters defined in QcMatrix library that will be used in Fortran APIs
#include "api/qcmatrix_f_basic.h90"
#include "api/qcmatrix_f_mat_symmetry.h90"
#if defined(QCMATRIX_STORAGE_MODE)
#include "api/qcmatrix_f_mat_storage.h90"
#endif
#include "api/qcmatrix_f_mat_duplicate.h90"
#include "api/qcmatrix_f_mat_operations.h90"
#if defined(QCMATRIX_ENABLE_VIEW)
#include "api/qcmatrix_f_mat_view.h90"
#endif

#if defined(ADAPTER_REAL_MAT)
#define RealMat Matrix_
#define RealMatCreate Matrix_Create
#define RealMatSetSymType Matrix_SetSymType
#define RealMatSetDimMat Matrix_SetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define RealMatSetStorageMode Matrix_SetStorageMode
#endif
#define RealMatAssemble Matrix_Assemble
#define RealMatGetSymType Matrix_GetSymType
#define RealMatGetDimMat Matrix_GetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define RealMatGetStorageMode Matrix_GetStorageMode
#endif
#define RealMatIsAssembled Matrix_IsAssembled
#define RealMatSetValues Matrix_SetValues
#define RealMatAddValues Matrix_AddValues
#define RealMatGetValues Matrix_GetValues
#define RealMatDuplicate Matrix_Duplicate
#define RealMatZeroEntries Matrix_ZeroEntries
#define RealMatGetTrace Matrix_GetTrace
#define RealMatGetMatProdTrace Matrix_GetMatProdTrace
#define RealMatDestroy Matrix_Destroy
#if defined(QCMATRIX_ENABLE_VIEW)
#define RealMatWrite Matrix_Write
#define RealMatRead Matrix_Read
#endif
#define RealMatScale Matrix_Scale
#define RealMatAXPY Matrix_AXPY
#define RealMatTranspose Matrix_Transpose
#define RealMatGEMM Matrix_GEMM
#endif

    integer(kind=QINT), save, private :: RealMat_ID = 0

    ! real matrix
    type, public :: RealMat
        private
        integer(kind=QINT) :: Mat_ID = 0
        integer(kind=QINT) :: sym_type = QNONSYMMAT                !symmetry of the matrix
        integer(kind=QINT) :: dim_mat = 0                          !dimension of the matrix
#if defined(QCMATRIX_STORAGE_MODE)
        integer(kind=QINT) :: storage_mode = UNKNOWN_STORAGE_MODE  !storage mode is not implemented
#endif
        real(kind=QREAL), pointer :: values(:)                     !numerical values, in column major order
    end type RealMat

    ! unit number of standard output
    integer(kind=4), parameter, private :: STDOUT = 6

    public :: RealMatCreate
    public :: RealMatSetSymType
    public :: RealMatSetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
    public :: RealMatSetStorageMode
#endif
    public :: RealMatAssemble
    public :: RealMatGetSymType
    public :: RealMatGetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
    public :: RealMatGetStorageMode
#endif
    public :: RealMatIsAssembled
    public :: RealMatSetValues
    public :: RealMatAddValues
    public :: RealMatGetValues
    public :: RealMatDuplicate
    public :: RealMatZeroEntries
    public :: RealMatGetTrace
    public :: RealMatGetMatProdTrace
    public :: RealMatDestroy
#if defined(QCMATRIX_ENABLE_VIEW)
    public :: RealMatWrite
    public :: RealMatRead
#endif
    public :: RealMatScale
    public :: RealMatAXPY
    public :: RealMatTranspose
    public :: RealMatGEMM

    contains

    subroutine RealMatCreate(A)
        type(RealMat), intent(out) :: A
        A%Mat_ID = RealMat_ID
        RealMat_ID = RealMat_ID+1
        A%sym_type = QNONSYMMAT
        A%dim_mat = 0
#if defined(QCMATRIX_STORAGE_MODE)
        A%storage_mode = 0
#endif
        nullify(A%values)
    end subroutine RealMatCreate

    subroutine RealMatSetSymType(A, sym_type)
        type(RealMat), intent(inout) :: A
        integer(kind=QINT), intent(in) :: sym_type
        select case (sym_type)
            case (QSYMMAT)
                A%sym_type = QSYMMAT
            case (QANTISYMMAT)
                A%sym_type = QANTISYMMAT
            case (QNONSYMMAT)
                A%sym_type = QNONSYMMAT
            case default
                write(STDOUT,100) sym_type
                call error_exit(STDOUT,   &
                                __LINE__, &
                                "RealMatSetSymType@tests/f90/adapter/src/real_matrix.F90")
        end select
100     format("RealMatSetSymType>> invalid symmetry type", I8)
    end subroutine RealMatSetSymType

    subroutine RealMatSetDimMat(A, num_row, num_col)
        type(RealMat), intent(inout) :: A
        integer(kind=QINT), intent(in) :: num_row
        integer(kind=QINT), intent(in) :: num_col
        if (num_row<1) then
            write(STDOUT,100) num_row
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatSetDimMate@tests/f90/adapter/src/real_matrix.F90")
        else if (num_col<1) then
            write(STDOUT,100) num_col
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatSetDimMate@tests/f90/adapter/src/real_matrix.F90")
!FIXME: to implemented non-square matrix
        else if (num_col/=num_row) then
            write(STDOUT,100) num_row
            write(STDOUT,100) num_col
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatSetDimMate@tests/f90/adapter/src/real_matrix.F90")
        else
            A%dim_mat = num_row
        end if
100     format("RealMatSetDimMat>> invalid dimension", I16)
    end subroutine RealMatSetDimMat

#if defined(QCMATRIX_STORAGE_MODE)
    subroutine RealMatSetStorageMode(A, storage_mode)
        type(RealMat), intent(inout) :: A
        integer(kind=QINT), intent(in) :: storage_mode
        A%storage_mode = storage_mode
    end subroutine RealMatSetStorageMode
#endif

    subroutine RealMatAssemble(A)
        type(RealMat), intent(inout) :: A
        if (associated(A%values)) then
            ! we re-allocate the values if its size is not correct
            if (size(A%values)/=A%dim_mat*A%dim_mat) then
                deallocate(A%values)
                nullify(A%values)
                allocate(A%values(A%dim_mat*A%dim_mat))
            end if
        else
            allocate(A%values(A%dim_mat*A%dim_mat))
        end if
        A%values = 0
    end subroutine RealMatAssemble

    subroutine RealMatGetSymType(A, sym_type)
        type(RealMat), intent(in) :: A
        integer(kind=QINT), intent(out) :: sym_type
        sym_type = A%sym_type
    end subroutine RealMatGetSymType

    subroutine RealMatGetDimMat(A, num_row, num_col)
        type(RealMat), intent(in) :: A
        integer(kind=QINT), intent(out) :: num_row
        integer(kind=QINT), intent(out) :: num_col
        num_row = A%dim_mat
        num_col = A%dim_mat
    end subroutine RealMatGetDimMat

#if defined(QCMATRIX_STORAGE_MODE)
    subroutine RealMatGetStorageMode(A, storage_mode)
        type(RealMat), intent(in) :: A
        integer(kind=QINT), intent(out) :: storage_mode
        storage_mode = A%storage_mode
    end subroutine RealMatGetStorageMode
#endif

    subroutine RealMatIsAssembled(A, assembled)
        type(RealMat), intent(in) :: A
        logical(kind=4), intent(out) :: assembled
        assembled = associated(A%values)
    end subroutine RealMatIsAssembled

    subroutine RealMatSetValues(A,             &
                                idx_first_row, &
                                num_row_set,   &
                                idx_first_col, &
                                num_col_set,   &
                                values)
        type(RealMat), intent(inout) :: A
        integer(kind=QINT), intent(in) :: idx_first_row
        integer(kind=QINT), intent(in) :: num_row_set
        integer(kind=QINT), intent(in) :: idx_first_col
        integer(kind=QINT), intent(in) :: num_col_set
        real(kind=QREAL), intent(in) :: values(num_row_set*num_col_set)
        integer(kind=QINT) idx_last_row            !index of the last row
        integer(kind=QINT) irow, icol, jcol, ival  !incremental recorders
        idx_last_row = idx_first_row+num_row_set-1
        ! checks the range of rows
        if (idx_first_row<1 .or. idx_first_row>A%dim_mat .or. &
            num_row_set<1 .or. idx_last_row>A%dim_mat) then
            write(STDOUT,100) "index of the first row", idx_first_row
            write(STDOUT,100) "number of rows to set", num_row_set
            write(STDOUT,100) "number of rows of the matrix", A%dim_mat
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatSetValues@tests/f90/adapter/src/real_matrix.F90")
        end if
        ! checks the range of columns
        if (idx_first_col<1 .or. idx_first_col>A%dim_mat .or. &
            num_col_set<1 .or. idx_first_col+num_col_set>A%dim_mat+1) then
            write(STDOUT,100) "index of the first column", idx_first_col
            write(STDOUT,100) "number of columns to set", num_col_set
            write(STDOUT,100) "number of columns of the matrix", A%dim_mat
            call error_exit(STDOUT,   &   
                            __LINE__, &
                            "RealMatSetValues@tests/f90/adapter/src/real_matrix.F90")
        end if
        if (.not.associated(A%values)) then
            call RealMatAssemble(A)
        end if
        ival = 0
        jcol = (idx_first_col-1)*A%dim_mat
        do icol = 1, num_col_set
            do irow = idx_first_row, idx_last_row
                ival = ival+1
                A%values(jcol+irow) = values(ival)
            end do
            jcol = jcol+A%dim_mat
        end do
100     format("RealMatSetValues>> ",A,2I12)
    end subroutine RealMatSetValues

    subroutine RealMatAddValues(A,             &
                                idx_first_row, &
                                num_row_add,   &
                                idx_first_col, &
                                num_col_add,   &
                                values)
        type(RealMat), intent(inout) :: A
        integer(kind=QINT), intent(in) :: idx_first_row
        integer(kind=QINT), intent(in) :: num_row_add
        integer(kind=QINT), intent(in) :: idx_first_col
        integer(kind=QINT), intent(in) :: num_col_add
        real(kind=QREAL), intent(in) :: values(num_row_add*num_col_add)
        integer(kind=QINT) idx_last_row                  !index of the last row
        integer(kind=QINT) irow, icol, jcol, ival, jval  !incremental recorders
        idx_last_row = idx_first_row+num_row_add-1
        ! checks the range of rows
        if (idx_first_row<1 .or. idx_first_row>A%dim_mat .or. &
            num_row_add<1 .or. idx_last_row>A%dim_mat) then
            write(STDOUT,100) "index of the first row", idx_first_row
            write(STDOUT,100) "number of rows to add", num_row_add
            write(STDOUT,100) "number of rows of the matrix", A%dim_mat
            call error_exit(STDOUT,   &   
                            __LINE__, &
                            "RealMatAddValues@tests/f90/adapter/src/real_matrix.F90")
        end if
        ! checks the range of columns
        if (idx_first_col<1 .or. idx_first_col>A%dim_mat .or. &
            num_col_add<1 .or. idx_first_col+num_col_add>A%dim_mat+1) then
            write(STDOUT,100) "index of the first column", idx_first_col
            write(STDOUT,100) "number of columns to add", num_col_add
            write(STDOUT,100) "number of columns of the matrix", A%dim_mat
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatAddValues@tests/f90/adapter/src/real_matrix.F90")
        end if
        if (.not.associated(A%values)) then
            call RealMatAssemble(A)
            call RealMatZeroEntries(A)
        end if
        ival = 0
        jcol = (idx_first_col-1)*A%dim_mat
        do icol = 1, num_col_add
            do irow = idx_first_row, idx_last_row
                ival = ival+1
                jval = jcol+irow
                A%values(jval) = A%values(jval)+values(ival)
            end do
            jcol = jcol+A%dim_mat
        end do
100     format("RealMatAddValues>> ",A,2I12)
    end subroutine RealMatAddValues

    subroutine RealMatGetValues(A,             &
                                idx_first_row, &
                                num_row_get,   &
                                idx_first_col, &
                                num_col_get,   &
                                values)
        type(RealMat), intent(inout) :: A
        integer(kind=QINT), intent(in) :: idx_first_row
        integer(kind=QINT), intent(in) :: num_row_get
        integer(kind=QINT), intent(in) :: idx_first_col
        integer(kind=QINT), intent(in) :: num_col_get
        real(kind=QREAL), intent(out) :: values(num_row_get*num_col_get)
        integer(kind=QINT) idx_last_row            !index of the last row
        integer(kind=QINT) irow, icol, jcol, ival  !incremental recorders
        idx_last_row = idx_first_row+num_row_get-1
        ! checks the range of rows
        if (idx_first_row<1 .or. idx_first_row>A%dim_mat .or. &
            num_row_get<1 .or. idx_first_row>A%dim_mat) then
            write(STDOUT,100) "index of the first row", idx_first_row
            write(STDOUT,100) "number of rows to get", num_row_get
            write(STDOUT,100) "number of rows of the matrix", A%dim_mat
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatGetValues@tests/f90/adapter/src/real_matrix.F90")
        end if
        ! checks the range of columns
        if (idx_first_col<1 .or. idx_first_col>A%dim_mat .or. &
            num_col_get<1 .or. idx_first_col+num_col_get>A%dim_mat+1) then
            write(STDOUT,100) "index of the first column", idx_first_col
            write(STDOUT,100) "number of columns to get", num_col_get
            write(STDOUT,100) "number of columns of the matrix", A%dim_mat
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatGetValues@tests/f90/adapter/src/real_matrix.F90")
        end if
        if (.not.associated(A%values)) then
            values = 0
        else
            ival = 0
            jcol = (idx_first_col-1)*A%dim_mat
            do icol = 1, num_col_get
                do irow = idx_first_row, idx_last_row
                    ival = ival+1
                    values(ival) = A%values(jcol+irow)
                end do
                jcol = jcol+A%dim_mat
            end do
        end if
100     format("RealMatGetValues>> ",A,2I12)
    end subroutine RealMatGetValues

    subroutine RealMatDuplicate(A, duplicate_option, B)
        type(RealMat), intent(in) :: A
        integer(kind=QINT), intent(in) :: duplicate_option
        type(RealMat), intent(inout) :: B
        if (A%Mat_ID==B%Mat_ID) return
        B%sym_type = A%sym_type
        B%dim_mat = A%dim_mat
#if defined(QCMATRIX_STORAGE_MODE)
        B%storage_mode = A%storage_mode
#endif
        if (associated(B%values)) then
            deallocate(B%values)
            nullify(B%values)
        end if
        if (associated(A%values)) then
            allocate(B%values(B%dim_mat*B%dim_mat))
            if (duplicate_option==COPY_PATTERN_AND_VALUE) then
                call Real_BLAS_COPY(A%dim_mat*A%dim_mat, &
                                    A%values,            &
                                    1_QINT,              &
                                    B%values,            &
                                    1_QINT)
            end if
        end if
    end subroutine RealMatDuplicate

    subroutine RealMatZeroEntries(A)
        type(RealMat), intent(inout) :: A
        if (.not.associated(A%values)) then
            call RealMatAssemble(A)
        end if
        call Real_BLAS_SCAL(A%dim_mat*A%dim_mat, 0.0_QREAL, A%values, 1_QINT)
        A%sym_type = QSYMMAT
    end subroutine RealMatZeroEntries

    subroutine RealMatGetTrace(A, trace)
        type(RealMat), intent(in) :: A
        real(kind=QREAL), intent(out) :: trace
        integer(kind=QINT) icol
        if (associated(A%values)) then
            trace = 0.0_QREAL
            do icol = 0, A%dim_mat-1
                trace = trace+A%values(icol*(A%dim_mat+1)+1)
            end do
        else
            write(STDOUT,"(A)") "RealMatGetTrace>> A is not associated"
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatGetTrace@tests/f90/adapter/src/real_matrix.F90")
        end if
    end subroutine RealMatGetTrace

    subroutine RealMatGetMatProdTrace(A, B, op_B, trace)
        type(RealMat), intent(in) :: A
        type(RealMat), intent(in) :: B
        integer(kind=QINT), intent(in) :: op_B
        real(kind=QREAL), intent(out) :: trace
        real(kind=QREAL) tmp_dot
        integer(kind=QINT) irow
        if (.not.associated(A%values)) then
            write(STDOUT,100) "A is not associated"
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatGetMatProdTrace@tests/f90/adapter/src/real_matrix.F90")
        end if
        if (.not.associated(B%values)) then
            write(STDOUT,100) "B is not associated"
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatGetMatProdTrace@tests/f90/adapter/src/real_matrix.F90")
        end if
        if (A%dim_mat/=B%dim_mat) then
            write(STDOUT,100) "dimension of matrix A", A%dim_mat
            write(STDOUT,100) "dimension of matrix B", B%dim_mat
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatGetMatProdTrace@tests/f90/adapter/src/real_matrix.F90")
        end if
        ! trace of anti-symmetric matrix is zero
        if (A%sym_type*B%sym_type==QANTISYMMAT) then
            trace = 0.0_QREAL
        else
            select case (op_B)
                ! trace = \sum_{ij}A_{ij}*B_{ji}
                case (MAT_NO_OPERATION)
                    trace = 0.0_QREAL
                    do irow = 0, A%dim_mat-1
                        ! gets \sum{j}A_{ij}*B_{ji}
                        call Real_BLAS_DOT(A%dim_mat,                  &
                                           A%values(irow+1),           &
                                           A%dim_mat,                  &
                                           B%values(irow*B%dim_mat+1), &
                                           1_QINT,                     &
                                           tmp_dot)
                        trace = trace+tmp_dot
                    end do
                ! trace = \sum_{ij}A_{ij}*B_{ij}
                case (MAT_TRANSPOSE)
                    call Real_BLAS_DOT(A%dim_mat*A%dim_mat, &
                                       A%values,            &
                                       1_QINT,              &
                                       B%values,            &
                                       1_QINT,              &
                                       trace)
                case default
                    write(STDOUT,100) "invalid matrix operation on B", op_B
                    call error_exit(STDOUT,   &
                                    __LINE__, &
                                    "RealMatGetMatProdTrace@tests/f90/adapter/src/real_matrix.F90")
            end select
        end if
100     format("RealMatGetMatProdTrace>> ",A,I12)
    end subroutine RealMatGetMatProdTrace

    subroutine RealMatDestroy(A)
        type(RealMat), intent(inout) :: A
        A%sym_type = QNONSYMMAT
        A%dim_mat = 0
#if defined(QCMATRIX_STORAGE_MODE)
        A%storage_mode = 0        
#endif  
        if (associated(A%values)) then
            deallocate(A%values)
            nullify(A%values)
        end if
    end subroutine RealMatDestroy

#if defined(QCMATRIX_ENABLE_VIEW)
    subroutine RealMatWrite(A, mat_label, view_option)
        type(RealMat), intent(in) :: A
        character*(*), intent(in) :: mat_label
        integer(kind=QINT), intent(in) :: view_option
        integer(kind=4) :: IO_REAL_MAT = 7
        if (.not.associated(A%values)) then
            write(STDOUT,120) "A is not associated"
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatWrite@tests/f90/adapter/src/real_matrix.F90")
        end if
        select case (view_option)
            case (BINARY_VIEW)
                open(IO_REAL_MAT,          &
                     file=trim(mat_label), &
                     status="unknown",     &
                     form="unformatted")
                !write(IO_REAL_MAT) A%Mat_ID 
                write(IO_REAL_MAT) A%sym_type 
                write(IO_REAL_MAT) A%dim_mat 
#if defined(QCMATRIX_STORAGE_MODE)
                write(IO_REAL_MAT) A%storage_mode 
#endif
                write(IO_REAL_MAT) A%values
                close(IO_REAL_MAT)
            case (ASCII_VIEW)
                open(IO_REAL_MAT, file=trim(mat_label), status="unknown")
                !write(IO_REAL_MAT,100) A%Mat_ID
                write(IO_REAL_MAT,100) A%sym_type
                write(IO_REAL_MAT,100) A%dim_mat
#if defined(QCMATRIX_STORAGE_MODE)
                write(IO_REAL_MAT,100) A%storage_mode
#endif
                write(IO_REAL_MAT,110) A%values
                close(IO_REAL_MAT)
            case default
                write(STDOUT,120) "invalid view option", view_option
                call error_exit(STDOUT,   &
                                __LINE__, &
                                "RealMatWrite@tests/f90/adapter/src/real_matrix.F90")
        end select
100     format(I12)
110     format(6Es24.12)
120     format("RealMatWrite>> ",A,I8)
    end subroutine RealMatWrite

    subroutine RealMatRead(A, mat_label, view_option)
        type(RealMat), intent(inout) :: A
        character*(*), intent(in) :: mat_label
        integer(kind=QINT), intent(in) :: view_option
        integer(kind=4) :: IO_REAL_MAT = 7
        select case (view_option)
            case (BINARY_VIEW)
                open(IO_REAL_MAT,          &
                     file=trim(mat_label), &
                     status="old",         &
                     form="unformatted")
                !read(IO_REAL_MAT) A%Mat_ID 
                read(IO_REAL_MAT) A%sym_type
                read(IO_REAL_MAT) A%dim_mat
#if defined(QCMATRIX_STORAGE_MODE)
                read(IO_REAL_MAT) A%storage_mode
#endif
                if (associated(A%values)) then
                    if (size(A%values)/=A%dim_mat*A%dim_mat) then
                        deallocate(A%values)
                        nullify(A%values)
                        call RealMatAssemble(A)
                    end if
                else
                    call RealMatAssemble(A)
                end if
                read(IO_REAL_MAT) A%values
                close(IO_REAL_MAT)
            case (ASCII_VIEW)
                open(IO_REAL_MAT, file=trim(mat_label), status="old")
                !read(IO_REAL_MAT,*) A%Mat_ID 
                read(IO_REAL_MAT,*) A%sym_type
                read(IO_REAL_MAT,*) A%dim_mat
#if defined(QCMATRIX_STORAGE_MODE)
                read(IO_REAL_MAT,*) A%storage_mode
#endif
                if (associated(A%values)) then
                    if (size(A%values)/=A%dim_mat*A%dim_mat) then
                        deallocate(A%values)
                        nullify(A%values)
                        call RealMatAssemble(A)
                    end if
                else
                    call RealMatAssemble(A)
                end if
                read(IO_REAL_MAT,*) A%values
                close(IO_REAL_MAT)
            case default
                write(STDOUT,100) "invalid view option", view_option
                call error_exit(STDOUT,   &
                                __LINE__, &
                                "RealMatRead@tests/f90/adapter/src/real_matrix.F90")
        end select
100     format("RealMatRead>> ",A,I8)
    end subroutine RealMatRead
#endif

    subroutine RealMatScale(scal_number, A)
        real(kind=QREAL), intent(in) :: scal_number
        type(RealMat), intent(inout) :: A
        if (associated(A%values)) then
            call Real_BLAS_SCAL(A%dim_mat*A%dim_mat, scal_number, A%values, 1_QINT)
        else
            write(STDOUT,"(A)") "RealMatScale>> A is not associated"
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatScale@tests/f90/adapter/src/real_matrix.F90")
        end if
    end subroutine RealMatScale

    subroutine RealMatAXPY(multiplier, X, Y)
        real(kind=QREAL), intent(in) :: multiplier
        type(RealMat), intent(in) :: X
        type(RealMat), intent(inout) :: Y
        ! Y = (a+1)*Y
        if (X%Mat_ID==Y%Mat_ID) then
            call RealMatScale(multiplier+1.0_QREAL, Y)
        else
            if (associated(Y%values)) then
                if (X%dim_mat/=Y%dim_mat) then
                    write(STDOUT,100) "dimension of matrix X", X%dim_mat
                    write(STDOUT,100) "dimension of matrix Y", Y%dim_mat
                    call error_exit(STDOUT,   &
                                    __LINE__, &
                                    "RealMatAXPY@tests/f90/adapter/src/real_matrix.F90")
                end if
                call Real_BLAS_AXPY(X%dim_mat*X%dim_mat, &
                                    multiplier,          &
                                    X%values,            &
                                    1_QINT,              &
                                    Y%values,            &
                                    1_QINT)
                if (Y%sym_type/=X%sym_type) Y%sym_type = QNONSYMMAT
            ! returns Y = a*X
            else
                call RealMatDuplicate(X, COPY_PATTERN_AND_VALUE, Y)
                call RealMatScale(multiplier, Y)
            end if
        end if
100     format("RealMatAXPY>> ",A,I12)
    end subroutine RealMatAXPY

    subroutine RealMatTranspose(op_A, A, B)
        integer(kind=QINT), intent(in) :: op_A
        type(RealMat), intent(in) :: A
        type(RealMat), intent(inout) :: B
        integer(kind=QINT) irow, icol, jcol, ival, jval
        real(kind=QREAL) tmp_val
        select case (op_A)
            case (MAT_NO_OPERATION)
                call RealMatDuplicate(A, COPY_PATTERN_AND_VALUE, B)
            case (MAT_TRANSPOSE)
                ! out-of-place transpose
                if (A%Mat_ID/=B%Mat_ID) then
                    call RealMatDuplicate(A, COPY_PATTERN_ONLY, B)
                    jcol = -A%dim_mat
                    do icol = 0, A%dim_mat-1
                        jcol = jcol+A%dim_mat
                        do irow = 0, A%dim_mat-1
                            ! B(JI) = A(IJ)
                            B%values(irow*A%dim_mat+icol+1) = A%values(jcol+irow+1)
                        end do
                    end do
                ! in-place transpose
                else
                    jcol = -B%dim_mat
                    do icol = 0, B%dim_mat-1
                        jcol = jcol+B%dim_mat
                        do irow = 0, icol-1
                            ! B(JI) = B(IJ)
                            ival = jcol+irow+1
                            jval = irow*B%dim_mat+icol+1
                            tmp_val = B%values(jval)
                            B%values(jval) = B%values(ival)
                            B%values(ival) = tmp_val
                        end do
                    end do
                end if
            case default
                write(STDOUT,100) "invalid matrix operation", op_A
                call error_exit(STDOUT,   &
                                __LINE__, &
                                "RealMatTranspose@tests/f90/adapter/src/real_matrix.F90")
        end select
100     format("RealMatTranspose>> ",A,I8)
    end subroutine RealMatTranspose

    subroutine RealMatGEMM(op_A, op_B, alpha, A, B, beta, C)
        integer(kind=QINT), intent(in) :: op_A
        integer(kind=QINT), intent(in) :: op_B
        real(kind=QREAL), intent(in) :: alpha
        type(RealMat), intent(in) :: A
        type(RealMat), intent(in) :: B
        real(kind=QREAL), intent(in) :: beta
        type(RealMat), intent(inout) :: C
        character trans_A
        character trans_B
        if (C%Mat_ID==A%Mat_ID) then
            write(STDOUT,100) "ID of matrix C", C%Mat_ID
            write(STDOUT,100) "ID of matrix A", A%Mat_ID
            call error_exit(STDOUT,   &
                            __LINE__, &
                            "RealMatGEMM@tests/f90/adapter/src/real_matrix.F90")
        end if
        if (C%Mat_ID==B%Mat_ID) then
            write(STDOUT,100) "ID of matrix C", C%Mat_ID
            write(STDOUT,100) "ID of matrix B", B%Mat_ID
            call error_exit(STDOUT,   &   
                            __LINE__, &
                            "RealMatGEMM@tests/f90/adapter/src/real_matrix.F90")
        end if
        ! C = beta*C
        if (abs(alpha)<QZEROTHRSH) then
            call RealMatScale(beta, C)
        else
            select case (op_A)
                case (MAT_NO_OPERATION)
                    trans_A = "N"
                case (MAT_TRANSPOSE)
                    trans_A = "T"
                case default
                    write(STDOUT,100) "invalid matrix operation on A", op_A
                    call error_exit(STDOUT,   &
                                    __LINE__, &
                                    "RealMatGEMM@tests/f90/adapter/src/real_matrix.F90")
            end select
            select case (op_B)
                case (MAT_NO_OPERATION)
                    trans_B = "N"
                case (MAT_TRANSPOSE)
                    trans_B = "T"
                case default
                    write(STDOUT,100) "invalid matrix operation on B", op_B
                    call error_exit(STDOUT,   &
                                    __LINE__, &
                                    "RealMatGEMM@tests/f90/adapter/src/real_matrix.F90")
            end select
            if (.not.associated(A%values)) then
                write(STDOUT,100) "A is not associated"
                call error_exit(STDOUT,   &
                                __LINE__, &
                                "RealMatGEMM@tests/f90/adapter/src/real_matrix.F90")
            end if
            if (.not.associated(B%values)) then
                write(STDOUT,100) "B is not associated"
                call error_exit(STDOUT,   &   
                                __LINE__, &
                                "RealMatGEMM@tests/f90/adapter/src/real_matrix.F90")
            end if
            if (A%dim_mat/=B%dim_mat) then
                write(STDOUT,100) "dimension of matrix A", A%dim_mat
                write(STDOUT,100) "dimension of matrix B", B%dim_mat
                call error_exit(STDOUT,   &
                                __LINE__, &
                                "RealMatGEMM@tests/f90/adapter/src/real_matrix.F90")
            end if
            ! C = alpha*op(A)*op(B)
            if (abs(beta)<QZEROTHRSH .or. .not.associated(C%values)) then
                call RealMatDuplicate(A, COPY_PATTERN_ONLY, C)
            else
                if (A%dim_mat/=C%dim_mat) then
                    write(STDOUT,100) "dimension of matrix A", A%dim_mat
                    write(STDOUT,100) "dimension of matrix C", C%dim_mat
                    call error_exit(STDOUT,   &
                                    __LINE__, &
                                    "RealMatGEMM@tests/f90/adapter/src/real_matrix.F90")
                end if
            end if
            call Real_BLAS_GEMM(trans_A,   &
                                trans_B,   &
                                A%dim_mat, &
                                B%dim_mat, &
                                A%dim_mat, &
                                alpha,     &
                                A%values,  &
                                A%dim_mat, &
                                B%values,  &
                                B%dim_mat, &
                                beta,      &
                                C%values,  &
                                C%dim_mat)
            C%sym_type = QNONSYMMAT
        end if
100     format("RealMatGEMM>> ",A,I12)
    end subroutine RealMatGEMM

end module real_matrix
