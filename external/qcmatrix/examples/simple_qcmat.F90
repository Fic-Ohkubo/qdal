!!  QcMatrix: an abstract matrix library
!!  Copyright 2012-2014 Bin Gao
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
!!  This file is a simple example of using QcMatrix in Fortran code.
!!
!!  2014-03-13, Bin Gao:
!!  * first version

    program simple_qcmat
        use qcmatrix_f, only: QcMat,                  &
                              QcMatCreate_f,          &
                              QcMatBlockCreate_f,     &
                              QcMatSetDimMat_f,       &
                              QcMatSetDataType_f,     &
                              QcMatAssemble_f,        &
                              QcMatSetValues_f,       &
#if defined(QCMATRIX_ENABLE_VIEW)
                              QcMatWrite_f,           &
                              ASCII_VIEW,             &
#endif
                              QcMatDuplicate_f,       &
                              QcMatGEMM_f,            &
                              QcMatDestroy_f,         &
                              QSUCCESS,               &
                              QREALMAT,               &
                              QIMAGMAT,               &
                              QCMPLXMAT,              &
                              QINT,                   &
                              QREAL,                  &
                              COPY_PATTERN_AND_VALUE, &
                              MAT_NO_OPERATION,       &
                              MAT_HERM_TRANSPOSE
        implicit none
        type(QcMat) A, B, C                               !matrices
        integer(kind=QINT), parameter :: DIM_BLOCK = 4    !dimension of blocks
        integer(kind=QINT), parameter :: DIM_MAT = 4      !dimension of each block
        integer(kind=QINT), parameter :: NUM_BLOCKS = 10  !number of blocks to set data types
        ! let us have a matrix A like:
        !     | QREALMAT   QNULLMAT   QIMAGMAT   QNULLMAT  |
        ! A = | QNULLMAT   QIMAGMAT   QCMPLXMAT  QNULLMAT  |,
        !     | QCMPLXMAT  QCMPLXMAT  QNULLMAT   QNULLMAT  |
        !     | QCMPLXMAT  QREALMAT   QIMAGMAT   QCMPLXMAT |
        ! which could be set by calling QcMatSetDataType()
#if defined(QCMATRIX_ZERO_BASED)
        ! zero-based numbering
        integer(kind=QINT), parameter :: &                !indices of the block rows to set data types
            IDX_BLOCK_ROW(NUM_BLOCKS) = (/0, 0, 1, 1, 2, 2, 3, 3, 3, 3/)
        integer(kind=QINT), parameter :: &                !indices of the block columns to set data types
            IDX_BLOCK_COL(NUM_BLOCKS) = (/0, 2, 1, 2, 0, 1, 0, 1, 2, 3/)
#else
        ! one-based numbering
        integer(kind=QINT), parameter :: &                !indices of the block rows to set data types
            IDX_BLOCK_ROW(NUM_BLOCKS) = (/1, 1, 2, 2, 3, 3, 4, 4, 4, 4/)                 
        integer(kind=QINT), parameter :: &                !indices of the block columns to set data types
            IDX_BLOCK_COL(NUM_BLOCKS) = (/1, 3, 2, 3, 1, 2, 1, 2, 3, 4/)
#endif
        integer(kind=QINT), parameter :: &                !data types of the blocks
            DATA_TYPE(NUM_BLOCKS) = (/QREALMAT, QIMAGMAT,   &
                                      QIMAGMAT, QCMPLXMAT,  &
                                      QCMPLXMAT, QCMPLXMAT, &
                                      QCMPLXMAT, QREALMAT, QIMAGMAT, QCMPLXMAT/)
        real(kind=QREAL) values(DIM_MAT*DIM_MAT)          !values of each block, must be rank-1 array
        integer(kind=QINT) iblk                           !incremental recorder over blocks
        integer(kind=4) ierr                              !error information
        ! creates the matrix A
        ierr = QcMatCreate_f(A=A)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatCreate_f(A)"
        end if
        ! creates the blocks of the matrix A
        ierr = QcMatBlockCreate_f(A=A, dim_block=DIM_BLOCK)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatBlockCreate_f(A)"
        end if
        ! sets the dimension of each block
        ierr = QcMatSetDimMat_f(A=A, num_row=DIM_MAT, num_col=DIM_MAT)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatSetDimMat_f(A)"
        end if
        ! sets the data type of each block of the matrix A
        ierr = QcMatSetDataType_f(A=A ,                        &
                                  num_blocks=NUM_BLOCKS,       &
                                  idx_block_row=IDX_BLOCK_ROW, &
                                  idx_block_col=IDX_BLOCK_COL, &
                                  data_type=DATA_TYPE)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatSetDataType_f(A)"
        end if
        ! assembles the matrix A
        ierr = QcMatAssemble_f(A=A)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatAssemble_f(A)"
        end if
        ! sets the values of the matrix A
        values = 1.0_QREAL
        do iblk = 1, NUM_BLOCKS
            select case (DATA_TYPE(iblk))
                case (QREALMAT)
                    ierr = QcMatSetValues_f(A=A,                               &
                                            idx_block_row=IDX_BLOCK_ROW(iblk), &
                                            idx_block_col=IDX_BLOCK_COL(iblk), &
#if defined(QCMATRIX_ZERO_BASED)
                                            idx_first_row=0,                   &
                                            num_row_set=DIM_MAT,               &
                                            idx_first_col=0,                   &
                                            num_col_set=DIM_MAT,               &
#else
                                            idx_first_row=1,                   &
                                            num_row_set=DIM_MAT,               &
                                            idx_first_col=1,                   &
                                            num_col_set=DIM_MAT,               &
#endif
                                            values_real=values)
                case (QIMAGMAT)
                    ierr = QcMatSetValues_f(A=A,                               &
                                            idx_block_row=IDX_BLOCK_ROW(iblk), &
                                            idx_block_col=IDX_BLOCK_COL(iblk), &
#if defined(QCMATRIX_ZERO_BASED)
                                            idx_first_row=0,                   &
                                            num_row_set=DIM_MAT,               &
                                            idx_first_col=0,                   &
                                            num_col_set=DIM_MAT,               &
#else
                                            idx_first_row=1,                   &
                                            num_row_set=DIM_MAT,               &
                                            idx_first_col=1,                   &
                                            num_col_set=DIM_MAT,               &
#endif
                                            values_imag=values)
                case (QCMPLXMAT)
                    ierr = QcMatSetValues_f(A=A,                               &
                                            idx_block_row=IDX_BLOCK_ROW(iblk), &
                                            idx_block_col=IDX_BLOCK_COL(iblk), &
#if defined(QCMATRIX_ZERO_BASED)
                                            idx_first_row=0,                   &
                                            num_row_set=DIM_MAT,               &
                                            idx_first_col=0,                   &
                                            num_col_set=DIM_MAT,               &
#else
                                            idx_first_row=1,                   &
                                            num_row_set=DIM_MAT,               &
                                            idx_first_col=1,                   &
                                            num_col_set=DIM_MAT,               &
#endif
                                            values_real=values,                &
                                            values_imag=values)
                case default
                    write(6,*) "block (", IDX_BLOCK_ROW(iblk), ",", &
                               IDX_BLOCK_COL(iblk), ")"
                    write(6,*) "data type: ", DATA_TYPE(iblk)
                    stop "unknown data type"
            end select
            if (ierr/=QSUCCESS) then
                write(6,*) "block (", IDX_BLOCK_ROW(iblk), ",", IDX_BLOCK_COL(iblk), ")"
                stop "failed to call QcMatSetValues_f(A)"
            end if
        end do
#if defined(QCMATRIX_ENABLE_VIEW)
        ierr = QcMatWrite_f(A=A, mat_label="matrix_A", view_option=ASCII_VIEW)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatWrite_f(A)"
        end if
#endif
        ! creates the matrix B
        ierr = QcMatCreate_f(A=B)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatCreate_f(B)"
        end if
        ! copies the matrix A (including the numerical values) to B
        ierr = QcMatDuplicate_f(A=A, duplicate_option=COPY_PATTERN_AND_VALUE, B=B)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatDuplicate(A, COPY_PATTERN_AND_VALUE, B)"
        end if
#if defined(QCMATRIX_ENABLE_VIEW)
        ierr = QcMatWrite_f(A=B, mat_label="matrix_B", view_option=ASCII_VIEW)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatWrite_f(B)"
        end if
#endif
        ! creates the matrix C
        ierr = QcMatCreate_f(A=C)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatCreate_f(C)"
        end if
        ! performs matrix-matrix multiplication: C = A*B^{\dagger}
        ierr = QcMatGEMM_f(op_A=MAT_NO_OPERATION,         &
                           op_B=MAT_HERM_TRANSPOSE,       &
                           alpha=(/1.0_QREAL,0.0_QREAL/), &
                           A=A,                           &
                           B=B,                           &
                           beta=(/0.0_QREAL,0.0_QREAL/),  &
                           C=C)
        if (ierr/=QSUCCESS) then
            stop "failed to calculate C = A*B^{\dagger}"
        end if
#if defined(QCMATRIX_ENABLE_VIEW)
        ierr = QcMatWrite_f(A=C, mat_label="matrix_C1", view_option=ASCII_VIEW)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatWrite_f(C1)"
        end if
#endif
        ! performs matrix-matrix multiplication: C = A*B^{\dagger}+i*C
        ierr = QcMatGEMM_f(op_A=MAT_NO_OPERATION,         &
                           op_B=MAT_HERM_TRANSPOSE,       &
                           alpha=(/1.0_QREAL,0.0_QREAL/), &
                           A=A,                           &
                           B=B,                           &
                           beta=(/0.0_QREAL,1.0_QREAL/),  &
                           C=C)
        if (ierr/=QSUCCESS) then
            stop "failed to calculate C = A*B^{\dagger}+i*C"
        end if
#if defined(QCMATRIX_ENABLE_VIEW)
        ierr = QcMatWrite_f(A=C, mat_label="matrix_C2", view_option=ASCII_VIEW)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatWrite_f(C2)"
        end if          
#endif
        ! cleans
        ierr = QcMatDestroy_f(A=A)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatDestroy_f(A)"
        end if
        ierr = QcMatDestroy_f(A=B)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatDestroy_f(B)"
        end if
        ierr = QcMatDestroy_f(A=C)
        if (ierr/=QSUCCESS) then
            stop "failed to call QcMatDestroy_f(C)"
        end if
    end program simple_qcmat
