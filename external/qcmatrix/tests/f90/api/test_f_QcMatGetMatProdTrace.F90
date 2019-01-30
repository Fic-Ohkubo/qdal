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
!!  This file tests the function QcMatGetMatProdTrace_f().
!!
!!  2014-06-19, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatGetMatProdTrace.F90"

    !% \brief tests the function QcMatGetMatProdTrace_f()
    !  \author Bin Gao
    !  \date 2014-06-19
    !  \param[QcMat:type]{in} A the matrix
    subroutine test_f_QcMatGetMatProdTrace(A, B, io_log)
        use qcmatrix_f, only: QINT,                   &
                              QREAL,                  &
                              QcMat,                  &
                              QcMatIsAssembled_f,     &
                              QcMatCreate_f,          &
                              QcMatGetDimBlock_f,     &
                              QcMatGetTrace_f,        &
                              QcMatGetMatProdTrace_f, &
                              QcMatGEMM_f,            &
                              QcMatDestroy_f,         &
                              MAT_NO_OPERATION,       &
                              MAT_TRANSPOSE,          &
                              MAT_HERM_TRANSPOSE,     &
                              MAT_COMPLEX_CONJUGATE,  &
                              QSUCCESS,               &
                              QZEROTHRSH
        implicit none
        type(QcMat), intent(in) :: A
        type(QcMat), intent(in) :: B
        integer(kind=4), intent(in) :: io_log
        logical(kind=4) assembled                     !indicates if the matrix is assembled or not
        integer(kind=QINT) dim_block                  !dimension of blocks
        type(QcMat) C                                 !product matrix
        real(kind=QREAL), allocatable :: cf_trace(:)  !trace to compare with
        real(kind=QREAL), allocatable :: trace(:)     !trace from QcMatGetMatProdTrace()
        integer(kind=QINT), parameter :: all_mat_operations(4) = (/MAT_NO_OPERATION,   &  !all matrix operations
                                                                   MAT_TRANSPOSE,      &
                                                                   MAT_HERM_TRANSPOSE, &
                                                                   MAT_COMPLEX_CONJUGATE/)
        integer(kind=QINT) iop, iblk, jblk            !incremental recorders
        integer(kind=4) ierr                          !error information
        ! checks if the matrices A and B are assembled
        ierr = QcMatIsAssembled_f(A, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix A is not assembled ..."
                write(io_log,100) "QcMatGetMatProdTrace_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ierr = QcMatIsAssembled_f(B, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix B is not assembled ..."
                write(io_log,100) "QcMatGetMatProdTrace_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! creates the matrix C
        ierr = QcMatCreate_f(C)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! gets the dimension of blocks
        ierr = QcMatGetDimBlock_f(A, dim_block)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! allocates memory for the traces
        allocate(cf_trace(2*dim_block), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(trace(2*dim_block), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        do iop = 1, 4
            write(io_log,100) "operation on B", all_mat_operations(iop)
            ! gets the trace by calling QcMatGetMatProdTrace()
            ierr = QcMatGetMatProdTrace_f(A, B, all_mat_operations(iop), dim_block, trace)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ! gets the trace by calling QcMatGEMM() and QcMatGetTrace()
            ierr = QcMatGEMM_f(MAT_NO_OPERATION,        &
                               all_mat_operations(iop), &
                               (/1.0_QREAL,0.0_QREAL/), &
                               A,                       &
                               B,                       &
                               (/0.0_QREAL,0.0_QREAL/), &
                               C)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ierr = QcMatGetTrace_f(C, dim_block, cf_trace)
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
            write(io_log,100) "QcMatGetMatProdTrace_f(A,B) passed ..."
        end do
        ierr = QcMatDestroy_f(C)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        deallocate(cf_trace)
        deallocate(trace)
        return
100     format("test_f_QcMatGetMatProdTrace>> ",A,I4,A,I4,A,3Es16.8)
    end subroutine test_f_QcMatGetMatProdTrace

#undef QCMATRIX_F_TEST_SRC
