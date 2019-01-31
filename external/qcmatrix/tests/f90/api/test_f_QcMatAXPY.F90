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
!!  This file tests the function QcMatAXPY_f().
!!
!!  2014-03-28, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatAXPY.F90"

    !% \brief tests the function QcMatAXPY_f()
    !  \author Bin Gao
    !  \date 2014-03-28
    !  \param[QcMat:type]{in} X the matrix
    !  \param[QcMat:type]{inout} Y the matrix
    !% \param[integer]{in} io_log IO of the logfile
    subroutine test_f_QcMatAXPY(X, Y, io_log)
        use qcmatrix_f, only: QINT,                &
                              QREAL,               &
                              QcMat,               &
                              QcMatIsAssembled_f,  &
                              QcMatGetDimBlock_f,  &
                              QcMatGetDimMat_f,    &
                              QcMatGetAllValues_f, &
                              QcMatAXPY_f,         &
                              QcMatCfArray_f,      &
                              QcMatDestroy_f,      &
#if defined(QCMATRIX_ENABLE_VIEW)
                              QcMatWrite_f,        &
                              ASCII_VIEW,          &
#endif
                              QSUCCESS
        implicit none
        type(QcMat), intent(in) :: X
        type(QcMat), intent(inout) :: Y
        integer(kind=4), intent(in) :: io_log
        logical(kind=4) assembled                   !indicates if the matrix is assembled or not
        integer(kind=QINT) dim_block                !dimension of blocks
        integer(kind=QINT) dim_mat                  !dimension of each block
        integer(kind=QINT) size_values              !number of elements in the matrix
        real(kind=QREAL), allocatable :: X_real(:)  !values of the real part of the matrix X
        real(kind=QREAL), allocatable :: X_imag(:)  !values of the imaginary part of the matrix X
        real(kind=QREAL), allocatable :: Y_real(:)  !values of the real part of the matrix Y
        real(kind=QREAL), allocatable :: Y_imag(:)  !values of the imaginary part of the matrix Y
        real(kind=QREAL) multiplier(2,4)            !the multiplier
        logical(kind=4) :: row_major = .false.      !values in column major order
        logical(kind=4) is_equal                    !indicates if the matrix and array have the same values
        integer(kind=QINT) idat                     !incremental recorder over the data types of the scaling number
        integer(kind=4) ierr                        !error information
        ! checks if the matrix is assembled
        ierr = QcMatIsAssembled_f(X, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix X is not assembled ..."
                write(io_log,100) "QcMatAXPY_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! gets the dimension of blocks
        ierr = QcMatGetDimBlock_f(X, dim_block)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimBlock_f(X) passed ..."
        ! gets the dimension of each block
        ierr = QcMatGetDimMat_f(X, dim_mat, dim_mat)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimMat_f(X) passed ..."
        size_values = dim_block*dim_block*dim_mat*dim_mat
        ! allocates memory for the elements of the matrices
        allocate(X_real(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(X_imag(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(Y_real(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(Y_imag(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! gets all the values of the matrices
        ierr = QcMatGetAllValues_f(X, row_major, size_values, X_real, X_imag)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ierr = QcMatGetAllValues_f(Y, row_major, size_values, Y_real, Y_imag)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! loops over different data types of the multiplier
        multiplier(:,1) = (/0.5_QREAL,0.0_QREAL/)  !real number
        multiplier(:,2) = (/0.0_QREAL,0.5_QREAL/)  !imaginary number
        multiplier(:,3) = (/0.5_QREAL,0.5_QREAL/)  !complex number
        multiplier(:,4) = (/0.0_QREAL,0.0_QREAL/)  !zero
        do idat = 1, 4
            ! performs Y = a*X+Y by BLAS routine
            ! (a_{R}+i*a_{I})*(X_{R}+i*X_{I})+(Y_{R}+i*Y_{I})
            ! = a_{R}*X_{R}-a_{I}*X_{I}+Y_{R}
            ! + i*(a_{R}*X_{I}+a_{I}*X_{R}+Y_{I})
            !
            ! a_{R}*X_{R}+Y_{R}
            call BLAS_AXPY(size_values, multiplier(1,idat), X_real, 1_QINT, Y_real, 1_QINT)
            ! -a_{I}*X_{I}+(a_{R}*X_{R}+Y_{R})
            call BLAS_AXPY(size_values, -multiplier(2,idat), X_imag, 1_QINT, Y_real, 1_QINT)
            ! a_{R}*X_{I}+Y_{I}
            call BLAS_AXPY(size_values, multiplier(1,idat), X_imag, 1_QINT, Y_imag, 1_QINT)
            ! a_{I}*X_{R}+(a_{R}*X_{I}+Y_{I})
            call BLAS_AXPY(size_values, multiplier(2,idat), X_real, 1_QINT, Y_imag, 1_QINT)
            ! calls QcMatAXPY()
            ierr = QcMatAXPY_f(multiplier(:,idat), X, Y)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ierr = QcMatCfArray_f(Y, row_major, size_values, Y_real, Y_imag, is_equal)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            if (is_equal) then
                write(io_log,100) "QcMatAXPY_f(X, Y) passed ..."
            else
                ! dumps results to check
                write(io_log,100) "multiplier", multiplier(:,idat)
#if defined(QCMATRIX_ENABLE_VIEW)
                ierr = QcMatWrite_f(X, "QcMatAXPY_X", ASCII_VIEW)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
                ierr = QcMatWrite_f(Y, "QcMatAXPY_Y", ASCII_VIEW)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
                write(io_log,100) "real part of Y from BLAS"
                write(io_log,"(4Es24.12)") Y_real
                write(io_log,100) "imaginary part of Y from BLAS"
                write(io_log,"(4Es24.12)") Y_imag
                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end if
        end do
        ! cleans
        deallocate(X_real)
        deallocate(X_imag)
        deallocate(Y_real)
        deallocate(Y_imag)
        return
100     format("test_f_QcMatAXPY>> ",A,2Es24.12)
    end subroutine test_f_QcMatAXPY

#undef QCMATRIX_F_TEST_SRC
