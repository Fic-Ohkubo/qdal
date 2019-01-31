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
!!  This file tests the functions QcMatWrite_f() and QcMatRead_f().
!!
!!  2014-03-25, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatView.F90"

    !% \brief tests the functions QcMatWrite_f() and QcMatRead_f()
    !  \author Bin Gao
    !  \date 2014-03-25
    !  \param[QcMat:type]{in} A the matrix
    !  \param[character]{in} A_label the label of the matrix
    !% \param[integer]{in} io_log IO of the logfile
    subroutine test_f_QcMatView(A, A_label, io_log)
        use qcmatrix_f, only: QINT,               &
                              QcMat,              &
                              QcMatIsAssembled_f, &
                              ASCII_VIEW,         &
                              BINARY_VIEW,        &
                              QcMatWrite_f,       &
                              QcMatRead_f,        &
                              QcMatCreate_f,      &
                              QcMatDestroy_f,     &
                              QcMatIsEqual_f,     &
                              QSUCCESS
        implicit none
        type(QcMat), intent(in) :: A
        character*(*), intent(in) :: A_label
        integer(kind=4), intent(in) :: io_log
        logical(kind=4) assembled              !indicates if the matrix is assembled or not
        type(QcMat) B                          !matrix for reading
        logical(kind=4) :: cf_values = .true.  !comparing values
        logical(kind=4) is_equal               !indicates if two matrices are equal (pattern and values)
        integer(kind=4) ierr                   !error information
        ierr = QcMatIsAssembled_f(A, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix A is not assembled ..."
                write(io_log,100) "QcMatWrite_f() and QcMatRead_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if

        !the values written using option ASCII_VIEW are normally not accurate,
        !also we normally save matrix using BINARY_VIEW, so we skip the test
        !using ASCII_VIEW
        !-! tests QcMatWrite_f() with option ASCII_VIEW
        !-ierr = QcMatWrite_f(A, A_label, ASCII_VIEW)
        !-call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        !-write(io_log,100) "QcMatWrite_f(A, ASCII_VIEW) passed ..."
        !-! creates the matrix B for reading
        !-ierr = QcMatCreate_f(B)
        !-call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        !-write(io_log,100) "QcMatCreate_f(B) passed ..."
        !-! tests QcMatRead_f() with option ASCII_VIEW
        !-ierr = QcMatRead_f(B,  A_label, ASCII_VIEW)
        !-call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        !-! tests if we read in exactly the same matrix A
        !-ierr = QcMatIsEqual_f(A, B, cf_values, is_equal)
        !-call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        !-if (is_equal) then
        !-    write(io_log,100) "QcMatRead_f(B, ASCII_VIEW) passed ..."
        !-else
        !-    ! dumps results to check
        !-    ierr = QcMatWrite_f(A, "QcMatView_A", ASCII_VIEW)
        !-    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        !-    ierr = QcMatWrite_f(B, "QcMatView_B", ASCII_VIEW)
        !-    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        !-    call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        !-end if
        !-! frees the space taken by the matrix B
        !-ierr = QcMatDestroy_f(B)
        !-call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        !-write(io_log,100) "QcMatDestroy_f(B) passed ..."

        ! tests QcMatWrite_f() with option BINARY_VIEW
        ierr = QcMatWrite_f(A, A_label, BINARY_VIEW)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatWrite_f(A, BINARY_VIEW) passed ..."
        ! creates the matrix B for reading
        ierr = QcMatCreate_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatCreate_f(B) passed ..."
        ! tests QcMatRead_f() with option BINARY_VIEW
        ierr = QcMatRead_f(B,  A_label, BINARY_VIEW)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! tests if we read in exactly the same matrix A
        ierr = QcMatIsEqual_f(A, B, cf_values, is_equal)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        if (is_equal) then
            write(io_log,100) "QcMatRead_f(B, BINARY_VIEW) passed ..."
        else
            ! dumps results to check
            ierr = QcMatWrite_f(A, "QcMatView_A", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ierr = QcMatWrite_f(B, "QcMatView_B", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! frees the space taken by the matrix B
        ierr = QcMatDestroy_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatDestroy_f(B) passed ..."
        return
100     format("test_f_QcMatView>> ",A)
    end subroutine test_f_QcMatView

#undef QCMATRIX_F_TEST_SRC
