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
!!  This file tests the function QcMatDuplicate_f().
!!
!!  2014-03-25, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatDuplicate.F90"

    !% \brief tests the function QcMatDuplicate_f()
    !  \author Bin Gao
    !  \date 2014-03-25
    !  \param[QcMat:type]{in} A the matrix
    !% \param[integer]{in} io_log IO of the logfile
    subroutine test_f_QcMatDuplicate(A, io_log)
        use qcmatrix_f, only: QINT,              &
                              QcMat,             &
                              QcMatCreate_f,     &
                              QcMatDuplicate_f,  &
                              QcMatIsEqual_f,    &
                              QcMatDestroy_f,    &
#if defined(QCMATRIX_ENABLE_VIEW)
                              QcMatWrite_f,      &
                              ASCII_VIEW,        &
#endif
                              COPY_PATTERN_ONLY, &
                              COPY_PATTERN_AND_VALUE
        implicit none
        type(QcMat), intent(in) :: A
        integer(kind=4), intent(in) :: io_log
        type(QcMat) B             !duplicated matrix
        logical(kind=4) cf_values !indicates if comparing values
        logical(kind=4) is_equal  !indicates if two matrices are equal (pattern and values)
        integer(kind=4) ierr      !error information
        ! creates the duplicated matrix first
        ierr = QcMatCreate_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! tests QcMatDuplicate_f() with the option COPY_PATTERN_ONLY
        ierr = QcMatDuplicate_f(A, COPY_PATTERN_ONLY, B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        cf_values = .false.
        ierr = QcMatIsEqual_f(A, B, cf_values, is_equal)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        if (is_equal) then
            write(io_log,100) "QcMatDuplicate_f(A, COPY_PATTERN_ONLY) passed ..." 
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! tests QcMatDuplicate_f() with the option COPY_PATTERN_AND_VALUE
        ierr = QcMatDuplicate_f(A, COPY_PATTERN_AND_VALUE, B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        cf_values = .true.
        ierr = QcMatIsEqual_f(A, B, cf_values, is_equal)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        if (is_equal) then
            write(io_log,100) "QcMatDuplicate_f(A, COPY_PATTERN_AND_VALUE) passed ..."
        else
            ! dumps results to check
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite_f(A, "QcMatDuplicate_A", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
            ierr = QcMatWrite_f(B, "QcMatDuplicate_B", ASCII_VIEW)
            call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#endif
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! frees the space taken by the matrix B
        ierr = QcMatDestroy_f(B)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        return
100     format("test_f_QcMatDuplicate>> ",A)
    end subroutine test_f_QcMatDuplicate

#undef QCMATRIX_F_TEST_SRC
