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
!!  This file tests the function QcMatGEMM_f().
!!
!!  2014-03-28, Bin Gao:
!!  * first version

! parameters for test suite
#include "tests/qcmatrix_test_param.h"

#define QCMATRIX_F_TEST_SRC "tests/f90/api/test_f_QcMatGEMM.F90"

    !% \brief tests the function QcMatGEMM_f()
    !  \author Bin Gao
    !  \date 2014-03-28
    !  \param[QcMat:type]{in} A the matrix
    !  \param[QcMat:type]{in} B the matrix
    !% \param[integer]{in} io_log IO of the logfile
    subroutine test_f_QcMatGEMM(A, B, io_log)
        use qcmatrix_f, only: QINT,                   &
                              QREAL,                  &
                              QcMat,                  &
                              QSYMMAT,                &
                              QANTISYMMAT,            &
                              QNONSYMMAT,             &
                              QREALMAT,               &
                              QIMAGMAT,               &
                              QCMPLXMAT,              &
                              QcMatIsAssembled_f,     &
                              QcMatGetDimBlock_f,     &
                              QcMatGetDimMat_f,       &
                              QcMatGetAllValues_f,    &
                              QcMatGEMM_f,            &
                              QcMatCreate_f,          &
                              MAT_NO_OPERATION,       &
                              MAT_TRANSPOSE,          &
                              MAT_HERM_TRANSPOSE,     &
                              MAT_COMPLEX_CONJUGATE,  &
                              QcMatCfArray_f,         &
                              QcMatDestroy_f,         &
                              QcMatSetRandMat_f,      &
                              QcMatDuplicate_f,       &
                              COPY_PATTERN_AND_VALUE, &
                              QcMatScale_f,           &
#if defined(QCMATRIX_ENABLE_VIEW)
                              QcMatWrite_f,           &
                              ASCII_VIEW,             &
#endif
                              QSUCCESS
        implicit none
        type(QcMat), intent(in) :: A
        type(QcMat), intent(in) :: B
        integer(kind=4), intent(in) :: io_log
        type(QcMat) C                                  !product matrix
        type(QcMat) T                                  !temporary matrix
        logical(kind=4) assembled                      !indicates if the matrix is assembled or not
        integer(kind=QINT) dim_block                   !dimension of blocks
        integer(kind=QINT) dim_mat                     !dimension of each block
        integer(kind=QINT) size_values                 !number of elements in the matrix
        real(kind=QREAL), allocatable :: A_real(:)     !values of the real part of the matrix A
        real(kind=QREAL), allocatable :: A_imag(:)     !values of the imaginary part of the matrix A
        real(kind=QREAL), allocatable :: B_real(:)     !values of the real part of the matrix B
        real(kind=QREAL), allocatable :: B_imag(:)     !values of the imaginary part of the matrix B
        real(kind=QREAL), allocatable :: C_real(:)     !values of the real part of the matrix C
        real(kind=QREAL), allocatable :: C_imag(:)     !values of the imaginary part of the matrix C
        real(kind=QREAL), allocatable :: T_real(:)     !temporary values
        real(kind=QREAL), allocatable :: T_imag(:)     !temporary values
        integer(kind=QINT), parameter :: sym_type(3) = (/QSYMMAT,QANTISYMMAT,QNONSYMMAT/)  !all symmetry types
        integer(kind=QINT), parameter :: data_type(3) = (/QREALMAT,QIMAGMAT,QCMPLXMAT/)    !all data types
        integer(kind=QINT), parameter :: all_mat_operations(4) = (/MAT_NO_OPERATION,   &   !all matrix operations
                                                                   MAT_TRANSPOSE,      &
                                                                   MAT_HERM_TRANSPOSE, &
                                                                   MAT_COMPLEX_CONJUGATE/)
        character trans_A                              !operation on A, for BLAS routine
        character trans_B                              !operation on B, for BLAS routine
        real(kind=QREAL) alpha(2,2)                    !the scalar number
        real(kind=QREAL) beta(2)                       !the scalar number
        logical(kind=4) :: row_major = .false.         !values in column major order
        logical(kind=4) is_equal                       !indicates if the matrix and array have the same values
        integer(kind=QINT) iop, jop, isym, idat, jdat  !incremental recorders
        integer(kind=4) ierr                           !error information
        ! checks if the matrices A and B are assembled
        ierr = QcMatIsAssembled_f(A, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix A is not assembled ..."
                write(io_log,100) "QcMatGEMM_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ierr = QcMatIsAssembled_f(B, assembled)
        if (ierr==QSUCCESS) then
            if (.not. assembled) then
                write(io_log,100) "matrix B is not assembled ..."
                write(io_log,100) "QcMatGEMM_f() will not be tested ..."
                return
            end if
        else
            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
        end if
        ! gets the dimension of blocks
        ierr = QcMatGetDimBlock_f(A, dim_block)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimBlock_f(A) passed ..."
        ! gets the dimension of each block
        ierr = QcMatGetDimMat_f(A, dim_mat, dim_mat)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        write(io_log,100) "QcMatGetDimMat_f(A) passed ..."
        size_values = dim_block*dim_block*dim_mat*dim_mat
        ! allocates memory for the elements of the matrices
        allocate(A_real(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(A_imag(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(B_real(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(B_imag(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(C_real(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(C_imag(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        allocate(T_real(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! gets all the values of the matrices
        ierr = QcMatGetAllValues_f(A, row_major, size_values, A_real, A_imag)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ierr = QcMatGetAllValues_f(B, row_major, size_values, B_real, B_imag)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! set the alpha as complex number, since other types have been tested in CmplxMatScale() or CmplxMatAXPY()
        alpha(:,1) = (/0.01_QREAL,0.01_QREAL/)  !complex number
        alpha(:,2) = (/0.0_QREAL,0.0_QREAL/)    !zero
        ! we first test the matrix C which is not assembled
        beta = 0.0_QREAL
        do iop = 1, 4
            select case (all_mat_operations(iop))
                case (MAT_NO_OPERATION)
                    trans_A = "N"
                case (MAT_TRANSPOSE)
                    trans_A = "T"
                case (MAT_HERM_TRANSPOSE)
                    trans_A = "T"
                    A_imag = -A_imag
                case (MAT_COMPLEX_CONJUGATE)
                    trans_A = "N"
                    A_imag = -A_imag
                case default
                    call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
            end select
            do jop = 1, 4
                select case (all_mat_operations(jop))
                    case (MAT_NO_OPERATION)
                        trans_B = "N"
                    case (MAT_TRANSPOSE)
                        trans_B = "T"
                    case (MAT_HERM_TRANSPOSE)
                        trans_B = "T"
                        B_imag = -B_imag
                    case (MAT_COMPLEX_CONJUGATE)
                        trans_B = "N"
                        B_imag = -B_imag
                    case default
                        call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
                end select
                ! tests different alpha's
                do jdat = 1, 1  !we can not test alpha=0 here, since the matrix C is not assembled
                    ierr = QcMatCreate_f(C)
                    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                    write(io_log,100) "QcMatCreate_f(C) passed ..."
#endif
                    ierr = QcMatGEMM_f(all_mat_operations(iop), &
                                       all_mat_operations(jop), &
                                       alpha(:,jdat),           &
                                       A,                       &
                                       B,                       &
                                       beta,                    &
                                       C)
                    if (ierr==QSUCCESS) then
                        ! C = (a_{R}+i*a_{I})*(A_{R}+i*A_{I})*(B_{R}+i*B_{I})
                        !   = (a_{R}+i*a_{I})*(A_{R}*B_{R}-A_{I}*B_{I}+i*(A_{R}*B_{I}+A_{I}*B_{R}))
                        !   =   a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})
                        !   -   a_{I}*(A_{R}*B_{I}+A_{I}*B_{R})
                        !   + i*a_{R}*(A_{R}*B_{I}+A_{I}*B_{R})
                        !   + i*a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})
                        call BLAS_GEMM(trans_A,           &  !A_{R}*B_{R}
                                       trans_B,           &
                                       dim_block*dim_mat, &
                                       dim_block*dim_mat, &
                                       dim_block*dim_mat, &
                                       1.0_QREAL,         &
                                       A_real,            &
                                       dim_block*dim_mat, &
                                       B_real,            &
                                       dim_block*dim_mat, &
                                       0.0_QREAL,         &
                                       C_real,            &
                                       dim_block*dim_mat)
                        call BLAS_GEMM(trans_A,           &  !A_{I}*B_{I}
                                       trans_B,           &
                                       dim_block*dim_mat, &
                                       dim_block*dim_mat, &
                                       dim_block*dim_mat, &
                                       1.0_QREAL,         &
                                       A_imag,            &
                                       dim_block*dim_mat, &
                                       B_imag,            &
                                       dim_block*dim_mat, &
                                       0.0_QREAL,         &
                                       T_real,            &
                                       dim_block*dim_mat)
                        call BLAS_AXPY(size_values, &  !-A_{I}*B_{I}+A_{R}*B_{R}
                                       -1.0_QREAL,  &
                                       T_real,      &
                                       1_QINT,      &
                                       C_real,      &
                                       1_QINT)
                        call BLAS_GEMM(trans_A,           &  !A_{R}*B_{I}
                                       trans_B,           &
                                       dim_block*dim_mat, &
                                       dim_block*dim_mat, &
                                       dim_block*dim_mat, &
                                       1.0_QREAL,         &
                                       A_real,            &
                                       dim_block*dim_mat, &
                                       B_imag,            &
                                       dim_block*dim_mat, &
                                       0.0_QREAL,         &
                                       C_imag,            &
                                       dim_block*dim_mat)
                        call BLAS_GEMM(trans_A,           &  !A_{I}*B_{R}
                                       trans_B,           &
                                       dim_block*dim_mat, &
                                       dim_block*dim_mat, &
                                       dim_block*dim_mat, &
                                       1.0_QREAL,         &
                                       A_imag,            &
                                       dim_block*dim_mat, &
                                       B_real,            &
                                       dim_block*dim_mat, &
                                       0.0_QREAL,         &
                                       T_real,            &
                                       dim_block*dim_mat)
                        ! A_{I}*B_{R}+A_{R}*B_{I}
                        call BLAS_AXPY(size_values, &
                                       1.0_QREAL,   &
                                       T_real,      &
                                       1_QINT,      &
                                       C_imag,      &
                                       1_QINT)
                        call BLAS_COPY(size_values, C_real, 1_QINT, T_real, 1_QINT)
                        ! a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})
                        call BLAS_SCAL(size_values, alpha(1,jdat), C_real, 1_QINT)
                        ! -a_{I}*(A_{R}*B_{I}+A_{I}*B_{R})+a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})
                        call BLAS_AXPY(size_values,    &
                                       -alpha(2,jdat), &
                                       C_imag,         &
                                       1_QINT,         &
                                       C_real,         &
                                       1_QINT)
                        ! a_{R}*(A_{R}*B_{I}+A_{I}*B_{R})
                        call BLAS_SCAL(size_values, alpha(1,jdat), C_imag, 1_QINT)
                        ! a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})+a_{R}*(A_{R}*B_{I}+A_{I}*B_{R})
                        call BLAS_AXPY(size_values,   &
                                       alpha(2,jdat), &
                                       T_real,        &
                                       1_QINT,        &
                                       C_imag,        &
                                       1_QINT)
                        ierr = QcMatCfArray_f(C,           &
                                              row_major,   &
                                              size_values, &
                                              C_real,      &
                                              C_imag,      &
                                              is_equal)
                        call QErrorCheckCode(io_log,   &
                                             ierr,     &
                                             __LINE__, &
                                             QCMATRIX_F_TEST_SRC)
                        if (.not. is_equal) then
                            ! dumps results to check
                            write(io_log,100) "parameters",            &
                                              all_mat_operations(iop), &
                                              all_mat_operations(jop), &
                                              alpha(:,jdat)
#if defined(QCMATRIX_ENABLE_VIEW)
                            ierr = QcMatWrite_f(A, "QcMatGEMM_A", ASCII_VIEW)
                            call QErrorCheckCode(io_log,   &
                                                 ierr,     &
                                                 __LINE__, &
                                                 QCMATRIX_F_TEST_SRC)
                            ierr = QcMatWrite_f(B, "QcMatGEMM_B", ASCII_VIEW)
                            call QErrorCheckCode(io_log,   &
                                                 ierr,     &
                                                 __LINE__, &
                                                 QCMATRIX_F_TEST_SRC)
                            ierr = QcMatWrite_f(C, "QcMatGEMM_C", ASCII_VIEW)
                            call QErrorCheckCode(io_log,   &
                                                 ierr,     &
                                                 __LINE__, &
                                                 QCMATRIX_F_TEST_SRC)
#endif
                            write(io_log,100) "real part of C from BLAS"
                            write(io_log,"(4Es24.12)") C_real
                            write(io_log,100) "imaginary part of C from BLAS"
                            write(io_log,"(4Es24.12)") C_imag
                            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                        else
                            write(io_log,100) "QcMatGEMM_f(beta=0) passed", &
                                              all_mat_operations(iop),      &
                                              all_mat_operations(jop),      &
                                              alpha(:,jdat)
#endif
                        end if
                    else
                        write(io_log,100) "parameters",            &
                                          all_mat_operations(iop), &
                                          all_mat_operations(jop), &
                                          alpha(:,jdat)
#if defined(QCMATRIX_ENABLE_VIEW)
                        ierr = QcMatWrite_f(A, "QcMatGEMM_A", ASCII_VIEW)
                        call QErrorCheckCode(io_log,   &
                                             ierr,     &
                                             __LINE__, &
                                             QCMATRIX_F_TEST_SRC)
                        ierr = QcMatWrite_f(B, "QcMatGEMM_B", ASCII_VIEW)
                        call QErrorCheckCode(io_log,   &   
                                             ierr,     &   
                                             __LINE__, &
                                             QCMATRIX_F_TEST_SRC)
#endif
                        call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
                    end if
                    ! cleans the matrix C
                    ierr = QcMatDestroy_f(C)
                    call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                    write(io_log,100) "QcMatDestroy_f(C) passed ..."
#endif
                end do
                if (all_mat_operations(jop)==MAT_HERM_TRANSPOSE .or. &
                    all_mat_operations(jop)==MAT_COMPLEX_CONJUGATE) then
                    B_imag = -B_imag
                end if
            end do
            if (all_mat_operations(iop)==MAT_HERM_TRANSPOSE .or. &
                all_mat_operations(iop)==MAT_COMPLEX_CONJUGATE) then
                A_imag = -A_imag
            end if
        end do
#if !defined(DEBUG)
        write(io_log,100) "QcMatGEMM_f(beta=0) passed ..."
#endif
        ! allocates memory for temporary values
        allocate(T_imag(size_values), stat=ierr)
        call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
        ! set the beta as complex number, since it is only used in QcMatScale() and/or CmplxMatScale(),
        ! different beta's have been tested before
        beta = (/0.1_QREAL,0.1_QREAL/)
        ! tests different symmetry types (symmetric, anti-symmetric, non-symmetric) for matrix C
        do isym = 1, 3
            ! tests different data types (real, imaginary, complex) for matrix C
            do idat = 1, 3
                ! generates a random matrix C according to its symmetry and data types
                ierr = QcMatCreate_f(C)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                write(io_log,100) "QcMatCreate_f(C) passed ..."
#endif
                ierr = QcMatSetRandMat_f(C,               &
                                         sym_type(isym),  &
                                         data_type(idat), &
                                         dim_block,       &
                                         dim_mat,         &
                                         dim_mat)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                write(io_log,100) "QcMatSetRandMat_f(C) passed ..."
#endif
                ! temporary matrix as beta*C
                ierr = QcMatCreate_f(T)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                write(io_log,100) "QcMatCreate_f(T) passed ..."
#endif
                ! we first test the matrix C which is not assembled
                do iop = 1, 4
                    select case (all_mat_operations(iop))
                        case (MAT_NO_OPERATION)
                            trans_A = "N"
                        case (MAT_TRANSPOSE)
                            trans_A = "T"
                        case (MAT_HERM_TRANSPOSE)
                            trans_A = "T"
                            A_imag = -A_imag
                        case (MAT_COMPLEX_CONJUGATE)
                            trans_A = "N"
                            A_imag = -A_imag
                        case default
                            call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
                    end select
                    do jop = 1, 4
                        select case (all_mat_operations(jop))
                            case (MAT_NO_OPERATION)
                                trans_B = "N"
                            case (MAT_TRANSPOSE)
                                trans_B = "T"
                            case (MAT_HERM_TRANSPOSE)
                                trans_B = "T"
                                B_imag = -B_imag
                            case (MAT_COMPLEX_CONJUGATE)
                                trans_B = "N"
                                B_imag = -B_imag
                            case default
                                call QErrorExit(io_log, __LINE__, QCMATRIX_F_TEST_SRC)
                        end select
                        ! tests different alpha's
                        do jdat = 1, 2
                            ! T = beta*C
                            ierr = QcMatDuplicate_f(C, COPY_PATTERN_AND_VALUE, T)
                            call QErrorCheckCode(io_log,   &
                                                 ierr,     &
                                                 __LINE__, &
                                                 QCMATRIX_F_TEST_SRC)
                            ierr = QcMatScale_f(beta, T)
                            call QErrorCheckCode(io_log,   &
                                                 ierr,     &
                                                 __LINE__, &
                                                 QCMATRIX_F_TEST_SRC)
                            ! calcualtes C = alpha*A*B+beta*C using QcMatGEMM()
                            ierr = QcMatGEMM_f(all_mat_operations(iop), &
                                               all_mat_operations(jop), &
                                               alpha(:,jdat),           &
                                               A,                       &
                                               B,                       &
                                               beta,                    &
                                               C)
                            if (ierr==QSUCCESS) then
                                ! C = (a_{R}+i*a_{I})*(A_{R}+i*A_{I})*(B_{R}+i*B_{I})+beta*C
                                !   = (a_{R}+i*a_{I})*(A_{R}*B_{R}-A_{I}*B_{I}+i*(A_{R}*B_{I}+A_{I}*B_{R}))+beta*C
                                !   =   a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})
                                !   -   a_{I}*(A_{R}*B_{I}+A_{I}*B_{R})
                                !   +   (beta*C)_{R}
                                !   + i*a_{R}*(A_{R}*B_{I}+A_{I}*B_{R})
                                !   + i*a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})
                                !   + i*(beta*C)_{I}
                                call BLAS_GEMM(trans_A,           &  !A_{R}*B_{R}
                                               trans_B,           &
                                               dim_block*dim_mat, &
                                               dim_block*dim_mat, &
                                               dim_block*dim_mat, &
                                               1.0_QREAL,         &
                                               A_real,            &
                                               dim_block*dim_mat, &
                                               B_real,            &
                                               dim_block*dim_mat, &
                                               0.0_QREAL,         &
                                               T_real,            &
                                               dim_block*dim_mat)
                                call BLAS_GEMM(trans_A,           &  !A_{I}*B_{I}
                                               trans_B,           &
                                               dim_block*dim_mat, &
                                               dim_block*dim_mat, &
                                               dim_block*dim_mat, &
                                               1.0_QREAL,         &
                                               A_imag,            &
                                               dim_block*dim_mat, &
                                               B_imag,            &
                                               dim_block*dim_mat, &
                                               0.0_QREAL,         &
                                               C_real,            &
                                               dim_block*dim_mat)
                                call BLAS_AXPY(size_values, &  !-A_{I}*B_{I}+A_{R}*B_{R}
                                               -1.0_QREAL,  &
                                               C_real,      &
                                               1_QINT,      &
                                               T_real,      &
                                               1_QINT)
                                call BLAS_GEMM(trans_A,           &  !A_{R}*B_{I}
                                               trans_B,           &
                                               dim_block*dim_mat, &
                                               dim_block*dim_mat, &
                                               dim_block*dim_mat, &
                                               1.0_QREAL,         &
                                               A_real,            &
                                               dim_block*dim_mat, &
                                               B_imag,            &
                                               dim_block*dim_mat, &
                                               0.0_QREAL,         &
                                               T_imag,            &
                                               dim_block*dim_mat)
                                call BLAS_GEMM(trans_A,           &  !A_{I}*B_{R}
                                               trans_B,           &
                                               dim_block*dim_mat, &
                                               dim_block*dim_mat, &
                                               dim_block*dim_mat, &
                                               1.0_QREAL,         &
                                               A_imag,            &
                                               dim_block*dim_mat, &
                                               B_real,            &
                                               dim_block*dim_mat, &
                                               0.0_QREAL,         &
                                               C_imag,            &
                                               dim_block*dim_mat)
                                call BLAS_AXPY(size_values, &  !A_{I}*B_{R}+A_{R}*B_{I}
                                               1.0_QREAL,   &
                                               C_imag,      &
                                               1_QINT,      &
                                               T_imag,      &
                                               1_QINT)
                                ! gets the values of original C matrix scaled by beta
                                ierr = QcMatGetAllValues_f(T,           &
                                                           row_major,   &
                                                           size_values, &
                                                           C_real,      &
                                                           C_imag)
                                call QErrorCheckCode(io_log,   &
                                                     ierr,     &
                                                     __LINE__, &
                                                     QCMATRIX_F_TEST_SRC)
                                ! a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})+(beta*C)_{R}
                                call BLAS_AXPY(size_values,   &
                                               alpha(1,jdat), &
                                               T_real,        &
                                               1_QINT,        &
                                               C_real,        &
                                               1_QINT)
                                !-a_{I}*(A_{R}*B_{I}+A_{I}*B_{R})+a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})+(beta*C)_{R}
                                call BLAS_AXPY(size_values,    &
                                               -alpha(2,jdat), &
                                               T_imag,         &
                                               1_QINT,         &
                                               C_real,         &
                                               1_QINT)
                                ! a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})+(beta*C)_{I}
                                call BLAS_AXPY(size_values,   &
                                               alpha(2,jdat), &
                                               T_real,        &
                                               1_QINT,        &
                                               C_imag,        &
                                               1_QINT)
                                ! a_{R}*(A_{R}*B_{I}+A_{I}*B_{R})+a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})+(beta*C)_{I}
                                call BLAS_AXPY(size_values,   &
                                               alpha(1,jdat), &
                                               T_imag,        &
                                               1_QINT,        &
                                               C_imag,        &
                                               1_QINT)
                                ierr = QcMatCfArray_f(C,           &
                                                      row_major,   &
                                                      size_values, &
                                                      C_real,      &
                                                      C_imag,      &
                                                      is_equal)
                                call QErrorCheckCode(io_log,   &
                                                     ierr,     &
                                                     __LINE__, &
                                                     QCMATRIX_F_TEST_SRC)
                                if (.not. is_equal) then
                                    ! dumps results to check
                                    write(io_log,100) "parameters",            &
                                                      all_mat_operations(iop), &
                                                      all_mat_operations(jop), &
                                                      alpha(:,jdat),           &
                                                      beta,                    &
                                                      data_type(idat),         &
                                                      sym_type(isym)
#if defined(QCMATRIX_ENABLE_VIEW)
                                    ierr = QcMatWrite_f(A, "QcMatGEMM_A", ASCII_VIEW)
                                    call QErrorCheckCode(io_log,   &
                                                         ierr,     &
                                                         __LINE__, &
                                                         QCMATRIX_F_TEST_SRC)
                                    ierr = QcMatWrite_f(B, "QcMatGEMM_B", ASCII_VIEW)
                                    call QErrorCheckCode(io_log,   &
                                                         ierr,     &
                                                         __LINE__, &
                                                         QCMATRIX_F_TEST_SRC)
                                    ierr = QcMatWrite_f(C, "QcMatGEMM_C", ASCII_VIEW)
                                    call QErrorCheckCode(io_log,   &
                                                         ierr,     &
                                                         __LINE__, &
                                                         QCMATRIX_F_TEST_SRC)
                                    ierr = QcMatWrite_f(T, "QcMatGEMM_T", ASCII_VIEW)
                                    call QErrorCheckCode(io_log,   &
                                                         ierr,     &
                                                         __LINE__, &
                                                         QCMATRIX_F_TEST_SRC)
#endif
                                    write(io_log,100) "real part of C from BLAS"
                                    write(io_log,"(4Es24.12)") C_real
                                    write(io_log,100) "imaginary part of C from BLAS"
                                    write(io_log,"(4Es24.12)") C_imag
                                    call QErrorExit(io_log,   &
                                                    __LINE__, &
                                                    QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                                else
                                    write(io_log,100) "QcMatGEMM_f(beta!=0) passed", &
                                                      all_mat_operations(iop),       &
                                                      all_mat_operations(jop),       &
                                                      alpha(:,jdat),                 &
                                                      beta,                          &
                                                      data_type(idat),               &
                                                      sym_type(isym)
#endif
                                end if
                            else
                                write(io_log,100) "parameters",            &
                                                  all_mat_operations(iop), &
                                                  all_mat_operations(jop), &
                                                  alpha(:,jdat),           &
                                                  beta,                    &
                                                  data_type(idat),         &
                                                  sym_type(isym)
#if defined(QCMATRIX_ENABLE_VIEW)
                                ierr = QcMatWrite_f(A, "QcMatGEMM_A", ASCII_VIEW)
                                call QErrorCheckCode(io_log,   &
                                                     ierr,     &
                                                     __LINE__, &
                                                     QCMATRIX_F_TEST_SRC)
                                ierr = QcMatWrite_f(B, "QcMatGEMM_B", ASCII_VIEW)
                                call QErrorCheckCode(io_log,   &   
                                                     ierr,     &   
                                                     __LINE__, &
                                                     QCMATRIX_F_TEST_SRC)
                                ierr = QcMatWrite_f(C, "QcMatGEMM_C", ASCII_VIEW)
                                call QErrorCheckCode(io_log,   &   
                                                     ierr,     &   
                                                     __LINE__, &
                                                     QCMATRIX_F_TEST_SRC)
#endif
                                call QErrorExit(io_log,   &
                                                __LINE__, &
                                                QCMATRIX_F_TEST_SRC)
                            end if
                        end do
                        if (all_mat_operations(jop)==MAT_HERM_TRANSPOSE .or. &
                            all_mat_operations(jop)==MAT_COMPLEX_CONJUGATE) then
                            B_imag = -B_imag
                        end if
                    end do
                    if (all_mat_operations(iop)==MAT_HERM_TRANSPOSE .or. &
                        all_mat_operations(iop)==MAT_COMPLEX_CONJUGATE) then
                        A_imag = -A_imag
                    end if
                end do
                ! cleans the matrix
                ierr = QcMatDestroy_f(T)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                write(io_log,100) "QcMatDestroy_f(T) passed ..."
#endif
                ierr = QcMatDestroy_f(C)
                call QErrorCheckCode(io_log, ierr, __LINE__, QCMATRIX_F_TEST_SRC)
#if defined(DEBUG)
                write(io_log,100) "QcMatDestroy_f(C) passed ..."
#endif
            end do
        end do
#if !defined(DEBUG)
        write(io_log,100) "QcMatGEMM_f(beta!=0) passed ..."
#endif
        ! cleans
        deallocate(A_real)
        deallocate(A_imag)
        deallocate(B_real)
        deallocate(B_imag)
        deallocate(C_real)
        deallocate(C_imag)
        deallocate(T_real)
        deallocate(T_imag)
        return
100     format("test_f_QcMatGEMM>> ",A,2I4,4F12.6,2I4)
    end subroutine test_f_QcMatGEMM

#undef QCMATRIX_F_TEST_SRC
