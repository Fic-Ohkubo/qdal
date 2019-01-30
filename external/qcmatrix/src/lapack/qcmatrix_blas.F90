!!  QcMatrix: square block complex matrix for quantum chemistry calculations
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
!!  This file provides the interface of calling BLAS library.
!!
!!  2014-09-19, Bin Gao:
!!  * first version

    subroutine BLAS_SCAL(N, A, X, INCX)
        implicit none
#include "api/qcmatrix_f_basic.h90"
        integer(kind=QINT), intent(in) :: N
        real(kind=QREAL), intent(in) :: A
        real(kind=QREAL), intent(inout) :: X(N)
        integer(kind=QINT), intent(in) :: INCX
#include "lapack/qcmatrix_blas.h90"
        integer(kind=4), parameter :: STDOUT = 6
        integer(kind=QBLAS_INT) BLAS_INT_N
        integer(kind=QBLAS_INT) BLAS_INT_INCX
#if defined(QCMATRIX_64BIT_INTEGER) && !defined(QCMATRIX_BLAS_64BIT)
        if (QBLAS_INT_MAX<N .or. QBLAS_INT_MAX<INCX) then
            write(STDOUT,"(A,3I6)") "BLAS_SCAL>> N, INCX, QBLAS_INT_MAX", &
                                    N, INCX, QBLAS_INT_MAX
            call QErrorExit(STDOUT, __LINE__, "BLAS_SCAL@src/lapack/qcmatrix_blas.F90")
        end if
#endif
        BLAS_INT_N = int(N, QBLAS_INT)
        BLAS_INT_INCX = int(INCX, QBLAS_INT)
#if defined(QCMATRIX_SINGLE_PRECISION)
        call SSCAL(BLAS_INT_N, A, X, BLAS_INT_INCX)
#else
        call DSCAL(BLAS_INT_N, A, X, BLAS_INT_INCX)
#endif
        return
    end subroutine BLAS_SCAL
    
    subroutine BLAS_COPY(N, X, INCX, Y, INCY)
        implicit none
#include "api/qcmatrix_f_basic.h90"
        integer(kind=QINT), intent(in) :: N
        real(kind=QREAL), intent(in) :: X(N)
        integer(kind=QINT), intent(in) :: INCX
        real(kind=QREAL), intent(inout) :: Y(N)
        integer(kind=QINT), intent(in) :: INCY
#include "lapack/qcmatrix_blas.h90"
        integer(kind=4), parameter :: STDOUT = 6
        integer(kind=QBLAS_INT) BLAS_INT_N
        integer(kind=QBLAS_INT) BLAS_INT_INCX
        integer(kind=QBLAS_INT) BLAS_INT_INCY
#if defined(QCMATRIX_64BIT_INTEGER) && !defined(QCMATRIX_BLAS_64BIT)
        if (QBLAS_INT_MAX<N .or. QBLAS_INT_MAX<INCX .or. QBLAS_INT_MAX<INCY) then
            write(STDOUT,"(A,4I16)") "BLAS_COPY>> N, INCX, INCY, QBLAS_INT_MAX", &
                                     N, INCX, INCY, QBLAS_INT_MAX
            call QErrorExit(STDOUT, __LINE__, "BLAS_COPY@src/lapack/qcmatrix_blas.F90")
        end if
#endif
        BLAS_INT_N = int(N, QBLAS_INT)
        BLAS_INT_INCX = int(INCX, QBLAS_INT)
        BLAS_INT_INCY = int(INCY, QBLAS_INT)
#if defined(QCMATRIX_SINGLE_PRECISION)
        call SCOPY(BLAS_INT_N, X, BLAS_INT_INCX, Y, BLAS_INT_INCY)
#else
        call DCOPY(BLAS_INT_N, X, BLAS_INT_INCX, Y, BLAS_INT_INCY)
#endif
        return
    end subroutine BLAS_COPY
    
    subroutine BLAS_AXPY(N, A, X, INCX, Y, INCY)
        implicit none
#include "api/qcmatrix_f_basic.h90"
        integer(kind=QINT), intent(in) :: N
        real(kind=QREAL), intent(in) :: A
        real(kind=QREAL), intent(in) :: X(N)
        integer(kind=QINT), intent(in) :: INCX
        real(kind=QREAL), intent(inout) :: Y(N)
        integer(kind=QINT), intent(in) :: INCY
#include "lapack/qcmatrix_blas.h90"
        integer(kind=4), parameter :: STDOUT = 6
        integer(kind=QBLAS_INT) BLAS_INT_N
        integer(kind=QBLAS_INT) BLAS_INT_INCX
        integer(kind=QBLAS_INT) BLAS_INT_INCY
#if defined(QCMATRIX_64BIT_INTEGER) && !defined(QCMATRIX_BLAS_64BIT)
        if (QBLAS_INT_MAX<N .or. QBLAS_INT_MAX<INCX .or. QBLAS_INT_MAX<INCY) then
            write(STDOUT,"(A,4I16)") "BLAS_AXPY>> N, INCX, INCY, QBLAS_INT_MAX", &
                                     N, INCX, INCY, QBLAS_INT_MAX
            call QErrorExit(STDOUT, __LINE__, "BLAS_AXPY@src/lapack/qcmatrix_blas.F90")
        end if
#endif
        BLAS_INT_N = int(N, QBLAS_INT)
        BLAS_INT_INCX = int(INCX, QBLAS_INT)
        BLAS_INT_INCY = int(INCY, QBLAS_INT)
#if defined(QCMATRIX_SINGLE_PRECISION)
        call SAXPY(BLAS_INT_N, A, X, BLAS_INT_INCX, Y, BLAS_INT_INCY)
#else
        call DAXPY(BLAS_INT_N, A, X, BLAS_INT_INCX, Y, BLAS_INT_INCY)
#endif
        return
    end subroutine BLAS_AXPY
    
    subroutine BLAS_DOT(N, X, INCX, Y, INCY, DOT_XY)
        implicit none
#include "api/qcmatrix_f_basic.h90"
        integer(kind=QINT), intent(in) :: N
        real(kind=QREAL), intent(in) :: X(N)
        integer(kind=QINT), intent(in) :: INCX
        real(kind=QREAL), intent(inout) :: Y(N)
        integer(kind=QINT), intent(in) :: INCY
        real(kind=QREAL), intent(out) :: DOT_XY
#if defined(QCMATRIX_SINGLE_PRECISION)
        real(kind=QREAL) SDOT
#else
        real(kind=QREAL) DDOT
#endif
#include "lapack/qcmatrix_blas.h90"
        integer(kind=4), parameter :: STDOUT = 6
        integer(kind=QBLAS_INT) BLAS_INT_N
        integer(kind=QBLAS_INT) BLAS_INT_INCX
        integer(kind=QBLAS_INT) BLAS_INT_INCY
#if defined(QCMATRIX_64BIT_INTEGER) && !defined(QCMATRIX_BLAS_64BIT)
        if (QBLAS_INT_MAX<N .or. QBLAS_INT_MAX<INCX .or. QBLAS_INT_MAX<INCY) then
            write(STDOUT,"(A,4I16)") "BLAS_DOT>> N, INCX, INCY, QBLAS_INT_MAX", &
                                     N, INCX, INCY, QBLAS_INT_MAX
            call QErrorExit(STDOUT, __LINE__, "BLAS_DOT@src/lapack/qcmatrix_blas.F90")
        end if
#endif
        BLAS_INT_N = int(N, QBLAS_INT)
        BLAS_INT_INCX = int(INCX, QBLAS_INT)
        BLAS_INT_INCY = int(INCY, QBLAS_INT)
#if defined(QCMATRIX_SINGLE_PRECISION)
        DOT_XY = SDOT(BLAS_INT_N, X, BLAS_INT_INCX, Y, BLAS_INT_INCY)
#else
        DOT_XY = DDOT(BLAS_INT_N, X, BLAS_INT_INCX, Y, BLAS_INT_INCY)
#endif
        return
    end subroutine BLAS_DOT
    
    subroutine BLAS_GEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
        implicit none
#include "api/qcmatrix_f_basic.h90"
        character*1, intent(in) :: TRANSA
        character*1, intent(in) :: TRANSB
        integer(kind=QINT), intent(in) :: M
        integer(kind=QINT), intent(in) :: N
        integer(kind=QINT), intent(in) :: K
        real(kind=QREAL), intent(in) :: ALPHA
        real(kind=QREAL), intent(in) :: A(M*K)
        integer(kind=QINT), intent(in) :: LDA
        real(kind=QREAL), intent(in) :: B(K*N)
        integer(kind=QINT), intent(in) :: LDB
        real(kind=QREAL), intent(in) :: BETA
        real(kind=QREAL), intent(inout) :: C(M*N)
        integer(kind=QINT), intent(in) :: LDC
#include "lapack/qcmatrix_blas.h90"
        integer(kind=4), parameter :: STDOUT = 6
        integer(kind=QBLAS_INT) BLAS_INT_M
        integer(kind=QBLAS_INT) BLAS_INT_N
        integer(kind=QBLAS_INT) BLAS_INT_K
        integer(kind=QBLAS_INT) BLAS_INT_LDA
        integer(kind=QBLAS_INT) BLAS_INT_LDB
        integer(kind=QBLAS_INT) BLAS_INT_LDC
#if defined(QCMATRIX_64BIT_INTEGER) && !defined(QCMATRIX_BLAS_64BIT)
        if (QBLAS_INT_MAX<M .or. QBLAS_INT_MAX<N .or. QBLAS_INT_MAX<K .or. &
            QBLAS_INT_MAX<LDA .or. QBLAS_INT_MAX<LDB .or. QBLAS_INT_MAX<LDC) then
            write(STDOUT,"(A,7I16)") "BLAS_GEMM>> M, N, K, LDA, LDB, LDC, QBLAS_INT_MAX", &
                                     M, N, K, LDA, LDB, LDC, QBLAS_INT_MAX
            call QErrorExit(STDOUT, __LINE__, "BLAS_GEMM@src/lapack/qcmatrix_blas.F90")
        end if
#endif
        BLAS_INT_M = int(M, QBLAS_INT)
        BLAS_INT_N = int(N, QBLAS_INT)
        BLAS_INT_K = int(K, QBLAS_INT)
        BLAS_INT_LDA = int(LDA, QBLAS_INT)
        BLAS_INT_LDB = int(LDB, QBLAS_INT)
        BLAS_INT_LDC = int(LDC, QBLAS_INT)
#if defined(QCMATRIX_SINGLE_PRECISION)
        call SGEMM(TRANSA,       &
                   TRANSB,       &
                   BLAS_INT_M,   &
                   BLAS_INT_N,   &
                   BLAS_INT_K,   &
                   ALPHA,        &
                   A,            &
                   BLAS_INT_LDA, &
                   B,            &
                   BLAS_INT_LDB, &
                   BETA,         &
                   C,            &
                   BLAS_INT_LDC)
#else
        call DGEMM(TRANSA,       &
                   TRANSB,       &
                   BLAS_INT_M,   &
                   BLAS_INT_N,   &
                   BLAS_INT_K,   &
                   ALPHA,        &
                   A,            &
                   BLAS_INT_LDA, &
                   B,            &
                   BLAS_INT_LDB, &
                   BETA,         &
                   C,            &
                   BLAS_INT_LDC)
#endif
        return
    end subroutine BLAS_GEMM
