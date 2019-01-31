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
!!  This file provides the interface of BLAS for C codes.
!!
!!  2014-04-04, Bin Gao:
!!  * first version

! data types between C/Fortran
#include "api/qcmatrix_c_type.h"

    subroutine C_BLAS_SCAL(N, A, X, INCX) bind(C, name="C_BLAS_SCAL")
        use, intrinsic :: iso_c_binding, only: C_QINT,C_QREAL
        implicit none
        integer(kind=C_QINT), value, intent(in) :: N
        real(kind=C_QREAL), value, intent(in) :: A
        real(kind=C_QREAL), intent(inout) :: X(N)
        integer(kind=C_QINT), value, intent(in) :: INCX
        call BLAS_SCAL(N, A, X, INCX)
        return
    end subroutine C_BLAS_SCAL
    
    subroutine C_BLAS_COPY(N, X, INCX, Y, INCY) bind(C, name="C_BLAS_COPY")
        use, intrinsic :: iso_c_binding, only: C_QINT,C_QREAL
        implicit none
        integer(kind=C_QINT), value, intent(in) :: N
        real(kind=C_QREAL), intent(in) :: X(N)
        integer(kind=C_QINT), value, intent(in) :: INCX
        real(kind=C_QREAL), intent(inout) :: Y(N)
        integer(kind=C_QINT), value, intent(in) :: INCY
        call BLAS_COPY(N, X, INCX, Y, INCY)
        return
    end subroutine C_BLAS_COPY
    
    subroutine C_BLAS_AXPY(N, A, X, INCX, Y, INCY) bind(C, name="C_BLAS_AXPY")
        use, intrinsic :: iso_c_binding, only: C_QINT,C_QREAL
        implicit none
        integer(kind=C_QINT), value, intent(in) :: N
        real(kind=C_QREAL), value, intent(in) :: A
        real(kind=C_QREAL), intent(in) :: X(N)
        integer(kind=C_QINT), value, intent(in) :: INCX
        real(kind=C_QREAL), intent(inout) :: Y(N)
        integer(kind=C_QINT), value, intent(in) :: INCY
        call BLAS_AXPY(N, A, X, INCX, Y, INCY)
        return
    end subroutine C_BLAS_AXPY
    
    subroutine C_BLAS_DOT(N, X, INCX, Y, INCY, DOT_XY) bind(C, name="C_BLAS_DOT")
        use, intrinsic :: iso_c_binding, only: C_QINT,C_QREAL
        implicit none
        integer(kind=C_QINT), value, intent(in) :: N
        real(kind=C_QREAL), intent(in) :: X(N)
        integer(kind=C_QINT), value, intent(in) :: INCX
        real(kind=C_QREAL), intent(inout) :: Y(N)
        integer(kind=C_QINT), value, intent(in) :: INCY
        real(kind=C_QREAL), intent(out) :: DOT_XY
        call BLAS_DOT(N, X, INCX, Y, INCY, DOT_XY)
        return
    end subroutine C_BLAS_DOT
    
    subroutine C_BLAS_GEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) &
        bind(C, name="C_BLAS_GEMM")
        use, intrinsic :: iso_c_binding, only: C_CHAR,C_QINT,C_QREAL
        implicit none
        character(C_CHAR), value, intent(in) :: TRANSA
        character(C_CHAR), value, intent(in) :: TRANSB
        integer(kind=C_QINT), value, intent(in) :: M
        integer(kind=C_QINT), value, intent(in) :: N
        integer(kind=C_QINT), value, intent(in) :: K
        real(kind=C_QREAL), value, intent(in) :: ALPHA
        real(kind=C_QREAL), intent(in) :: A(M*K)
        integer(kind=C_QINT), value, intent(in) :: LDA
        real(kind=C_QREAL), intent(in) :: B(K*N)
        integer(kind=C_QINT), value, intent(in) :: LDB
        real(kind=C_QREAL), value, intent(in) :: BETA
        real(kind=C_QREAL), intent(inout) :: C(M*N)
        integer(kind=C_QINT), value, intent(in) :: LDC
        call BLAS_GEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
        return
    end subroutine C_BLAS_GEMM
