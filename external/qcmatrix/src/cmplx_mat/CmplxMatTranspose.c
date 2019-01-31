/* QcMatrix: an abstract matrix library
   Copyright 2012-2015 Bin Gao

   QcMatrix is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   QcMatrix is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with QcMatrix. If not, see <http://www.gnu.org/licenses/>.

   This file implements the function CmplxMatTranspose().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatTranspose QcMatTranspose
#define CmplxMatDuplicate QcMatDuplicate
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief performs an in-place or out-of-place matrix operation B = op(A)
    \author Bin Gao
    \date 2012-04-04
    \param[QcMatOperation:int]{in} op_A the operation on the matrix A, see file
         include/types/mat_operations.h
    \param[CmplxMat:struct]{in} A the matrix to perform matrix operation, should be
        at least assembled by CmplxMatAssemble()
    \param[CmplxMat:struct]{inout} B the result matrix; it could be A, or it should
        be at least created by CmplxMatCreate()
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatTranspose(const QcMatOperation op_A, CmplxMat *A, CmplxMat *B)
{
    QReal negative_one=-1;
    QErrorCode err_code;
    /* checks if the matrix A is assembled */
    if (A->data_type==QNULLMAT) {
        printf("CmplxMatTranspose>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    switch (op_A) {
    /* no matrix operation will be performed */
    case MAT_NO_OPERATION:
        /* we have B = A */
        if (A!=B) {
            err_code = CmplxMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
        }
        break;
    /* B = A^{*} */
    case MAT_COMPLEX_CONJUGATE:
        /* out-of-place complex conjugate */
        if (A!=B) {
            err_code = CmplxMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
        }
        if (B->data_type==QIMAGMAT || B->data_type==QCMPLXMAT) {
            err_code = RealMatScale(negative_one, &B->cmplx_mat[B->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
        }
        break;
    /* B = A^{T} */
    case MAT_TRANSPOSE:
        /* A is a non-Hermitian matrix */
        if (A->sym_type==QNONSYMMAT) {
            /* out-of-place transpose */
            if (A!=B) {
                err_code = CmplxMatDuplicate(A, COPY_PATTERN_ONLY, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
            }
            /* real part, B_{R} = A_{R}^{T} */
            if (A->data_type==QREALMAT || A->data_type==QCMPLXMAT) {
                err_code = RealMatTranspose(MAT_TRANSPOSE,
                                            &A->cmplx_mat[A->real_part],
                                            &B->cmplx_mat[B->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatTranspose");
            }
            /* imaginary part, B_{I} = A_{I}^{T} */
            if (A->data_type==QIMAGMAT || A->data_type==QCMPLXMAT) {
                err_code = RealMatTranspose(MAT_TRANSPOSE,
                                            &A->cmplx_mat[A->imag_part],
                                            &B->cmplx_mat[B->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatTranspose");
            }
        }
        /* A is a Hermitian or anti-Hermitian matrix */
        else {
            /* out-of-place transpose */
            if (A!=B) {
                err_code = CmplxMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
            }
            /* for Hermitian matrix: A_{R}^{T} = A_{R}, A_{I}^{T} = -A_{I}, so
               B_{R} = A_{R}^{T} = A_{R}, B_{I} = A_{I}^{T} = -A_{I};
               for anti-Hermitian matrix: A_{R}^{T} = -A_{R}, A_{I}^{T} = A_{I}, so
               B_{R} = A_{R}^{T} = -A_{R}, B_{I} = A_{I}^{T} = A_{I} */
            if (B->sym_type==QSYMMAT) {
                if (B->data_type==QIMAGMAT || B->data_type==QCMPLXMAT) {
                    err_code = RealMatScale(negative_one, &B->cmplx_mat[B->imag_part]);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                }
            }
            else {
                if (B->data_type==QREALMAT || B->data_type==QCMPLXMAT) {
                    err_code = RealMatScale(negative_one, &B->cmplx_mat[B->real_part]);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                }
            }
        }
        break;
    /* B = A^{\dagger} */
    case MAT_HERM_TRANSPOSE:
        /* A is a non-Hermitian matrix */
        if (A->sym_type==QNONSYMMAT) {
            /* out-of-place Hermitian transpose */
            if (A!=B) {
                err_code = CmplxMatDuplicate(A, COPY_PATTERN_ONLY, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
            }
            /* real part, B_{R} = A_{R}^{T} */
            if (A->data_type==QREALMAT || A->data_type==QCMPLXMAT) {
                err_code = RealMatTranspose(MAT_TRANSPOSE,
                                            &A->cmplx_mat[A->real_part],
                                            &B->cmplx_mat[B->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatTranspose");
            }
            /* imaginary part, B_{I} = -A_{I}^{T} */
            if (A->data_type==QIMAGMAT || A->data_type==QCMPLXMAT) {
                err_code = RealMatTranspose(MAT_TRANSPOSE,
                                            &A->cmplx_mat[A->imag_part],
                                            &B->cmplx_mat[B->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatTranspose");
                err_code = RealMatScale(negative_one, &B->cmplx_mat[B->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
            }
        }
        /* A is a Hermitian or anti-Hermitian matrix */
        else {
            /* out-of-place Hermitian transpose */
            if (A!=B) {
                err_code = CmplxMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
            }
            /* Hermitian transpose of an anti-Hermitian matrix */
            if (B->sym_type==QANTISYMMAT) {
                if (B->data_type==QREALMAT || B->data_type==QCMPLXMAT) {
                    err_code = RealMatScale(negative_one, &B->cmplx_mat[B->real_part]);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                }
                if (B->data_type==QIMAGMAT || B->data_type==QCMPLXMAT) {
                    err_code = RealMatScale(negative_one, &B->cmplx_mat[B->imag_part]);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                }
            }
        }
        break;
    default:
        printf("CmplxMatTranspose>> operation on matrix A: %d\n", op_A);
        QErrorExit(FILE_AND_LINE, "invalid operation on matrix");
    }
    return QSUCCESS;
}
