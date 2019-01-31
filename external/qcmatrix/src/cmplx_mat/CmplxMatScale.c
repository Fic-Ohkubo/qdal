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

   This file implements the function CmplxMatScale().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatScale QcMatScale
#define CmplxMatZeroEntries QcMatZeroEntries
#else
#include "impls/cmplx_mat.h"
#endif

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief scales all elements of a matrix by a given (complex) number
    \author Bin Gao
    \date 2012-04-04
    \param[QReal:real]{in} scal_number the scaling number with scal_number[0] being the
        real part and scal_number[1] the imaginary part
    \param[CmplxMat:struct]{inout} A the matrix to be scaled, should be at least assembled
        by CmplxMatAssemble()
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatScale(const QReal scal_number[], CmplxMat *A)
{
    QReal ratio_imag_real;  /* the ratio of imaginary part to the real part of the scaling number */
    QErrorCode err_code;
    /* the imaginary part of the scaling number is not 0 */
    if (QAbs(scal_number[1])>QZEROTHRSH) {
        /* a complex scaling number
           A_{R} = -scal_number[1]*A_{I}+scal_number[0]*A_{R}
           A_{I} = scal_number[1]*A_{R}+scal_number[0]*A_{I} */
        if (QAbs(scal_number[0])>QZEROTHRSH) {
            switch (A->data_type) {
            case QREALMAT:
                /* A_{I} = scal_number[1]*A_{R}
                   A'_{R} = scal_number[0]*A_{R} */
                err_code = RealMatDuplicate(&A->cmplx_mat[A->real_part],
                                            COPY_PATTERN_AND_VALUE,
                                            &A->cmplx_mat[A->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                err_code = RealMatScale(scal_number[1], &A->cmplx_mat[A->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                err_code = RealMatScale(scal_number[0], &A->cmplx_mat[A->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                break;
            case QIMAGMAT:
                /* A_{R} = -scal_number[1]*A_{I}
                   A'_{I} = scal_number[0]*A_{I} */
                err_code = RealMatDuplicate(&A->cmplx_mat[A->imag_part],
                                            COPY_PATTERN_AND_VALUE,
                                            &A->cmplx_mat[A->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                err_code = RealMatScale(-scal_number[1], &A->cmplx_mat[A->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                err_code = RealMatScale(scal_number[0], &A->cmplx_mat[A->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                break;
            case QCMPLXMAT:
                /* A'_{R} = scal_number[0]*A_{R} */
                err_code = RealMatScale(scal_number[0], &A->cmplx_mat[A->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                /* A''_{R} = -scal_number[1]*A_{I}+A'_{R} */
                err_code = RealMatAXPY(-scal_number[1],
                                       &A->cmplx_mat[A->imag_part],
                                       &A->cmplx_mat[A->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                /* A''_{I} = scal_number[1]*A_{R}+scal_number[0]*A_{I}
                           = scal_number[1]/scal_number[0]*A''_{R}
                           + (scal_number[0]^2+scal_number[1]^1)/scal_number[0]*A_{I} */
                ratio_imag_real = scal_number[1]/scal_number[0];
                err_code = RealMatScale(scal_number[0]+ratio_imag_real*scal_number[1],
                                        &A->cmplx_mat[A->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                err_code = RealMatAXPY(ratio_imag_real,
                                       &A->cmplx_mat[A->real_part],
                                       &A->cmplx_mat[A->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                break;
            default:
                printf("CmplxMatScale>> data type of matrix A: %d\n", A->data_type);
                QErrorExit(FILE_AND_LINE, "invalid data type");
            }
            /* changes the symmetry and data types of the matrix */
            A->sym_type = QNONSYMMAT;
            A->data_type = QCMPLXMAT;
        }
        /* an imaginary scaling number
           A_{R} = -scal_number[1]*A_{I}
           A_{I} = scal_number[1]*A_{R} */
        else {
            switch (A->data_type) {
            case QREALMAT:
                err_code = RealMatScale(scal_number[1], &A->cmplx_mat[A->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                break;
            case QIMAGMAT:
                err_code = RealMatScale(-scal_number[1], &A->cmplx_mat[A->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                break;
            case QCMPLXMAT:
                err_code = RealMatScale(scal_number[1], &A->cmplx_mat[A->real_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                err_code = RealMatScale(-scal_number[1], &A->cmplx_mat[A->imag_part]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                break;
            default:
                printf("CmplxMatScale>> data type of matrix A: %d\n", A->data_type);
                QErrorExit(FILE_AND_LINE, "invalid data type");
            }
            /* swaps the real and imaginary parts of the matrix */
            A->real_part = A->imag_part;
            A->imag_part = 1-A->imag_part;
            /* changes the symmetry and data types of the matrix */
            A->sym_type = -A->sym_type;
            A->data_type = -A->data_type;
        }
    }
    /* a real scaling number, A = scal_number[0]*A */
    else if (QAbs(scal_number[0])>QZEROTHRSH) {
        if (A->data_type==QREALMAT || A->data_type==QCMPLXMAT) {
            err_code = RealMatScale(scal_number[0], &A->cmplx_mat[A->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
        }
        if (A->data_type==QIMAGMAT || A->data_type==QCMPLXMAT) {
            err_code = RealMatScale(scal_number[0], &A->cmplx_mat[A->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
        }
    }
    /* the scaling number is zero, we have A = 0 */
    else {
        err_code = CmplxMatZeroEntries(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatZeroEntries");
    }
    return QSUCCESS;
}
