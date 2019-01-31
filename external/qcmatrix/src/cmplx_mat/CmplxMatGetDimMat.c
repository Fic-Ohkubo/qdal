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

   This file implements the function CmplxMatGetDimMat().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatGetDimMat QcMatGetDimMat
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief gets the dimension of a matrix
    \param[CmplxMat:struct]{in} A the matrix, should be at least created by CmplxMatCreate()
    \param[QInt:int]{out} num_row number of rows of the matrix
    \param[QInt:int]{out} num_col number of columns of the matrix
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatGetDimMat(CmplxMat *A, QInt *num_row, QInt *num_col)
{
    QInt imag_nrow, imag_ncol;
    QErrorCode err_code;
    /* the dimensions of the real and imaginary parts should be the same */
    err_code = RealMatGetDimMat(&A->cmplx_mat[A->real_part], num_row, num_col);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetDimMat");
    err_code = RealMatGetDimMat(&A->cmplx_mat[A->imag_part], &imag_nrow, &imag_ncol);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetDimMat");
    if (*num_row<=0) {
        *num_row = imag_nrow;
    }
    else if ((*num_row!=imag_nrow) && (imag_nrow>0)) {
        printf("number of rows of the real part %"QINT_FMT"\n", *num_row);
        printf("number of rows of the imaginary part %"QINT_FMT"\n", imag_nrow);
        QErrorExit(FILE_AND_LINE, "invalid number of rows");
    }
    if (*num_col<=0) {
        *num_col = imag_ncol;
    }
    else if ((*num_col!=imag_ncol) && (imag_ncol>0)) {
        printf("number of columns of the real part %"QINT_FMT"\n", *num_col);
        printf("number of columns of the imaginary part %"QINT_FMT"\n", imag_ncol);
        QErrorExit(FILE_AND_LINE, "invalid number of columns");
    }
    return QSUCCESS;
}
