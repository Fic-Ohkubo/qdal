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

   This file implements the function CmplxMatSetDimMat().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatSetDimMat QcMatSetDimMat
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief sets the dimension of a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{inout} A the matrix, should be created by CmplxMatCreate()
    \param[QInt:int]{in} num_row number of rows of the matrix
    \param[QInt:int]{in} num_col number of columns of the matrix
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatSetDimMat(CmplxMat *A, const QInt num_row, const QInt num_col)
{
    QInt which_part;
    QErrorCode err_code;
    if (num_row<1) {
        printf("CmplxMatSetDimMat>> input number of rows %"QINT_FMT"\n", num_row);
        QErrorExit(FILE_AND_LINE, "invalid number of rows");
    }
    if (num_col<1) {
        printf("CmplxMatSetDimMat>> input number of columns %"QINT_FMT"\n", num_col);
        QErrorExit(FILE_AND_LINE, "invalid number of columns");
    }
    /* FIXME: to remove */
    if (num_row!=num_col) {
        QErrorExit(FILE_AND_LINE, "non-square matrix not implemented yet");
    }
    for (which_part=0; which_part<2; which_part++) {
        err_code = RealMatSetDimMat(&A->cmplx_mat[which_part], num_row, num_col);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetDimMat");
    }
    return QSUCCESS;
}
