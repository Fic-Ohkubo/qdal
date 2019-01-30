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

   This file implements the function CmplxMatZeroEntries().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatZeroEntries QcMatZeroEntries
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief zeros all entries of a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{inout} A the matrix, should be at least assembled by CmplxMatAssemble()
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatZeroEntries(CmplxMat *A)
{
    QErrorCode err_code;
    switch (A->data_type) {
    case QREALMAT:
        err_code = RealMatZeroEntries(&A->cmplx_mat[A->real_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatZeroEntries");
        break;
    case QIMAGMAT:
        err_code = RealMatZeroEntries(&A->cmplx_mat[A->imag_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatZeroEntries");
        break;
    case QCMPLXMAT:
        err_code = RealMatZeroEntries(&A->cmplx_mat[A->real_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatZeroEntries");
        err_code = RealMatZeroEntries(&A->cmplx_mat[A->imag_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatZeroEntries");
        break;
    default:
        printf("CmplxMatZeroEntries>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    A->sym_type = QSYMMAT;
    return QSUCCESS;
}
