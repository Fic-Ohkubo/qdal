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

   This file implements the function CmplxMatAssemble().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatAssemble QcMatAssemble
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief assembles a matrix (e.g. allocating memory) so that it could be used
        in further matrix calculations, this function should be invoked after
        CmplxMatCreate() and CmplxMatSet...()
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{inout} A the matrix to be assembled
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatAssemble(CmplxMat *A)
{
    QErrorCode err_code;
    switch (A->data_type) {
    case QREALMAT:
        err_code = RealMatAssemble(&A->cmplx_mat[A->real_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
        break;
    case QIMAGMAT:
        err_code = RealMatAssemble(&A->cmplx_mat[A->imag_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
        break;
    case QCMPLXMAT:
        err_code = RealMatAssemble(&A->cmplx_mat[A->real_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
        err_code = RealMatAssemble(&A->cmplx_mat[A->imag_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
        break;
    default:
        break;
    }
    return QSUCCESS;
}
