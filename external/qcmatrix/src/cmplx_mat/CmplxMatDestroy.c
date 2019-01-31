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

   This file implements the function CmplxMatDestroy().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatDestroy QcMatDestroy
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief frees space taken by a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{delete} A the matrix, should be at least created by CmplxMatCreate()
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatDestroy(CmplxMat *A)
{
    QInt which_part;
    QErrorCode err_code;
    /* prevents repeatedly frees the space */
    if (A->cmplx_mat!=NULL) {
        for (which_part=0; which_part<2; which_part++) {
            err_code = RealMatDestroy(&A->cmplx_mat[which_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDestroy");
        }
        free(A->cmplx_mat);
        A->cmplx_mat = NULL;
    }
    A->real_part = 0;
    A->imag_part = 1;
    A->data_type = QNULLMAT;
    A->sym_type = QNONSYMMAT;
    return QSUCCESS;
}
