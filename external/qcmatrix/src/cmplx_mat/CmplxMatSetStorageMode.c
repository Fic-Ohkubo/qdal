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

   This file implements the function CmplxMatSetStorageMode().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatSetStorageMode QcMatSetStorageMode
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief sets the matrix storage mode
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{inout} A the matrix, should be created by CmplxMatCreate()
    \param[QcStorageMode:int]{in} storage_mode given matrix storage mode,
        should be defined and implemented in external library
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatSetStorageMode(CmplxMat *A, const QcStorageMode storage_mode)
{
    QInt which_part;
    QErrorCode err_code;
    for (which_part=0; which_part<2; which_part++) {
        err_code = RealMatSetStorageMode(&A->cmplx_mat[which_part], storage_mode);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetStorageMode");
    }
    return QSUCCESS;
}
