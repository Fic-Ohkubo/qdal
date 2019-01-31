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

   This file implements the function CmplxMatIsAssembled().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatIsAssembled QcMatIsAssembled
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief checks if a matrix is assembled or not
    \param[CmplxMat:struct]{in} A the matrix, should be at least created by CmplxMatCreate()
    \param[QBool:int]{out} assembled indicates if the matrix is assembled or not
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatIsAssembled(CmplxMat *A, QBool *assembled)
{
    if (A->data_type==QIMAGMAT || A->data_type==QREALMAT || A->data_type==QCMPLXMAT) {
        *assembled = QTRUE;
    }
    else {
        *assembled = QFALSE;
    }
    return QSUCCESS;
}
