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

   This file implements the function CmplxMatGetAdapterMat().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatGetAdapterMat QcMatGetAdapterMat
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief gets the adapter matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{in} A the matrix, should be created by CmplxMatCreate()
    \param[QcDataType:int]{in} data_type which part to extract, see file
        include/types/mat_data.h
    \param[RealMat:struct]{out} A_adapter the real or imaginary part
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatGetAdapterMat(CmplxMat *A, const QcDataType data_type, RealMat **A_adapter)
{
    switch (data_type) {
    case QREALMAT:
        *A_adapter = &A->cmplx_mat[A->real_part];
        break;
    case QIMAGMAT:
        *A_adapter = &A->cmplx_mat[A->imag_part];
        break;
    default:
        printf("CmplxMatGetAdapterMat>> data type: %d\n", data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    return QSUCCESS;
}
