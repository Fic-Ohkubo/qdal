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

   This file implements the function QcMatIsDiagonal().

   2014-04-03, Bin Gao:
   * first version
*/

#if !defined(ADAPTER_BLOCK_CMPLX) && !defined(ADAPTER_BLOCK_REAL)

#include "qcmatrix.h"

/*@% \brief checks if a matrix is block diagonal or not
     \author Bin Gao
     \date 2014-04-03
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate() and
         QcMatBlockCreate()
     \param[QBool:int]{out} is_diag indicates if the matrix is block diagonal or not
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatIsDiagonal(QcMat *A, QBool *is_diag)
{
    QInt irow, icol;
    *is_diag = QTRUE;
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<irow-1; icol++) {
            if (A->assembled[irow][icol]==QTRUE) {
                *is_diag = QFALSE;
                break;
            }
        }
        for (icol=irow+1; icol<A->dim_block; icol++) {
            if (A->assembled[irow][icol]==QTRUE) {
                *is_diag = QFALSE;
                break;
            }
        }
        if (*is_diag==QFALSE) break;
    }
    return QSUCCESS;
}

#endif
