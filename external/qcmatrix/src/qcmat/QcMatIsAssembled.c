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

   This file implements the function QcMatIsAssembled().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief checks if a matrix is assembled or not
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate() and
         QcMatBlockCreate()
     \param[QBool:int]{out} assembled indicates if the matrix is assembled or not
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatIsAssembled(QcMat *A, QBool *assembled)
{
    QInt irow, icol;
    *assembled = QFALSE;
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            if (A->assembled[irow][icol]==QTRUE) {
                *assembled = QTRUE;
                break;
            }
        }
        if (*assembled==QTRUE) break;
    }
    return QSUCCESS;
}
