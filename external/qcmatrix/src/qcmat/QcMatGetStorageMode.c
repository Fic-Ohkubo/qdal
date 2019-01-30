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

   This file implements the function QcMatGetStorageMode().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief gets the matrix storage mode of a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate() and
         QcMatBlockCreate()
     \param[QcStorageMode:int]{out} storage_mode return matrix storage mode,
         should be defined and implemented in external library
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetStorageMode(QcMat *A, QcStorageMode *storage_mode)
{
    QInt irow, icol;
    QErrorCode err_code;
    if (A->blocks==NULL) {
        QErrorExit(FILE_AND_LINE, "A->blocks is not allocated");
    }
    /* gets the storage mode of the first assembled block */
    *storage_mode = UNKNOWN_STORAGE_MODE;
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            if (A->assembled[irow][icol]==QTRUE) {
                err_code = CmplxMatGetStorageMode(&A->blocks[irow][icol], storage_mode);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetStorageMode");
                break;
            }
        }
        if (*storage_mode!=UNKNOWN_STORAGE_MODE) break;
    }
    return QSUCCESS;
}
