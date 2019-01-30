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

   This file implements the function QcMatSetStorageMode().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief sets the matrix storage mode of a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate() and
         QcMatBlockCreate()
     \param[QcStorageMode:int]{in} storage_mode given matrix storage mode,
         should be defined and implemented in external library
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatSetStorageMode(QcMat *A, const QcStorageMode storage_mode)
{
    QInt irow, icol;
    QErrorCode err_code;
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            err_code = CmplxMatSetStorageMode(&A->blocks[irow][icol], storage_mode);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatSetStorageMode");
        }
    }
    return QSUCCESS;
}
