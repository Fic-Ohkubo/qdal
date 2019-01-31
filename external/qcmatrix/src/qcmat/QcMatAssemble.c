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

   This file implements the function QcMatAssemble().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief assembles a matrix (e.g. allocating memory) so that it could be used
         in further matrix calculations, this function should be invoked after
         QcMatCreate(), QcMatBlockCreate() and QcMatSet...()
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{inout} A the matrix to be assembled
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatAssemble(QcMat *A)
{
    QInt irow, icol;
    QErrorCode err_code;
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            if (A->assembled[irow][icol]==QTRUE) {
                err_code = CmplxMatAssemble(&A->blocks[irow][icol]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatAssemble");
            }
        }
    }
    return QSUCCESS;
}
