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

   This file implements the function QcMatGetDimMat().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief gets the dimension of each block of a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate() and
         QcMatBlockCreate()
     \param[QInt:int]{out} num_row number of rows of each block
     \param[QInt:int]{out} num_col number of columns of each block
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetDimMat(QcMat *A, QInt *num_row, QInt *num_col)
{
    QInt irow, icol;
    QErrorCode err_code;
    if (A->blocks==NULL) {
        QErrorExit(FILE_AND_LINE, "A->blocks is not allocated");
    }
    /* gets the dimension of the first assembled block */
    *num_row = 0;
    *num_col = 0;
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            if (A->assembled[irow][icol]==QTRUE) {
                err_code = CmplxMatGetDimMat(&A->blocks[irow][icol], num_row, num_col);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetDimMat");
                break;
            }
        }
        if (*num_row!=0 && *num_col!=0) break;
    }
    return QSUCCESS;
}
