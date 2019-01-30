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

   This file implements the function QcMatSetDimMat().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief sets the dimension of each block of a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate() and
         QcMatBlockCreate()
     \param[QInt:int]{in} num_row number of rows of each block
     \param[QInt:int]{in} num_col number of columns of each block
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatSetDimMat(QcMat *A, const QInt num_row, const QInt num_col)
{
    QInt irow, icol;
    QErrorCode err_code;
    if (num_row<1) {
        printf("QMatSetDimMat>> input number of rows %"QINT_FMT"\n", num_row);
        QErrorExit(FILE_AND_LINE, "invalid number of rows");
    }
    if (num_col<1) {
        printf("QMatSetDimMat>> input number of columns %"QINT_FMT"\n", num_col);
        QErrorExit(FILE_AND_LINE, "invalid number of columns");
    }
    /* FIXME: to remove */
    if (num_row!=num_col) {
        QErrorExit(FILE_AND_LINE, "non-square matrix not implemented yet");
    }
    /* sets the dimension of each block */
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            err_code = CmplxMatSetDimMat(&A->blocks[irow][icol], num_row, num_col);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatSetDimMat");
        }
    }
    return QSUCCESS;
}
