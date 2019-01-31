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

   This file implements the function RealMatSetDimMat().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief sets the dimension of a matrix
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{inout} A the matrix, should be created by RealMatCreate()
    \param[QInt:int]{in} num_row number of rows of the matrix
    \param[QInt:int]{in} num_col number of columns of the matrix
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatSetDimMat(RealMat *A, const QInt num_row, const QInt num_col)
{
    /* FIXME: to remove */
    if (num_row!=num_col) {
        QErrorExit(FILE_AND_LINE, "non-square matrix not implemented yet");
    }
    A->num_row = num_row;
    A->num_col = num_col;
    return QSUCCESS;
}
