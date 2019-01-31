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

   This file implements the function RealMatGetTrace().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief gets the trace of a matrix
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{in} A the matrix, should be at least assembled by RealMatAssemble()
    \param[QReal:real]{out} trace the trace
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatGetTrace(RealMat *A, QReal *trace)
{
    QInt icol;
    if (A->values==NULL) {
        QErrorExit(FILE_AND_LINE, "A is not assembled");
    }
    else if (A->num_row!=A->num_col) {
        printf("RealMatGetTrace>> number of rows %"QINT_FMT"\n", A->num_row);
        printf("RealMatGetTrace>> number of columns %"QINT_FMT"\n", A->num_col);
        QErrorExit(FILE_AND_LINE, "A is not a square matrix");
    }
    else {
        *trace = 0;
        for (icol=0; icol<A->num_col; icol++) {
            *trace += A->values[icol*(A->num_row+1)];
        }
    }
    return QSUCCESS;
}
