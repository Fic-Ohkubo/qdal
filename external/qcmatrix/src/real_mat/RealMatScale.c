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

   This file implements the function RealMatScale().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief scales all elements of a matrix by a given number
    \author Bin Gao
    \date 2014-06-16
    \param[QReal:real]{in} scal_number the scaling number
    \param[RealMat:struct]{inout} A the matrix to be scaled, should be at least assembled
        by RealMatAssemble()
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatScale(const QReal scal_number, RealMat *A)
{
    QInt inc_A=1;
    if (A->values==NULL) {
        QErrorExit(FILE_AND_LINE, "A is not assembled");
    }
    else {
        C_BLAS_SCAL(A->num_row*A->num_col, scal_number, A->values, inc_A);
    }
    return QSUCCESS;
}
