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

   This file implements the function RealMatZeroEntries().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief zeros all entries of a matrix
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{inout} A the matrix, should be at least assembled by RealMatAssemble()
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatZeroEntries(RealMat *A)
{
    //QInt inc_A=1;
    //QReal real_zero=0;
    QInt ival;
    QErrorCode err_code;
    /* assembles the matrix A if it was not */
    if (A->values==NULL) {
        err_code = RealMatAssemble(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
    }
    //C_BLAS_SCAL(A->num_row*A->num_col, real_zero, A->values, inc_A);
    for (ival=0; ival<A->num_row*A->num_col; ival++) {
        A->values[ival] = 0;
    }
    A->sym_type = QSYMMAT;
    return QSUCCESS;
}
