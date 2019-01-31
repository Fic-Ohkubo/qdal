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

   This file implements the function RealMatAXPY().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief computes Y = a*X+Y
    \author Bin Gao
    \date 2014-06-16
    \param[QReal:real]{in} multiplier the multiplier a
    \param[RealMat:struct]{in} X the first matrix, should be at least assembled
        by RealMatAssemble()
    \param[RealMat:struct]{inout} Y the second matrix, should be at least created
        by RealMatCreate()
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatAXPY(const QReal multiplier, RealMat *X, RealMat *Y)
{
    QInt inc_X=1;
    QInt inc_Y=1;
    QReal scal_number;
    QErrorCode err_code;  
    /* the multiplier is zero, we have Y = Y */
    if (QAbs(multiplier)<QZEROTHRSH) {
        return QSUCCESS;
    }
    /* Y = (a+1)*Y */
    if (X==Y) {
        scal_number = multiplier+1;
        err_code = RealMatScale(scal_number, Y);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
    }
    else {
        /* Y is not assembled, we have Y = a*X */
        if (Y->values==NULL) {
            err_code = RealMatDuplicate(X, COPY_PATTERN_AND_VALUE, Y);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
            err_code = RealMatScale(multiplier, Y);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
        }
        /* Y is assembled */
        else {
            if (X->num_row!=Y->num_row) {
                printf("RealMatAXPY>> number of rows %"QINT_FMT", %"QINT_FMT"\n",
                       X->num_row,
                       Y->num_row);
                QErrorExit(FILE_AND_LINE, "invalid number of rows");
            }
            if (X->num_col!=Y->num_col) {
                printf("RealMatAXPY>> number of columns %"QINT_FMT", %"QINT_FMT"\n",
                       X->num_col,
                       Y->num_col);
                QErrorExit(FILE_AND_LINE, "invalid number of columns");
            }
            C_BLAS_AXPY(X->num_row*X->num_col,
                        multiplier,
                        X->values,
                        inc_X,
                        Y->values,
                        inc_Y);
            if (Y->sym_type!=X->sym_type) {
                Y->sym_type = QNONSYMMAT;
            }
        }
    }
    return QSUCCESS;
}
