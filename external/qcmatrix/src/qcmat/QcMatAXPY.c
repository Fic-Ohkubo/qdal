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

   This file implements the function QcMatAXPY().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief computes Y = a*X+Y
     \author Bin Gao
     \date 2012-04-04
     \param[QReal:real]{in} multiplier the complex multiplier a with multiplier[0]
         being the real part and multiplier[1] the imaginary part
     \param[QcMat:struct]{in} X the first matrix, should be at least assembled
         by QcMatAssemble()
     \param[QcMat:struct]{inout} Y the second matrix, should be at least created
         by QcMatCreate()
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatAXPY(const QReal multiplier[], QcMat *X, QcMat *Y)
{
    QReal scal_number[2];
    QInt irow, icol;
    QErrorCode err_code;
    /* Y = (a+1)*Y */
    if (X==Y) {
        scal_number[0] = multiplier[0]+1;
        scal_number[1] = multiplier[1];
        err_code = QcMatScale(scal_number, Y);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatScale");
    }
    else {
        /* the blocks of matrix Y is not created, we have Y = a*X */
        if (Y->blocks==NULL) {
            err_code = QcMatDuplicate(X, COPY_PATTERN_AND_VALUE, Y);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDuplicate");
            err_code = QcMatScale(multiplier, Y);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatScale");
        }
        /* the blocks of matrix Y is created */
        else {
            if (Y->dim_block!=X->dim_block) {
                printf("QcMatAXPY>> dimension of blocks (X) %"QINT_FMT"\n", X->dim_block);
                printf("QcMatAXPY>> dimension of blocks (Y) %"QINT_FMT"\n", Y->dim_block);
                QErrorExit(FILE_AND_LINE, "invalid dimension of blocks");
            }
            for (irow=0; irow<X->dim_block; irow++) {
                for (icol=0; icol<X->dim_block; icol++) {
                    if (X->assembled[irow][icol]==QTRUE) {
                        err_code = CmplxMatAXPY(multiplier,
                                                &X->blocks[irow][icol],
                                                &Y->blocks[irow][icol]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatAXPY");
                        Y->assembled[irow][icol] = QTRUE;
                    }
                }
            }
        }
    }
    return QSUCCESS;
}
