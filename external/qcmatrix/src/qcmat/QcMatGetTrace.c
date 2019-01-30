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

   This file implements the function QcMatGetTrace().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief gets the traces of the first few diagonal blocks of a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate()
         and QcMatBlockCreate()
     \param[QInt:int]{in} num_blocks the number of diagonal blocks
     \param[QReal:real]{out} trace the traces, size is 2*\var{num_blocks}
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetTrace(QcMat *A, const QInt num_blocks, QReal *trace)
{
    QInt irow, jrow;
    QErrorCode err_code;
    if (num_blocks>A->dim_block) {
        printf("QcMatGetTrace>> input number of diagonal blocks %"QINT_FMT"\n",
               num_blocks);
        printf("QcMatGetTrace>> number of diagonal blocks of the matrix %"QINT_FMT"\n",
               A->dim_block);
        QErrorExit(FILE_AND_LINE, "invalid input number of diagonal blocks");
    }
    for (irow=0,jrow=0; irow<num_blocks; irow++) {
        if (A->assembled[irow][irow]==QTRUE) {
            err_code = CmplxMatGetTrace(&A->blocks[irow][irow], &trace[jrow]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetTrace");
            jrow += 2;
        }
        else {
            trace[jrow++] = 0;
            trace[jrow++] = 0;
        }
    }
    return QSUCCESS;
}
