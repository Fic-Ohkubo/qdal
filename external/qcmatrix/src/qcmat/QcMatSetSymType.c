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

   This file implements the function QcMatSetSymType().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief sets the symmetry type of a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate()
         and QcMatBlockCreate()
     \param[QcSymType:int]{in} sym_type given symmetry type, see file
         include/types/mat_symmetry.h
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatSetSymType(QcMat *A, const QcSymType sym_type)
{
    QInt irow, icol;
    QErrorCode err_code;
    if (sym_type!=QSYMMAT && sym_type!=QANTISYMMAT && sym_type!=QNONSYMMAT) {
        printf("QcMatSetSymType>> input symmetry type %d\n", sym_type);
        QErrorExit(FILE_AND_LINE, "invalid symmetry type");
    }
    /* sets the symmetry type of blocks */
    for (irow=0; irow<A->dim_block; irow++) {
        /* lower triangular blocks */
        for (icol=0; icol<irow; icol++) {
            err_code = CmplxMatSetSymType(&A->blocks[irow][icol], QNONSYMMAT);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatSetSymType");
        }
        /* normally, only the diagonal blocks take the same symmetry type as the matrix */
        err_code = CmplxMatSetSymType(&A->blocks[irow][irow], sym_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatSetSymType");
        /* upper triangular blocks */
        for (icol=irow+1; icol<A->dim_block; icol++) {
            err_code = CmplxMatSetSymType(&A->blocks[irow][icol], QNONSYMMAT);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatSetSymType");
        }   
    }
/* FIXME: we could avoid of assembling upper triangular blocks for matrices with QSYMMAT or QANTISYMMAT */
    A->sym_type = sym_type;
    return QSUCCESS;
}
