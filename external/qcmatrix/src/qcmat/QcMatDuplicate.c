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

   This file implements the function QcMatDuplicate().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief duplicates a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate()
         and QcMatBlockCreate()
     \param[QcDuplicateOption:int]{in} duplicate_option duplicate option, see file
         include/types/mat_duplicate.h
     \param[QcMat:struct]{inout} B the new matrix, should be at least created by
         QcMatCreate(), and all its previous information will be destroyed
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatDuplicate(QcMat *A, const QcDuplicateOption duplicate_option, QcMat *B)
{
    QInt irow, icol;
    QErrorCode err_code;
    if (A==B) return QSUCCESS;
    /* erase all previous information of the matrix B */
    if (B->blocks!=NULL) {
        err_code = QcMatDestroy(B);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDestroy");
        err_code = QcMatCreate(B);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatCreate");
    }
    /* sets the dimension of blocks and creates the blocks of B */
    err_code = QcMatBlockCreate(B, A->dim_block);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatBlockCreate");
    B->sym_type = A->sym_type;
    /* duplicates all the blocks */
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            B->assembled[irow][icol] = A->assembled[irow][icol];
            /* we require that CmplxMatMatDuplicate() is able to erase the information of
               the matrix B, which for instance could return B as QNULLMAT if the
               matrix A is QNULLMAT */
            err_code = CmplxMatDuplicate(&A->blocks[irow][icol],
                                         duplicate_option,
                                         &B->blocks[irow][icol]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
        }
    }
    return QSUCCESS;
}
