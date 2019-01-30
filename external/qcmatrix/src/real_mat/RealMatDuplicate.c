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

   This file implements the function RealMatDuplicate().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief duplicates a matrix
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{in} A the matrix, should be at least
        created by RealMatCreate()
    \param[QcDuplicateOption:int]{in} duplicate_option duplicate option,
        see file include/types/mat_duplicate.h
    \param[RealMat:struct]{inout} B the new matrix, should be at least
        created by RealMatCreate(), and all its previous information
        will be destroyed
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatDuplicate(RealMat *A,
                            const QcDuplicateOption duplicate_option,
                            RealMat *B)
{
    QInt inc_A=1;
    QInt inc_B=1;
    QErrorCode err_code;
    if (A!=B) {
#if defined(QCMATRIX_STORAGE_MODE)
        B->storage_mode = A->storage_mode;
#endif
        B->sym_type = A->sym_type;
        B->num_row = A->num_row;
        B->num_col = A->num_col;
        if (B->values!=NULL) {
            free(B->values);
            B->values = NULL;
        }
        if (A->values!=NULL) {
            err_code = RealMatAssemble(B);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            if (duplicate_option==COPY_PATTERN_AND_VALUE) {
                C_BLAS_COPY(A->num_row*A->num_col, A->values, inc_A, B->values, inc_B);
            }
        }
    }
    return QSUCCESS;
}
