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

   This file implements the function CmplxMatDuplicate().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatDuplicate QcMatDuplicate
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief duplicates a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{in} A the matrix, should be at least
        created by CmplxMatCreate()
    \param[QcDuplicateOption:int]{in} duplicate_option duplicate option,
        see file include/types/mat_duplicate.h
    \param[CmplxMat:struct]{inout} B the new matrix, should be at least
        created by CmplxMatCreate(), and all its previous information
        will be destroyed
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatDuplicate(CmplxMat *A,
                             const QcDuplicateOption duplicate_option,
                             CmplxMat *B)
{
    QInt which_part;
    QErrorCode err_code;
    if (A==B) return QSUCCESS;
    B->sym_type = A->sym_type;
    B->data_type = A->data_type;
    B->real_part = A->real_part;
    B->imag_part = A->imag_part;
    for (which_part=0; which_part<2; which_part++) {
        /* we require that MatDuplicate() is able to erase the information of
           the matrix B, which for instance could return B as QNULLMAT if the
           matrix A is QNULLMAT */
        err_code = RealMatDuplicate(&A->cmplx_mat[which_part],
                                    duplicate_option,
                                    &B->cmplx_mat[which_part]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
    }
    return QSUCCESS;
}
