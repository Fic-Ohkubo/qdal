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

   This file implements the function RealMatGetMatProdTrace().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief gets the trace of a matrix-matrix product A*op(B)
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{in} A the left matrix, should be at least assembled
        by RealMatAssemble()
    \param[RealMat:struct]{in} B the right matrix, should be at least assembled
        by RealMatAssemble()
    \param[QcMatOperation:int]{in} op_B the operation on the matrix B, see file
        include/types/mat_operations.h
    \param[QReal:real]{out} trace the trace
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatGetMatProdTrace(RealMat *A,
                                  RealMat *B,
                                  const QcMatOperation op_B,
                                  QReal *trace)
{
    QInt inc_A;
    QInt inc_B;
    QReal tmp_dot;
    QInt irow;
    if (A->values==NULL) {
        QErrorExit(FILE_AND_LINE, "A is not assembled");
    }
    else if (B->values==NULL) {
        QErrorExit(FILE_AND_LINE, "B is not assembled");
    }
    else {
        switch (op_B) {
        /* trace = \sum_{ij}A_{ij}*B_{ji} */
        case MAT_NO_OPERATION:
            /* ensures matrix-matrix multiplication */
            if (A->num_col!=B->num_row) {
                printf("RealMatGetMatProdTrace>> number of columns (A) %"QINT_FMT"\n",
                       A->num_col);
                printf("RealMatGetMatProdTrace>> number of rows (B) %"QINT_FMT"\n",
                       B->num_row);
                QErrorExit(FILE_AND_LINE, "invalid matrix dimensions");
            }
            /* ensures square product matrix */
            if (A->num_row!=B->num_col) {
                printf("RealMatGetMatProdTrace>> number of rows (A) %"QINT_FMT"\n",
                       A->num_row);
                printf("RealMatGetMatProdTrace>> number of columns (B) %"QINT_FMT"\n",
                       B->num_col);
                QErrorExit(FILE_AND_LINE, "invalid matrix dimensions");
            }
            *trace = 0;
            /* trace of anti-symmetric matrix is zero */
            if (A->sym_type*B->sym_type!=QANTISYMMAT) {
#if defined(QCMATRIX_ROW_MAJOR)
                inc_A = 1;
                for (irow=0; irow<A->num_row; irow++) {
                    /* gets \sum{j}A_{ij}*B_{ji} */
                    C_BLAS_DOT(A->num_col,
                               &A->values[irow*A->num_col],
                               inc_A,
                               &B->values[irow],
                               B->num_col,
                               &tmp_dot);
                    *trace += tmp_dot;
                }
#else
                inc_B = 1;
                for (irow=0; irow<A->num_row; irow++) {
                    /* gets \sum{j}A_{ij}*B_{ji} */
                    C_BLAS_DOT(A->num_col,
                               &A->values[irow],
                               A->num_row,
                               &B->values[irow*B->num_row],
                               inc_B,
                               &tmp_dot);
                    *trace += tmp_dot;
                }
#endif
            }
            break;
        /* trace = \sum_{ij}A_{ij}*B_{ij} */
        case MAT_TRANSPOSE:
            /* ensures matrix-matrix multiplication */
            if (A->num_col!=B->num_col) {
                printf("RealMatGetMatProdTrace>> number of columns (A) %"QINT_FMT"\n",
                       A->num_col);
                printf("RealMatGetMatProdTrace>> number of columns (B) %"QINT_FMT"\n",
                       B->num_col);
                QErrorExit(FILE_AND_LINE, "invalid matrix dimensions");
            }
            /* ensures square product matrix */
            if (A->num_row!=B->num_row) {
                printf("RealMatGetMatProdTrace>> number of rows (A) %"QINT_FMT"\n",
                       A->num_row);
                printf("RealMatGetMatProdTrace>> number of rows (B) %"QINT_FMT"\n",
                       B->num_row);
                QErrorExit(FILE_AND_LINE, "invalid matrix dimensions");
            }
            /* trace of anti-symmetric matrix is zero */
            if (A->sym_type*B->sym_type==QANTISYMMAT) {
                *trace = 0;
            }
            else {
                inc_A = 1;
                inc_B = 1;
                C_BLAS_DOT(A->num_row*A->num_col,
                           A->values,
                           inc_A,
                           B->values,
                           inc_B,
                           trace);
            }
            break;
        default:
            printf("RealMatGetMatProdTrace>> operation on matrix B: %d\n", op_B);
            QErrorExit(FILE_AND_LINE, "invalid matrix operation");
        }
    }
    return QSUCCESS;
}
