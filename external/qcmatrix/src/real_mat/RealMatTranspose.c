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

   This file implements the function RealMatTranspose().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief performs an in-place or out-of-place matrix operation B = op(A)
    \author Bin Gao
    \date 2014-06-16
    \param[QcMatOperation:int]{in} op_A the operation on the matrix A, see file
         include/types/mat_operations.h
    \param[RealMat:struct]{in} A the matrix to perform matrix operation, should be
        at least assembled by RealMatAssemble()
    \param[RealMat:struct]{inout} B the result matrix; it could be A, or it should
        be at least created by RealMatCreate()
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatTranspose(const QcMatOperation op_A, RealMat *A, RealMat *B)
{
    QInt irow, icol, offset_row;
    QInt size_val;       /* size of elements */
    QReal val_replaced;  /* holds element to be replaced, eventually becomes next element to move */
    QReal tmp_val;       /* temporary value for swapping */
    QInt new_location;   /* new location */
    QInt cycleBegin;     /* holds start of cycle */
    QInt ival;           /* iterator */
    QBool *val_moved;    /* hash to mark moved elements */
    QReal negative_one=-1;
    QErrorCode err_code;
    /* checks if the matrix A is assembled */
    if (A->values==NULL) {
        QErrorExit(FILE_AND_LINE, "A is not assembled");
    }
    switch (op_A) {
    /* no matrix operation will be performed */
    case MAT_NO_OPERATION:
        /* we have B = A */
        if (A!=B) {
            err_code = RealMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
        }
        break;
    /* B = A^{T} */
    case MAT_TRANSPOSE:
        /* A is a non-symmetric matrix */
        if (A->sym_type==QNONSYMMAT) {
            /* out-of-place transpose */
            if (A!=B) {
#if defined(QCMATRIX_STORAGE_MODE)
                B->storage_mode = A->storage_mode;
#endif
                B->num_row = A->num_col;
                B->num_col = A->num_row;
                err_code = RealMatAssemble(B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
                /* sets the values of the matrix B */
                offset_row = -A->num_row;
                for (icol=0; icol<A->num_col; icol++) {
                    offset_row += A->num_row;
                    for (irow=0; irow<A->num_row; irow++) {
                        /* B(JI) = A(IJ) */
                        B->values[irow*A->num_col+icol] = A->values[offset_row+irow];
                    }
                }
            }
            /* in place transpose A[R][C] --> A[C][R], modified from
               http://www.geeksforgeeks.org/inplace-m-x-n-size-matrix-transpose

               A[or][oc] --> A[nr][nc], in column major order:

               ol = oc x R + or
               nl = nc x C + nr
               nr = oc
               nc = or ==> nl = or x C + oc
               N  = R x C

               ol x C = oc x R x C + or x C
                      = oc x N     + or x C
                      = oc x N     + (nl - oc)
                      = oc x (N-1) + nl

               or, nl = ol x C - oc x (N-1) ==> nl = (ol x C) mod (N-1)
             */
            else {
                size_val = A->num_row*A->num_col-1;
                val_moved = (QBool *)malloc(size_val*sizeof(QBool));
                if (val_moved==NULL) {
                    QErrorExit(FILE_AND_LINE, "failed to allocate memory for val_moved");
                }
                /* the first and last elements won't move */
                val_moved[0] = QTRUE;
                val_moved[size_val-1] = QTRUE;
                for (ival=1; ival<size_val-1; ival++) {
                    val_moved[ival] = QFALSE;
                }
                ival = 1;
                while (ival<size_val) {
                    cycleBegin = ival;
                    val_replaced = A->values[ival];
                    do {
#if defined(QCMATRIX_ROW_MAJOR)
                        new_location = (ival*A->num_row)%size_val;
#else
                        new_location = (ival*A->num_col)%size_val;
#endif
                        tmp_val = val_replaced;
                        val_replaced = A->values[new_location];
                        A->values[new_location] = tmp_val;
                        val_moved[ival] = QTRUE;
                        ival = new_location;
                    } while (ival!=cycleBegin);
                    /* gets next move */
                    for (ival=1; ival<size_val && val_moved[ival]==QTRUE; ival++);
                }
                free(val_moved);
                val_moved = NULL;
            }
        }
        /* A is a symmetric or anti-symmetric matrix */
        else {
            /* out-of-place transpose */
            if (A!=B) {
                err_code = RealMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
            }
            if (B->sym_type==QANTISYMMAT) {
                err_code = RealMatScale(negative_one, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
            }
        }
        break;
    default:
        printf("RealMatTranspose>> operation on matrix A: %d\n", op_A);
        QErrorExit(FILE_AND_LINE, "invalid operation on matrix");
    }
    return QSUCCESS;
}
