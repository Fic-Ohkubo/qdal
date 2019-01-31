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

   This file implements the function RealMatAddValues().

   2014-08-22, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief adds the values to a matrix
    \author Bin Gao
    \date 2014-08-22
    \param[RealMat:struct]{inout} A the matrix, should be at least created
        by RealMatCreate()
    \param[QInt:int]{in} idx_first_row index of the first row from which
        the values are added
    \param[QInt:int]{in} num_row_add number of rows to add
    \param[QInt:int]{in} idx_first_col index of the first column from which
        the values are added
    \param[QInt:int]{in} num_col_add number of columns to add
    \param[QReal:real]{in} *values values of the matrix
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatAddValues(RealMat *A,
                            const QInt idx_first_row,
                            const QInt num_row_add,
                            const QInt idx_first_col,
                            const QInt num_col_add,
                            const QReal *values)
{
    QInt p_idx_first_row;  /* index of the first row */
    QInt p_idx_last_row;   /* index of the last row */
    QInt p_idx_first_col;  /* index of the first column */
    QInt p_idx_last_col;   /* index of the last column */
    QInt offset_val;       /* offset of values in the memory */
    QInt irow,icol,ival;   /* incremental recorders */
    QErrorCode err_code;   /* error information */
    if (values!=NULL) {
        if (A->num_row<0 || A->num_col<0) {
            printf("RealMatAddValues>> number of rows %"QINT_FMT"\n", A->num_row);
            printf("RealMatAddValues>> number of columns %"QINT_FMT"\n", A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid number of rows and/or columns");
        }
        p_idx_first_row = idx_first_row;
        p_idx_last_row = idx_first_row+num_row_add-1;
        p_idx_first_col = idx_first_col;
        p_idx_last_col = idx_first_col+num_col_add-1;
#if defined(QCMATRIX_ZERO_BASED)
        /* checks the index of the first row */
        if (p_idx_first_row<0 || p_idx_first_row>=A->num_row) {
            printf("RealMatAddValues>> index of the first row %"QINT_FMT"\n",
                   p_idx_first_row);
            printf("RealMatAddValues>> number of rows of the matrix %"QINT_FMT"\n",
                   A->num_row);
            QErrorExit(FILE_AND_LINE, "invalid index of the first row");
        }
        /* checks the number of rows to add */
        if (num_row_add<1 || p_idx_last_row>=A->num_row) {
            printf("RealMatAddValues>> index of the first row %"QINT_FMT"\n",
                   p_idx_first_row);
            printf("RealMatAddValues>> number of rows to add %"QINT_FMT"\n",
                   num_row_add);
            printf("RealMatAddValues>> number of rows of the matrix %"QINT_FMT"\n",
                   A->num_row);
            QErrorExit(FILE_AND_LINE, "invalid number of rows to add");
        }
        /* checks the index of the first column */
        if (p_idx_first_col<0 || p_idx_first_col>=A->num_col) {
            printf("RealMatAddValues>> index of the first column %"QINT_FMT"\n",
                   p_idx_first_col);
            printf("RealMatAddValues>> number of columns of the matrix %"QINT_FMT"\n",
                   A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid index of the first column");
        }
        /* checks the number of columns to add */
        if (num_col_add<1 || p_idx_last_col>=A->num_col) {
            printf("RealMatAddValues>> index of the first column %"QINT_FMT"\n",
                   p_idx_first_col);
            printf("RealMatAddValues>> number of columns to add %"QINT_FMT"\n",
                   num_col_add);
            printf("RealMatAddValues>> number of columns of the matrix %"QINT_FMT"\n",
                   A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid number of columns to add");
        }
#else
        if (p_idx_first_row<1 || p_idx_first_row>A->num_row) {
            printf("RealMatAddValues>> index of the first row %"QINT_FMT"\n",
                   p_idx_first_row);
            printf("RealMatAddValues>> number of rows of the matrix %"QINT_FMT"\n",
                   A->num_row);
            QErrorExit(FILE_AND_LINE, "invalid index of the first row");
        }
        if (num_row_add<1 || p_idx_last_row>A->num_row) {
            printf("RealMatAddValues>> index of the first row %"QINT_FMT"\n",
                   p_idx_first_row);
            printf("RealMatAddValues>> number of rows to add %"QINT_FMT"\n",
                   num_row_add);
            printf("RealMatAddValues>> number of rows of the matrix %"QINT_FMT"\n",
                   A->num_row);
            QErrorExit(FILE_AND_LINE, "invalid number of rows to add");
        }
        if (p_idx_first_col<1 || p_idx_first_col>A->num_col) {
            printf("RealMatAddValues>> index of the first column %"QINT_FMT"\n",
                   p_idx_first_col);
            printf("RealMatAddValues>> number of columns of the matrix %"QINT_FMT"\n",
                   A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid index of the first column");
        }
        if (num_col_add<1 || p_idx_last_col>A->num_col) {
            printf("RealMatAddValues>> index of the first column %"QINT_FMT"\n",
                   p_idx_first_col);
            printf("RealMatAddValues>> number of columns to add %"QINT_FMT"\n",
                   num_col_add);
            printf("RealMatAddValues>> number of columns of the matrix %"QINT_FMT"\n",
                   A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid number of columns to add");
        }
#endif
        /* assembles and zeros the matrix A if it was not */
        if (A->values==NULL) {
            err_code = RealMatAssemble(A);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            err_code = RealMatZeroEntries(A);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatZeroEntries");
        }
#if defined(QCMATRIX_ROW_MAJOR)
        /* adds the values, in row major order */
#if defined(QCMATRIX_ZERO_BASED)
        offset_val = p_idx_first_row*A->num_col;
#else
        offset_val = (p_idx_first_row-1)*A->num_col;
        p_idx_first_col--;
#endif
        for (irow=0,ival=0; irow<num_row_add; irow++) {
#if defined(QCMATRIX_ZERO_BASED)
            for (icol=p_idx_first_col; icol<=p_idx_last_col; icol++,ival++) {
#else
            for (icol=p_idx_first_col; icol<p_idx_last_col; icol++,ival++) {
#endif
                A->values[offset_val+icol] += values[ival];
            }
            offset_val += A->num_col;
        }
#else
        /* adds the values, in column major order */
#if defined(QCMATRIX_ZERO_BASED)
        offset_val = p_idx_first_col*A->num_row;
#else
        offset_val = (p_idx_first_col-1)*A->num_row;
        p_idx_first_row--;
#endif
        for (icol=0,ival=0; icol<num_col_add; icol++) {
#if defined(QCMATRIX_ZERO_BASED)
            for (irow=p_idx_first_row; irow<=p_idx_last_row; irow++,ival++) {
#else
            for (irow=p_idx_first_row; irow<p_idx_last_row; irow++,ival++) {
#endif
                A->values[offset_val+irow] += values[ival];
            }
            offset_val += A->num_row;
        }
#endif
    }
    return QSUCCESS;
}
