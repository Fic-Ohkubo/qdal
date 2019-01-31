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

   This file implements the function RealMatGetValues().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief gets the values of a matrix
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{in} A the matrix, should be at least created by RealMatCreate()
    \param[QInt:int]{in} idx_first_row index of the first row from which
        the values are got
    \param[QInt:int]{in} num_row_get number of rows that the values are got
    \param[QInt:int]{in} idx_first_col index of the first column from which
        the values are got
    \param[QInt:int]{in} num_col_get number of columns that the values are got
    \param[QReal:real]{out} *values values of the matrix
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatGetValues(RealMat *A,
                            const QInt idx_first_row,
                            const QInt num_row_get,
                            const QInt idx_first_col,
                            const QInt num_col_get,
                            QReal *values)
{
    QInt p_idx_first_row;  /* index of the first row */
    QInt p_idx_last_row;   /* index of the last row */
    QInt p_idx_first_col;  /* index of the first column */
    QInt p_idx_last_col;   /* index of the last column */
    QInt offset_val;       /* offset of values in the memory */
    QInt irow,icol,ival;   /* incremental recorders */
    if (values!=NULL) {
        if (A->num_row<0 || A->num_col<0) {
            printf("RealMatGetValues>> number of rows %"QINT_FMT"\n", A->num_row);
            printf("RealMatGetValues>> number of columns %"QINT_FMT"\n", A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid number of rows and/or columns");
        }
        p_idx_first_row = idx_first_row;
        p_idx_last_row = idx_first_row+num_row_get-1;
        p_idx_first_col = idx_first_col;
        p_idx_last_col = idx_first_col+num_col_get-1;
#if defined(QCMATRIX_ZERO_BASED)
        /* checks the index of the first row */
        if (p_idx_first_row<0 || p_idx_first_row>=A->num_row) {
            printf("RealMatGetValues>> index of the first row %"QINT_FMT"\n",
                   p_idx_first_row);
            printf("RealMatGetValues>> number of rows of the matrix %"QINT_FMT"\n",
                   A->num_row);
            QErrorExit(FILE_AND_LINE, "invalid index of the first row");
        }
        /* checks the number of rows to get */
        if (num_row_get<1 || p_idx_last_row>=A->num_row) {
            printf("RealMatGetValues>> index of the first row %"QINT_FMT"\n",
                   p_idx_first_row);
            printf("RealMatGetValues>> number of rows to get %"QINT_FMT"\n",
                   num_row_get);
            printf("RealMatGetValues>> number of rows of the matrix %"QINT_FMT"\n",
                   A->num_row);
            QErrorExit(FILE_AND_LINE, "invalid number of rows to get");
        }
        /* checks the index of the first column */
        if (p_idx_first_col<0 || p_idx_first_col>=A->num_col) {
            printf("RealMatGetValues>> index of the first column %"QINT_FMT"\n",
                   p_idx_first_col);
            printf("RealMatGetValues>> number of columns of the matrix %"QINT_FMT"\n",
                   A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid index of the first column");
        }
        /* checks the number of columns to get */
        if (num_col_get<1 || p_idx_last_col>=A->num_col) {
            printf("RealMatGetValues>> index of the first column %"QINT_FMT"\n",
                   p_idx_first_col);
            printf("RealMatGetValues>> number of columns to get %"QINT_FMT"\n",
                   num_col_get);
            printf("RealMatGetValues>> number of columns of the matrix %"QINT_FMT"\n",
                   A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid number of columns to get");
        }
#else
        if (p_idx_first_row<1 || p_idx_first_row>A->num_row) {
            printf("RealMatGetValues>> index of the first row %"QINT_FMT"\n",
                   p_idx_first_row);
            printf("RealMatGetValues>> number of rows of the matrix %"QINT_FMT"\n",
                   A->num_row);
            QErrorExit(FILE_AND_LINE, "invalid index of the first row");
        }
        if (num_row_get<1 || p_idx_last_row>A->num_row) {
            printf("RealMatGetValues>> index of the first row %"QINT_FMT"\n",
                   p_idx_first_row);
            printf("RealMatGetValues>> number of rows to get %"QINT_FMT"\n",
                   num_row_get);
            printf("RealMatGetValues>> number of rows of the matrix %"QINT_FMT"\n",
                   A->num_row);
            QErrorExit(FILE_AND_LINE, "invalid number of rows to get");
        }
        if (p_idx_first_col<1 || p_idx_first_col>A->num_col) {
            printf("RealMatGetValues>> index of the first column %"QINT_FMT"\n",
                   p_idx_first_col);
            printf("RealMatGetValues>> number of columns of the matrix %"QINT_FMT"\n",
                   A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid index of the first column");
        }
        if (num_col_get<1 || p_idx_last_col>A->num_col) {
            printf("RealMatGetValues>> index of the first column %"QINT_FMT"\n",
                   p_idx_first_col);
            printf("RealMatGetValues>> number of columns to get %"QINT_FMT"\n",
                   num_col_get);
            printf("RealMatGetValues>> number of columns of the matrix %"QINT_FMT"\n",
                   A->num_col);
            QErrorExit(FILE_AND_LINE, "invalid number of columns to get");
        }
#endif
        /* returns zero if the matrix A is not assembled */
        if (A->values==NULL) {
            irow = num_row_get*num_col_get;
            for (ival=0; ival<irow; ival++) {
                values[ival] = 0;
            }
        }
        else {
#if defined(QCMATRIX_ROW_MAJOR)
        /* gets the values, in row major order */
#if defined(QCMATRIX_ZERO_BASED) 
            offset_val = p_idx_first_row*A->num_col;
#else
            offset_val = (p_idx_first_row-1)*A->num_col;
            p_idx_first_col--;
#endif 
            for (irow=0,ival=0; irow<num_row_get; irow++) {
#if defined(QCMATRIX_ZERO_BASED)
                for (icol=p_idx_first_col; icol<=p_idx_last_col; icol++,ival++) {
#else
                for (icol=p_idx_first_col; icol<p_idx_last_col; icol++,ival++) {
#endif
                    values[ival] = A->values[offset_val+icol];
                }
                offset_val += A->num_col;
            }
#else
        /* gets the values, in column major order */
#if defined(QCMATRIX_ZERO_BASED)
            offset_val = p_idx_first_col*A->num_row;
#else
            offset_val = (p_idx_first_col-1)*A->num_row;
            p_idx_first_row--;
#endif
            for (icol=0,ival=0; icol<num_col_get; icol++) {
#if defined(QCMATRIX_ZERO_BASED)
                for (irow=p_idx_first_row; irow<=p_idx_last_row; irow++,ival++) {
#else
                for (irow=p_idx_first_row; irow<p_idx_last_row; irow++,ival++) {
#endif
                    values[ival] = A->values[offset_val+irow];
                }
                offset_val += A->num_row;
            }
#endif
        }
    }
    return QSUCCESS;
}
