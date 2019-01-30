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

   This file implements the function QcMatSetValues().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief sets the values of a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{inout} A the matrix, should be at least created by QcMatCreate() and
         QcMatBlockCreate()
     \param[QInt:int]{in} idx_block_row index of the block row
     \param[QInt:int]{in} idx_block_col index of the block column
     \param[QInt:int]{in} idx_first_row index of the first row from which
         the values are set
     \param[QInt:int]{in} num_row_set number of rows that the values are set
     \param[QInt:int]{in} idx_first_col index of the first column from which
         the values are set
     \param[QInt:int]{in} num_col_set number of columns that the values are set
     \param[QReal:real]{in} *values_real values of the real part
     \param[QReal:real]{in} *values_imag values of the imaginary part
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatSetValues(QcMat *A,
                          const QInt idx_block_row,
                          const QInt idx_block_col,
                          const QInt idx_first_row,
                          const QInt num_row_set,
                          const QInt idx_first_col,
                          const QInt num_col_set,
                          const QReal *values_real,
                          const QReal *values_imag)
{
    QErrorCode err_code;
    /* checks the indices of the block row and column */
#if defined(QCMATRIX_ZERO_BASED)
    if (idx_block_row<0 || idx_block_row>=A->dim_block) {
#else
    if (idx_block_row<1 || idx_block_row>A->dim_block) {
#endif
        printf("QcMatSetValues>> input index of block row %"QINT_FMT"\n",
               idx_block_row);
        QErrorExit(FILE_AND_LINE, "invalid index of the block row");
    }
#if defined(QCMATRIX_ZERO_BASED)
    if (idx_block_col<0 || idx_block_col>=A->dim_block) {
#else
    if (idx_block_col<1 || idx_block_col>A->dim_block) {
#endif
        printf("QcMatSetValues>> input index of block column %"QINT_FMT"\n",
               idx_block_col);
        QErrorExit(FILE_AND_LINE, "invalid index of the block column");
    }
    /* assembles the block if it was not previously */
#if defined(QCMATRIX_ZERO_BASED)
    if (A->assembled[idx_block_row][idx_block_col]==QFALSE) {
        err_code = CmplxMatAssemble(&A->blocks[idx_block_row][idx_block_col]);
#else
    if (A->assembled[idx_block_row-1][idx_block_col-1]==QFALSE) {
        err_code = CmplxMatAssemble(&A->blocks[idx_block_row-1][idx_block_col-1]);
#endif
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatAssemble");
#if defined(QCMATRIX_ZERO_BASED)
        A->assembled[idx_block_row][idx_block_col] = QTRUE;
#else
        A->assembled[idx_block_row-1][idx_block_col-1] = QTRUE;
#endif
    }
#if defined(QCMATRIX_ZERO_BASED)
    err_code = CmplxMatSetValues(&A->blocks[idx_block_row][idx_block_col],
                                 idx_first_row,
                                 num_row_set,
                                 idx_first_col,
                                 num_col_set,
                                 values_real,
                                 values_imag);
#else
    err_code = CmplxMatSetValues(&A->blocks[idx_block_row-1][idx_block_col-1],
                                 idx_first_row,
                                 num_row_set,
                                 idx_first_col,
                                 num_col_set,
                                 values_real,
                                 values_imag);
#endif
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatSetValues");
    return QSUCCESS;
}
