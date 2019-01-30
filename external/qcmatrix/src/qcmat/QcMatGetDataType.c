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

   This file implements the function QcMatGetDataType().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief gets the data types of matrix elements of some blocks
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate()
         and QcMatBlockCreate()
     \param[QInt:int]{in} num_blocks number of blocks to get the data types
     \param[QInt:int]{in} idx_block_row row indices of the blocks
     \param[QInt:int]{in} idx_block_col column indices of the blocks
     \param[QcDataType:int]{out} data_type data types of the blocks, see file
         include/types/mat_data.h
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetDataType(QcMat *A,
                            const QInt num_blocks,
                            const QInt idx_block_row[],
                            const QInt idx_block_col[],
                            QcDataType *data_type)
{
    QInt iblk;
    QErrorCode err_code;
    for (iblk=0; iblk<num_blocks; iblk++) {
        /* checks the row and column indices */
#if defined(QCMATRIX_ZERO_BASED)
        if (idx_block_row[iblk]<0 || idx_block_row[iblk]>=A->dim_block) {
#else
        if (idx_block_row[iblk]<1 || idx_block_row[iblk]>A->dim_block) {
#endif
            printf("QcMatGetDataType>> input row index %"QINT_FMT" of block %"QINT_FMT"\n",
                   idx_block_row[iblk],
                   iblk);
            QErrorExit(FILE_AND_LINE, "invalid row index");
        }
#if defined(QCMATRIX_ZERO_BASED)
        if (idx_block_col[iblk]<0 || idx_block_col[iblk]>=A->dim_block) {
#else
        if (idx_block_col[iblk]<1 || idx_block_col[iblk]>A->dim_block) {
#endif
            printf("QcMatGetDataType>> input column index %"QINT_FMT" of block %"QINT_FMT"\n",
                   idx_block_col[iblk],
                   iblk);
            QErrorExit(FILE_AND_LINE, "invalid column index");
        }
#if defined(QCMATRIX_ZERO_BASED)
        if (A->assembled[idx_block_row[iblk]][idx_block_col[iblk]]==QTRUE) {
            err_code = CmplxMatGetDataType(&A->blocks[idx_block_row[iblk]][idx_block_col[iblk]],
                                           &data_type[iblk]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetDataType");
        }
#else
        if (A->assembled[idx_block_row[iblk]-1][idx_block_col[iblk]-1]==QTRUE) {
            err_code = CmplxMatGetDataType(&A->blocks[idx_block_row[iblk]-1][idx_block_col[iblk]-1],
                                           &data_type[iblk]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetDataType");
        }
#endif
        else {
            data_type[iblk] = QNULLMAT;
        }
    }
    return QSUCCESS;
}
