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

   This file implements the function QcMatSetDataType().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief sets the data types of matrix elements of some blocks
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate()
         and QcMatBlockCreate()
     \param[QInt:int]{in} num_blocks number of blocks to set the data types
     \param[QInt:int]{in} idx_block_row row indices of the blocks
     \param[QInt:int]{in} idx_block_col column indices of the blocks
     \param[QcDataType:int]{in} block_data_types given data types of the blocks,
         see file include/types/mat_data.h
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatSetDataType(QcMat *A,
                            const QInt num_blocks,
                            const QInt idx_block_row[],
                            const QInt idx_block_col[],
                            const QcDataType block_data_types[])
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
            printf("QcMatSetDataType>> input row index %"QINT_FMT" of block %"QINT_FMT"\n",
                   idx_block_row[iblk],
                   iblk);
            QErrorExit(FILE_AND_LINE, "invalid row index");
        }
#if defined(QCMATRIX_ZERO_BASED)
        if (idx_block_col[iblk]<0 || idx_block_col[iblk]>=A->dim_block) {
#else
        if (idx_block_col[iblk]<1 || idx_block_col[iblk]>A->dim_block) {
#endif
            printf("QcMatSetDataType>> input column index %"QINT_FMT" of block %"QINT_FMT"\n",
                   idx_block_col[iblk],
                   iblk);
            QErrorExit(FILE_AND_LINE, "invalid column index");
        }
        /* sets the data type */
        if (block_data_types[iblk]==QREALMAT ||
            block_data_types[iblk]==QIMAGMAT ||
            block_data_types[iblk]==QCMPLXMAT) {
#if defined(QCMATRIX_ZERO_BASED)
            err_code = CmplxMatSetDataType(&A->blocks[idx_block_row[iblk]][idx_block_col[iblk]],
                                           block_data_types[iblk]);
#else
            err_code = CmplxMatSetDataType(&A->blocks[idx_block_row[iblk]-1][idx_block_col[iblk]-1],
                                           block_data_types[iblk]);
#endif
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatSetDataType");
            /* this block will be assembled later */
#if defined(QCMATRIX_ZERO_BASED)
            A->assembled[idx_block_row[iblk]][idx_block_col[iblk]] = QTRUE;
#else
            A->assembled[idx_block_row[iblk]-1][idx_block_col[iblk]-1] = QTRUE;
#endif
        }
    }
    return QSUCCESS;
}
