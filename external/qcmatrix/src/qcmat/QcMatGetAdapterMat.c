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

   This file implements the function QcMatGetAdapterMat().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/* external library has implemented complex matrix */
#if defined(ADAPTER_CMPLX_MAT)
/*@% \brief gets the adapter matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate() and
         QcMatBlockCreate()
     \param[QInt:int]{in} idx_block_row index of the block row
     \param[QInt:int]{in} idx_block_col index of the block column
     \param[CmplxMat:struct]{out} *A_adapter the adapter matrix
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetAdapterMat(QcMat *A,
                              const QInt idx_block_row,
                              const QInt idx_block_col,
                              CmplxMat **A_adapter)
{
    /* checks the indices of the block row and column */
#if defined(QCMATRIX_ZERO_BASED)
    if (idx_block_row<0 || idx_block_row>=A->dim_block) {
#else
    if (idx_block_row<1 || idx_block_row>A->dim_block) {
#endif
        printf("QcMatGetAdapterMat>> input index of block row %"QINT_FMT"\n",
               idx_block_row);
        QErrorExit(FILE_AND_LINE, "invalid index of the block row");
    }
#if defined(QCMATRIX_ZERO_BASED)
    if (idx_block_col<0 || idx_block_col>=A->dim_block) {
#else
    if (idx_block_col<1 || idx_block_col>A->dim_block) {
#endif
        printf("QcMatGetAdapterMat>> input index of block column %"QINT_FMT"\n",
               idx_block_col);
        QErrorExit(FILE_AND_LINE, "invalid index of the block column");
    }
#if defined(QCMATRIX_ZERO_BASED)
    *A_adapter = &A->blocks[idx_block_row][idx_block_col];
#else
    *A_adapter = &A->blocks[idx_block_row-1][idx_block_col-1];
#endif
    return QSUCCESS;
}
/* external library has implemented real matrix */
#elif defined(ADAPTER_REAL_MAT)
/*@% \brief gets the adapter matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate()
         and QcMatBlockCreate()
     \param[QInt:int]{in} idx_block_row index of the block row
     \param[QInt:int]{in} idx_block_col index of the block column
     \param[QcDataType:int]{in} data_type which part to get, see file
         include/types/mat_data.h
     \param[RealMat:struct]{out} A_adapter the real or imaginary part
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetAdapterMat(QcMat *A,
                              const QInt idx_block_row,
                              const QInt idx_block_col,
                              const QcDataType data_type,
                              RealMat **A_adapter)
{
    QErrorCode err_code;
    /* checks the indices of the block row and column */
#if defined(QCMATRIX_ZERO_BASED)
    if (idx_block_row<0 || idx_block_row>=A->dim_block) {
#else
    if (idx_block_row<1 || idx_block_row>A->dim_block) {
#endif
        printf("QcMatGetAdapterMat>> input index of block row %"QINT_FMT"\n",
               idx_block_row);
        QErrorExit(FILE_AND_LINE, "invalid index of the block row");
    }
#if defined(QCMATRIX_ZERO_BASED)
    if (idx_block_col<0 || idx_block_col>=A->dim_block) {
#else
    if (idx_block_col<1 || idx_block_col>A->dim_block) {
#endif
        printf("QcMatGetAdapterMat>> input index of block column %"QINT_FMT"\n",
               idx_block_col);
        QErrorExit(FILE_AND_LINE, "invalid index of the block column");
    }
#if defined(QCMATRIX_ZERO_BASED)
    err_code = CmplxMatGetAdapterMat(&A->blocks[idx_block_row][idx_block_col],
                                     data_type,
                                     A_adapter);
#else
    err_code = CmplxMatGetAdapterMat(&A->blocks[idx_block_row-1][idx_block_col-1],
                                     data_type,
                                     A_adapter);
#endif
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetAdapterMat");
    return QSUCCESS;
}
#endif
