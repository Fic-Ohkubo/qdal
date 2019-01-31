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

   This file implements the function CmplxMatSetDataType().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"

/*% \brief sets the data types of matrix elements of some blocks
    \author Bin Gao
    \date 2012-04-04
    \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate()
        and QcMatBlockCreate()
    \param[QInt:int]{in} num_blocks number of blocks to set the data types
    \param[QInt:int]{in} idx_block_row row indices of the blocks
    \param[QInt:int]{in} idx_block_col column indices of the blocks
    \param[QcDataType:int]{in} data_type given data types of the blocks, see file
        include/types/mat_data.h
    \return[QErrorCode:int] error information
*/
QErrorCode QcMatSetDataType(QcMat *A,
                            const QInt num_blocks,
                            const QInt idx_block_row[],
                            const QInt idx_block_col[],
                            const QcDataType data_type[])
{
    QInt real_num_blocks;
    QInt *real_idx_block_row;
    QInt *real_idx_block_col;
    QInt imag_num_blocks;
    QInt *imag_idx_block_row;
    QInt *imag_idx_block_col;
    QInt iblk;
    QErrorCode err_code;
    /* allocates memory for the non-zero blocks of the real and imaginary parts */
    real_idx_block_row = (QInt *)malloc(sizeof(QInt)*num_blocks);
    if (real_idx_block_row==NULL) {
        printf("QcMatSetDataType>> input number of blocks %"QINT_FMT"\n", num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for real_idx_block_row");
    }
    real_idx_block_col = (QInt *)malloc(sizeof(QInt)*num_blocks);
    if (real_idx_block_col==NULL) {
        printf("QcMatSetDataType>> input number of blocks %"QINT_FMT"\n", num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for real_idx_block_col");
    }
    imag_idx_block_row = (QInt *)malloc(sizeof(QInt)*num_blocks);
    if (imag_idx_block_row==NULL) {
        printf("QcMatSetDataType>> input number of blocks %"QINT_FMT"\n", num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for imag_idx_block_row");
    }
    imag_idx_block_col = (QInt *)malloc(sizeof(QInt)*num_blocks);
    if (imag_idx_block_col==NULL) {
        printf("QcMatSetDataType>> input number of blocks %"QINT_FMT"\n", num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for imag_idx_block_col");
    }
    /* gets the number of non-zero blocks and their indices according to the given data types */
    real_num_blocks = -1;
    imag_num_blocks = -1;
    for (iblk=0; iblk<num_blocks; iblk++) {
        switch (data_type[iblk]) {
        case QREALMAT:
            real_idx_block_row[++real_num_blocks] = idx_block_row[iblk];
            real_idx_block_col[real_num_blocks] = idx_block_col[iblk];
            break;
        case QIMAGMAT:
            imag_idx_block_row[++imag_num_blocks] = idx_block_row[iblk];
            imag_idx_block_col[imag_num_blocks] = idx_block_col[iblk];
            break;
        case QCMPLXMAT:
            real_idx_block_row[++real_num_blocks] = idx_block_row[iblk];
            real_idx_block_col[real_num_blocks] = idx_block_col[iblk];
            imag_idx_block_row[++imag_num_blocks] = idx_block_row[iblk];
            imag_idx_block_col[imag_num_blocks] = idx_block_col[iblk];
            break;
        default:
            break;
        }
    }
    real_num_blocks++;
    imag_num_blocks++;
    /* sets the non-zero blocks and data type */
    if (real_num_blocks>0) {
        err_code = RealMatSetNonZeroBlocks(&A->cmplx_mat[A->real_part],
                                           real_num_blocks,
                                           real_idx_block_row,
                                           real_idx_block_col);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetNonZeroBlocks");
        if (imag_num_blocks>0) {
            err_code = RealMatSetNonZeroBlocks(&A->cmplx_mat[A->imag_part],
                                               imag_num_blocks,
                                               imag_idx_block_row,
                                               imag_idx_block_col);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetNonZeroBlocks");
            A->data_type = QCMPLXMAT;
        }
        else {
            A->data_type = QREALMAT;
        }
    }
    else {
        if (imag_num_blocks>0) {
            err_code = RealMatSetNonZeroBlocks(&A->cmplx_mat[A->imag_part],
                                               imag_num_blocks,
                                               imag_idx_block_row,
                                               imag_idx_block_col);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetNonZeroBlocks");
            A->data_type = QIMAGMAT;
        }
        else {
            A->data_type = QNULLMAT;
        }
    }
    /* cleans */
    free(real_idx_block_row);
    real_idx_block_row = NULL;
    free(real_idx_block_col);
    real_idx_block_col = NULL;
    free(imag_idx_block_row);
    imag_idx_block_row = NULL;
    free(imag_idx_block_col);
    imag_idx_block_col = NULL;
    return QSUCCESS;
}

#else
#include "impls/cmplx_mat.h"

/*% \brief sets the data type of matrix elements
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{inout} A the matrix, should be created
        by CmplxMatCreate()
    \param[QcDataType:int]{in} data_type given data type, see file
        include/types/mat_data.h
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatSetDataType(CmplxMat *A, const QcDataType data_type)
{
    switch (data_type) {
    case QREALMAT:
        A->data_type = QREALMAT;
        break;
    case QIMAGMAT:
        A->data_type = QIMAGMAT;
        break;
    case QCMPLXMAT:
        A->data_type = QCMPLXMAT;
        break;
    default:
        A->data_type = QNULLMAT;
    }
    return QSUCCESS;
}
#endif /* defined(ADAPTER_BLOCK_REAL) */
