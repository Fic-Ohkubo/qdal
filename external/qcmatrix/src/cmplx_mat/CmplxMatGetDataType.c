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

   This file implements the function CmplxMatGetDataType().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"

/*% \brief gets the data types of matrix elements of some blocks
    \author Bin Gao
    \date 2012-04-04
    \param[QcMat:struct]{in} A the matrix, should be created by QcMatCreate() and
        QcMatBlockCreate()
    \param[QInt:int]{in} num_blocks number of blocks to set the data types
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
    QBool *real_non_zero;
    QBool *imag_non_zero;
    QInt iblk;
    QErrorCode err_code;
    switch case (A->data_type) {
    case QREALMAT:
        /* allocates memory for the non-zero blocks of the real part */
        real_non_zero = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (real_non_zero==NULL) {
            printf("QcMatGetDataType>> input number of blocks %"QINT_FMT"\n",
                   num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for real_non_zero");
        }
        /* gets the zero blocks of the real part */
        err_code = RealMatGetNonZeroBlocks(&A->cmplx_mat[A->real_part],
                                           num_blocks,
                                           idx_block_row,
                                           idx_block_col,
                                           real_non_zero);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetNonZeroBlocks");
        for (iblk=0; iblk<num_blocks; iblk++) {
            if (real_non_zero[iblk]==QTRUE) {
                data_type[iblk] = QREALMAT;
            }
            else {
                data_type[iblk] = QNULLMAT;
            }
        }
        free(real_non_zero);
        real_non_zero = NULL;
        break;
    case QIMAGMAT:
        /* allocates memory for the non-zero blocks of the imaginary part */
        imag_non_zero = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (imag_non_zero==NULL) {
            printf("QcMatGetDataType>> input number of blocks %"QINT_FMT"\n",
                   num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for imag_non_zero");
        }
        /* gets the zero blocks of the imaginary part */
        err_code = RealMatGetNonZeroBlocks(&A->cmplx_mat[A->imag_part],
                                           num_blocks,
                                           idx_block_row,
                                           idx_block_col,
                                           imag_non_zero);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetNonZeroBlocks");
        for (iblk=0; iblk<num_blocks; iblk++) {
            if (imag_non_zero[iblk]==QTRUE) {
                data_type[iblk] = QREALMAT;
            }
            else {
                data_type[iblk] = QNULLMAT;
            }
        }
        free(imag_non_zero);
        imag_non_zero = NULL;
        break;
    case QCMPLXMAT:
        /* allocates memory for the non-zero blocks of the real and imaginary parts */
        real_non_zero = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (real_non_zero==NULL) {
            printf("QcMatGetDataType>> input number of blocks %"QINT_FMT"\n",
                   num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for real_non_zero");
        }
        imag_non_zero = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (imag_non_zero==NULL) {
            printf("QcMatGetDataType>> input number of blocks %"QINT_FMT"\n",
                   num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for imag_non_zero");
        }
        /* gets the zero blocks of the real and imaginary parts */
        err_code = RealMatGetNonZeroBlocks(&A->cmplx_mat[A->real_part],
                                           num_blocks,
                                           idx_block_row,
                                           idx_block_col,
                                           real_non_zero);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetNonZeroBlocks");
        err_code = RealMatGetNonZeroBlocks(&A->cmplx_mat[A->imag_part],
                                           num_blocks,
                                           idx_block_row,
                                           idx_block_col,
                                           imag_non_zero);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetNonZeroBlocks");
        for (iblk=0; iblk<num_blocks; iblk++) {
            if (real_non_zero[iblk]==QTRUE) {
                if (imag_non_zero[iblk]==QTRUE) {
                    data_type[iblk] = QCMPLXMAT;
                }
                else {
                    data_type[iblk] = QREALMAT;
                }
            }
            else {
                if (imag_non_zero[iblk]==QTRUE) {
                    data_type[iblk] = QIMAGMAT;
                }
                else {
                    data_type[iblk] = QNULLMAT;
                }
            }
        }
        free(real_non_zero);
        real_non_zero = NULL;
        free(imag_non_zero);
        imag_non_zero = NULL;
        break;
    default:
        for (iblk=0; iblk<num_blocks; iblk++) {
            data_type[iblk] = QNULLMAT;
        }
    }
    return QSUCCESS;
}

#else
#include "impls/cmplx_mat.h"

/*% \brief gets the data type of matrix elements
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{in} A the matrix, should be at least created
        by CmplxMatCreate()
    \param[QcDataType:int]{out} data_type data type of the matrix, see file
        include/types/mat_data.h
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatGetDataType(CmplxMat *A, QcDataType *data_type)
{
    *data_type = A->data_type;
    return QSUCCESS;
}
#endif /* defined(ADAPTER_BLOCK_REAL) */
