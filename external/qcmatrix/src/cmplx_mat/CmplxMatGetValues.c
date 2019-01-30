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

   This file implements the function CmplxMatGetValues().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"

/*% \brief gets the values of a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[QcMat:struct]{in} A the matrix, should be at least assembled by QcMatAssemble()
    \param[QInt:int]{in} idx_block_row index of the block row
    \param[QInt:int]{in} idx_block_col index of the block column
    \param[QInt:int]{in} idx_first_row index of the first row from which
        the values are got
    \param[QInt:int]{in} num_row_get number of rows that the values are got
    \param[QInt:int]{in} idx_first_col index of the first column from which
        the values are got
    \param[QInt:int]{in} num_col_get number of columns that the values are got
    \param[QReal:real]{out} *values_real values of the real part
    \param[QReal:real]{out} *values_imag values of the imaginary part
    \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetValues(QcMat *A,
                          const QInt idx_block_row,
                          const QInt idx_block_col,
                          const QInt idx_first_row,
                          const QInt num_row_get,
                          const QInt idx_first_col,
                          const QInt num_col_get,
                          QReal *values_real,
                          QReal *values_imag)
{
    QInt size_values;
    QInt ival;
    QErrorCode err_code;
    /* gets the values of the real part */
    if (values_real!=NULL) {
        switch (A->data_type) {
        case QREALMAT:
            err_code = RealMatGetValues(&A->cmplx_mat[A->real_part],
                                        idx_block_row,
                                        idx_block_col,
                                        idx_first_row,
                                        num_row_get,
                                        idx_first_col,
                                        num_col_get,
                                        values_real);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetValues");
            break;
        case QCMPLXMAT:
            err_code = RealMatGetValues(&A->cmplx_mat[A->real_part],
                                        idx_block_row,
                                        idx_block_col,
                                        idx_first_row,
                                        num_row_get,
                                        idx_first_col,
                                        num_col_get,
                                        values_real);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetValues");
            break;
        /* returns zero if the real part is not assembled */
        default:
            size_values = num_row_get*num_col_get;
            for (ival=0; ival<size_values; ival++) {
                values_real[ival] = 0;
            }
        }
    }
    /* gets the values of the imaginary part */
    if (values_imag!=NULL) {
        switch (A->data_type) {
        case QIMAGMAT:
            err_code = RealMatGetValues(&A->cmplx_mat[A->imag_part],
                                        idx_block_row,
                                        idx_block_col,
                                        idx_first_row,
                                        num_row_get,
                                        idx_first_col,
                                        num_col_get,
                                        values_imag);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetValues");
            break;
        case QCMPLXMAT:
            err_code = RealMatGetValues(&A->cmplx_mat[A->imag_part],
                                        idx_block_row,
                                        idx_block_col,
                                        idx_first_row,
                                        num_row_get,
                                        idx_first_col,
                                        num_col_get,
                                        values_imag);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetValues");
            break;
        /* returns zero if the imaginary part is not assembled */
        default:
            size_values = num_row_get*num_col_get;
            for (ival=0; ival<size_values; ival++) {
                values_imag[ival] = 0;
            }
        }
    }
    return QSUCCESS;
}
#else
#include "impls/cmplx_mat.h"

/*% \brief gets the values of a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{in} A the matrix, should be at least assembled by CmplxMatAssemble()
    \param[QInt:int]{in} idx_first_row index of the first row from which
        the values are got
    \param[QInt:int]{in} num_row_get number of rows that the values are got
    \param[QInt:int]{in} idx_first_col index of the first column from which
        the values are got
    \param[QInt:int]{in} num_col_get number of columns that the values are got
    \param[QReal:real]{out} *values_real values of the real part
    \param[QReal:real]{out} *values_imag values of the imaginary part
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatGetValues(CmplxMat *A,
                             const QInt idx_first_row,
                             const QInt num_row_get,
                             const QInt idx_first_col,
                             const QInt num_col_get,
                             QReal *values_real,
                             QReal *values_imag)
{
    QInt size_values;
    QInt ival;
    QErrorCode err_code;
    /* gets the values of the real part */
    if (values_real!=NULL) {
        switch (A->data_type) {
        case QREALMAT:
            err_code = RealMatGetValues(&A->cmplx_mat[A->real_part],
                                        idx_first_row,
                                        num_row_get,
                                        idx_first_col,
                                        num_col_get,
                                        values_real);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetValues");
            break;
        case QCMPLXMAT:
            err_code = RealMatGetValues(&A->cmplx_mat[A->real_part],
                                        idx_first_row,
                                        num_row_get,
                                        idx_first_col,
                                        num_col_get,
                                        values_real);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetValues");
            break;
        /* returns zero if the real part is not assembled */
        default:
            size_values = num_row_get*num_col_get;
            for (ival=0; ival<size_values; ival++) {
                values_real[ival] = 0;
            }
        }
    }
    /* gets the values of the imaginary part */
    if (values_imag!=NULL) {
        switch (A->data_type) {
        case QIMAGMAT:
            err_code = RealMatGetValues(&A->cmplx_mat[A->imag_part],
                                        idx_first_row,
                                        num_row_get,
                                        idx_first_col,
                                        num_col_get,
                                        values_imag);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetValues");
            break;
        case QCMPLXMAT:
            err_code = RealMatGetValues(&A->cmplx_mat[A->imag_part],
                                        idx_first_row,
                                        num_row_get,
                                        idx_first_col,
                                        num_col_get,
                                        values_imag);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetValues");
            break;
        /* returns zero if the imaginary part is not assembled */
        default:
            size_values = num_row_get*num_col_get;
            for (ival=0; ival<size_values; ival++) {
                values_imag[ival] = 0;
            }
        }
    }
    return QSUCCESS;
}
#endif  /* defined(ADAPTER_BLOCK_REAL) */
