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

   This file implements the function CmplxMatSetValues().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"

/*% \brief sets the values of a matrix
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
    /* sets the values of the real part */
    if (values_real!=NULL) {
        if (A->data_type==QNULLMAT) {
            err_code = RealMatAssemble(&A->cmplx_mat[A->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            A->data_type = QREALMAT;
        }
        else if (A->data_type==QIMAGMAT) {
            err_code = RealMatAssemble(&A->cmplx_mat[A->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            A->data_type = QCMPLXMAT;
        }
        err_code = RealMatSetValues(&A->cmplx_mat[A->real_part],
                                    idx_block_row,
                                    idx_block_col,
                                    idx_first_row,
                                    num_row_set,
                                    idx_first_col,
                                    num_col_set,
                                    values_real);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetValues");
    }
    /* sets the values of the imaginary part */
    if (values_imag!=NULL) {
        if (A->data_type==QNULLMAT) {
            err_code = RealMatAssemble(&A->cmplx_mat[A->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            A->data_type = QIMAGMAT;
        }
        else if (A->data_type==QREALMAT) {
            err_code = RealMatAssemble(&A->cmplx_mat[A->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            A->data_type = QCMPLXMAT;
        }
        err_code = RealMatSetValues(&A->cmplx_mat[A->imag_part],
                                    idx_block_row,
                                    idx_block_col,
                                    idx_first_row,
                                    num_row_set,
                                    idx_first_col,
                                    num_col_set,
                                    values_imag);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetValues");
    }
    return QSUCCESS;
}
#else
#include "impls/cmplx_mat.h"

/*% \brief sets the values of a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{inout} A the matrix, should be at least created by CmplxMatCreate()
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
QErrorCode CmplxMatSetValues(CmplxMat *A,
                             const QInt idx_first_row,
                             const QInt num_row_set,
                             const QInt idx_first_col,
                             const QInt num_col_set,
                             const QReal *values_real,
                             const QReal *values_imag)
{
    QErrorCode err_code;
    /* sets the values of the real part */
    if (values_real!=NULL) {
        if (A->data_type==QNULLMAT) {
            err_code = RealMatAssemble(&A->cmplx_mat[A->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            A->data_type = QREALMAT;
        }
        else if (A->data_type==QIMAGMAT) {
            err_code = RealMatAssemble(&A->cmplx_mat[A->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            A->data_type = QCMPLXMAT;
        }
        err_code = RealMatSetValues(&A->cmplx_mat[A->real_part],
                                    idx_first_row,
                                    num_row_set,
                                    idx_first_col,
                                    num_col_set,
                                    values_real);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetValues");
    }
    /* sets the values of the imaginary part */
    if (values_imag!=NULL) {
        if (A->data_type==QNULLMAT) {
            err_code = RealMatAssemble(&A->cmplx_mat[A->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            A->data_type = QIMAGMAT;
        }
        else if (A->data_type==QREALMAT) {
            err_code = RealMatAssemble(&A->cmplx_mat[A->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
            A->data_type = QCMPLXMAT;
        }
        err_code = RealMatSetValues(&A->cmplx_mat[A->imag_part],
                                    idx_first_row,
                                    num_row_set,
                                    idx_first_col,
                                    num_col_set,
                                    values_imag);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetValues");
    }
    return QSUCCESS;
}
#endif  /* defined(ADAPTER_BLOCK_REAL) */
