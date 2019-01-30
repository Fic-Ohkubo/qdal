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

   This file implements the function QcMatGetAllValues().

   2013-03-29, Bin Gao:
   * first version
*/

#include "qcmatrix.h"
#include "utilities/qcmatrix_algebra.h"
#include "tests/qcmatrix_test_param.h"
#include "tests/qcmatrix_check_dim.h"

/*@% \brief gets all values of a matrix, may be only used for test suite
     \author Bin Gao
     \date 2013-03-29
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate()
     \param[QBool:int]{in} row_major if returning values in row major order
     \param[QInt:int]{in} size_values the size of values of the real and imaginary parts
     \param[QReal:real]{out} values_real values of the real part
     \param[QReal:real]{out} values_imag values of the imaginary part
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetAllValues(QcMat *A,
                             const QBool row_major,
                             const QInt size_values,
                             QReal *values_real,
                             QReal *values_imag)
{
    QBool assembled;
    QInt dim_block;
    QInt num_row;
    QInt num_col;
    QInt dim_qcmat;
    QInt size_mat;
    QInt idx_first_row;
    QInt num_row_get;
    QInt idx_first_col;
    QInt num_col_get;
    QReal *block_real;
    QReal *block_imag;
    QInt offset_block_row,offset_block_col,offset_block;
    QInt offset_val;
    QInt iblk,jblk;
    QInt irow,icol;
    QInt ival;
    QErrorCode err_code;
    err_code = QcMatIsAssembled(A, &assembled);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatIsAssembled(A)");
    if (assembled==QTRUE) {
        /* gets the dimensions */
        err_code = QcMatGetDimBlock(A, &dim_block);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDimBlock(A)");
        err_code = QcMatGetDimMat(A, &num_row, &num_col);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDimMat(A)");
        /* checks the dimensions */
        QCheckDimension(dim_block, num_row, num_col, FILE_AND_LINE);
        size_mat = num_row*num_col;
        if (size_values!=dim_block*dim_block*size_mat) {
            printf("QcMatGetAllValues>> total number of elements %"QINT_FMT"\n",
                   dim_block*dim_block*size_mat);
            printf("QcMatGetAllValues>> input size of the elements %"QINT_FMT"\n",
                   size_values);
            QErrorExit(FILE_AND_LINE, "invalid size of elements");
        }
/*FIXME: related to num_col for row major? */
        dim_qcmat = dim_block*num_row;
        /* gets the values block by block */
#if defined(QCMATRIX_ZERO_BASED)
        idx_first_row = 0;
        idx_first_col = 0;
#else
        idx_first_row = 1;
        idx_first_col = 1;
#endif
        num_row_get = num_row;
        num_col_get = num_col;
        block_real = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (block_real==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for block_real");
        }
        block_imag = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (block_imag==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for block_imag");
        }
        /* returns values in row major order */
        if (row_major==QTRUE) {
            for (iblk=0; iblk<dim_block; iblk++) {
                offset_block_row = iblk*dim_block*size_mat;
                for (jblk=0; jblk<dim_block; jblk++) {
/*FIXME: offset_block_col */
                    offset_block_col = jblk*num_row;
                    offset_block = offset_block_row+offset_block_col;
                    err_code = QcMatGetValues(A,
#if defined(QCMATRIX_ZERO_BASED)
                                              iblk,
                                              jblk,
#else
                                              iblk+1,
                                              jblk+1,
#endif
                                              idx_first_row,
                                              num_row_get,
                                              idx_first_col,
                                              num_col_get,
                                              block_real,
                                              block_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(A)");
                    if (values_real!=NULL) {
#if defined(QCMATRIX_ROW_MAJOR)
                        for (irow=0,ival=0; irow<num_row; irow++) {
                            offset_val = offset_block+irow*dim_qcmat;
                            for (icol=0; icol<num_col; icol++) {
                                values_real[offset_val+icol] = block_real[ival++];
                            }
                        }
#else
                        for (icol=0,ival=0; icol<num_col; icol++) {
                            offset_val = offset_block+icol;
                            for (irow=0; irow<num_row; irow++) {
                                values_real[offset_val+irow*dim_qcmat] = block_real[ival++];
                            }
                        }
#endif
                    }
                    if (values_imag!=NULL) {
#if defined(QCMATRIX_ROW_MAJOR)
                        for (irow=0,ival=0; irow<num_row; irow++) {
                            offset_val = offset_block+irow*dim_qcmat;
                            for (icol=0; icol<num_col; icol++) {
                                values_imag[offset_val+icol] = block_imag[ival++];
                            }
                        }
#else
                        for (icol=0,ival=0; icol<num_col; icol++) {
                            offset_val = offset_block+icol;
                            for (irow=0; irow<num_row; irow++) {
                                values_imag[offset_val+irow*dim_qcmat] = block_imag[ival++];
                            }
                        }
#endif
                    }
                }
            }
        }
        /* returns values in column major order */
        else {
            for (jblk=0; jblk<dim_block; jblk++) {
                offset_block_col = jblk*dim_block*size_mat;
                for (iblk=0; iblk<dim_block; iblk++) {
/*FIXME: offset_block_row? */
                    offset_block_row = iblk*num_col;
                    offset_block = offset_block_col+offset_block_row;
                    err_code = QcMatGetValues(A,
#if defined(QCMATRIX_ZERO_BASED)
                                              iblk,
                                              jblk,
#else
                                              iblk+1,
                                              jblk+1,
#endif
                                              idx_first_row,
                                              num_row_get,
                                              idx_first_col,
                                              num_col_get,
                                              block_real,
                                              block_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(A)");
                    if (values_real!=NULL) {
#if defined(QCMATRIX_ROW_MAJOR)
                        for (irow=0,ival=0; irow<num_row; irow++) {
                            offset_val = offset_block+irow;
                            for (icol=0; icol<num_col; icol++) {
                                values_real[offset_val+icol*dim_qcmat] = block_real[ival++];
                            }
                        }
#else
                        for (icol=0,ival=0; icol<num_col; icol++) {
                            offset_val = offset_block+icol*dim_qcmat;
                            for (irow=0; irow<num_row; irow++) {
                                values_real[offset_val+irow] = block_real[ival++];
                            }
                        }
#endif
                    }
                    if (values_imag!=NULL) {
#if defined(QCMATRIX_ROW_MAJOR)
                        for (irow=0,ival=0; irow<num_row; irow++) {
                            offset_val = offset_block+irow;
                            for (icol=0; icol<num_col; icol++) {
                                values_imag[offset_val+icol*dim_qcmat] = block_imag[ival++];
                            }
                        }
#else
                        for (icol=0,ival=0; icol<num_col; icol++) {
                            offset_val = offset_block+icol*dim_qcmat;
                            for (irow=0; irow<num_row; irow++) {
                                values_imag[offset_val+irow] = block_imag[ival++];
                            }
                        }
#endif
                    }
                }
            }
        }
        free(block_real);
        block_real = NULL;
        free(block_imag);
        block_imag = NULL;
    }
    /* returns zero if the matrix is not assembled */
    else {
        if (values_real!=NULL) {
            for (ival=0; ival<size_values; ival++) {
                values_real[ival] = 0;
            }
        }
        if (values_imag!=NULL) {
            for (ival=0; ival<size_values; ival++) {
                values_imag[ival] = 0;
            }
        }
    }
    return QSUCCESS;
}
