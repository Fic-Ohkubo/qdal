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

   This file implements the function QcMatCfArray().

   2013-03-23, Bin Gao:
   * first version
*/

#include "qcmatrix.h"
#include "utilities/qcmatrix_algebra.h"
#include "tests/qcmatrix_test_param.h"
#include "tests/qcmatrix_check_dim.h"

/*@% \brief compares if the values of a matrix and two arrays (real and imaginary parts)
         are equal, may be only used for test suite
     \author Bin Gao
     \date 2013-03-23
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate()
     \param[QBool:int]{in} row_major if given values in row major order
     \param[QInt:int]{in} size_values the size of values of the real and imaginary parts
     \param[QReal:real]{in} values_real the values of real part
     \param[QReal:real]{in} values_imag the values of imaginary part
     \param[QBool:int]{out} is_equal indicates if the values of the matrix and the arrays are equal
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatCfArray(QcMat *A,
                        const QBool row_major,
                        const QInt size_values,
                        const QReal *values_real,
                        const QReal *values_imag,
                        QBool *is_equal)
{
    QBool assembled;
    QInt dim_block;
    QInt num_row;
    QInt num_col;
    QInt dim_qcmat;
    QInt num_blocks;
    QInt *idx_block_row;
    QInt *idx_block_col;
    QcDataType *data_type;
    QInt size_mat;
    QInt idx_first_row;
    QInt num_row_get;
    QInt idx_first_col;
    QInt num_col_get;
    QReal *A_real,*A_imag;
    QInt offset_block_row,offset_block_col,offset_block;
    QInt offset_val;
    QInt iblk,jblk,kblk;
    QInt irow,icol;
    QInt ival,jval;
    QErrorCode err_code;
    /* checks if the matrix is assembled */
    err_code = QcMatIsAssembled(A, &assembled);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatIsAssembled(A)");
    if (assembled==QFALSE) {
        *is_equal = QTRUE;
        for (ival=0; ival<size_values; ival++) {
            if (QAbs(values_real[ival])>CF_THRESHOLD ||
                QAbs(values_imag[ival])>CF_THRESHOLD) {
                printf("QcMatCfArray>> %"QINT_FMT": %f, %f\n",
                       ival,
                       values_real[ival],
                       values_imag[ival]);
                *is_equal = QFALSE;
                break;
            }
        }
    }
    else {
        /* gets the dimensions */
        err_code = QcMatGetDimBlock(A, &dim_block);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDimBlock(A)");
        err_code = QcMatGetDimMat(A, &num_row, &num_col);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDimMat(A)");
        /* checks the dimensions */
        QCheckDimension(dim_block, num_row, num_col, FILE_AND_LINE);
        num_blocks = dim_block*dim_block;
        size_mat = num_row*num_col;
        if (size_values!=num_blocks*size_mat) {
            printf("QcMatCfArray>> sizes %"QINT_FMT" %"QINT_FMT" %"QINT_FMT"\n",
                   size_values,
                   num_blocks,
                   size_mat);
            *is_equal = QFALSE;
            return QSUCCESS;
        }
/*FIXME: num_col for row major? */
        dim_qcmat = dim_block*num_row;
        /* gets the data types of the blocks */
        idx_block_row = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (idx_block_row==NULL) {
            printf("QcMatCfArray>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for idx_block_row");
        }
        idx_block_col = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (idx_block_col==NULL) {
            printf("QcMatCfArray>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for idx_block_col");
        }
        data_type = (QcDataType *)malloc(sizeof(QcDataType)*num_blocks);
        if (data_type==NULL) {
            printf("QcMatCfArray>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for data_type");
        }
        /* generates the indices of the blocks */
#if defined(QCMATRIX_ZERO_BASED)
        for (iblk=0,kblk=0; iblk<dim_block; iblk++) {
            for (jblk=0; jblk<dim_block; jblk++) {
#else
        for (iblk=1,kblk=0; iblk<=dim_block; iblk++) {
            for (jblk=1; jblk<=dim_block; jblk++) {
#endif
                if (row_major==QTRUE) {
                    idx_block_row[kblk] = iblk;
                    idx_block_col[kblk++] = jblk;
                }
                else {
                    idx_block_row[kblk] = jblk;
                    idx_block_col[kblk++] = iblk;
                }
            }
        }
        err_code = QcMatGetDataType(A, num_blocks, idx_block_row, idx_block_col, data_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDataType(A)");
        /* allocates memory for the values of real and imaginary parts */
#if defined(QCMATRIX_ZERO_BASED)
        idx_first_row = 0;
        idx_first_col = 0;
#else
        idx_first_row = 1;
        idx_first_col = 1;
#endif
        num_row_get = num_row;
        num_col_get = num_col;
        A_real = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (A_real==NULL) {
            printf("QcMatCfArray>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for A_real");
        }
        A_imag = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (A_imag==NULL) {
            printf("QcMatCfArray>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for A_imag");
        }
        /* compare the values block by block */
        *is_equal = QTRUE;
        for (iblk=0; iblk<num_blocks; iblk++) {
            if (row_major==QTRUE) {
/*FIXME: offset related to num_col? */
#if defined(QCMATRIX_ZERO_BASED)
                offset_block_row = idx_block_row[iblk]*dim_block*size_mat;
                offset_block_col = idx_block_col[iblk]*num_row;
#else
                offset_block_row = (idx_block_row[iblk]-1)*dim_block*size_mat;
                offset_block_col = (idx_block_col[iblk]-1)*num_row;
#endif
            }
            else {
#if defined(QCMATRIX_ZERO_BASED)
                offset_block_row = idx_block_row[iblk]*num_row;
                offset_block_col = idx_block_col[iblk]*dim_block*size_mat;
#else
                offset_block_row = (idx_block_row[iblk]-1)*num_row;
                offset_block_col = (idx_block_col[iblk]-1)*dim_block*size_mat;
#endif
            }
            offset_block = offset_block_row+offset_block_col;
            /* get the values of a block according to its data type */
            switch (data_type[iblk]) {
            case QREALMAT:
                err_code = QcMatGetValues(A,
                                          idx_block_row[iblk],
                                          idx_block_col[iblk],
                                          idx_first_row,
                                          num_row_get,
                                          idx_first_col,
                                          num_col_get,
                                          A_real,
                                          NULL);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(A)");
                if (row_major==QTRUE) {
#if defined(QCMATRIX_ROW_MAJOR)
                    for (irow=0,ival=0; irow<num_row; irow++) {
                        offset_val = offset_block+irow*dim_qcmat;
                        for (icol=0; icol<num_col; icol++,ival++) {
                            jval = offset_val+icol;
                            if (QAbs(values_real[jval]-A_real[ival])>CF_THRESHOLD ||
                                QAbs(values_imag[jval])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> real block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f %f, %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       A_real[ival],
                                       values_imag[jval]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#else
                    for (icol=0,ival=0; icol<num_col; icol++) {
                        offset_val = offset_block+icol;
                        for (irow=0; irow<num_row; irow++,ival++) {
                            jval = offset_val+irow*dim_qcmat;
                            if (QAbs(values_real[jval]-A_real[ival])>CF_THRESHOLD ||
                                QAbs(values_imag[jval])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> real block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f %f, %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       A_real[ival],
                                       values_imag[jval]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#endif
                }
                else {
#if defined(QCMATRIX_ROW_MAJOR)
                    for (irow=0,ival=0; irow<num_row; irow++) {
                        offset_val = offset_block+irow;
                        for (icol=0; icol<num_col; icol++,ival++) {
                            jval = offset_val+icol*dim_qcmat;
                            if (QAbs(values_real[jval]-A_real[ival])>CF_THRESHOLD ||
                                QAbs(values_imag[jval])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> real block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f %f, %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       A_real[ival],
                                       values_imag[jval]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#else
                    for (icol=0,ival=0; icol<num_col; icol++) {
                        offset_val = offset_block+icol*dim_qcmat;
                        for (irow=0; irow<num_row; irow++,ival++) {
                            jval = offset_val+irow;
                            if (QAbs(values_real[jval]-A_real[ival])>CF_THRESHOLD ||
                                QAbs(values_imag[jval])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> real block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f %f, %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       A_real[ival],
                                       values_imag[jval]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#endif
                }
                break;
            case QIMAGMAT:
                err_code = QcMatGetValues(A,
                                          idx_block_row[iblk],
                                          idx_block_col[iblk],
                                          idx_first_row,
                                          num_row_get,
                                          idx_first_col,
                                          num_col_get,
                                          NULL,
                                          A_imag);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(A)");
                if (row_major==QTRUE) {
#if defined(QCMATRIX_ROW_MAJOR)
                    for (irow=0,ival=0; irow<num_row; irow++) {
                        offset_val = offset_block+irow*dim_qcmat;
                        for (icol=0; icol<num_col; icol++,ival++) {
                            jval = offset_val+icol;
                            if (QAbs(values_real[jval])>CF_THRESHOLD ||
                                QAbs(values_imag[jval]-A_imag[ival])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> imaginary block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f, %f %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       values_imag[jval],
                                       A_imag[ival]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#else
                    for (icol=0,ival=0; icol<num_col; icol++) {
                        offset_val = offset_block+icol;
                        for (irow=0; irow<num_row; irow++,ival++) {
                            jval = offset_val+irow*dim_qcmat;
                            if (QAbs(values_real[jval])>CF_THRESHOLD ||
                                QAbs(values_imag[jval]-A_imag[ival])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> imaginary block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f, %f %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       values_imag[jval],
                                       A_imag[ival]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#endif
                }
                else {
#if defined(QCMATRIX_ROW_MAJOR)
                    for (irow=0,ival=0; irow<num_row; irow++) {
                        offset_val = offset_block+irow;
                        for (icol=0; icol<num_col; icol++,ival++) {
                            jval = offset_val+icol*dim_qcmat;
                            if (QAbs(values_real[jval])>CF_THRESHOLD ||
                                QAbs(values_imag[jval]-A_imag[ival])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> imaginary block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f, %f %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       values_imag[jval],
                                       A_imag[ival]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#else
                    for (icol=0,ival=0; icol<num_col; icol++) {
                        offset_val = offset_block+icol*dim_qcmat;
                        for (irow=0; irow<num_row; irow++,ival++) {
                            jval = offset_val+irow;
                            if (QAbs(values_real[jval])>CF_THRESHOLD ||
                                QAbs(values_imag[jval]-A_imag[ival])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> imaginary block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f, %f %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       values_imag[jval],
                                       A_imag[ival]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#endif
                }
                break;
            case QCMPLXMAT:
                err_code = QcMatGetValues(A,
                                          idx_block_row[iblk],
                                          idx_block_col[iblk],
                                          idx_first_row,
                                          num_row_get,
                                          idx_first_col,
                                          num_col_get,
                                          A_real,
                                          A_imag);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(A)");
                if (row_major==QTRUE) {
#if defined(QCMATRIX_ROW_MAJOR)
                    for (irow=0,ival=0; irow<num_row; irow++) {
                        offset_val = offset_block+irow*dim_qcmat;
                        for (icol=0; icol<num_col; icol++,ival++) {
                            jval = offset_val+icol;
                            if (QAbs(values_real[jval]-A_real[ival])>CF_THRESHOLD ||
                                QAbs(values_imag[jval]-A_imag[ival])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> complex block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f %f, %f %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       A_real[ival],
                                       values_imag[jval],
                                       A_imag[ival]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#else
                    for (icol=0,ival=0; icol<num_col; icol++) {
                        offset_val = offset_block+icol;
                        for (irow=0; irow<num_row; irow++,ival++) {
                            jval = offset_val+irow*dim_qcmat;
                            if (QAbs(values_real[jval]-A_real[ival])>CF_THRESHOLD ||
                                QAbs(values_imag[jval]-A_imag[ival])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> complex block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f %f, %f %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       A_real[ival],
                                       values_imag[jval],
                                       A_imag[ival]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#endif
                }
                else {
#if defined(QCMATRIX_ROW_MAJOR)
                    for (irow=0,ival=0; irow<num_row; irow++) {
                        offset_val = offset_block+irow;
                        for (icol=0; icol<num_col; icol++,ival++) {
                            jval = offset_val+icol*dim_qcmat;
                            if (QAbs(values_real[jval]-A_real[ival])>CF_THRESHOLD ||
                                QAbs(values_imag[jval]-A_imag[ival])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> complex block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f %f, %f %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       A_real[ival],
                                       values_imag[jval],
                                       A_imag[ival]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#else
                    for (icol=0,ival=0; icol<num_col; icol++) {
                        offset_val = offset_block+icol*dim_qcmat;
                        for (irow=0; irow<num_row; irow++,ival++) {
                            jval = offset_val+irow;
                            if (QAbs(values_real[jval]-A_real[ival])>CF_THRESHOLD ||
                                QAbs(values_imag[jval]-A_imag[ival])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> complex block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f %f, %f %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       A_real[ival],
                                       values_imag[jval],
                                       A_imag[ival]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#endif
                }
                break;
            default:
                if (row_major==QTRUE) {
#if defined(QCMATRIX_ROW_MAJOR)
                    for (irow=0,ival=0; irow<num_row; irow++) {
                        offset_val = offset_block+irow*dim_qcmat;
                        for (icol=0; icol<num_col; icol++,ival++) {
                            jval = offset_val+icol;
                            if (QAbs(values_real[jval])>CF_THRESHOLD ||
                                QAbs(values_imag[jval])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> zero block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f, %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       values_imag[jval]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#else
                    for (icol=0,ival=0; icol<num_col; icol++) {
                        offset_val = offset_block+icol;
                        for (irow=0; irow<num_row; irow++,ival++) {
                            jval = offset_val+irow*dim_qcmat;
                            if (QAbs(values_real[jval])>CF_THRESHOLD ||
                                QAbs(values_imag[jval])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> zero block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f, %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       values_imag[jval]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#endif
                }
                else {
#if defined(QCMATRIX_ROW_MAJOR)
                    for (irow=0,ival=0; irow<num_row; irow++) {
                        offset_val = offset_block+irow;
                        for (icol=0; icol<num_col; icol++,ival++) {
                            jval = offset_val+icol*dim_qcmat;
                            if (QAbs(values_real[jval])>CF_THRESHOLD ||
                                QAbs(values_imag[jval])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> zero block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f, %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       values_imag[jval]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#else
                    for (icol=0,ival=0; icol<num_col; icol++) {
                        offset_val = offset_block+icol*dim_qcmat;
                        for (irow=0; irow<num_row; irow++,ival++) {
                            jval = offset_val+irow;
                            if (QAbs(values_real[jval])>CF_THRESHOLD ||
                                QAbs(values_imag[jval])>CF_THRESHOLD) {
                                printf("QcMatCfArray>> zero block (%"QINT_FMT",%"QINT_FMT")/(%"QINT_FMT",%"QINT_FMT"): %f, %f\n",
                                       irow,
                                       icol,
                                       idx_block_row[iblk],
                                       idx_block_col[iblk],
                                       values_real[jval],
                                       values_imag[jval]);
                                *is_equal = QFALSE;
                                break;
                            }
                        }
                        if (*is_equal==QFALSE) break;
                    }
#endif
                }
            }
            if (*is_equal==QFALSE) break;
        }
        /* cleans */
        free(A_real);
        A_real = NULL;
        free(A_imag);
        A_imag = NULL;
        free(idx_block_row);
        idx_block_row = NULL;
        free(idx_block_col);
        idx_block_col = NULL;
        free(data_type);
        data_type = NULL;
    }
    return QSUCCESS;
}
