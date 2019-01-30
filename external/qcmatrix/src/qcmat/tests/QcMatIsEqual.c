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

   This file implements the function QcMatIsEqual().

   2013-03-23, Bin Gao:
   * first version
*/

#include "qcmatrix.h"
#include "utilities/qcmatrix_algebra.h"
#include "tests/qcmatrix_test_param.h"
#include "tests/qcmatrix_check_dim.h"

/*@% \brief compares if two matrices are equal, may be only used for test suite
     \author Bin Gao
     \date 2013-03-23
     \param[QcMat:struct]{in} A the matrix, should be at least created by QcMatCreate()
     \param[QcMat:struct]{in} B the matrix, should be at least created by QcMatCreate()
     \param[QBool:int]{in} cf_values indicates if comparing values
     \param[QBool:int]{out} is_equal indicates if two matrices are equal (pattern and/or values)
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatIsEqual(QcMat *A, QcMat *B, const QBool cf_values, QBool *is_equal)
{
    QInt A_dim_block,B_dim_block;
    QInt A_nrow,A_ncol,B_nrow,B_ncol;
    QInt num_blocks;
    QInt *idx_block_row;
    QInt *idx_block_col;
    QcDataType *A_data_type,*B_data_type;
    QInt size_mat;
    QInt idx_first_row;
    QInt num_row_get;
    QInt idx_first_col;
    QInt num_col_get;
    QReal *A_real,*A_imag,*B_real,*B_imag;
    QInt iblk,jblk,kblk;
    QInt ival;
    QErrorCode err_code;
    /* compares the symmetry type */
    if (A->sym_type!=B->sym_type) {
        *is_equal = QFALSE;
        printf("QcMatIsEqual>> symmetry not equal %d %d\n", A->sym_type, B->sym_type);
        return QSUCCESS;
    }
    /* gets the dimensions */
    err_code = QcMatGetDimBlock(A, &A_dim_block);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDimBlock(A)");
    err_code = QcMatGetDimBlock(B, &B_dim_block);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDimBlock(B)");
    err_code = QcMatGetDimMat(A, &A_nrow, &A_ncol);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDimMat(A)");
    err_code = QcMatGetDimMat(B, &B_nrow, &B_ncol);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDimMat(B)");
    /* checks the dimensions */
    QCheckDimension(A_dim_block, A_nrow, A_ncol, FILE_AND_LINE);
    QCheckDimension(B_dim_block, B_nrow, B_ncol, FILE_AND_LINE);
    /* compares the dimensions */
    if (A_dim_block!=B_dim_block) {
        *is_equal = QFALSE;
        printf("QcMatIsEqual>> dimension of blocks not equal %"QINT_FMT" %"QINT_FMT"\n",
               A_dim_block,
               B_dim_block);
        return QSUCCESS;
    }
    if (A_nrow!=B_nrow) {
        *is_equal = QFALSE;
        printf("QcMatIsEqual>> number of rows not equal %"QINT_FMT" %"QINT_FMT"\n",
               A_nrow,
               B_nrow);
        return QSUCCESS;
    }
    if (A_ncol!=B_ncol) {
        *is_equal = QFALSE;
        printf("QcMatIsEqual>> number of columns not equal %"QINT_FMT" %"QINT_FMT"\n",
               A_ncol,
               B_ncol);
        return QSUCCESS;
    }
    /* gets the data types of the blocks */
    num_blocks = A_dim_block*A_dim_block;
    idx_block_row = (QInt *)malloc(sizeof(QInt)*num_blocks);
    if (idx_block_row==NULL) {
        printf("QcMatIsEqual>> number of blocks %"QINT_FMT"\n", num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for idx_block_row");
    }
    idx_block_col = (QInt *)malloc(sizeof(QInt)*num_blocks);
    if (idx_block_col==NULL) {
        printf("QcMatIsEqual>> number of blocks %"QINT_FMT"\n", num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for idx_block_col");
    }
    A_data_type = (QcDataType *)malloc(sizeof(QcDataType)*num_blocks);
    if (A_data_type==NULL) {
        printf("QcMatIsEqual>> number of blocks %"QINT_FMT"\n", num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for A_data_type");
    }
    B_data_type = (QcDataType *)malloc(sizeof(QcDataType)*num_blocks);
    if (B_data_type==NULL) {
        printf("QcMatIsEqual>> number of blocks %"QINT_FMT"\n", num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for B_data_type");
    }
    /* generates the indices of the blocks */
#if defined(QCMATRIX_ZERO_BASED)
    for (iblk=0,kblk=0; iblk<A_dim_block; iblk++) {
        for (jblk=0; jblk<A_dim_block; jblk++) {
#else
    for (iblk=1,kblk=0; iblk<=A_dim_block; iblk++) {
        for (jblk=1; jblk<=A_dim_block; jblk++) {
#endif
            idx_block_row[kblk] = iblk;
            idx_block_col[kblk++] = jblk;
        }
    }
    err_code = QcMatGetDataType(A, num_blocks, idx_block_row, idx_block_col, A_data_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDataType(A)");
    err_code = QcMatGetDataType(B, num_blocks, idx_block_row, idx_block_col, B_data_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetDataType(B)");
    /* compares the data types of the blocks */
    *is_equal = QTRUE;
    for (iblk=0; iblk<num_blocks; iblk++) {
        if (A_data_type[iblk]!=B_data_type[iblk]) {
            *is_equal = QFALSE;
            printf("QcMatIsEqual>> data type of block (%"QINT_FMT",%"QINT_FMT") not equal %d %d\n",
                   idx_block_row[iblk],
                   idx_block_col[iblk],
                   A_data_type[iblk],
                   B_data_type[iblk]);
            break;
        }
    }
    free(idx_block_row);
    idx_block_row = NULL;
    free(idx_block_col);
    idx_block_col = NULL;
    free(B_data_type);
    B_data_type = NULL;
    /* compare the values of two matrices */
    if (*is_equal==QTRUE && cf_values==QTRUE) {
        size_mat = A_nrow*A_ncol;
#if defined(QCMATRIX_ZERO_BASED)
        idx_first_row = 0;
        idx_first_col = 0;
#else
        idx_first_row = 1;
        idx_first_col = 1;
#endif
        num_row_get = A_nrow;
        num_col_get = A_ncol;
        /* allocates memory for the values of real and imaginary parts */
        A_real = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (A_real==NULL) {
            printf("QcMatIsEqual>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for A_real");
        }
        A_imag = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (A_imag==NULL) {
            printf("QcMatIsEqual>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for A_imag");
        }
        B_real = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (B_real==NULL) {
            printf("QcMatIsEqual>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for B_real");
        }
        B_imag = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (B_imag==NULL) {
            printf("QcMatIsEqual>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for B_imag");
        }
        /* gets the values of the blocks one by one and compare them */
#if defined(QCMATRIX_ZERO_BASED)
        for (iblk=0,kblk=0; iblk<A_dim_block; iblk++) {
            for (jblk=0; jblk<A_dim_block; jblk++) {
#else
        for (iblk=1,kblk=0; iblk<=A_dim_block; iblk++) {
            for (jblk=1; jblk<=A_dim_block; jblk++) {
#endif
                switch (A_data_type[kblk++]) {
                case QREALMAT:
                    err_code = QcMatGetValues(A,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_get,
                                              idx_first_col,
                                              num_col_get,
                                              A_real,
                                              NULL);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(A)");
                    err_code = QcMatGetValues(B,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_get,
                                              idx_first_col,
                                              num_col_get,
                                              B_real,
                                              NULL);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(B)");
                    for (ival=0; ival<size_mat; ival++) {
                        if (QAbs(A_real[ival]-B_real[ival])>CF_THRESHOLD) {
                            *is_equal = QFALSE;
                            printf("QcMatIsEqual>> real %"QINT_FMT"/(%"QINT_FMT",%"QINT_FMT") not equal %f %f\n",
                                   ival,
                                   iblk,
                                   jblk,
                                   A_real[ival],
                                   B_real[ival]);
                            break;
                        }
                    }
                    break;
                case QIMAGMAT:
                    err_code = QcMatGetValues(A,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_get,
                                              idx_first_col,
                                              num_col_get,
                                              NULL,
                                              A_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(A)");
                    err_code = QcMatGetValues(B,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_get,
                                              idx_first_col,
                                              num_col_get,
                                              NULL,
                                              B_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(B)");
                    for (ival=0; ival<size_mat; ival++) {
                        if (QAbs(A_imag[ival]-B_imag[ival])>CF_THRESHOLD) {
                            *is_equal = QFALSE;
                            printf("QcMatIsEqual>> imag %"QINT_FMT"/(%"QINT_FMT",%"QINT_FMT") not equal %f %f\n",
                                   ival,
                                   iblk,
                                   jblk,
                                   A_imag[ival],
                                   B_imag[ival]);
                            break;
                        }
                    }
                    break;
                case QCMPLXMAT:
                    err_code = QcMatGetValues(A,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_get,
                                              idx_first_col,
                                              num_col_get,
                                              A_real,
                                              A_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(A)");
                    err_code = QcMatGetValues(B,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_get,
                                              idx_first_col,
                                              num_col_get,
                                              B_real,
                                              B_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetValues(B)");
                    for (ival=0; ival<size_mat; ival++) {
                        if (QAbs(A_real[ival]-B_real[ival])>CF_THRESHOLD ||
                            QAbs(A_imag[ival]-B_imag[ival])>CF_THRESHOLD) {
                            *is_equal = QFALSE;
                            printf("QcMatIsEqual>> value %"QINT_FMT"/(%"QINT_FMT",%"QINT_FMT") not equal (%f,%f) (%f,%f)\n",
                                   ival,
                                   iblk,
                                   jblk,
                                   A_real[ival],
                                   A_imag[ival],
                                   B_real[ival],
                                   B_imag[ival]);
                            break;
                        }
                    }
                    break;
                default:
                    break;
                }
                if (*is_equal==QFALSE) break;
            }
            if (*is_equal==QFALSE) break;
        }
        /* cleans */
        free(A_real);
        A_real = NULL;
        free(A_imag);
        A_imag = NULL;
        free(B_real);
        B_real = NULL;
        free(B_imag);
        B_imag = NULL;
    }
    free(A_data_type);
    A_data_type = NULL;
    return QSUCCESS;
}
