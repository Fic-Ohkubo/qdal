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

   This file implements the function QcMatSetRandMat().

   2013-03-23, Bin Gao:
   * first version
*/

#include <time.h>

#include "qcmatrix.h"
#include "utilities/qcmatrix_algebra.h"
#include "tests/qcmatrix_test_param.h"
#include "tests/qcmatrix_check_dim.h"

/*@% \brief sets the data types and values of a matrix randomly according to its
         symmetry and data types, may be only for test suite
     \author Bin Gao
     \date 2013-03-23
     \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate()
     \param[QcSymType:int]{in} sym_type given symmetry type, see file
         include/types/mat_symmetry.h
     \param[QcDataType:int]{in} data_type given data type of the matrix,
         see file include/types/mat_data.h
     \param[QInt:int]{in} dim_block the dimension of blocks
     \param[QInt:int]{in} num_row number of rows of each block
     \param[QInt:int]{in} num_col number of columns of each block
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatSetRandMat(QcMat *A,
                           const QcSymType sym_type,
                           const QcDataType data_type,
                           const QInt dim_block,
                           const QInt num_row,
                           const QInt num_col)
{
    QInt num_data_types;
    QcDataType *all_data_types;
    QInt num_blocks;
    QInt *idx_block_row;
    QInt *idx_block_col;
    QcDataType *block_data_type;
    QInt *diag_ind;
    QcDataType *diag_data_type;
    QBool assembled;
    QInt size_mat;
    QInt idx_first_row;
    QInt num_row_set;
    QInt idx_first_col;
    QInt num_col_set;
    QReal *values_real;
    QReal *values_imag;
    QInt iblk, jblk, kblk;
    QInt irow, icol;
    QInt ival, jval, kval;
    QErrorCode err_code;
    /* checks the dimensions */
    QCheckDimension(dim_block, num_row, num_col, FILE_AND_LINE);
    /* creates the blocks */
    err_code = QcMatBlockCreate(A, dim_block);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatBlockCreate");
    /* sets the symmetry type, possible data types of blocks and the dimension of each block */
    err_code = QcMatSetSymType(A, sym_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetSymType");
    switch (data_type) {
    case QREALMAT:
        num_data_types = 1;
        all_data_types = (QcDataType *)malloc(sizeof(QInt)*(num_data_types+1));
        if (all_data_types==NULL) {
            printf("QcMatSetRandMat>> number of data types %"QINT_FMT"\n",
                   num_data_types+1);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for all_data_types");
        }
        all_data_types[0] = QNULLMAT;
        all_data_types[1] = QREALMAT;
        break;
    case QIMAGMAT:
        num_data_types = 1;
        all_data_types = (QcDataType *)malloc(sizeof(QInt)*(num_data_types+1));
        if (all_data_types==NULL) {
            printf("QcMatSetRandMat>> number of data types %"QINT_FMT"\n",
                   num_data_types+1);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for all_data_types");
        }
        all_data_types[0] = QNULLMAT;
        all_data_types[1] = QIMAGMAT;
        break;
    case QCMPLXMAT:
        num_data_types = 3;
        all_data_types = (QcDataType *)malloc(sizeof(QInt)*(num_data_types+1));
        if (all_data_types==NULL) {
            printf("QcMatSetRandMat>> number of data types %"QINT_FMT"\n",
                   num_data_types+1);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for all_data_types");
        }
        all_data_types[0] = QNULLMAT;
        all_data_types[1] = QREALMAT;
        all_data_types[2] = QIMAGMAT;
        all_data_types[3] = QCMPLXMAT;
        break;
    default:
        printf("QcMatSetRandMat>> input data type %d\n", data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    err_code = QcMatSetDimMat(A, num_row, num_col);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetDimMat");
    size_mat = num_row*num_col;
#if defined(QCMATRIX_ZERO_BASED)
    idx_first_row = 0;
    idx_first_col = 0;
#else
    idx_first_row = 1;
    idx_first_col = 1;
#endif
    num_row_set = num_row;
    num_col_set = num_col;
    /* randomizes seed */
    srand(time(NULL));
    /* Hermitian or anti-Hermitian matrix */
    if (sym_type==QSYMMAT || sym_type==QANTISYMMAT) {
        /* allocates memory for row and column indices and data types */
        num_blocks = dim_block*(dim_block-1)/2;
        idx_block_row = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (idx_block_row==NULL) {
            printf("QcMatSetRandMat>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for idx_block_row");
        }
        idx_block_col = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (idx_block_col==NULL) {
            printf("QcMatSetRandMat>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for idx_block_col");
        }
        block_data_type = (QcDataType *)malloc(sizeof(QcDataType)*num_blocks);
        if (block_data_type==NULL) {
            printf("QcMatSetRandMat>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for block_data_type");
        }
        /* generates the data types of the upper blocks randomly */
#if defined(QCMATRIX_ZERO_BASED)
        for (iblk=0,kblk=0; iblk<dim_block; iblk++) {
            for (jblk=iblk+1; jblk<dim_block; jblk++) {
#else
        for (iblk=1,kblk=0; iblk<=dim_block; iblk++) {
            for (jblk=iblk+1; jblk<=dim_block; jblk++) {
#endif
                idx_block_row[kblk] = iblk;
                idx_block_col[kblk] = jblk;
                block_data_type[kblk++] = all_data_types[QRandInt(0,num_data_types)];
            }
        }
        err_code = QcMatSetDataType(A,
                                    num_blocks,
                                    idx_block_row,
                                    idx_block_col,
                                    block_data_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetDataType");
        /* the lower blocks take exactly the same data types as their corresponding blocks */
        err_code = QcMatSetDataType(A,
                                    num_blocks,
                                    idx_block_col,
                                    idx_block_row,
                                    block_data_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetDataType");
        free(idx_block_row);
        idx_block_row = NULL;
        free(idx_block_col);
        idx_block_col = NULL;
        /* generates the data types of the diagonal blocks randomly */
        diag_ind = (QInt *)malloc(sizeof(QInt)*dim_block);
        if (diag_ind==NULL) {
            printf("QcMatSetRandMat>> number of diagonal blocks %"QINT_FMT"\n",
                   dim_block);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for diag_ind");
        }
        diag_data_type = (QcDataType *)malloc(sizeof(QcDataType)*dim_block);
        if (diag_data_type==NULL) {
            printf("QcMatSetRandMat>> number of diagonal blocks %"QINT_FMT"\n",
                   dim_block);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for diag_data_type");
        }
        for (iblk=0; iblk<dim_block; iblk++) {
#if defined(QCMATRIX_ZERO_BASED)
            diag_ind[iblk] = iblk;
#else
            diag_ind[iblk] = iblk+1;
#endif
            diag_data_type[iblk] = all_data_types[QRandInt(0,num_data_types)];
        }
        err_code = QcMatSetDataType(A, dim_block, diag_ind, diag_ind, diag_data_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetDataType");
        /* we will not generate QNULLMAT */
        err_code = QcMatIsAssembled(A, &assembled);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatIsAssembled");
        if (assembled==QFALSE) {
            diag_data_type[0] = QCMPLXMAT;
            err_code = QcMatSetDataType(A, 1, diag_ind, diag_ind, diag_data_type);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetDataType");
        }
        free(diag_ind);
        diag_ind = NULL;
        err_code = QcMatAssemble(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatAssemble");
        /* allocates memory for the random values of real and imaginary parts */
        values_real = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (values_real==NULL) {
            printf("QcMatSetRandMat>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for values_real");
        }
        values_imag = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (values_imag==NULL) {
            printf("QcMatSetRandMat>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for values_imag");
        }
        /* generates the values of each off-diagonal block randomly */
#if defined(QCMATRIX_ZERO_BASED)
        for (iblk=0,kblk=0; iblk<dim_block; iblk++) {
            for (jblk=iblk+1; jblk<dim_block; jblk++) {
#else
        for (iblk=1,kblk=0; iblk<=dim_block; iblk++) {
            for (jblk=iblk+1; jblk<=dim_block; jblk++) {
#endif
                /* real part */
                if (block_data_type[kblk]==QREALMAT || block_data_type[kblk]==QCMPLXMAT) {
                    if (sym_type==QSYMMAT) {
/*FIXME: row or column major? */
                        ival = -num_row;
                        for (irow=0; irow<num_row; irow++) {
                            ival += num_row;
                            kval = -num_row+irow;
                            for (icol=0; icol<num_col; icol++) {
                                jval = ival+icol;
                                values_real[jval] = QRandReal(0,1);
                                kval += num_row;
                                values_imag[kval] = values_real[jval];
                            }
                        }
                    }
                    else {
                        ival = -num_col;
                        for (irow=0; irow<num_row; irow++) {
                            ival += num_col;
                            kval = -num_row+irow;
                            for (icol=0; icol<num_col; icol++) {
                                jval = ival+icol;
                                values_real[jval] = QRandReal(0,1);
                                kval += num_row;
                                values_imag[kval] = -values_real[jval];
                            }
                        }
                    }
                    err_code = QcMatSetValues(A,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_set,
                                              idx_first_col,
                                              num_col_set,
                                              values_real,
                                              NULL);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                    err_code = QcMatSetValues(A,
                                              jblk,
                                              iblk,
                                              idx_first_row,
                                              num_row_set,
                                              idx_first_col,
                                              num_col_set,
                                              values_imag,
                                              NULL);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                }
                /* imaginary part */
                if (block_data_type[kblk]==QIMAGMAT || block_data_type[kblk]==QCMPLXMAT) {
                    if (sym_type==QSYMMAT) {
/*FIXME: row or column major? */
                        ival = -num_row;
                        for (irow=0; irow<num_row; irow++) {
                            ival += num_row;
                            kval = -num_row+irow;
                            for (icol=0; icol<num_col; icol++) {
                                jval = ival+icol;
                                values_real[jval] = QRandReal(0,1);
                                kval += num_row;
                                values_imag[kval] = -values_real[jval];
                            }
                        }
                    }
                    else {
                        ival = -num_col;
                        for (irow=0; irow<num_row; irow++) {
                            ival += num_col;
                            kval = -num_row+irow;
                            for (icol=0; icol<num_col; icol++) {
                                jval = ival+icol;
                                values_real[jval] = QRandReal(0,1);
                                kval += num_row;
                                values_imag[kval] = values_real[jval];
                            }
                        }
                    }
                    err_code = QcMatSetValues(A,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_set,
                                              idx_first_col,
                                              num_col_set,
                                              NULL,
                                              values_real);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                    err_code = QcMatSetValues(A,
                                              jblk,
                                              iblk,
                                              idx_first_row,
                                              num_row_set,
                                              idx_first_col,
                                              num_col_set,
                                              NULL,
                                              values_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                }
                kblk++;
            }
        }
        /* generates the values of each diagonal block randomly */
        for (iblk=0; iblk<dim_block; iblk++) {
            /* real part */
            if (diag_data_type[iblk]==QREALMAT || diag_data_type[iblk]==QCMPLXMAT) {
                if (sym_type==QSYMMAT) {
/*FIXME: row or column major? */
                    ival = -num_row;
                    for (irow=0; irow<num_row; irow++) {
                        ival += num_row;
                        values_real[ival+irow] = QRandReal(0,1);
                        kval = irow*(num_row+1);
                        for (icol=irow+1; icol<num_col; icol++) {
                            jval = ival+icol;
                            values_real[jval] = QRandReal(0,1);
                            kval += num_row;
                            values_real[kval] = values_real[jval];
                        }
                    }
                }
                else {
                    ival = -num_col;
                    for (irow=0; irow<num_row; irow++) {
                        ival += num_col;
                        values_real[ival+irow] = 0;
                        kval = irow*(num_col+1);
                        for (icol=irow+1; icol<num_col; icol++) {
                            jval = ival+icol;
                            values_real[jval] = QRandReal(0,1);
/*FIXME: row or column major? */
                            kval += num_row;
                            values_real[kval] = -values_real[jval];
                        }
                    }
                }
            }
            /* imaginary part */
            if (diag_data_type[iblk]==QIMAGMAT || diag_data_type[iblk]==QCMPLXMAT) {
                if (sym_type==QSYMMAT) {
                    ival = -num_row;
                    for (irow=0; irow<num_row; irow++) {
                        ival += num_row;
                        values_imag[ival+irow] = 0;
                        kval = irow*(num_row+1);
                        for (icol=irow+1; icol<num_col; icol++) {
                            jval = ival+icol;
                            values_imag[jval] = QRandReal(0,1);
                            kval += num_row;
                            values_imag[kval] = -values_imag[jval];
                        }
                    }
                }
                else {
                    ival = -num_col;
                    for (irow=0; irow<num_row; irow++) {
                        ival += num_col;
                        values_imag[ival+irow] = QRandReal(0,1);
                        kval = irow*(num_col+1);
                        for (icol=irow+1; icol<num_col; icol++) {
                            jval = ival+icol;
                            values_imag[jval] = QRandReal(0,1);
/*FIXME: row or column major? */
                            kval += num_row;
                            values_imag[kval] = values_imag[jval];
                        }
                    }
                }
            }
            switch (diag_data_type[iblk]) {
            case QNULLMAT:
                break;
            case QREALMAT:
                err_code = QcMatSetValues(A,
#if defined(QCMATRIX_ZERO_BASED)
                                          iblk,
                                          iblk,
#else
                                          iblk+1,
                                          iblk+1,
#endif
                                          idx_first_row,
                                          num_row_set,
                                          idx_first_col,
                                          num_col_set,
                                          values_real,
                                          NULL);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                break;
            case QIMAGMAT:
                err_code = QcMatSetValues(A,
#if defined(QCMATRIX_ZERO_BASED)
                                          iblk,
                                          iblk,
#else
                                          iblk+1,
                                          iblk+1,
#endif
                                          idx_first_row,
                                          num_row_set,
                                          idx_first_col,
                                          num_col_set,
                                          NULL,
                                          values_imag);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                break;
            case QCMPLXMAT:
                err_code = QcMatSetValues(A,
#if defined(QCMATRIX_ZERO_BASED)
                                          iblk,
                                          iblk,
#else
                                          iblk+1,
                                          iblk+1,
#endif
                                          idx_first_row,
                                          num_row_set,
                                          idx_first_col,
                                          num_col_set,
                                          values_real,
                                          values_imag);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
            }
        }
        free(diag_data_type);
        diag_data_type = NULL;
    }
    /* non-Hermitian matrix */
    else if (sym_type==QNONSYMMAT) {
        /* allocates memory for row and column indices and data types */
        num_blocks = dim_block*dim_block;
        idx_block_row = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (idx_block_row==NULL) {
            printf("QcMatSetRandMat>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for idx_block_row");
        }
        idx_block_col = (QInt *)malloc(sizeof(QInt)*num_blocks);
        if (idx_block_col==NULL) {
            printf("QcMatSetRandMat>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for idx_block_col");
        }
        block_data_type = (QcDataType *)malloc(sizeof(QcDataType)*num_blocks);
        if (block_data_type==NULL) {
            printf("QcMatSetRandMat>> number of blocks %"QINT_FMT"\n", num_blocks);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for block_data_type");
        }
        /* generates the data types of the blocks randomly */
#if defined(QCMATRIX_ZERO_BASED)
        for (iblk=0,kblk=0; iblk<dim_block; iblk++) {
            for (jblk=0; jblk<dim_block; jblk++) {
#else
        for (iblk=1,kblk=0; iblk<=dim_block; iblk++) {
            for (jblk=1; jblk<=dim_block; jblk++) {
#endif
                idx_block_row[kblk] = iblk;
                idx_block_col[kblk] = jblk;
                block_data_type[kblk++] = all_data_types[QRandInt(0,num_data_types)];
            }
        }
        err_code = QcMatSetDataType(A,
                                    num_blocks,
                                    idx_block_row,
                                    idx_block_col,
                                    block_data_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetDataType");
        /* we will not generate QNULLMAT */
        err_code = QcMatIsAssembled(A, &assembled);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatIsAssembled");
        if (assembled==QFALSE) {
            block_data_type[0] = QCMPLXMAT;
            err_code = QcMatSetDataType(A,
                                        1,
                                        idx_block_row,
                                        idx_block_col,
                                        block_data_type);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetDataType");
        }
        free(idx_block_row);
        idx_block_row = NULL;
        free(idx_block_col);
        idx_block_col = NULL;
        err_code = QcMatAssemble(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatAssemble");
        /* allocates memory for the random values of real and imaginary parts */
        values_real = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (values_real==NULL) {
            printf("QcMatSetRandMat>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for values_real");
        }
        values_imag = (QReal *)malloc(sizeof(QReal)*size_mat);
        if (values_imag==NULL) {
            printf("QcMatSetRandMat>> size of each block %"QINT_FMT"\n", size_mat);
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for values_imag");
        }
        /* generates the values of each block randomly */
#if defined(QCMATRIX_ZERO_BASED)
        for (iblk=0,kblk=0; iblk<dim_block; iblk++) {
            for (jblk=0; jblk<dim_block; jblk++) {
#else
        for (iblk=1,kblk=0; iblk<=dim_block; iblk++) {
            for (jblk=1; jblk<=dim_block; jblk++) {
#endif
                switch (block_data_type[kblk++]) {
                case QNULLMAT:
                    break;
                case QREALMAT:
                    for (ival=0; ival<size_mat; ival++) {
                        values_real[ival] = QRandReal(0,1);
                    }
                    err_code = QcMatSetValues(A,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_set,
                                              idx_first_col,
                                              num_col_set,
                                              values_real,
                                              NULL);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                    break;
                case QIMAGMAT:
                    for (ival=0; ival<size_mat; ival++) {
                        values_imag[ival] = QRandReal(0,1);
                    }
                    err_code = QcMatSetValues(A,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_set,
                                              idx_first_col,
                                              num_col_set,
                                              NULL,
                                              values_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                    break;
                case QCMPLXMAT:
                    for (ival=0; ival<size_mat; ival++) {
                        values_real[ival] = QRandReal(0,1);
                        values_imag[ival] = QRandReal(0,1);
                    }
                    err_code = QcMatSetValues(A,
                                              iblk,
                                              jblk,
                                              idx_first_row,
                                              num_row_set,
                                              idx_first_col,
                                              num_col_set,
                                              values_real,
                                              values_imag);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatSetValues");
                }
            }
        }
    }
    else {
        printf("QcMatSetRandMat>> input symmetry type: %d\n", sym_type);
        QErrorExit(FILE_AND_LINE, "invalid symmetry type");
    }
    free(all_data_types);
    all_data_types = NULL;
    free(block_data_type);
    block_data_type = NULL;
    free(values_real);
    values_real = NULL;
    free(values_imag);
    values_imag = NULL;
    return QSUCCESS;
}
