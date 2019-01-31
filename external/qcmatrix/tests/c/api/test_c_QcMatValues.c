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

   This file tests the functions QcMatSetValues(), QcMatAddValues() and QcMatGetValues().

   2014-08-22, Bin Gao:
   * adds the test of QcMatAddValues()

   2014-03-25, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"

/*% \brief tests the functions QcMatSetValues(), QcMatAddValues() and QcMatGetValues()
    \author Bin Gao
    \date 2014-03-25
    \param[QInt]{in} dim_block the dimension of blocks
*/
QVoid test_c_QcMatValues(const QInt dim_block)
{
    QcMat A;                                  /* matrix for test */
    QInt num_blocks;                          /* number of blocks, as \var{dim_block}*\var{dim_block} */
    QInt *idx_block_row;                      /* row indices of the blocks */
    QInt *idx_block_col;                      /* column indices of the blocks */
    QcDataType *data_type;                    /* data types of the blocks */
    const QInt dim_mat=6;                     /* dimension of each block */
    QInt size_mat;                            /* number of elements in each block */
    QInt idx_first_row;                       /* index of the first row from which the values are set/added */
    QInt num_row_set=dim_mat;                 /* number of rows to set/add */
    QInt idx_first_col;                       /* index of the first column from which the values are set/added */
    QInt num_col_set=dim_mat;                 /* number of columns to set/add */
    QReal *values_real;                       /* values of the real part */
    QBool is_equal;                           /* indicates if the matrix and array have the same values */
    QInt iblk, jblk, kblk;                    /* incremental recorders for blocks */
    QErrorCode ierr;                          /* error information */
    ierr = QcMatCreate(&A);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatValues>> QcMatCreate(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatValues>> failed to call QcMatCreate(A)\n");
        exit(ierr);
    }
    ierr = QcMatBlockCreate(&A, dim_block);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatValues>> QcMatBlockCreate(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatValues>> failed to call QcMatBlockCreate(A)\n");
        exit(ierr);
    }
    /* sets the data types of the blocks */
    num_blocks = dim_block*dim_block;
    idx_block_row = (QInt *)malloc(sizeof(QInt)*num_blocks);
    if (idx_block_row==NULL) {
        printf("test_c_QcMatValues>> failed to allocate idx_block_row\n");
        exit(QFAILURE);
    }
    idx_block_col = (QInt *)malloc(sizeof(QInt)*num_blocks);
    if (idx_block_col==NULL) {
        printf("test_c_QcMatValues>> failed to allocate idx_block_col\n");
        exit(QFAILURE);
    }
    data_type = (QcDataType *)malloc(sizeof(QcDataType)*num_blocks);
    if (data_type==NULL) {
        printf("test_c_QcMatValues>> failed to allocate data_type\n");
        exit(QFAILURE);
    }
    kblk = 0;
#if defined(QCMATRIX_ZERO_BASED)
    for (iblk=0; iblk<dim_block; iblk++) {
        for (jblk=0; jblk<dim_block; jblk++) {
#else
    for (iblk=1; iblk<=dim_block; iblk++) {
        for (jblk=1; jblk<=dim_block; jblk++) {
#endif
            /* QcMatrix uses row major order for the blocks */
            idx_block_row[kblk] = iblk;
            idx_block_col[kblk] = jblk;
            data_type[kblk++] = QCMPLXMAT;
        }
    }
    ierr = QcMatSetDataType(&A, num_blocks, idx_block_row, idx_block_col, data_type);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatValues>> failed to call QcMatSetDataType(A)\n");
        exit(ierr);
    }
    free(data_type);
    data_type = NULL;
    /* sets the dimension of each block */
    ierr = QcMatSetDimMat(&A, dim_mat, dim_mat);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatValues>> failed to call QcMatSetDimMat(A)\n");
        exit(ierr);
    }
    /* assembles the matrix */
    ierr = QcMatAssemble(&A);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatValues>> failed to call QcMatAssemble(A)\n");
        exit(ierr);
    }
    /* allocates memory for setting the elements of each block */
    size_mat = dim_mat*dim_mat;
    values_real = (QReal *)malloc(sizeof(QReal)*num_blocks*size_mat);
    if (values_real==NULL) {
        printf("test_c_QcMatValues>> failed to allocate values_real\n");
        exit(QFAILURE);
    }
    /* we sets all the elements as 1+i */
    for (iblk=0; iblk<num_blocks*size_mat; iblk++) {
        values_real[iblk] = 1;
    }
    /* sets the elements block by block */
#if defined(QCMATRIX_ZERO_BASED)
    idx_first_row = 0;
    idx_first_col = 0;
#else
    idx_first_row = 1;
    idx_first_col = 1;
#endif
    for (iblk=0; iblk<num_blocks; iblk++) {
        ierr = QcMatSetValues(&A,
                              idx_block_row[iblk],
                              idx_block_col[iblk],
                              idx_first_row,
                              num_row_set,
                              idx_first_col,
                              num_col_set,
                              values_real,
                              values_real);
        if (ierr!=QSUCCESS) {
            printf("test_c_QcMatValues>> block (%"QINT_FMT", %"QINT_FMT")\n",
                   idx_block_row[iblk],
                   idx_block_col[iblk]);
            printf("test_c_QcMatValues>> failed to call QcMatSetValues(A)\n");
            exit(ierr);
        }
    }
    /* checks the elements set by QcMatSetValues() */
    ierr = QcMatCfArray(&A,
                        QFALSE,
                        num_blocks*size_mat,
                        values_real,
                        values_real,
                        &is_equal);
    if (ierr==QSUCCESS) {
        if (is_equal==QTRUE) {
            printf("test_c_QcMatValues>> QcMatSetValues(A) passed ...\n");
        }
        else {
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite(&A, "QcMatValues_A", ASCII_VIEW);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatValues>> failed to call QcMatWrite(A)\n");
                exit(ierr);
            }
#endif
            printf("test_c_QcMatValues>> QcMatSetValues(A) failed\n");
            exit(is_equal);
        }
    }
    else {
        printf("test_c_QcMatValues>> failed to call QcMatCfArray(A)\n");
        exit(ierr);
    }
    /* adds values to the matrix */
    for (iblk=0; iblk<num_blocks; iblk++) {
        ierr = QcMatAddValues(&A,
                              idx_block_row[iblk],
                              idx_block_col[iblk],
                              idx_first_row,
                              num_row_set,
                              idx_first_col,
                              num_col_set,
                              values_real,
                              values_real);
        if (ierr!=QSUCCESS) {
            printf("test_c_QcMatValues>> block (%"QINT_FMT", %"QINT_FMT")\n",
                   idx_block_row[iblk],
                   idx_block_col[iblk]);
            printf("test_c_QcMatValues>> failed to call QcMatAddValues(A)\n");
            exit(ierr);
        }
    }
    /* all the elements become 2+2*i */
    for (iblk=0; iblk<num_blocks*size_mat; iblk++) {
        values_real[iblk] = 2;
    }
    /* checks the elements set by QcMatAddValues() */
    ierr = QcMatCfArray(&A,
                        QFALSE,
                        num_blocks*size_mat,
                        values_real,
                        values_real,
                        &is_equal);
    if (ierr==QSUCCESS) {
        if (is_equal==QTRUE) {
            printf("test_c_QcMatValues>> QcMatAddValues(A) passed ...\n");
        }
        else {
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite(&A, "QcMatValues_A", ASCII_VIEW);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatValues>> failed to call QcMatWrite(A)\n");
                exit(ierr);
            }
#endif
            printf("test_c_QcMatValues>> QcMatAddValues(A) failed\n");
            exit(is_equal);
        }
    }
    else {
        printf("test_c_QcMatValues>> failed to call QcMatCfArray(A)\n");
        exit(ierr);
    }
    /* cleans up */
    free(idx_block_row);
    idx_block_row = NULL;
    free(idx_block_col);
    idx_block_col = NULL;
    free(values_real);
    values_real = NULL;
    ierr = QcMatDestroy(&A);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatValues>> failed to call QcMatDestroy(A)\n");
        exit(ierr);
    }
}
