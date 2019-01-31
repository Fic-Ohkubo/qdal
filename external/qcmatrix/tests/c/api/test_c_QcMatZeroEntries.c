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

   This file tests the function QcMatZeroEntries().

   2014-03-28, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"

/*% \brief tests the function QcMatZeroEntries()
    \author Bin Gao
    \date 2014-03-28
    \param[QcMat:type]{in} A the matrix
*/
QVoid test_c_QcMatZeroEntries(QcMat *A)
{
    QcMat B;             /* duplication of the matrix A */
    QBool assembled;     /* indicates if the matrix is assembled or not */
    QInt dim_block;      /* dimension of blocks */
    QInt dim_mat;        /* dimension of each block */
    QInt size_values;    /* number of elements in the matrix */
    QReal *values_real;  /* values of the real part */
    QBool is_equal;      /* indicates if the matrix and array have the same values */
    QInt ival;           /* incremental recorder over values */
    QErrorCode ierr;     /* error information */
    /* checks if the matrix is assembled */
    ierr = QcMatIsAssembled(A, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatZeroEntries>> matrix A is not assembled ...\n");
            printf("test_c_QcMatZeroEntries>> QcMatZeroEntries() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatZeroEntries>> failed to call QcMatIsAssembled(A)\n");
        exit(ierr);
    }
    /* duplicates the matrix, and uses duplication for the test */
    ierr = QcMatCreate(&B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatZeroEntries>> failed to call QcMatCreate(B)\n");
        exit(ierr);
    }
    ierr = QcMatDuplicate(A, COPY_PATTERN_AND_VALUE, &B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatZeroEntries>> failed to call QcMatDuplicate(A, COPY_PATTERN_AND_VALUE)\n");
        exit(ierr);
    }
    /* gets the dimension of blocks */
    ierr = QcMatGetDimBlock(&B, &dim_block);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatZeroEntries>> QcMatGetDimBlock(B) passed ...\n");
    }
    else {
        printf("test_c_QcMatZeroEntries>> failed to call QcMatGetDimBlock(B)\n");
        exit(ierr);
    }
    /* gets the dimension of each block */
    ierr = QcMatGetDimMat(&B, &dim_mat, &dim_mat);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatZeroEntries>> QcMatGetDimMat(B) passed ...\n");
    }
    else {
        printf("test_c_QcMatZeroEntries>> failed to call QcMatGetDimMat(B)\n");
        exit(ierr);
    }
    size_values = dim_block*dim_block*dim_mat*dim_mat;
    /* allocates memory for the elements of the matrix */
    values_real = (QReal *)malloc(sizeof(QReal)*size_values);
    if (values_real==NULL) {
        printf("test_c_QcMatZeroEntries>> failed to allocate values_real\n");
        exit(QFAILURE);
    }
    for (ival=0; ival<size_values; ival++) {
        values_real[ival] = 0;
    }
    /* zeros all entries of the matrix by QcMatZeroEntries() */
    ierr = QcMatZeroEntries(&B);
    if (ierr==QSUCCESS) {
        ierr = QcMatCfArray(&B,
                            QFALSE,
                            size_values,
                            values_real,
                            values_real,
                            &is_equal);
        if (ierr==QSUCCESS) {
            if (is_equal==QTRUE) {
                printf("test_c_QcMatZeroEntries>> QcMatZeroEntries(B) passed ...\n");
            }
            else {
                /* dumps results to check */
#if defined(QCMATRIX_ENABLE_VIEW)
                ierr = QcMatWrite(&B, "QcMatZeroEntries_B", ASCII_VIEW);
                if (ierr!=QSUCCESS) {
                    printf("test_c_QcMatZeroEntries>> failed to call QcMatWrite(B)\n");
                    exit(ierr);
                }
#endif
                printf("test_c_QcMatZeroEntries>> QcMatZeroEntries(B) failed\n");
                exit(is_equal);
            }
        }
        else {
            printf("test_c_QcMatZeroEntries>> failed to call QcMatCfArray(B)\n");
            exit(ierr);
        }
    }
    else {
        printf("test_c_QcMatZeroEntries>> failed to call QcMatZeroEntries(B)\n");
        exit(ierr);
    }
    /* cleans */
    free(values_real);
    values_real = NULL;
    ierr = QcMatDestroy(&B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatZeroEntries>> QcMatDestroy(B) failed\n");
        exit(ierr);
    }
}
