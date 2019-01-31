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

   This file tests the function QcMatTranspose().

   2014-03-28, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"

/*% \brief tests the function QcMatTranspose()
    \author Bin Gao
    \date 2014-03-28
    \param[QcMat:type]{in} A the matrix
*/
QVoid test_c_QcMatTranspose(QcMat *A)
{
    QcMat B;                           /* duplication of the matrix A */
    QBool assembled;                   /* indicates if the matrix is assembled or not */
    QInt dim_block;                    /* dimension of blocks */
    QInt dim_mat;                      /* dimension of each block */
    QInt size_values;                  /* number of elements in the matrix */
    QReal *values_real;                /* values of the real part */
    QReal *values_imag;                /* values of the imaginary part */
    QReal value_tmp;                   /* temporary value */
    QcMatOperation all_mat_operations[3] = {MAT_TRANSPOSE,
                                            MAT_HERM_TRANSPOSE,
                                            MAT_COMPLEX_CONJUGATE};
    QBool is_equal;                    /* indicates if the matrix and array have the same values */
    QInt iop, irow, icol, ival, jval;  /* incremental recorders */
    QErrorCode ierr;                   /* error information */
    /* checks if the matrix is assembled */
    ierr = QcMatIsAssembled(A, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatTranspose>> matrix A is not assembled ...\n");
            printf("test_c_QcMatTranspose>> QcMatTranspose() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatTranspose>> failed to call QcMatIsAssembled(A)\n");
        exit(ierr);
    }
    /* duplicates the matrix A, and uses the duplication for the test */
    ierr = QcMatCreate(&B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatTranspose>> failed to call QcMatCreate(B)\n");
        exit(ierr);
    }
    ierr = QcMatTranspose(MAT_NO_OPERATION, A, &B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatTranspose>> failed to call QcMatTranspose(MAT_NO_OPERATION, A)\n");
        exit(ierr);
    }
    /* checks if B = A */
    ierr = QcMatIsEqual(A, &B, QTRUE, &is_equal);
    if (ierr==QSUCCESS) {
        if (is_equal==QTRUE) {
            printf("test_c_QcMatTranspose>> QcMatTranspose(MAT_NO_OPERATION, A) passed ...\n");
        }
        else {
            /* dumps results to check */
#if defined(QCMATRIX_ENABLE_VIEW)
            ierr = QcMatWrite(A, "QcMatTranspose_A", ASCII_VIEW);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatTranspose>> failed to call QcMatWrite(A)\n");
                exit(ierr);
            }
            ierr = QcMatWrite(&B, "QcMatTranspose_B", ASCII_VIEW);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatTranspose>> failed to call QcMatWrite(B)\n");
                exit(ierr);
            }
#endif
            printf("test_c_QcMatTranspose>> QcMatTranspose(MAT_NO_OPERATION, A) failed\n");
            exit(is_equal);
        }
    }
    else {
        printf("test_c_QcMatTranspose>> failed to call QcMatIsEqual(A, B)\n");
        exit(ierr);
    }
    /* gets the dimension of blocks */
    ierr = QcMatGetDimBlock(A, &dim_block);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatTranspose>> QcMatGetDimBlock(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatTranspose>> failed to call QcMatGetDimBlock(A)\n");
        exit(ierr);
    }
    /* gets the dimension of each block */
    ierr = QcMatGetDimMat(A, &dim_mat, &dim_mat);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatTranspose>> QcMatGetDimMat(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatTranspose>> failed to call QcMatGetDimMat(A)\n");
        exit(ierr);
    }
    size_values = dim_block*dim_block*dim_mat*dim_mat;
    /* allocates memory for the elements of the matrix */
    values_real = (QReal *)malloc(sizeof(QReal)*size_values);
    if (values_real==NULL) {
        printf("test_c_QcMatTranspose>> failed to allocate values_real\n");
        exit(QFAILURE);
    }
    values_imag = (QReal *)malloc(sizeof(QReal)*size_values);
    if (values_imag==NULL) {
        printf("test_c_QcMatTranspose>> failed to allocate values_imag\n");
        exit(QFAILURE);
    }
    for (iop=0; iop<3; iop++) {
        /* tests QcMatTranspose() out-of-place with option all_mat_operations[iop] */
        ierr = QcMatTranspose(all_mat_operations[iop], A, &B);
        if (ierr==QSUCCESS) {
            /* gets all the values of the matrix */
            ierr = QcMatGetAllValues(A,
                                     QTRUE,
                                     size_values,
                                     values_real,
                                     values_imag);
            if (ierr==QSUCCESS) {
                printf("test_c_QcMatTranspose>> QcMatGetAllValues(A) passed ...\n");
            }
            else {
                printf("test_c_QcMatTranspose>> failed to call QcMatGetAllValues(A)\n");
                exit(ierr);
            }
            /* manually transposes */
            switch (all_mat_operations[iop]) {
            case MAT_TRANSPOSE:
                ival = -dim_block*dim_mat;
                for (irow=0; irow<dim_block*dim_mat; irow++) {
                    ival += dim_block*dim_mat;
                    for (icol=0; icol<irow; icol++) {
                        jval = icol*dim_block*dim_mat+irow;
                        /* real part */
                        value_tmp = values_real[ival+icol];
                        values_real[ival+icol] = values_real[jval];
                        values_real[jval] = value_tmp;
                        /* imaginary part */
                        value_tmp = values_imag[ival+icol];
                        values_imag[ival+icol] = values_imag[jval];
                        values_imag[jval] = value_tmp;
                    }
                }
                break;
            case MAT_HERM_TRANSPOSE:
                ival = -dim_block*dim_mat;
                for (irow=0; irow<dim_block*dim_mat; irow++) {
                    ival += dim_block*dim_mat;
                    for (icol=0; icol<=irow; icol++) {
                        jval = icol*dim_block*dim_mat+irow;
                        /* real part */
                        value_tmp = values_real[ival+icol];
                        values_real[ival+icol] = values_real[jval];
                        values_real[jval] = value_tmp;
                        /* imaginary part */
                        value_tmp = values_imag[ival+icol];
                        values_imag[ival+icol] = -values_imag[jval];
                        values_imag[jval] = -value_tmp;
                    }
                }
                break;
            case MAT_COMPLEX_CONJUGATE:
                for (ival=0; ival<size_values; ival++) {
                    values_imag[ival] = -values_imag[ival];
                }
                break;
            default:
                printf("test_c_QcMatTranspose>> invalid operation on matrix\n");
                exit(all_mat_operations[iop]);
            }
            ierr = QcMatCfArray(&B,
                                QTRUE,
                                size_values,
                                values_real,
                                values_imag,
                                &is_equal);
            if (ierr==QSUCCESS) {
                if (is_equal==QTRUE) {
                    printf("test_c_QcMatTranspose>> QcMatTranspose(A, %d) out-of-place passed ...\n",
                           all_mat_operations[iop]);
                }
                else {
                    /* dumps results to check */
#if defined(QCMATRIX_ENABLE_VIEW)
                    ierr = QcMatWrite(A, "QcMatTranspose_A", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatTranspose>> failed to call QcMatWrite(A)\n");
                        exit(ierr);
                    }
                    ierr = QcMatWrite(&B, "QcMatTranspose_B", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatTranspose>> failed to call QcMatWrite(B)\n");
                        exit(ierr);
                    }
#endif
                    printf("test_c_QcMatTranspose>> real part of B from hand coding\n");
                    for (ival=0; ival<size_values; ival++) {
                        if (ival%4==3) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                            printf("%20.12lf\n", values_real[ival]);
#else
                            printf("%20.12f\n", values_real[ival]);
#endif
                        }
                        else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                            printf("%20.12lf  ", values_real[ival]);
#else
                            printf("%20.12f  ", values_real[ival]);
#endif
                        }
                    }
                    printf("test_c_QcMatTranspose>> imaginary part of B from hand coding\n");
                    for (ival=0; ival<size_values; ival++) {
                        if (ival%4==3) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                            printf("%20.12lf\n", values_imag[ival]);
#else
                            printf("%20.12f\n", values_imag[ival]);
#endif
                        }
                        else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                            printf("%20.12lf  ", values_imag[ival]);
#else
                            printf("%20.12f  ", values_imag[ival]);
#endif
                        }
                    }
                    printf("test_c_QcMatTranspose>> QcMatTranspose(A) out-of-place failed\n");
                    exit(all_mat_operations[iop]);
                }
            }
            else {
                printf("test_c_QcMatTranspose>> failed to call QcMatCfArray(B)\n");
                exit(ierr);
            }
        }
        else {
            printf("test_c_QcMatTranspose>> failed to call QcMatTranspose(A) out-of-place\n");
            exit(all_mat_operations[iop]);
        }
        /* tests QcMatTranspose() in-place with option all_mat_operations[iop] */
        ierr = QcMatTranspose(all_mat_operations[iop], &B, &B);
        if (ierr==QSUCCESS) {
            /* manually transposes */
            switch (all_mat_operations[iop]) {
            case MAT_TRANSPOSE:
                ival = -dim_block*dim_mat;
                for (irow=0; irow<dim_block*dim_mat; irow++) {
                    ival += dim_block*dim_mat;
                    for (icol=0; icol<irow; icol++) {
                        jval = icol*dim_block*dim_mat+irow;
                        /* real part */
                        value_tmp = values_real[ival+icol];
                        values_real[ival+icol] = values_real[jval];
                        values_real[jval] = value_tmp;
                        /* imaginary part */
                        value_tmp = values_imag[ival+icol];
                        values_imag[ival+icol] = values_imag[jval];
                        values_imag[jval] = value_tmp;
                    }
                }
                break;
            case MAT_HERM_TRANSPOSE:
                ival = -dim_block*dim_mat;
                for (irow=0; irow<dim_block*dim_mat; irow++) {
                    ival += dim_block*dim_mat;
                    for (icol=0; icol<=irow; icol++) {
                        jval = icol*dim_block*dim_mat+irow;
                        /* real part */
                        value_tmp = values_real[ival+icol];
                        values_real[ival+icol] = values_real[jval];
                        values_real[jval] = value_tmp;
                        /* imaginary part */
                        value_tmp = values_imag[ival+icol];
                        values_imag[ival+icol] = -values_imag[jval];
                        values_imag[jval] = -value_tmp;
                    }
                }
                break;
            case MAT_COMPLEX_CONJUGATE:
                for (ival=0; ival<size_values; ival++) {
                    values_imag[ival] = -values_imag[ival];
                }
                break;
            default:
                printf("test_c_QcMatTranspose>> invalid operation on matrix\n");
                exit(all_mat_operations[iop]);
            }
            ierr = QcMatCfArray(&B,
                                QTRUE,
                                size_values,
                                values_real,
                                values_imag,
                                &is_equal);
            if (ierr==QSUCCESS) {
                if (is_equal==QTRUE) {
                    printf("test_c_QcMatTranspose>> QcMatTranspose(B, %d) in-place passed ...\n",
                           all_mat_operations[iop]);
                }
                else {
                    /* dumps results to check */
#if defined(QCMATRIX_ENABLE_VIEW)
                    ierr = QcMatWrite(&B, "QcMatTranspose_B", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatTranspose>> failed to call QcMatWrite(B)\n");
                        exit(ierr);
                    }
#endif
                    printf("test_c_QcMatTranspose>> real part of B from hand coding\n");
                    for (ival=0; ival<size_values; ival++) {
                        printf("%e\n",values_real[ival]);
                    }
                    printf("test_c_QcMatTranspose>> imaginary part of B from hand coding\n");
                    for (ival=0; ival<size_values; ival++) {
                        printf("%e\n",values_imag[ival]);
                    }
                    printf("test_c_QcMatTranspose>> QcMatTranspose(B) in-place failed\n");
                    exit(all_mat_operations[iop]);
                }
            }
            else {
                printf("test_c_QcMatTranspose>> failed to call QcMatCfArray(B)\n");
                exit(ierr);
            }
        }
        else {
            printf("test_c_QcMatTranspose>> failed to call QcMatTranspose(B) in-place\n");
            exit(all_mat_operations[iop]);
        }
    }
    /* cleans */
    free(values_real);
    values_real = NULL;
    free(values_imag);
    values_imag = NULL;
    ierr = QcMatDestroy(&B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatTranspose>> QcMatDestroy(B) failed\n");
        exit(ierr);
    }
}
