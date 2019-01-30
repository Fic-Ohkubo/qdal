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

   This file tests the function QcMatScale().

   2014-03-28, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"
/* BLAS routines */
#include "lapack/qcmatrix_c_blas.h"

/*% \brief tests the function QcMatScale()
    \author Bin Gao
    \date 2014-03-28
    \param[QcMat:type]{in} A the matrix
*/
QVoid test_c_QcMatScale(QcMat *A)
{
    QcMat B;                  /* duplication of the matrix A */
    QBool assembled;          /* indicates if the matrix is assembled or not */
    QInt dim_block;           /* dimension of blocks */
    QInt dim_mat;             /* dimension of each block */
    QInt size_values;         /* number of elements in the matrix */
    QReal *values_real;       /* values of the real part */
    QReal *values_imag;       /* values of the imaginary part */
    QReal *values_tmp;        /* temporary values */
    QReal scal_number[4][2];  /* scaling number */
    QBool is_equal;           /* indicates if the matrix and array have the same values */
    QInt idat;                /* incremental recorder over the data types of the scaling number */
    QInt ival;                /* incremental recorder over values */
    QErrorCode ierr;          /* error information */
    /* checks if the matrix is assembled */
    ierr = QcMatIsAssembled(A, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatScale>> matrix A is not assembled ...\n");
            printf("test_c_QcMatScale>> QcMatScale() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatScale>> failed to call QcMatIsAssembled(A)\n");
        exit(ierr);
    }
    /* duplicates the matrix A, and uses the duplication for the test */
    ierr = QcMatCreate(&B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatScale>> failed to call QcMatCreate(B)\n");
        exit(ierr);
    }
    ierr = QcMatDuplicate(A, COPY_PATTERN_AND_VALUE, &B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatScale>> failed to call QcMatDuplicate(A, COPY_PATTERN_AND_VALUE)\n");
        exit(ierr);
    }
    /* gets the dimension of blocks */
    ierr = QcMatGetDimBlock(&B, &dim_block);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatScale>> QcMatGetDimBlock(B) passed ...\n");
    }
    else {
        printf("test_c_QcMatScale>> failed to call QcMatGetDimBlock(B)\n");
        exit(ierr);
    }
    /* gets the dimension of each block */
    ierr = QcMatGetDimMat(&B, &dim_mat, &dim_mat);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatScale>> QcMatGetDimMat(B) passed ...\n");
    }
    else {
        printf("test_c_QcMatScale>> failed to call QcMatGetDimMat(B)\n");
        exit(ierr);
    }
    size_values = dim_block*dim_block*dim_mat*dim_mat;
    /* allocates memory for the elements of the matrix */
    values_real = (QReal *)malloc(sizeof(QReal)*size_values);
    if (values_real==NULL) {
        printf("test_c_QcMatScale>> failed to allocate values_real\n");
        exit(QFAILURE);
    }
    values_imag = (QReal *)malloc(sizeof(QReal)*size_values);
    if (values_imag==NULL) {
        printf("test_c_QcMatScale>> failed to allocate values_imag\n");
        exit(QFAILURE);
    }
    values_tmp = (QReal *)malloc(sizeof(QReal)*size_values);
    if (values_tmp==NULL) {
        printf("test_c_QcMatScale>> failed to allocate values_tmp\n");
        exit(QFAILURE);
    }
    /* gets all the values of the matrix */
    ierr = QcMatGetAllValues(&B, QFALSE, size_values, values_real, values_imag);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatScale>> failed to call QcMatGetAllValues(B)\n");
        exit(ierr);
    }
    /* loops over different data types of the scaling number */
    scal_number[0][0] = 0.5; scal_number[0][1] = 0.0;  /* real number */
    scal_number[1][0] = 0.0; scal_number[1][1] = 0.5;  /* imaginary number */
    scal_number[2][0] = 0.5; scal_number[2][1] = 0.5;  /* complex number */
    scal_number[3][0] = 0.0; scal_number[3][1] = 0.0;  /* zero */
    for (idat=0; idat<4; idat++) {
        /* scales the matrix by BLAS routine
           (a_{R}+i*a_{I})*(A_{R}+i*A_{I})
           = a_{R}*A_{R}-a_{I}*A_{I}
           + i*(a_{R}*A_{I}+a_{I}*A_{R}) */
        C_BLAS_COPY(size_values, values_real, 1, values_tmp, 1);
        /* a_{R}*A_{R} */
        C_BLAS_SCAL(size_values, scal_number[idat][0], values_real, 1);
        /* -a_{I}*A_{I}+(a_{R}*A_{R}) */
        C_BLAS_AXPY(size_values,
                    -scal_number[idat][1],
                    values_imag,
                    1,
                    values_real,
                    1);
        /* a_{R}*A_{I} */
        C_BLAS_SCAL(size_values, scal_number[idat][0], values_imag, 1);
        /* a_{I}*A_{R}+(a_{R}*A_{I}) */
        C_BLAS_AXPY(size_values,
                    scal_number[idat][1],
                    values_tmp,
                    1,
                    values_imag,
                    1);
        /* scales the matrix by QcMatScale() */
        ierr = QcMatScale(scal_number[idat], &B);
        if (ierr==QSUCCESS) {
            ierr = QcMatCfArray(&B,
                                QFALSE,
                                size_values,
                                values_real,
                                values_imag,
                                &is_equal);
            if (ierr==QSUCCESS) {
                if (is_equal==QTRUE) {
                    printf("test_c_QcMatScale>> QcMatScale(B) passed ...\n");
                }
                else {
                    /* dumps results to check */
                    printf("test_c_QcMatScale>> scaling number (%f, %f)\n",
                           scal_number[idat][0],
                           scal_number[idat][1]);
#if defined(QCMATRIX_ENABLE_VIEW)
                    ierr = QcMatWrite(A, "QcMatScale_A", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatScale>> failed to call QcMatWrite(A)\n");
                        exit(ierr);
                    }
                    ierr = QcMatWrite(&B, "QcMatScale_B", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatScale>> failed to call QcMatWrite(B)\n");
                        exit(ierr);
                    }
#endif
                    printf("test_c_QcMatScale>> real part of B from BLAS\n");
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
                    printf("test_c_QcMatScale>> imaginary part of B from BLAS\n");
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
                    printf("test_c_QcMatScale>> QcMatScale(B) failed\n");
                    exit(is_equal);
                }
            }
            else {
                printf("test_c_QcMatScale>> failed to call QcMatCfArray(B)\n");
                exit(ierr);
            }
        }
        else {
            printf("test_c_QcMatScale>> failed to call QcMatScale(B)\n");
            exit(ierr);
        }
    }
    /* cleans */
    free(values_real);
    values_real = NULL;
    free(values_imag);
    values_imag = NULL;
    free(values_tmp);
    values_tmp = NULL;
    ierr = QcMatDestroy(&B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatScale>> QcMatDestroy(B) failed\n");
        exit(ierr);
    }
}
