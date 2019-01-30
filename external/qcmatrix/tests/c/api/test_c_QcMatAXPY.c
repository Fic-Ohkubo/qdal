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

   This file tests the function QcMatAXPY().

   2014-03-28, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"
/* BLAS routines */
#include "lapack/qcmatrix_c_blas.h"

/*% \brief tests the function QcMatAXPY()
    \author Bin Gao
    \date 2014-03-28
    \param[QcMat:type]{in} X the matrix
    \param[QcMat:type]{inout} Y the matrix
*/
QVoid test_c_QcMatAXPY(QcMat *X, QcMat *Y)
{
    QBool assembled;         /* indicates if the matrix is assembled or not */
    QInt dim_block;          /* dimension of blocks */
    QInt dim_mat;            /* dimension of each block */
    QInt size_values;        /* number of elements in the matrix */
    QReal *X_real;           /* values of the real part of the matrix X */
    QReal *X_imag;           /* values of the imaginary part of the matrix X */
    QReal *Y_real;           /* values of the real part of the matrix Y */
    QReal *Y_imag;           /* values of the imaginary part of the matrix Y */
    QReal multiplier[4][2];  /* the multiplier */
    QBool is_equal;          /* indicates if the matrix and array have the same values */
    QInt idat;               /* incremental recorder over the data types of the scaling number */
    QInt ival;               /* incremental recorder over values */
    QErrorCode ierr;         /* error information */
    /* checks if the matrix is assembled */
    ierr = QcMatIsAssembled(X, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatAXPY>> matrix X is not assembled ...\n");
            printf("test_c_QcMatAXPY>> QcMatAXPY() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatAXPY>> failed to call QcMatIsAssembled(X)\n");
        exit(ierr);
    }
    /* gets the dimension of blocks */
    ierr = QcMatGetDimBlock(X, &dim_block);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatAXPY>> QcMatGetDimBlock(X) passed ...\n");
    }
    else {
        printf("test_c_QcMatAXPY>> failed to call QcMatGetDimBlock(X)\n");
        exit(ierr);
    }
    /* gets the dimension of each block */
    ierr = QcMatGetDimMat(X, &dim_mat, &dim_mat);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatAXPY>> QcMatGetDimMat(X) passed ...\n");
    }
    else {
        printf("test_c_QcMatAXPY>> failed to call QcMatGetDimMat(X)\n");
        exit(ierr);
    }
    size_values = dim_block*dim_block*dim_mat*dim_mat;
    /* allocates memory for the elements of the matrices */
    X_real = (QReal *)malloc(sizeof(QReal)*size_values);
    if (X_real==NULL) {
        printf("test_c_QcMatAXPY>> failed to allocate X_real\n");
        exit(QFAILURE);
    }
    X_imag = (QReal *)malloc(sizeof(QReal)*size_values);
    if (X_imag==NULL) {
        printf("test_c_QcMatAXPY>> failed to allocate X_imag\n");
        exit(QFAILURE);
    }
    Y_real = (QReal *)malloc(sizeof(QReal)*size_values);
    if (Y_real==NULL) {
        printf("test_c_QcMatAXPY>> failed to allocate Y_real\n");
        exit(QFAILURE);
    }
    Y_imag = (QReal *)malloc(sizeof(QReal)*size_values);
    if (Y_imag==NULL) {
        printf("test_c_QcMatAXPY>> failed to allocate Y_imag\n");
        exit(QFAILURE);
    }
    /* gets all the values of the matrices */
    ierr = QcMatGetAllValues(X, QFALSE, size_values, X_real, X_imag);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatAXPY>> failed to call QcMatGetAllValues(X)\n");
        exit(ierr);
    }
    ierr = QcMatGetAllValues(Y, QFALSE, size_values, Y_real, Y_imag);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatAXPY>> failed to call QcMatGetAllValues(Y)\n");
        exit(ierr);
    }
    /* loops over different data types of the multiplier */
    multiplier[0][0] = 0.5; multiplier[0][1] = 0.0;  /* real number */
    multiplier[1][0] = 0.0; multiplier[1][1] = 0.5;  /* imaginary number */
    multiplier[2][0] = 0.5; multiplier[2][1] = 0.5;  /* complex number */
    multiplier[3][0] = 0.0; multiplier[3][1] = 0.0;  /* zero */
    for (idat=0; idat<4; idat++) {
        /* performs Y = a*X+Y by BLAS routine
           (a_{R}+i*a_{I})*(X_{R}+i*X_{I})+(Y_{R}+i*Y_{I})
           = a_{R}*X_{R}-a_{I}*X_{I}+Y_{R}
           + i*(a_{R}*X_{I}+a_{I}*X_{R}+Y_{I}) */
        /* a_{R}*X_{R}+Y_{R} */
        C_BLAS_AXPY(size_values, multiplier[idat][0], X_real, 1, Y_real, 1);
        /* -a_{I}*X_{I}+(a_{R}*X_{R}+Y_{R}) */
        C_BLAS_AXPY(size_values, -multiplier[idat][1], X_imag, 1, Y_real, 1);
        /* a_{R}*X_{I}+Y_{I} */
        C_BLAS_AXPY(size_values, multiplier[idat][0], X_imag, 1, Y_imag, 1);
        /* a_{I}*X_{R}+(a_{R}*X_{I}+Y_{I}) */
        C_BLAS_AXPY(size_values, multiplier[idat][1], X_real, 1, Y_imag, 1);
        /* calls QcMatAXPY() */
        ierr = QcMatAXPY(multiplier[idat], X, Y);
        if (ierr==QSUCCESS) {
            ierr = QcMatCfArray(Y, QFALSE, size_values, Y_real, Y_imag, &is_equal);
            if (ierr==QSUCCESS) {
                if (is_equal==QTRUE) {
                    printf("test_c_QcMatAXPY>> QcMatAXPY(X, Y) passed ...\n");
                }
                else {
                    /* dumps results to check */
                    printf("test_c_QcMatAXPY>> multiplier (%f, %f)\n",
                           multiplier[idat][0],
                           multiplier[idat][1]);
#if defined(QCMATRIX_ENABLE_VIEW)
                    ierr = QcMatWrite(X, "QcMatAXPY_X", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatAXPY>> failed to call QcMatWrite(X)\n");
                        exit(ierr);
                    }
                    ierr = QcMatWrite(Y, "QcMatAXPY_Y", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatAXPY>> failed to call QcMatWrite(Y)\n");
                        exit(ierr);
                    }
#endif
                    printf("test_c_QcMatAXPY>> real part of Y from BLAS\n");
                    for (ival=0; ival<size_values; ival++) {
                        if (ival%4==3) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                            printf("%20.12lf\n", Y_real[ival]);
#else
                            printf("%20.12f\n", Y_real[ival]);
#endif
                        }
                        else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                            printf("%20.12lf  ", Y_real[ival]);
#else
                            printf("%20.12f  ", Y_real[ival]);
#endif
                        }
                    }
                    printf("test_c_QcMatAXPY>> imaginary part of Y from BLAS\n");
                    for (ival=0; ival<size_values; ival++) {
                        if (ival%4==3) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                            printf("%20.12lf\n", Y_imag[ival]);
#else
                            printf("%20.12f\n", Y_imag[ival]);
#endif
                        }
                        else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                            printf("%20.12lf  ", Y_imag[ival]);
#else
                            printf("%20.12f  ", Y_imag[ival]);
#endif
                        }
                    }
                    printf("test_c_QcMatAXPY>> QcMatAXPY(X, Y) failed\n");
                    exit(is_equal);
                }
            }
            else {
                printf("test_c_QcMatAXPY>> failed to call QcMatCfArray(Y)\n");
                exit(ierr);
            }
        }
        else {
            printf("test_c_QcMatAXPY>> failed to call QcMatAXPY(X, Y)\n");
            exit(ierr);
        }
    }
    /* cleans */
    free(X_real);
    X_real = NULL;
    free(X_imag);
    X_imag = NULL;
    free(Y_real);
    Y_real = NULL;
    free(Y_imag);
    Y_imag = NULL;
}
