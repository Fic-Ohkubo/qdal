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

   This file is the test suite of C APIs.

   2014-03-13, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"

/* declaration of different test routines */
extern QVoid test_c_QcMatValues(const QInt dim_block);
extern QVoid test_c_QcMatDuplicate(QcMat *A);
extern QVoid test_c_QcMatZeroEntries(QcMat *A);
extern QVoid test_c_QcMatGetTrace(QcMat *A);
extern QVoid test_c_QcMatScale(QcMat *A);
extern QVoid test_c_QcMatAXPY(QcMat *X, QcMat *Y);
extern QVoid test_c_QcMatTranspose(QcMat *A);
extern QVoid test_c_QcMatView(QcMat *A, const QChar *A_label);
extern QVoid test_c_QcMatGEMM(QcMat *A, QcMat *B);
extern QVoid test_c_QcMatGetMatProdTrace(QcMat *A, QcMat *B);

#if defined(QCMATRIX_TEST_EXECUTABLE)
QErrorCode main()
#else
QErrorCode test_c_QcMatrix()
#endif
{
    QcMat A, B;        /* matrices for tests */
    QInt dim_block;    /* dimension of blocks */
    QInt dim_mat = 6;  /* dimension of each block */
    QcSymType sym_type[3] = {QSYMMAT,QANTISYMMAT,QNONSYMMAT};  /* all symmetry types */
    QcDataType data_type[3] = {QREALMAT,QIMAGMAT,QCMPLXMAT};   /* all data types */
#if defined(QCMATRIX_ENABLE_VIEW)
    QChar A_label[6] = "B S D ";
#endif
    QInt isym, jsym;   /* incremental recorders for symmetry types */
    QInt idat, jdat;   /* incremental recorders for data types */
    QInt ierr;         /* error information */
    /* tests different numbers of blocks */
    for (dim_block=1; dim_block<=MAX_DIM_BLOCK; dim_block++) {
        /* tests QcMatSetValues() and QcMatGetValues() */
        test_c_QcMatValues(dim_block);
        /* tests different symmetry types (symmetric, anti-symmetric, non-symmetric) for matrix A */
        for (isym=0; isym<3; isym++) {
            /* tests different data types (real, imaginary, complex) for matrix A */
            for (idat=0; idat<3; idat++) {
                ierr = QcMatCreate(&A);
                if (ierr==QSUCCESS) {
                    printf("test_c_QcMatrix>> QcMatCreate(A) passed ...\n");
                }
                else {
                    printf("test_c_QcMatrix>> failed to call QcMatCreate(A)\n");
                    exit(ierr);
                }
                /* generates a random matrix A according to its symmetry and data types */
                ierr = QcMatSetRandMat(&A,
                                       sym_type[isym],
                                       data_type[idat],
                                       dim_block,
                                       dim_mat,
                                       dim_mat);
                if (ierr==QSUCCESS) {
                    printf("test_c_QcMatrix>> QcMatSetRandMat(A) passed ...\n");
                }
                else {
                    printf("test_c_QcMatrix>> failed to call QcMatSetRandMat(A)\n"); 
                    exit(ierr);
                }
                /* tests QcMatDuplicate() */
                test_c_QcMatDuplicate(&A);
                /* tests QcMatZeroEntries() */
                test_c_QcMatZeroEntries(&A);
                /* tests QcMatGetTrace() */
                test_c_QcMatGetTrace(&A);
                /* tests QcMatScale() */
                test_c_QcMatScale(&A);
                /* tests QcMatAXPY() */
                ierr = QcMatCreate(&B);
                if (ierr==QSUCCESS) {
                    test_c_QcMatAXPY(&A, &B);
                    ierr = QcMatDestroy(&B);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatrix>> failed to call QcMatDestroy(B)\n");
                        exit(ierr);
                    }
                }
                else {
                    printf("test_c_QcMatrix>> failed to call QcMatCreate(B)\n");
                    exit(ierr);
                }
                /* tests QcMatTranspose() in-place and out-of-place */
                test_c_QcMatTranspose(&A);
#if defined(QCMATRIX_ENABLE_VIEW)
                /* tests QcMatWrite() and QcMatRead() */
                A_label[1] = (QChar)(dim_block+'0');
                A_label[3] = (QChar)(isym+1+'0');
                A_label[5] = (QChar)(idat+1+'0');
                test_c_QcMatView(&A, A_label);
#endif
#if defined(ADAPTER_C_LANG)
                /* tests QcMatGetExternalMat() */
                /*test_c_QcMatGetExternalMat(&A);*/
#endif
                /* tests different symmetry types (symmetric, anti-symmetric, non-symmetric) for matrix B */
                for (jsym=0; jsym<3; jsym++) {
                    /* tests different data types (real, imaginary, complex) for matrix B */
                    for (jdat=0; jdat<3; jdat++) {
                        ierr = QcMatCreate(&B);
                        if (ierr==QSUCCESS) {
                            printf("test_c_QcMatrix>> QcMatCreate(B) passed ...\n");
                        }
                        else {
                            printf("test_c_QcMatrix>> failed to call QcMatCreate(B)\n");
                            exit(ierr);
                        }
                        /* generates a random matrix B according to its symmetry and data types */
                        ierr = QcMatSetRandMat(&B,
                                               sym_type[jsym],
                                               data_type[jdat],
                                               dim_block,
                                               dim_mat,
                                               dim_mat);
                        if (ierr==QSUCCESS) {
                            printf("test_c_QcMatrix>> QcMatSetRandMat(B) passed ...\n");
                        }
                        else {
                            printf("test_c_QcMatrix>> failed to call QcMatSetRandMat(B)\n");
                            exit(ierr);
                        }
                        /* tests QcMatAXPY() */
                        test_c_QcMatAXPY(&A, &B);
                        /* tests QcMatGEMM() */
                        test_c_QcMatGEMM(&A, &B);
                        /* tests QcMatGetMatProdTrace() */
                        test_c_QcMatGetMatProdTrace(&A, &B);
                        /* cleans */
                        ierr = QcMatDestroy(&B);
                        if (ierr==QSUCCESS) {
                            printf("test_c_QcMatrix>> QcMatDestroy(B) passed ...\n");
                        }
                        else {
                            printf("test_c_QcMatrix>> failed to call QcMatDestroy(B)\n");
                            exit(ierr);
                        }
                    }
                }
                /* cleans */
                ierr = QcMatDestroy(&A);
                if (ierr==QSUCCESS) {
                    printf("test_c_QcMatrix>> QcMatDestroy(A) passed ...\n");
                }
                else {
                    printf("test_c_QcMatrix>> failed to call QcMatDestroy(A)\n");
                    exit(ierr);
                }
            }
        }
    }
    return QSUCCESS;
}
