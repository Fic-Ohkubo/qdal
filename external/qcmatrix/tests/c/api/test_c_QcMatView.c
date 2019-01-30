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

   This file tests the functions QcMatWrite() and QcMatRead().

   2014-03-25, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"

/*% \brief tests the functions QcMatWrite() and QcMatRead()
    \author Bin Gao
    \date 2014-03-25
    \param[QcMat:type]{in} A the matrix
    \param[character]{in} A_label the label of the matrix
*/
QVoid test_c_QcMatView(QcMat *A, const QChar *A_label)
{
    QBool assembled;  /* indicates if the matrix is assembled or not */
    QcMat B;          /* matrix for reading */
    QBool is_equal;   /* indicates if two matrices are equal (pattern and values) */
    QErrorCode ierr;  /* error information */
    ierr = QcMatIsAssembled(A, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatView>> matrix A is not assembled ...\n");
            printf("test_c_QcMatView>> QcMatWrite() and QcMatRead() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatView>> failed to call QcMatIsAssembled(A)\n");
        exit(ierr);
    }
    /* the values written using option ASCII_VIEW are normally not accurate,
       also we normally save matrix using BINARY_VIEW, so we skip the test
       using ASCII_VIEW */
    //-/* tests QcMatWrite() with option ASCII_VIEW */
    //-ierr = QcMatWrite(A, A_label, ASCII_VIEW);
    //-if (ierr==QSUCCESS) {
    //-    printf("test_c_QcMatView>> QcMatWrite(A, ASCII_VIEW) passed ...\n");
    //-}
    //-else {
    //-    printf("test_c_QcMatView>> QcMatWrite(A, ASCII_VIEW) failed\n");
    //-    exit(ierr);
    //-}
    //-/* creates the matrix B for reading */
    //-ierr = QcMatCreate(&B);
    //-if (ierr==QSUCCESS) {
    //-    printf("test_c_QcMatView>> QcMatCreate(B) passed ...\n");
    //-}
    //-else {
    //-    printf("test_c_QcMatView>> QcMatCreate(B) failed\n");
    //-    exit(ierr);
    //-}
    //-/* tests QcMatRead() with option ASCII_VIEW */
    //-ierr = QcMatRead(&B,  A_label, ASCII_VIEW);
    //-if (ierr==QSUCCESS) {
    //-    /* tests if we read in exactly the same matrix A */
    //-    ierr = QcMatIsEqual(A, &B, QTRUE, &is_equal);
    //-    if (ierr==QSUCCESS) {
    //-        if (is_equal==QTRUE) {
    //-            printf("test_c_QcMatView>> QcMatRead(B, ASCII_VIEW) passed ...\n");
    //-        }
    //-        else {
    //-            /* dumps results to check */
    //-            ierr = QcMatWrite(A, "QcMatView_A", ASCII_VIEW);
    //-            if (ierr!=QSUCCESS) {
    //-                printf("test_c_QcMatView>> failed to call QcMatWrite(A)\n");
    //-                exit(ierr);
    //-            }
    //-            ierr = QcMatWrite(&B, "QcMatView_B", ASCII_VIEW);
    //-            if (ierr!=QSUCCESS) {
    //-                printf("test_c_QcMatView>> failed to call QcMatWrite(B)\n");
    //-                exit(ierr);
    //-            }
    //-            printf("test_c_QcMatView>> QcMatRead(B, ASCII_VIEW) failed\n");
    //-            exit(is_equal);
    //-        }
    //-    }
    //-    else {
    //-        printf("test_c_QcMatView>> QcMatIsEqual(ASCII_VIEW) failed\n");
    //-        exit(ierr);
    //-    }
    //-}
    //-else {
    //-    printf("test_c_QcMatView>> failed to call QcMatRead(B, ASCII_VIEW)\n");
    //-    exit(ierr);
    //-}
    //-/* frees the space taken by the matrix B */
    //-ierr = QcMatDestroy(&B);
    //-if (ierr==QSUCCESS) {
    //-    printf("test_c_QcMatView>> QcMatDestroy(B) passed ...\n");
    //-}
    //-else {
    //-    printf("test_c_QcMatView>> QcMatDestroy(B) failed\n");
    //-    exit(ierr);
    //-}
    /* tests QcMatWrite() with option BINARY_VIEW */
    ierr = QcMatWrite(A, A_label, BINARY_VIEW);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatView>> QcMatWrite(A, BINARY_VIEW) passed ...\n");
    }
    else {
        printf("test_c_QcMatView>> QcMatWrite(A, BINARY_VIEW) failed\n");
        exit(ierr);
    }
    /* creates the matrix B for reading */
    ierr = QcMatCreate(&B);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatView>> QcMatCreate(B) passed ...\n");
    }
    else {
        printf("test_c_QcMatView>> QcMatCreate(B) failed\n");
        exit(ierr);
    }
    /* tests QcMatRead() with option BINARY_VIEW */
    ierr = QcMatRead(&B,  A_label, BINARY_VIEW);
    if (ierr==QSUCCESS) {
        /* tests if we read in exactly the same matrix A */
        ierr = QcMatIsEqual(A, &B, QTRUE, &is_equal);
        if (ierr==QSUCCESS) {
            if (is_equal==QTRUE) {
                printf("test_c_QcMatView>> QcMatRead(B, BINARY_VIEW) passed ...\n");
            }
            else {
                /* dumps results to check */
                ierr = QcMatWrite(A, "QcMatView_A", ASCII_VIEW);
                if (ierr!=QSUCCESS) {
                    printf("test_c_QcMatView>> failed to call QcMatWrite(A)\n");
                    exit(ierr);
                }
                ierr = QcMatWrite(&B, "QcMatView_B", ASCII_VIEW);
                if (ierr!=QSUCCESS) {
                    printf("test_c_QcMatView>> failed to call QcMatWrite(B)\n");
                    exit(ierr);
                }
                printf("test_c_QcMatView>> QcMatRead(B, BINARY_VIEW) failed\n");
                exit(is_equal);
            }
        }
        else {
            printf("test_c_QcMatView>> QcMatIsEqual(BINARY_VIEW) failed\n");
            exit(ierr);
        }
    }
    else {
        printf("test_c_QcMatView>> failed to call QcMatRead(B, BINARY_VIEW)\n");
        exit(ierr);
    }
    /* frees the space taken by the matrix B */
    ierr = QcMatDestroy(&B);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatView>> QcMatDestroy(B) passed ...\n");
    }
    else {
        printf("test_c_QcMatView>> QcMatDestroy(B) failed\n");
        exit(ierr);
    }
}
