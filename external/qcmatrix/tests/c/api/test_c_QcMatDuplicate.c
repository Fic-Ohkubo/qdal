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

   This file tests the function QcMatDuplicate().

   2014-03-25, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"

/*% \brief tests the function QcMatDuplicate()
    \author Bin Gao
    \date 2014-03-25
    \param[QcMat:type]{in} A the matrix
*/
QVoid test_c_QcMatDuplicate(QcMat *A)
{
    QcMat B;          /* duplicated matrix */
    QBool is_equal;   /* indicates if two matrices are equal (pattern and values) */
    QErrorCode ierr;  /* error information */
    /* creates the duplicated matrix first */
    ierr = QcMatCreate(&B);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatDuplicate>> QcMatCreate(B) passed ...\n");
    }
    else {
        printf("test_c_QcMatDuplicate>> QcMatCreate(B) failed\n");
        exit(ierr);
    }
    /* tests QcMatDuplicate() with the option COPY_PATTERN_ONLY */
    ierr = QcMatDuplicate(A, COPY_PATTERN_ONLY, &B);
    if (ierr==QSUCCESS) {
        ierr = QcMatIsEqual(A, &B, QFALSE, &is_equal);
        if (ierr==QSUCCESS) {
            if (is_equal==QTRUE) {
                printf("test_c_QcMatDuplicate>> QcMatDuplicate(A, COPY_PATTERN_ONLY) passed ...\n");
            }
            else {
                printf("test_c_QcMatDuplicate>> QcMatDuplicate(A, COPY_PATTERN_ONLY) failed\n");
                exit(is_equal);
            }
        }
        else {
            printf("test_c_QcMatDuplicate>> failed to call QcMatIsEqual(A, B)\n");
            exit(ierr);
        }
    }
    else {
        printf("test_c_QcMatDuplicate>> failed to call QcMatDuplicate(A, COPY_PATTERN_ONLY)\n");
        exit(ierr);
    }
    /* tests QcMatDuplicate() with the option COPY_PATTERN_AND_VALUE */
    ierr = QcMatDuplicate(A, COPY_PATTERN_AND_VALUE, &B);
    if (ierr==QSUCCESS) {
        ierr = QcMatIsEqual(A, &B, QTRUE, &is_equal);
        if (ierr==QSUCCESS) { 
            if (is_equal==QTRUE) {
                printf("test_c_QcMatDuplicate>> QcMatDuplicate(A, COPY_PATTERN_AND_VALUE) passed ...\n");
            }
            else {
                /* dumps results to check */
#if defined(QCMATRIX_ENABLE_VIEW)
                ierr = QcMatWrite(A, "QcMatDuplicate_A", ASCII_VIEW);
                if (ierr!=QSUCCESS) {
                    printf("test_c_QcMatDuplicate>> failed to call QcMatWrite(A)\n");
                    exit(ierr);
                }
                ierr = QcMatWrite(&B, "QcMatDuplicate_B", ASCII_VIEW);
                if (ierr!=QSUCCESS) {
                    printf("test_c_QcMatDuplicate>> failed to call QcMatWrite(B)\n");
                    exit(ierr);
                }
#endif
                printf("test_c_QcMatDuplicate>> QcMatDuplicate(A, COPY_PATTERN_AND_VALUE) failed\n");
                exit(is_equal);
            }
        }
        else {
            printf("test_c_QcMatDuplicate>> failed to call QcMatIsEqual(A, B)\n");
            exit(ierr);
        }
    }
    else {
        printf("test_c_QcMatDuplicate>> failed to call QcMatDuplicate(A, COPY_PATTERN_AND_VALUE)\n");
        exit(ierr);
    }
    /* frees the space taken by the matrix B */
    ierr = QcMatDestroy(&B);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatDuplicate>> QcMatDestroy(B) failed\n");
        exit(ierr);
    }
}
