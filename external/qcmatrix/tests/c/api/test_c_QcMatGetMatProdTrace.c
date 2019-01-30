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

   This file tests the function QcMatGetMatProdTrace().

   2014-06-19, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"

/*% \brief tests the function QcMatGetMatProdTrace()
    \author Bin Gao
    \date 2014-06-19
    \param[QcMat:type]{in} A the matrix
*/
QVoid test_c_QcMatGetMatProdTrace(QcMat *A, QcMat *B)
{
    QBool assembled;                                          /* indicates if the matrix is assembled or not */
    QInt dim_block;                                           /* dimension of blocks */
    QcMat C;                                                  /* product matrix */
    QReal *cf_trace;                                          /* trace to compare with */
    QReal *trace;                                             /* trace from QcMatGetMatProdTrace() */
    QcMatOperation all_mat_operations[4] = {MAT_NO_OPERATION, /* all matrix operations */
                                            MAT_TRANSPOSE,
                                            MAT_HERM_TRANSPOSE,
                                            MAT_COMPLEX_CONJUGATE};
    QReal positive_one[2] = {1.0,0.0};                        /* positive one */
    QReal real_zero[2] = {0.0,0.0};                           /* zero */
    QInt iop, iblk, pos_real, pos_imag;                       /* incremental recorders */
    QErrorCode ierr;                                          /* error information */
    /* checks if the matrices A and B are assembled */
    ierr = QcMatIsAssembled(A, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatGetMatProdTrace>> matrix A is not assembled ...\n");
            printf("test_c_QcMatGetMatProdTrace>> QcMatGetMatProdTrace() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatIsAssembled(A)\n");
        exit(ierr);
    }
    ierr = QcMatIsAssembled(B, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatGetMatProdTrace>> matrix B is not assembled ...\n");
            printf("test_c_QcMatGetMatProdTrace>> QcMatGetMatProdTrace() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatIsAssembled(B)\n");
        exit(ierr);
    }
    ierr = QcMatCreate(&C);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatCreate(C)\n");
        exit(ierr);
    }
    /* gets the dimension of blocks */
    ierr = QcMatGetDimBlock(A, &dim_block);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatGetMatProdTrace>> QcMatGetDimBlock(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatGetDimBlock(A)\n");
        exit(ierr);
    }
    /* allocates memory for the traces */
    cf_trace = (QReal *)malloc(sizeof(QReal)*2*dim_block);
    if (cf_trace==NULL) {
        printf("test_c_QcMatGetMatProdTrace>> failed to allocate cf_trace\n");
        exit(QFAILURE);
    }
    trace = (QReal *)malloc(sizeof(QReal)*2*dim_block);
    if (trace==NULL) {
        printf("test_c_QcMatGetMatProdTrace>> failed to allocate trace\n");
        exit(QFAILURE);
    }
    for (iop=0; iop<4; iop++) {
        printf("test_c_QcMatGetMatProdTrace>> operation on B %d\n",
               all_mat_operations[iop]);
        ierr = QcMatGetMatProdTrace(A, B, all_mat_operations[iop], dim_block, trace);
        if (ierr==QSUCCESS) {
            /* gets the trace by calling QcMatGEMM() and QcMatGetTrace() */
            ierr = QcMatGEMM(MAT_NO_OPERATION,
                             all_mat_operations[iop],
                             positive_one,
                             A,
                             B,
                             real_zero,
                             &C);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatGEMM(A,B)\n");
                exit(ierr);
            }
            ierr = QcMatGetTrace(&C, dim_block, cf_trace);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatGetTrace(C)\n");
                exit(ierr);
            }
            for (iblk=0,pos_real=0,pos_imag=1; iblk<dim_block; iblk++) {
                if (QAbs(trace[pos_real]-cf_trace[pos_real])>CF_THRESHOLD ||
                    QAbs(trace[pos_imag]-cf_trace[pos_imag])>CF_THRESHOLD) {
                    printf("test_c_QcMatGetMatProdTrace>> %"QINT_FMT": (%f,%f), (%f,%f)\n",
                           iblk,
                           trace[pos_real],
                           trace[pos_imag],
                           cf_trace[pos_real],
                           cf_trace[pos_imag]);
#if defined(QCMATRIX_ENABLE_VIEW)
                    ierr = QcMatWrite(A, "QcMatGetMatProdTrace_A", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatWrite(A)\n");
                        exit(ierr);
                    }
                    ierr = QcMatWrite(B, "QcMatGetMatProdTrace_B", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatWrite(B)\n");
                        exit(ierr);
                    }
                    ierr = QcMatWrite(&C, "QcMatGetMatProdTrace_C", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatWrite(C)\n");
                        exit(ierr);
                    }
#endif
                    printf("test_c_QcMatGetMatProdTrace>> QcMatGetMatProdTrace(A,B) failed\n");
                    exit(QFAILURE);
                }
                pos_real += 2;
                pos_imag += 2;
            }
            printf("test_c_QcMatGetMatProdTrace>> QcMatGetMatProdTrace(A,B) passed ...\n");
        }
        else {
            printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatGetMatProdTrace(A,B)\n");
            exit(ierr);
        }
    }
    ierr = QcMatDestroy(&C);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatGetMatProdTrace>> failed to call QcMatDestroy(C)\n");
        exit(ierr);
    }
    free(cf_trace);
    cf_trace = NULL;
    free(trace);
    trace = NULL;
}
