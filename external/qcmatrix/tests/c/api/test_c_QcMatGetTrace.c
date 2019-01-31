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

   This file tests the function QcMatGetTrace().

   2014-03-28, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"

/*% \brief tests the function QcMatGetTrace()
    \author Bin Gao
    \date 2014-03-28
    \param[QcMat:type]{in} A the matrix
*/
QVoid test_c_QcMatGetTrace(QcMat *A)
{
    QInt dim_block;         /* dimension of blocks */
    QInt *idx_block_row;    /* row indices of the blocks */
    QcDataType *data_type;  /* data types of the blocks */
    QInt dim_mat;           /* dimension of each block */
    QInt size_mat;          /* number of elements in each block */
    QInt idx_first_row;     /* index of the first row to get values */
    QInt num_row_get;       /* number of rows to get */
    QInt idx_first_col;     /* index of the first column to get values */
    QInt num_col_get;       /* number of columns to get */
    QReal *values_real;     /* values of the real part */
    QReal *values_imag;     /* values of the imaginary part */
    QReal *cf_trace;        /* trace to compare with */
    QReal *trace;           /* trace from QcMatGetTrace() */
    QInt iblk, jblk, kblk;  /* incremental recorders for blocks */
    QInt ival, jval;        /* incremental recorders for elements of each block */
    QErrorCode ierr;        /* error information */
    /* gets the dimension of blocks */
    ierr = QcMatGetDimBlock(A, &dim_block);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatGetTrace>> QcMatGetDimBlock(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatGetTrace>> failed to call QcMatGetDimBlock(A)\n");
        exit(ierr);
    }
    /* gets the data types of the diagonal blocks */
    idx_block_row = (QInt *)malloc(sizeof(QInt)*dim_block);
    if (idx_block_row==NULL) {
        printf("test_c_QcMatGetTrace>> failed to allocate idx_block_row\n");
        exit(QFAILURE);
    }
    data_type = (QcDataType *)malloc(sizeof(QcDataType)*dim_block);
    if (data_type==NULL) {
        printf("test_c_QcMatGetTrace>> failed to allocate data_type\n");
        exit(QFAILURE);
    }
    for (iblk=0; iblk<dim_block; iblk++) {
#if defined(QCMATRIX_ZERO_BASED)
        idx_block_row[iblk] = iblk;
#else
        idx_block_row[iblk] = iblk+1;
#endif
    }
    ierr = QcMatGetDataType(A, dim_block, idx_block_row, idx_block_row, data_type);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatGetTrace>> failed to call QcMatGetDataType(A)\n");
        exit(ierr);
    }
    /* gets the dimension of each block */
    ierr = QcMatGetDimMat(A, &dim_mat, &dim_mat);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatGetTrace>> QcMatGetDimMat(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatGetTrace>> failed to call QcMatGetDimMat(A)\n");
        exit(ierr);
    }
    size_mat = dim_mat*dim_mat;
    /* allocates memory for getting the elements of each diagonal block */
    values_real = (QReal *)malloc(sizeof(QReal)*size_mat);
    if (values_real==NULL) {
        printf("test_c_QcMatGetTrace>> failed to allocate values_real\n");
        exit(QFAILURE);
    }
    values_imag = (QReal *)malloc(sizeof(QReal)*size_mat);
    if (values_imag==NULL) {
        printf("test_c_QcMatGetTrace>> failed to allocate values_imag\n");
        exit(QFAILURE);
    }
    /* calculates the trace block by block */
    cf_trace = (QReal *)malloc(sizeof(QReal)*2*dim_block);
    if (cf_trace==NULL) {
        printf("test_c_QcMatGetTrace>> failed to allocate cf_trace\n");
        exit(QFAILURE);
    }
#if defined(QCMATRIX_ZERO_BASED)
    idx_first_row = 0;
    idx_first_col = 0;
#else
    idx_first_row = 1;
    idx_first_col = 1;
#endif
    num_row_get = dim_mat;
    num_col_get = dim_mat;
    for (iblk=0; iblk<dim_block; iblk++) {
        jblk = 2*iblk;
        kblk = 2*iblk+1;
        cf_trace[jblk] = 0;  /* real part */
        cf_trace[kblk] = 0;  /* imaginary part */
        switch (data_type[iblk]) {
        /* real block */
        case QREALMAT:
            ierr = QcMatGetValues(A,
                                  idx_block_row[iblk],
                                  idx_block_row[iblk],
                                  idx_first_row,
                                  num_row_get,
                                  idx_first_col,
                                  num_col_get,
                                  values_real,
                                  NULL);
            if (ierr==QSUCCESS) {
                jval = -dim_mat-1;
                for (ival=0; ival<dim_mat; ival++) {
                    jval += dim_mat+1;
                    cf_trace[jblk] += values_real[jval];
                }
            }
            else {
                printf("test_c_QcMatGetTrace>> real block (%"QINT_FMT", %"QINT_FMT")\n",
                       idx_block_row[iblk],
                       idx_block_row[iblk]);
                printf("test_c_QcMatGetTrace>> failed to call QcMatGetValues(A)\n");
                exit(ierr);
            }
            break;
        /* imaginary block */
        case QIMAGMAT:
            ierr = QcMatGetValues(A,
                                  idx_block_row[iblk],
                                  idx_block_row[iblk],
                                  idx_first_row,
                                  num_row_get,
                                  idx_first_col,
                                  num_col_get,
                                  NULL,
                                  values_imag);
            if (ierr==QSUCCESS) {
                jval = -dim_mat-1;
                for (ival=0; ival<dim_mat; ival++) {
                    jval += dim_mat+1;
                    cf_trace[kblk] += values_imag[jval];
                }
            }
            else {
                printf("test_c_QcMatGetTrace>> imaginary block (%"QINT_FMT", %"QINT_FMT")\n",
                       idx_block_row[iblk],
                       idx_block_row[iblk]);
                printf("test_c_QcMatGetTrace>> failed to call QcMatGetValues(A)\n");
                exit(ierr);
            }
            break;
        /* complex block */
        case QCMPLXMAT:
            ierr = QcMatGetValues(A,
                                  idx_block_row[iblk],
                                  idx_block_row[iblk],
                                  idx_first_row,
                                  num_row_get,
                                  idx_first_col,
                                  num_col_get,
                                  values_real,
                                  values_imag);
            if (ierr==QSUCCESS) {
                jval = -dim_mat-1;
                for (ival=0; ival<dim_mat; ival++) {
                    jval += dim_mat+1;
                    cf_trace[jblk] += values_real[jval];
                    cf_trace[kblk] += values_imag[jval];
                }
            }
            else {
                printf("test_c_QcMatGetTrace>> complex block (%"QINT_FMT", %"QINT_FMT")\n",
                       idx_block_row[iblk],
                       idx_block_row[iblk]);
                printf("test_c_QcMatGetTrace>> failed to call QcMatGetValues(A)\n");
                exit(ierr);
            }
            break;
        default:
            break;
        }
    }
    free(idx_block_row);
    idx_block_row = NULL;
    free(data_type);
    data_type = NULL;
    free(values_real);
    values_real = NULL;
    free(values_imag);
    values_imag = NULL;
    /* gets the trace by calling QcMatGetTrace() */
    trace = (QReal *)malloc(sizeof(QReal)*2*dim_block);
    if (trace==NULL) {
        printf("test_c_QcMatGetTrace>> failed to allocate trace\n");
        exit(QFAILURE);
    }
    ierr = QcMatGetTrace(A, dim_block, trace);
    if (ierr==QSUCCESS) {
        for (iblk=0,jblk=0,kblk=1; iblk<dim_block; iblk++) {
            if (QAbs(trace[jblk]-cf_trace[jblk])>CF_THRESHOLD ||
                QAbs(trace[kblk]-cf_trace[kblk])>CF_THRESHOLD) {
                printf("test_c_QcMatGetTrace>> %"QINT_FMT": (%f,%f), (%f,%f)\n",
                       iblk,
                       trace[jblk],
                       trace[kblk],
                       cf_trace[jblk],
                       cf_trace[kblk]);
                printf("test_c_QcMatGetTrace>> QcMatGetTrace(A) failed\n");
                exit(QFAILURE);
            }
            jblk += 2;
            kblk += 2;
        }
        printf("test_c_QcMatGetTrace>> QcMatGetTrace(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatGetTrace>> failed to call QcMatGetTrace(A)\n");
        exit(ierr);
    }
    free(cf_trace);
    cf_trace = NULL;
    free(trace);
    trace = NULL;
}
