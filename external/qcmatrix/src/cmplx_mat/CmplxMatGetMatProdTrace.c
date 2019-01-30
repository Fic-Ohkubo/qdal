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

   This file implements the function CmplxMatGetMatProdTrace().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"

/*% \brief gets the traces of the first few diagonal blocks of a matrix-matrix product A*op(B)
    \author Bin Gao
    \date 2012-04-04
    \param[QcMat:struct]{in} A the left matrix, should be at least assembled
        by QcMatAssemble()
    \param[QcMat:struct]{in} B the right matrix, should be at least assembled
        by QcMatAssemble()
    \param[QcMatOperation:int]{in} op_B the operation on the matrix B, see file
        include/types/mat_operations.h
    \param[QInt:int]{in} num_blocks is the number of diagonal blocks
    \param[QReal:real]{out} trace the traces, size is 2*\var{num_blocks}
    \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetMatProdTrace(QcMat *A,
                                QcMat *B,
                                const QcMatOperation op_B,
                                const QInt num_blocks,
                                QReal *trace)
{
    QcMatOperation op_B_part;  /* operation on the real and imaginary parts of the matrix B */
    QReal sign_B_imag;         /* sign of the imaginary part of the matrix B after operation */
    QReal *part_trace;
    QInt irow, pos_tr;
    QErrorCode err_code;
    if (num_blocks<1) {
        printf("QcMatGetMatProdTrace>> input number of diagonal blocks %"QINT_FMT"\n",
               num_blocks);
        QErrorExit(FILE_AND_LINE, "invalid input number of diagonal blocks");
    }
    /* checks if the matrix B is assembled */
    if (B->data_type==QNULLMAT) {
        printf("QcMatGetMatProdTrace>> data type of matrix B: %d\n", B->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    switch (op_B) {
    case MAT_NO_OPERATION:
        op_B_part = MAT_NO_OPERATION;
        sign_B_imag = 1;
        break;
    case MAT_TRANSPOSE:
        op_B_part = MAT_TRANSPOSE;
        sign_B_imag = 1;
        break;
    case MAT_HERM_TRANSPOSE:
        op_B_part = MAT_TRANSPOSE;
        sign_B_imag = -1;
        break;
    case MAT_COMPLEX_CONJUGATE:
        op_B_part = MAT_NO_OPERATION;
        sign_B_imag = -1;
        break;
    default:
        printf("QcMatGetMatProdTrace>> operation on matrix B: %d\n", op_B);
        QErrorExit(FILE_AND_LINE, "invalid matrix operation");
    }
    /* allocates memory for the traces of real or imaginary part */
    part_trace = (QReal *)malloc(num_blocks*sizeof(QReal));
    if (part_trace==NULL) {
        printf("QcMatGetMatProdTrace>> input number of diagonal blocks %"QINT_FMT"\n",
               num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for part_trace");
    }
    switch (A->data_type) {
    case QREALMAT:
        switch (B->data_type) {
        case QREALMAT:
            /* A_{R}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                trace[pos_tr++] = part_trace[irow];  /* real part */
                trace[pos_tr++] = 0;                 /* imaginary part */
            }
            break;
        case QIMAGMAT:
            /* A_{R}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            if (sign_B_imag==1) {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr++] = 0;                 /* real part */
                    trace[pos_tr++] = part_trace[irow];  /* imaginary part */
                }
            }
            else {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr++] = 0;                  /* real part */
                    trace[pos_tr++] = -part_trace[irow];  /* imaginary part */
                }
            }
            break;
        case QCMPLXMAT:
            /* A_{R}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                trace[pos_tr] = part_trace[irow];
                pos_tr += 2;
            }
            /* A_{R}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            if (sign_B_imag==1) {
                for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                    trace[pos_tr] = part_trace[irow];
                    pos_tr += 2;
                }
            }
            else {
                for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                    trace[pos_tr] = -part_trace[irow];
                    pos_tr += 2;
                }
            }
            break;
        default:
            printf("QcMatGetMatProdTrace>> data type of matrix B: %d\n", B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    case QIMAGMAT:
        switch (B->data_type) {
        case QREALMAT:
            /* A_{I}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                trace[pos_tr++] = 0;                 /* real part */
                trace[pos_tr++] = part_trace[irow];  /* imaginary part */
            }
            break;
        case QIMAGMAT:
            /* A_{I}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            if (sign_B_imag==1) {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr++] = -part_trace[irow];  /* real part */
                    trace[pos_tr++] = 0;                  /* imaginary part */
                }
            }
            else {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr++] = part_trace[irow];  /* real part */
                    trace[pos_tr++] = 0;                 /* imaginary part */
                }
            }
            break;
        case QCMPLXMAT:
            /* A_{I}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            if (sign_B_imag==1) {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr] = -part_trace[irow];
                    pos_tr += 2;
                }
            }
            else {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr] = part_trace[irow];
                    pos_tr += 2;
                }
            }
            /* A_{I}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                trace[pos_tr] = part_trace[irow];
                pos_tr += 2;
            }
            break;
        default:
            printf("QcMatGetMatProdTrace>> data type of matrix B: %d\n", B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    case QCMPLXMAT:
        switch (B->data_type) {
        case QREALMAT:
            /* A_{R}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                trace[pos_tr] = part_trace[irow];
                pos_tr += 2;
            }
            /* A_{I}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                trace[pos_tr] = part_trace[irow];
                pos_tr += 2;
            }
            break;
        case QIMAGMAT:
            /* A_{I}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            if (sign_B_imag==1) {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr] = -part_trace[irow];
                    pos_tr += 2;
                }
            }
            else {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr] = part_trace[irow];
                    pos_tr += 2;
                }
            }
            /* A_{R}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            if (sign_B_imag==1) {
                for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                    trace[pos_tr] = part_trace[irow];
                    pos_tr += 2;
                }
            }
            else {
                for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                    trace[pos_tr] = -part_trace[irow];
                    pos_tr += 2;
                }
            }
            break;
        case QCMPLXMAT:
            /* A_{R}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                trace[pos_tr] = part_trace[irow];
                pos_tr += 2;
            }
            /* A_{I}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                trace[pos_tr] = part_trace[irow];
                pos_tr += 2;
            }
            /* A_{I}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            if (sign_B_imag==1) {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr] += -part_trace[irow];
                    pos_tr += 2;
                }
            }
            else {
                for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
                    trace[pos_tr] += part_trace[irow];
                    pos_tr += 2;
                }
            }
            /* A_{R}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              num_blocks
                                              part_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            if (sign_B_imag==1) {
                for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                    trace[pos_tr] += part_trace[irow];
                    pos_tr += 2;
                }
            }
            else {
                for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
                    trace[pos_tr] += -part_trace[irow];
                    pos_tr += 2;
                }
            }
            break;
        default:
            printf("QcMatGetMatProdTrace>> data type of matrix B: %d\n", B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    default:
        free(part_trace);
        part_trace = NULL;
        printf("QcMatGetMatProdTrace>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    free(part_trace);
    part_trace = NULL;
    return QSUCCESS;
}
#else
#include "impls/cmplx_mat.h"

/*% \brief gets the trace of a matrix-matrix product A*op(B)
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{in} A the left matrix, should be at least assembled
        by CmplxMatAssemble()
    \param[CmplxMat:struct]{in} B the right matrix, should be at least assembled
        by CmplxMatAssemble()
    \param[QcMatOperation:int]{in} op_B the operation on the matrix B, see file
        include/types/mat_operations.h
    \param[QReal:real]{out} trace the trace
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatGetMatProdTrace(CmplxMat *A,
                                   CmplxMat *B,
                                   const QcMatOperation op_B,
                                   QReal *trace)
{
    QcMatOperation op_B_part;  /* operation on the real and imaginary parts of the matrix B */
    QReal sign_B_imag;        /* sign of the imaginary part of the matrix B after operation */
    QReal tmp_trace;
    QErrorCode err_code;
    switch (op_B) {
    case MAT_NO_OPERATION:
        op_B_part = MAT_NO_OPERATION;
        sign_B_imag = 1;
        break;
    case MAT_TRANSPOSE:
        op_B_part = MAT_TRANSPOSE;
        sign_B_imag = 1;
        break;
    case MAT_HERM_TRANSPOSE:
        op_B_part = MAT_TRANSPOSE;
        sign_B_imag = -1;
        break;
    case MAT_COMPLEX_CONJUGATE:
        op_B_part = MAT_NO_OPERATION;
        sign_B_imag = -1;
        break;
    default:
        printf("CmplxMatGetMatProdTrace>> operation on matrix B: %d\n", op_B);
        QErrorExit(FILE_AND_LINE, "invalid matrix operation");
    }
    switch (A->data_type) {
    case QREALMAT:
        switch (B->data_type) {
        case QREALMAT:
            /* A_{R}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              &trace[0]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[1] = 0;
            break;
        case QIMAGMAT:
            trace[0] = 0;
            /* A_{R}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              &tmp_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[1] = (sign_B_imag==1)?tmp_trace:-tmp_trace;
            break;
        case QCMPLXMAT:
            /* A_{R}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              &trace[0]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            /* A_{R}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              &tmp_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[1] = (sign_B_imag==1)?tmp_trace:-tmp_trace;
            break;
        default:
            printf("CmplxMatGetMatProdTrace>> data type of matrix B: %d\n",
                   B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    case QIMAGMAT:
        switch (B->data_type) {
        case QREALMAT:
            trace[0] = 0;
            /* A_{I}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              &trace[1]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            break;
        case QIMAGMAT:
            /* A_{I}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              &tmp_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[0] = (sign_B_imag==1)?-tmp_trace:tmp_trace;
            trace[1] = 0;
            break;
        case QCMPLXMAT:
            /* A_{I}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              &tmp_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[0] = (sign_B_imag==1)?-tmp_trace:tmp_trace;
            /* A_{I}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              &trace[1]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            break;
        default:
            printf("CmplxMatGetMatProdTrace>> data type of matrix B: %d\n",
                   B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    case QCMPLXMAT:
        switch (B->data_type) {
        case QREALMAT:
            /* A_{R}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              &trace[0]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            /* A_{I}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              &trace[1]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            break;
        case QIMAGMAT:
            /* A_{I}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              &tmp_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[0] = (sign_B_imag==1)?-tmp_trace:tmp_trace;
            /* A_{R}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              &tmp_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[1] = (sign_B_imag==1)?tmp_trace:-tmp_trace;
            break;
        case QCMPLXMAT:
            /* A_{R}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              &trace[0]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            /* A_{I}*op(B_{R}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->real_part],
                                              op_B_part,
                                              &trace[1]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            /* A_{I}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->imag_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              &tmp_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[0] += (sign_B_imag==1)?-tmp_trace:tmp_trace;
            /* A_{R}*op(B_{I}) */
            err_code = RealMatGetMatProdTrace(&A->cmplx_mat[A->real_part],
                                              &B->cmplx_mat[B->imag_part],
                                              op_B_part,
                                              &tmp_trace);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetMatProdTrace");
            trace[1] += (sign_B_imag==1)?tmp_trace:-tmp_trace;
            break;
        default:
            printf("CmplxMatGetMatProdTrace>> data type of matrix B: %d\n",
                   B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    default:
        printf("CmplxMatGetMatProdTrace>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    return QSUCCESS;
}
#endif  /* defined(ADAPTER_BLOCK_REAL) */
