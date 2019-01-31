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

   This file implements the function RealMatGEMM().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief performs matrix-matrix multiplication C = alpha*op(A)*op(B)+beta*C,
        where valid operations op(...) can be found in file include/types/mat_operations.h
    \author Bin Gao
    \date 2014-06-16
    \param[QcMatOperation:int]{in} op_A the operation on the matrix A, see file
        include/types/mat_operations.h
    \param[QcMatOperation:int]{in} op_B the operation on the matrix B, see file
        include/types/mat_operations.h
    \param[QReal:real]{in} alpha the scalar number
    \param[RealMat:struct]{in} A the left matrix, should be at least assembled
        by RealMatAssemble()
    \param[RealMat:struct]{in} B the right matrix, should be at least assembled
        by RealMatAssemble()
    \param[QReal:real]{in} beta the scalar number
    \param[RealMat:struct]{inout} C the product matrix, should be at least created
        by RealMatCreate(), so that we require function MatGEMM() could assemble
        the matrix C if it is not
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatGEMM(const QcMatOperation op_A,
                       const QcMatOperation op_B,
                       const QReal alpha,
                       RealMat *A,
                       RealMat *B,
                       const QReal beta,
                       RealMat *C)
{
#if defined(QCMATRIX_ROW_MAJOR)
    RealMat *T;
    QReal real_zero=0;
#endif
    QInt A_num_row, A_num_col;  /* numbers of rows and columns of op(A) */
    QInt B_num_row, B_num_col;  /* numbers of rows and columns of op(B) */
    QChar trans_A;              /* operation on A, for BLAS routine */
    QChar trans_B;              /* operation on B, for BLAS routine */
    QErrorCode err_code;
    if (C==A || C==B) {
        printf("QcMatGEMM>> matrix C should not be the same as the matrix A or B\n");
        QErrorExit(FILE_AND_LINE, "matrix C is the same as the matrix A or B");
    }
    /* if the scalar alpha is zero, we have C = beta*C */
    if (QAbs(alpha)<QZEROTHRSH) {
        err_code = RealMatScale(beta, C);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
        return QSUCCESS;
    }
    /* checks if the matrix A is assembled */
    if (A->values==NULL) {
        QErrorExit(FILE_AND_LINE, "A is not assembled");
    }
    /* checks if the matrix B is assembled */
    if (B->values==NULL) {
        QErrorExit(FILE_AND_LINE, "B is not assembled");
    }
    /* sets the operations for BLAS routine */
    switch (op_A) {
    case MAT_NO_OPERATION:
        A_num_row = A->num_row;
        A_num_col = A->num_col;
        trans_A = 'N';
        break;
    case MAT_TRANSPOSE:
        A_num_row = A->num_col;
        A_num_col = A->num_row;
        trans_A = 'T';
        break;
    default:
        printf("RealMatGEMM>> operation on matrix A: %d\n", op_A);
        QErrorExit(FILE_AND_LINE, "invalid matrix operation");
    }
    switch (op_B) {
    case MAT_NO_OPERATION:
        B_num_row = B->num_row;
        B_num_col = B->num_col;
        trans_B = 'N';
        break;
    case MAT_TRANSPOSE:
        B_num_row = B->num_col;
        B_num_col = B->num_row;
        trans_B = 'T';
        break;
    default:
        printf("RealMatGEMM>> operation on matrix B: %d\n", op_B);
        QErrorExit(FILE_AND_LINE, "invalid matrix operation");
    }
    /* checks the dimension */
    if (A_num_col!=B_num_row) {
        printf("RealMatGEMM>> number of columns of op(A) %"QINT_FMT"\n", A_num_col);
        printf("RealMatGEMM>> number of rows of op(B) %"QINT_FMT"\n", B_num_row);
        QErrorExit(FILE_AND_LINE, "invalid dimensions");
    }
    /* if the scalar beta is zero or C is not assembled, we have C = alpha*op(A)*op(B) */
    if (QAbs(beta)<QZEROTHRSH || C->values==NULL) {
#if defined(QCMATRIX_STORAGE_MODE)
        C->storage_mode = A->storage_mode;
#endif
        C->num_row = A_num_row;
        C->num_col = B_num_col;
        err_code = RealMatAssemble(C);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
    }
    else {
        /* checks the dimension */
        if (A_num_row!=C->num_row) {
            printf("RealMatGEMM>> number of rows of op(A) %"QINT_FMT"\n", A_num_row);
            printf("RealMatGEMM>> number of rows of C %"QINT_FMT"\n", C->num_row);
            QErrorExit(FILE_AND_LINE, "invalid dimensions");
        }
        if (B_num_col!=C->num_col) {
            printf("RealMatGEMM>> number of columns of op(B) %"QINT_FMT"\n", B_num_col);
            printf("RealMatGEMM>> number of columns of C %"QINT_FMT"\n", C->num_col);
            QErrorExit(FILE_AND_LINE, "invalid dimensions");
        }
    }
#if defined(QCMATRIX_ROW_MAJOR)
    if (QAbs(beta)>QZEROTHRSH) {
        /* T = C */
        T = (RealMat *)malloc(sizeof(RealMat));
        if (T==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for T");
        }
        err_code = RealMatCreate(T);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatCreate");
        err_code = RealMatDuplicate(C, COPY_PATTERN_AND_VALUE, T);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
        /* C = alpha*op(B)*op(A), BLAS (column major order) will then calculate
           alpha*op(B^{T})*op(A^{T}) = [alpha*op(A)*op(B)]^{T} */
        C_BLAS_GEMM(trans_B,
                    trans_A,
                    C->num_col,
                    C->num_row,
                    B_num_row,
                    alpha,
                    B->values,
                    B_num_col,
                    A->values,
                    A_num_col,
                    real_zero,
                    C->values,
                    C->num_col);
        /* C = beta*T + C */
        err_code = RealMatAXPY(beta, T, C);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
        /* cleans */
        err_code = RealMatDestroy(T);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDestroy");
        free(T);
        T = NULL;
    }
    else {
        /* C = alpha*op(B)*op(A) */
        C_BLAS_GEMM(trans_B,
                    trans_A,
                    C->num_col,
                    C->num_row,
                    B_num_row,
                    alpha,
                    B->values,
                    B_num_col,
                    A->values,
                    A_num_col,
                    real_zero,
                    C->values,
                    C->num_col);
    }
#else
    C_BLAS_GEMM(trans_A,
                trans_B,
                C->num_row,
                C->num_col,
                A_num_col,
                alpha,
                A->values,
                A_num_row,
                B->values,
                B_num_row,
                beta,
                C->values,
                C->num_row);
#endif
    C->sym_type = QNONSYMMAT;
    return QSUCCESS;
}
