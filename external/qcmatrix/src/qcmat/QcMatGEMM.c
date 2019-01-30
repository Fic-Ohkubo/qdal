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

   This file implements the function QcMatGEMM().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*@% \brief performs matrix-matrix multiplication C = alpha*op(A)*op(B)+beta*C,
         where valid operations op(...) can be found in file include/types/mat_operations.h
     \author Bin Gao
     \date 2012-04-04
     \param[QcMatOperation:int]{in} op_A the operation on the matrix A, see file
         include/types/mat_operations.h
     \param[QcMatOperation:int]{in} op_B the operation on the matrix B, see file
         include/types/mat_operations.h
     \param[QReal:real]{in} alpha the scalar number
     \param[QcMat:struct]{in} A the left matrix, should be at least assembled
         by QcMatAssemble()
     \param[QcMat:struct]{in} B the right matrix, should be at least assembled
         by QcMatAssemble()
     \param[QReal:real]{in} beta the scalar number
     \param[QcMat:struct]{inout} C the product matrix, should be at least created
         by QcMatCreate(), so that we require function CmplxMatGEMM() could assemble
         the matrix C if it is not
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGEMM(const QcMatOperation op_A,
                     const QcMatOperation op_B,
                     const QReal alpha[],
                     QcMat *A,
                     QcMat *B,
                     const QReal beta[],
                     QcMat *C)
{
    QBool A_assembled;
    QBool B_assembled;
    QBool C_assembled;
    QReal positive_one[2]={1,0};
    QInt irow, jcol, krow;
    QErrorCode err_code;
    if (C==A || C==B) {
        printf("QcMatGEMM>> matrix C should not be the same as the matrix A or B\n");
        QErrorExit(FILE_AND_LINE, "matrix C is the same as the matrix A or B");
    }
    /* if the scalar alpha is zero, we have C = beta*C */
    if (QAbs(alpha[0])<QZEROTHRSH && QAbs(alpha[1])<QZEROTHRSH) {
        err_code = QcMatScale(beta, C);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatScale");
        return QSUCCESS;
    }
    /* checks if some of the blocks of the matrix A is assembled */
    if (A->blocks==NULL) {
        QErrorExit(FILE_AND_LINE, "blocks of the matrix A is not created");
    }
    else {
        err_code = QcMatIsAssembled(A, &A_assembled);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatIsAssembled(A)");
        if (A_assembled==QFALSE) {
            QErrorExit(FILE_AND_LINE, "blocks of the matrix A is not assembled");
        }
    }
    /* checks if some of the blocks of the matrix B is assembled */
    if (B->blocks==NULL) {
        QErrorExit(FILE_AND_LINE, "blocks of the matrix B is not created");
    }
    else {
        err_code = QcMatIsAssembled(B, &B_assembled);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatIsAssembled(B)");
        if (B_assembled==QFALSE) {
            QErrorExit(FILE_AND_LINE, "blocks of the matrix B is not assembled");
        }
    }
    /* checks the dimension of blocks of A and B */
    if (A->dim_block!=B->dim_block) {
        printf("QcMatGEMM>> dimension of blocks (A) %"QINT_FMT"\n", A->dim_block);
        printf("QcMatGEMM>> dimension of blocks (B) %"QINT_FMT"\n", B->dim_block);
        QErrorExit(FILE_AND_LINE, "invalid dimension of blocks");
    }
    /* the scalar beta is not zero, we have C = alpha*op(A)*op(B)+beta*C */
    if (QAbs(beta[0])>QZEROTHRSH || QAbs(beta[1])>QZEROTHRSH) {
        /* checks the dimension of blocks of the matrix C */
        if (C->dim_block!=A->dim_block) {
            printf("QcMatGEMM>> dimension of blocks (A) %"QINT_FMT"\n", A->dim_block);
            printf("QcMatGEMM>> dimension of blocks (C) %"QINT_FMT"\n", C->dim_block);
            QErrorExit(FILE_AND_LINE, "invalid dimension of blocks");
        }
        err_code = QcMatIsAssembled(C, &C_assembled);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatIsAssembled(C)");
        /* scales C as C = beta*C if beta!=1 */
        if (C_assembled==QTRUE) {
            if (QAbs(beta[0]-1)>QZEROTHRSH || QAbs(beta[1])>QZEROTHRSH) {
                err_code = QcMatScale(beta, C);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatScale");
            }
        }
        else {
            QErrorExit(FILE_AND_LINE, "blocks of the matrix C is not assembled");
        }
    }
    /* the scalar beta is zero, so C = alpha*op(A)*op(B), and we will erase all previous information of C */
    else {
       if (C->blocks!=NULL) {
           err_code = QcMatDestroy(C);
           QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDestroy");
           err_code = QcMatCreate(C);
           QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatCreate");
       }
       /* sets the dimension of blocks and creates the blocks of C */
       err_code = QcMatBlockCreate(C, A->dim_block);
       QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatBlockCreate");
    }
    /* C_{IJ} = \sum_{K}alpha*op(A)_{IK}op(B)_{KJ}+beta*C_{IJ}, where
       op(A)_{IK} = A_{IK} (MAT_NO_OPERATION)
                  = op(A_{KI}) (MAT_TRANSPOSE)
                  = op(A_{KI}) (MAT_HERM_TRANSPOSE)
                  = op(A_{IK}) (MAT_COMPLEX_CONJUGATE) */
    if (op_A==MAT_NO_OPERATION || op_A==MAT_COMPLEX_CONJUGATE) {
        if (op_B==MAT_NO_OPERATION || op_B==MAT_COMPLEX_CONJUGATE) {
#if defined(QCMATRIX_STRASSEN_METHOD)
#endif
            for (irow=0; irow<A->dim_block; irow++) {
                for (krow=0; krow<A->dim_block; krow++) {
                    for (jcol=0; jcol<A->dim_block; jcol++) {
                        if (A->assembled[irow][krow]==QTRUE && B->assembled[krow][jcol]==QTRUE) {
                            err_code = CmplxMatGEMM(op_A,
                                                    op_B,
                                                    alpha,
                                                    &A->blocks[irow][krow],
                                                    &B->blocks[krow][jcol],
                                                    positive_one,
                                                    &C->blocks[irow][jcol]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGEMM");
                            C->assembled[irow][jcol] = QTRUE;
                        }
                    }
                }
            }
        }
        else if (op_B==MAT_TRANSPOSE || op_B==MAT_HERM_TRANSPOSE) {
#if defined(QCMATRIX_STRASSEN_METHOD)
#endif
            for (irow=0; irow<A->dim_block; irow++) {
                for (jcol=0; jcol<A->dim_block; jcol++) {
                    for (krow=0; krow<A->dim_block; krow++) {
                        if (A->assembled[irow][krow]==QTRUE && B->assembled[jcol][krow]==QTRUE) {
                            err_code = CmplxMatGEMM(op_A,
                                                    op_B,
                                                    alpha,
                                                    &A->blocks[irow][krow],
                                                    &B->blocks[jcol][krow],
                                                    positive_one,
                                                    &C->blocks[irow][jcol]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGEMM");
                            C->assembled[irow][jcol] = QTRUE;
                        }
                    }
                }
            }
        }
        else {
            printf("QcMatGEMM>> operation on matrix B: %d\n", op_B);
            QErrorExit(FILE_AND_LINE, "invalid matrix operation");
        }
    }
    else if (op_A==MAT_TRANSPOSE || op_A==MAT_HERM_TRANSPOSE) {
        if (op_B==MAT_NO_OPERATION || op_B==MAT_COMPLEX_CONJUGATE) {
#if defined(QCMATRIX_STRASSEN_METHOD)
#endif
            for (krow=0; krow<A->dim_block; krow++) {
                for (irow=0; irow<A->dim_block; irow++) {
                    for (jcol=0; jcol<A->dim_block; jcol++) {
                        if (A->assembled[krow][irow]==QTRUE && B->assembled[krow][jcol]==QTRUE) {
                            err_code = CmplxMatGEMM(op_A,
                                                    op_B,
                                                    alpha,
                                                    &A->blocks[krow][irow],
                                                    &B->blocks[krow][jcol],
                                                    positive_one,
                                                    &C->blocks[irow][jcol]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGEMM");
                            C->assembled[irow][jcol] = QTRUE;
                        }
                    }
                }
            }
        }
        else if (op_B==MAT_TRANSPOSE || op_B==MAT_HERM_TRANSPOSE) {
#if defined(QCMATRIX_STRASSEN_METHOD)
#endif
            for (irow=0; irow<A->dim_block; irow++) {
                for (jcol=0; jcol<A->dim_block; jcol++) {
                    for (krow=0; krow<A->dim_block; krow++) {
                        if (A->assembled[krow][irow]==QTRUE && B->assembled[jcol][krow]==QTRUE) {
                            err_code = CmplxMatGEMM(op_A,
                                                    op_B,
                                                    alpha,
                                                    &A->blocks[krow][irow],
                                                    &B->blocks[jcol][krow],
                                                    positive_one,
                                                    &C->blocks[irow][jcol]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGEMM");
                            C->assembled[irow][jcol] = QTRUE;
                        }
                    }
                }
            }
        }
        else {
            printf("QcMatGEMM>> operation on matrix B: %d\n", op_B);
            QErrorExit(FILE_AND_LINE, "invalid matrix operation");
        }
    }
    else {
        printf("QcMatGEMM>> operation on matrix A: %d\n", op_A);
        QErrorExit(FILE_AND_LINE, "invalid matrix operation");
    }
    /* usually the matrix C becomes a non-symmetric (non-Hermitian) matrix */
    C->sym_type = QNONSYMMAT;
    return QSUCCESS;
}
