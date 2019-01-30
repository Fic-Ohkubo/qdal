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

   This file implements the function QcMatGetMatProdTrace().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief gets the traces of the first few diagonal blocks of matrix-matrix product A*op(B)
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the left matrix, should be at least created by QcMatCreate()
         and QcMatBlockCreate()
     \param[QcMat:struct]{in} B the right matrix, should be at least created by QcMatCreate()
         and QcMatBlockCreate()
     \param[QcMatOperation:int]{in} op_B the operation on the matrix B, see file
         include/types/mat_operations.h
     \param[QInt:int]{in} num_blocks the number of diagonal blocks
     \param[QReal:real]{out} trace the traces, size is 2*\var{num_blocks}
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetMatProdTrace(QcMat *A,
                                QcMat *B,
                                const QcMatOperation op_B,
                                const QInt num_blocks,
                                QReal *trace)
{
    QBool A_assembled;
    QBool B_assembled;
    QReal tmp_trace[2];
    QInt pos_real, pos_imag;
    QInt irow, icol;
    QErrorCode err_code;
    /* checks the input number of diagonal blocks */
    if (num_blocks>A->dim_block) {
        printf("QcMatGetMatProdTrace>> input number of diagonal blocks %"QINT_FMT"\n",
               num_blocks);
        printf("QcMatGetMatProdTrace>> dimension of blocks (A) %"QINT_FMT"\n",
               A->dim_block);
        QErrorExit(FILE_AND_LINE, "invalid input number of diagonal blocks");
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
        printf("QcMatGetMatProdTrace>> dimension of blocks (A) %"QINT_FMT"\n",
               A->dim_block);
        printf("QcMatGetMatProdTrace>> dimension of blocks (B) %"QINT_FMT"\n",
               B->dim_block);
        QErrorExit(FILE_AND_LINE, "invalid dimension of blocks");
    }
    /* Tr(C_{II}) = Tr(\sum_{J}A_{IJ}*B_{JI}) = \sum_{J}Tr(A_{IJ}*B_{JI}) */
    if (op_B==MAT_NO_OPERATION || op_B==MAT_COMPLEX_CONJUGATE) {
        for (irow=0,pos_real=0,pos_imag=1; irow<num_blocks; irow++) {
            trace[pos_real] = 0;
            trace[pos_imag] = 0;
            for (icol=0; icol<A->dim_block; icol++) {
                if (A->assembled[irow][icol]==QTRUE && B->assembled[icol][irow]==QTRUE) {
                    err_code = CmplxMatGetMatProdTrace(&A->blocks[irow][icol],
                                                       &B->blocks[icol][irow],
                                                       op_B,
                                                       tmp_trace);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetMatProdTrace");
                    trace[pos_real] += tmp_trace[0];
                    trace[pos_imag] += tmp_trace[1];
                }
            }
            pos_real += 2;
            pos_imag += 2;
        }
    }
    /* Tr(C_{II}) = Tr(\sum_{J}A_{IJ}*B_{IJ}) = \sum_{J}Tr(A_{IJ}*B_{IJ}) */
    else if (op_B==MAT_TRANSPOSE || op_B==MAT_HERM_TRANSPOSE) {
        for (irow=0,pos_real=0,pos_imag=1; irow<num_blocks; irow++) {
            trace[pos_real] = 0;
            trace[pos_imag] = 0;
            for (icol=0; icol<A->dim_block; icol++) {
                if (A->assembled[irow][icol]==QTRUE && B->assembled[irow][icol]==QTRUE) {
                    err_code = CmplxMatGetMatProdTrace(&A->blocks[irow][icol],
                                                       &B->blocks[irow][icol],
                                                       op_B,
                                                       tmp_trace);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetMatProdTrace");
                    trace[pos_real] += tmp_trace[0];
                    trace[pos_imag] += tmp_trace[1];
                }
            }
            pos_real += 2;
            pos_imag += 2;
        }
    }
    else {
        printf("QcMatGetMatProdTrace>> operation on matrix B: %d\n", op_B);
        QErrorExit(FILE_AND_LINE, "invalid matrix operation");
    }
    return QSUCCESS;
}
