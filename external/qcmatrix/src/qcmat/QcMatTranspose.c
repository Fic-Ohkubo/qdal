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

   This file implements the function QcMatTranspose().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief performs an in-place or out-of-place matrix operation B = op(A)
     \author Bin Gao
     \date 2012-04-04
     \param[QcMatOperation:int]{in} op_A the operation on the matrix A, see file
         include/types/mat_operations.h
     \param[QcMat:struct]{in} A the matrix to perform matrix operation, should be
         at least assembled by QcMatAssemble()
     \param[QcMat:struct]{inout} B the result matrix; it could be A, or it should
         be at least created by QcMatCreate()
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatTranspose(const QcMatOperation op_A, QcMat *A, QcMat *B)
{
    QBool A_assembled;
    CmplxMat *T;
    QReal negative_one[2]={-1,0};
    QInt irow, icol;
    QErrorCode err_code;
    /* checks if some of the blocks of the matrix A is assembled */
    if (A->blocks==NULL) {
        QErrorExit(FILE_AND_LINE, "blocks of the matrix A is not created");
    }
    else {
        err_code = QcMatIsAssembled(A, &A_assembled);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatIsAssembled");
        if (A_assembled==QFALSE) {
            QErrorExit(FILE_AND_LINE, "blocks of the matrix A is not assembled");
        }
    }
    switch (op_A) {
    /* no matrix operation will be performed */
    case MAT_NO_OPERATION:
        /* we have B = A */
        if (A!=B) {
            err_code = QcMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDuplicate");
        }
        break;
    /* B = A^{*} */
    case MAT_COMPLEX_CONJUGATE:
        /* out-of-place complex conjugate */
        if (A!=B) {
            err_code = QcMatDuplicate(A, COPY_PATTERN_ONLY, B);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDuplicate");
        }
        for (irow=0; irow<A->dim_block; irow++) {
            for (icol=0; icol<A->dim_block; icol++) {
                if (A->assembled[irow][icol]==QTRUE) {
                    err_code = CmplxMatTranspose(MAT_COMPLEX_CONJUGATE,
                                                 &A->blocks[irow][icol],
                                                 &B->blocks[irow][icol]);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                }
            }
        }
        break;
    /* B = A^{T}, or B_{JI} = A_{IJ}^{T} */
    case MAT_TRANSPOSE:
        switch (A->sym_type) {
        /* A is a non-Hermitian matrix */
        case QNONSYMMAT:
            /* in-place transpose */
            if (A==B) {
                /* temporary matrix T to save the off-diagonal blocks */
                T = (CmplxMat *)malloc(sizeof(CmplxMat));
                if (T==NULL) {
                    QErrorExit(FILE_AND_LINE, "failed to allocate memory for T");
                }
                err_code = CmplxMatCreate(T);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatCreate");
                for (irow=0; irow<A->dim_block; irow++) {
                    /* does transpose of the off-diagonal blocks */
                    for (icol=0; icol<irow; icol++) {
                        if (A->assembled[irow][icol]==QTRUE) {
                            /* saves A_{JI} into the matrix T, which could be QNULLMAT */
                            err_code = CmplxMatDuplicate(&A->blocks[icol][irow],
                                                         COPY_PATTERN_AND_VALUE,
                                                         T);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
                            /* A_{JI} = A_{IJ}^{T} */
                            err_code = CmplxMatTranspose(MAT_TRANSPOSE,
                                                         &A->blocks[irow][icol],
                                                         &A->blocks[icol][irow]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                            /* A_{IJ} = A_{JI}^{T} */
                            if (A->assembled[icol][irow]==QTRUE) {
                                err_code = CmplxMatTranspose(MAT_TRANSPOSE,
                                                             T,
                                                             &A->blocks[irow][icol]);
                                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                            }
                            /* A_{IJ} = QNULLMAT */
                            else {
                                err_code = CmplxMatDuplicate(T,
                                                             COPY_PATTERN_AND_VALUE,
                                                             &A->blocks[irow][icol]);
                                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
                                A->assembled[irow][icol] = QFALSE;
                                A->assembled[icol][irow] = QTRUE;
                            }
                        }
                        else if (A->assembled[icol][irow]==QTRUE) {
                            /* saves A_{IJ} into the matrix T, which becomes QNULLMAT */
                            err_code = CmplxMatDuplicate(&A->blocks[irow][icol],
                                                         COPY_PATTERN_AND_VALUE,
                                                         T); 
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
                            /* A_{IJ} = A_{JI}^{T} */
                            err_code = CmplxMatTranspose(MAT_TRANSPOSE,
                                                         &A->blocks[icol][irow],
                                                         &A->blocks[irow][icol]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                            /* A_{JI} = QNULLMAT */
                            err_code = CmplxMatDuplicate(T,
                                                         COPY_PATTERN_AND_VALUE,
                                                         &A->blocks[icol][irow]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
                            A->assembled[irow][icol] = QTRUE;
                            A->assembled[icol][irow] = QFALSE;
                        }
                    }
                    /* does transpose of the diagonal block */
                    if (A->assembled[irow][irow]==QTRUE) {
                        err_code = CmplxMatTranspose(MAT_TRANSPOSE,
                                                     &A->blocks[irow][irow],
                                                     &A->blocks[irow][irow]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                    }
                }
                /* cleans */
                err_code = CmplxMatDestroy(T);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDestroy");
                free(T);
                T = NULL;
            }
            /* out-of-place transpose */
            else {
                /* erase all previous information of the matrix B */
                if (B->dim_block!=0) {
                    err_code = QcMatDestroy(B);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDestroy");
                    err_code = QcMatCreate(B);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatCreate");
                }
                /* sets the dimension of blocks and creates the blocks of B */
                err_code = QcMatBlockCreate(B, A->dim_block);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatBlockCreate");
                B->sym_type = QNONSYMMAT;
                /* does transpose block by block */
                for (irow=0; irow<A->dim_block; irow++) {
                    for (icol=0; icol<A->dim_block; icol++) {
                        if (A->assembled[irow][icol]==QTRUE) {
                            err_code = CmplxMatTranspose(MAT_TRANSPOSE,
                                                         &A->blocks[irow][icol],
                                                         &B->blocks[icol][irow]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                            /* this block of the matrix B is assembled */
                            B->assembled[icol][irow] = QTRUE;
                        }
                    }
                }
            }
            break;
        /* A is a Hermitian matrix */
        case QSYMMAT:
            /* out-of-place transpose */
            if (A!=B) {
                err_code = QcMatDuplicate(A, COPY_PATTERN_ONLY, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDuplicate");
            }
            /* we have B = A^{T} = A^{*} */
            for (irow=0; irow<A->dim_block; irow++) {
                for (icol=0; icol<A->dim_block; icol++) {
                    if (A->assembled[irow][icol]==QTRUE) {
                        err_code = CmplxMatTranspose(MAT_COMPLEX_CONJUGATE,
                                                     &A->blocks[irow][icol],
                                                     &B->blocks[irow][icol]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                    }
                }
            }
            break;
        /* A is an anti-Hermitian matrix */
        case QANTISYMMAT:
            /* out-of-place transpose */
            if (A!=B) {
                err_code = QcMatDuplicate(A, COPY_PATTERN_ONLY, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDuplicate");
            }
            /* we have B = A^{T} = -A^{*} */
            for (irow=0; irow<A->dim_block; irow++) {
                for (icol=0; icol<A->dim_block; icol++) {
                    if (A->assembled[irow][icol]==QTRUE) {
                        err_code = CmplxMatTranspose(MAT_COMPLEX_CONJUGATE,
                                                     &A->blocks[irow][icol],
                                                     &B->blocks[irow][icol]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                        err_code = CmplxMatScale(negative_one, &B->blocks[irow][icol]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatScale");
                    }
                }
            }
            break;
        default:
            printf("QcMatTranspose>> symmetry type of matrix A: %d\n", A->sym_type);
            QErrorExit(FILE_AND_LINE, "invalid symmetry type");
        }
        break;
    /* B = A^{\dagger}, or B_{JI} = A_{IJ}^{\dagger} */
    case MAT_HERM_TRANSPOSE:
        switch (A->sym_type) {
        /* A is a non-Hermitian matrix */
        case QNONSYMMAT:
            /* in-place Hermitian transpose */
            if (A==B) {
                /* temporary matrix T to save the off-diagonal blocks */
                T = (CmplxMat *)malloc(sizeof(CmplxMat));
                if (T==NULL) {
                    QErrorExit(FILE_AND_LINE, "failed to allocate memory for T");
                }
                err_code = CmplxMatCreate(T);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatCreate");
                for (irow=0; irow<A->dim_block; irow++) {
                    /* does Hermitian transpose of the off-diagonal blocks */
                    for (icol=0; icol<irow; icol++) {
                        if (A->assembled[irow][icol]==QTRUE) {
                            /* saves A_{JI} into the matrix T, which could be QNULLMAT */
                            err_code = CmplxMatDuplicate(&A->blocks[icol][irow],
                                                         COPY_PATTERN_AND_VALUE,
                                                         T);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
                            /* A_{JI} = A_{IJ}^{\dagger} */
                            err_code = CmplxMatTranspose(MAT_HERM_TRANSPOSE,
                                                         &A->blocks[irow][icol],
                                                         &A->blocks[icol][irow]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                            /* A_{IJ} = A_{JI}^{\dagger} */
                            if (A->assembled[icol][irow]==QTRUE) {
                                err_code = CmplxMatTranspose(MAT_HERM_TRANSPOSE,
                                                             T,
                                                             &A->blocks[irow][icol]);
                                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                            }
                            /* A_{IJ} = QNULLMAT */
                            else {
                                err_code = CmplxMatDuplicate(T,
                                                             COPY_PATTERN_AND_VALUE,
                                                             &A->blocks[irow][icol]);
                                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
                                A->assembled[irow][icol] = QFALSE;
                                A->assembled[icol][irow] = QTRUE;
                            }
                        }
                        else if (A->assembled[icol][irow]==QTRUE) {
                            /* saves A_{IJ} into the matrix T, which becomes QNULLMAT */
                            err_code = CmplxMatDuplicate(&A->blocks[irow][icol],
                                                         COPY_PATTERN_AND_VALUE,
                                                         T); 
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
                            /* A_{IJ} = A_{JI}^{\dagger} */
                            err_code = CmplxMatTranspose(MAT_HERM_TRANSPOSE,
                                                         &A->blocks[icol][irow],
                                                         &A->blocks[irow][icol]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                            /* A_{JI} = QNULLMAT */
                            err_code = CmplxMatDuplicate(T,
                                                         COPY_PATTERN_AND_VALUE,
                                                         &A->blocks[icol][irow]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
                            A->assembled[irow][icol] = QTRUE;
                            A->assembled[icol][irow] = QFALSE;
                        }
                    }
                    /* does Hermitian transpose of the diagonal block */
                    if (A->assembled[irow][irow]==QTRUE) {
                        err_code = CmplxMatTranspose(MAT_HERM_TRANSPOSE,
                                                     &A->blocks[irow][irow],
                                                     &A->blocks[irow][irow]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                    }
                }
                /* cleans */
                err_code = CmplxMatDestroy(T);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDestroy");
                free(T);
                T = NULL;
            }
            /* out-of-place Hermitian transpose */
            else {
                /* erase all previous information of the matrix B */
                if (B->dim_block!=0) {
                    err_code = QcMatDestroy(B);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDestroy");
                    err_code = QcMatCreate(B);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatCreate");
                }
                /* sets the dimension of blocks and creates the blocks of B */
                err_code = QcMatBlockCreate(B, A->dim_block);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatBlockCreate");
                B->sym_type = QNONSYMMAT;
                /* does Hermitian transpose block by block */
                for (irow=0; irow<A->dim_block; irow++) {
                    for (icol=0; icol<A->dim_block; icol++) {
                        if (A->assembled[irow][icol]==QTRUE) {
                            err_code = CmplxMatTranspose(MAT_HERM_TRANSPOSE,
                                                         &A->blocks[irow][icol],
                                                         &B->blocks[icol][irow]);
                            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatTranspose");
                            /* this block of the matrix B is assembled */
                            B->assembled[icol][irow] = QTRUE;
                        }
                    }
                }
            }
            break;
        /* A is a Hermitian matrix */
        case QSYMMAT:
            /* out-of-place Hermitian transpose, B = A^{\dagger} = A */
            if (A!=B) {
                err_code = QcMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDuplicate");
            }
            break;
        /* A is an anti-Hermitian matrix */
        case QANTISYMMAT:
            /* out-of-place Hermitian transpose */
            if (A!=B) {
                err_code = QcMatDuplicate(A, COPY_PATTERN_AND_VALUE, B);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDuplicate");
            }
            /* we have B = A^{\dagger} = -A */
            for (irow=0; irow<B->dim_block; irow++) {
                for (icol=0; icol<B->dim_block; icol++) {
                    if (B->assembled[irow][icol]==QTRUE) {
                        err_code = CmplxMatScale(negative_one,
                                                 &B->blocks[irow][icol]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatScale");
                    }
                }
            }
            break;
        default:
            printf("QcMatTranspose>> symmetry type of matrix A: %d\n", A->sym_type);
            QErrorExit(FILE_AND_LINE, "invalid symmetry type");
        }
        break;
    default:
        printf("QcMatTranspose>> operation on matrix A: %d\n", op_A);
        QErrorExit(FILE_AND_LINE, "invalid operation on matrix");
    }
    return QSUCCESS;
}
