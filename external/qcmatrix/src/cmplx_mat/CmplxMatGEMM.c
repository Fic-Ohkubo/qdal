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

   This file implements the function CmplxMatGEMM().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatGEMM QcMatGEMM
#define CmplxMatScale QcMatScale
#define CmplxMatCreate QcMatCreate
#define CmplxMatAXPY QcMatAXPY
#define CmplxMatDestroy QcMatDestroy
#else
#include "impls/cmplx_mat.h"
#endif

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief performs matrix-matrix multiplication C = alpha*op(A)*op(B)+beta*C,
        where valid operations op(...) can be found in file include/types/mat_operations.h
    \author Bin Gao
    \date 2012-04-04
    \param[QcMatOperation:int]{in} op_A the operation on the matrix A, see file
        include/types/mat_operations.h
    \param[QcMatOperation:int]{in} op_B the operation on the matrix B, see file
        include/types/mat_operations.h
    \param[QReal:real]{in} alpha the scalar number
    \param[CmplxMat:struct]{in} A the left matrix, should be at least assembled
        by CmplxMatAssemble()
    \param[CmplxMat:struct]{in} B the right matrix, should be at least assembled
        by CmplxMatAssemble()
    \param[QReal:real]{in} beta the scalar number
    \param[CmplxMat:struct]{inout} C the product matrix, should be at least created
        by CmplxMatCreate(), so that we require function MatGEMM() could assemble
        the matrix C if it is not
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatGEMM(const QcMatOperation op_A,
                        const QcMatOperation op_B,
                        const QReal alpha[],
                        CmplxMat *A,
                        CmplxMat *B,
                        const QReal beta[],
                        CmplxMat *C)
{
    QcMatOperation op_A_part;  /* operation on the real and imaginary parts of the matrix A */
    QReal sign_A_imag;        /* sign of the imaginary part of the matrix A after operation */
    QcMatOperation op_B_part;  /* operation on the real and imaginary parts of the matrix B */
    QReal sign_B_imag;        /* sign of the imaginary part of the matrix B after operation */
    QReal positive_one;
    QReal real_zero;
    CmplxMat *T;              /* T = op(A)*op(B) */
#if defined(QCMATRIX_3M_METHOD)
    RealMat *T1;
#endif
    QErrorCode err_code;
    if (C==A || C==B) {
        printf("QcMatGEMM>> matrix C should not be the same as the matrix A or B\n");
        QErrorExit(FILE_AND_LINE, "matrix C is the same as the matrix A or B");
    }
    /* if the scalar alpha is zero, we have C = beta*C */
    if (QAbs(alpha[0])<QZEROTHRSH && QAbs(alpha[1])<QZEROTHRSH) {
        err_code = CmplxMatScale(beta, C);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatScale");
        return QSUCCESS;
    }
    /* if the scalar beta is zero, we have C = alpha*op(A)*op(B) */
    if (QAbs(beta[0])<QZEROTHRSH && QAbs(beta[1])<QZEROTHRSH) {
        T = C;
    }
    /* if the matrix C is not assembled, we return C = alpha*op(A)*op(B) */
    else if (C->data_type==QNULLMAT) {
        T = C;
    }
    /* C_{R} = alpha_{R}*(A_{R}*B_{R}-A_{I}*B_{I})+beta_{R}*C_{R}
             - alpha_{I}*(A_{R}*B_{I}+A_{I}*B_{R})-beta_{I}*C_{I}
             = alpha_{R}*T_{R}+beta_{R}*C_{R}
             - alpha_{I}*T_{I}-beta_{I}*C_{I},
       C_{I} = alpha_{R}*(A_{R}*B_{I}+A_{I}*B_{R})+beta_{R}*C_{I}
             + alpha_{I}*(A_{R}*B_{R}-A_{I}*B_{I})+beta_{I}*C_{R}
             = alpha_{R}*T_{I}+beta_{R}*C_{I}
             + alpha_{I}*T_{R}+beta_{I}*C_{R},
       where
       T_{R} = A_{R}*B_{R}-A_{I}*B_{I},
       T_{I} = A_{R}*B_{I}+A_{I}*B_{R} */
    else {
        /* creates the matrix T */
        T = (CmplxMat *)malloc(sizeof(CmplxMat));
        if (T==NULL) {
            QErrorExit(FILE_AND_LINE, "failed to allocate memory for T");
        }
        err_code = CmplxMatCreate(T);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatCreate");
    }
    /* sets the operations on the real and imaginary parts, and the sign of the imaginary part */
    switch (op_A) {
    case MAT_NO_OPERATION:
        op_A_part = MAT_NO_OPERATION;
        sign_A_imag = 1;
        break;
    case MAT_TRANSPOSE:
        op_A_part = MAT_TRANSPOSE;
        sign_A_imag = 1;
        break;
    case MAT_HERM_TRANSPOSE:
        op_A_part = MAT_TRANSPOSE;
        sign_A_imag = -1;
        break;
    case MAT_COMPLEX_CONJUGATE:
        op_A_part = MAT_NO_OPERATION;
        sign_A_imag = -1;
        break;
    default:
        printf("CmplxMatGEMM>> operation on matrix A: %d\n", op_A);
        QErrorExit(FILE_AND_LINE, "invalid matrix operation");
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
        printf("CmplxMatGEMM>> operation on matrix B: %d\n", op_B);
        QErrorExit(FILE_AND_LINE, "invalid matrix operation");
    }
    /* performs matrix-matrix multiplication T = op(A)*op(B) */
    positive_one = 1;
    real_zero = 0;
    switch (A->data_type) {
    case QREALMAT:
        switch (B->data_type) {
        case QREALMAT:
            /* T_{R} = op(A_{R})*op(B_{R}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->real_part],
                                   real_zero,
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            break;
        case QIMAGMAT:
            /* T_{I} = op(A_{R})*op(B_{I}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->imag_part],
                                   real_zero,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            break;
        case QCMPLXMAT:
            /* T_{R} = op(A_{R})*op(B_{R}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->real_part],
                                   real_zero,
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            /* T_{I} = op(A_{R})*op(B_{I}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->imag_part],
                                   real_zero,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            break;
        default:
            printf("CmplxMatGEMM>> data type of matrix B: %d\n", B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    case QIMAGMAT:
        switch (B->data_type) {
        case QREALMAT:
            /* T_{I} = op(A_{I})*op(B_{R}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_A_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->real_part],
                                   real_zero,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            break;
        case QIMAGMAT:
            /* T_{R} = -op(A_{I})*op(B_{I}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   -sign_A_imag*sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->imag_part],
                                   real_zero,
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            break;
        case QCMPLXMAT:
            /* T_{I} = op(A_{I})*op(B_{R}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_A_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->real_part],
                                   real_zero,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            /* T_{R} = -op(A_{I})*op(B_{I}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   -sign_A_imag*sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->imag_part],
                                   real_zero,
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            break;
        default:
            printf("CmplxMatGEMM>> data type of matrix B: %d\n", B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    case QCMPLXMAT:
        switch (B->data_type) {
        case QREALMAT:
            /* T_{R} = op(A_{R})*op(B_{R}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->real_part],
                                   real_zero,
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            /* T_{I} = op(A_{I})*op(B_{R}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_A_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->real_part],
                                   real_zero,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            break;
        case QIMAGMAT:
            /* T_{R} = -op(A_{I})*op(B_{I}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   -sign_A_imag*sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->imag_part],
                                   real_zero,
                                   &T->cmplx_mat[T->real_part]);
            /* T_{I} = op(A_{R})*op(B_{I}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->imag_part],
                                   real_zero,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            break;
        case QCMPLXMAT:
#if defined(QCMATRIX_3M_METHOD)
            /* copies A_{R} to T_{R} */
            err_code = RealMatDuplicate(&A->cmplx_mat[A->real_part],
                                        COPY_PATTERN_AND_VALUE,
                                        &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
            /* calculates A_{R} +/- A_{I}, and saves it to T_{R} */
            err_code = RealMatAXPY(sign_A_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
            /* creates the temporary matrix T_{1}, and copies B_{R} to it */
            T1 = (RealMat *)malloc(sizeof(RealMat));
            if (T1==NULL) {
                QErrorExit(FILE_AND_LINE, "failed to allocate memory for T1");
            }
            err_code = RealMatCreate(T1);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatCreate");
            err_code = RealMatDuplicate(&B->cmplx_mat[B->real_part],
                                        COPY_PATTERN_AND_VALUE,
                                        T1);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
            /* calculates B_{R} +/- B_{I}, and saves it to T_{1} */
            err_code = RealMatAXPY(sign_B_imag*positive_one,
                                   &B->cmplx_mat[B->imag_part],
                                   T1);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
            /* calculates [op(A_{R})+op(A_{I})]*[op(B_{R})+op(B_{I})] and saves it to T_{I} */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   positive_one,
                                   &T->cmplx_mat[T->real_part],
                                   T1,
                                   real_zero,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            /* calculates op(A_{R})*op(B_{R}) and saves it to T_{R} */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->real_part],
                                   real_zero,
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            /* updates T_{I} = [op(A_{R})+op(A_{I})]*[op(B_{R})+op(B_{I})]-op(A_{R})*op(B_{R}) */
            err_code = RealMatAXPY(-positive_one,
                                   &T->cmplx_mat[T->real_part],
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
            /* calculates op(A_{I})*op(B_{I}) and saves it to T_{1} */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_A_imag*sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->imag_part],
                                   real_zero,
                                   T1);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            /* updates T_{R} = op(A_{R})*op(B_{R})-op(A_{I})*op(B_{I}) */
            err_code = RealMatAXPY(-positive_one, T1, &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
            /* updates T_{I} = [op(A_{R})+op(A_{I})]*[op(B_{R})+op(B_{I})]
                             - op(A_{R})*op(B_{R})-op(A_{I})*op(B_{I}) */
            err_code = RealMatAXPY(-positive_one, T1, &T->cmplx_mat[T->imag_part]);
            /* cleans */
            err_code = RealMatDestroy(T1);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDestroy");
            free(T1);
            T1 = NULL;
#else
            /* T_{R} = op(A_{R})*op(B_{R})-op(A_{I})*op(B_{I}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->real_part],
                                   real_zero,
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   -sign_A_imag*sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->imag_part],
                                   positive_one,
                                   &T->cmplx_mat[T->real_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            /* T_{I} = op(A_{R})*op(B_{I})+op(A_{I})*op(B_{R}) */
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_B_imag*positive_one,
                                   &A->cmplx_mat[A->real_part],
                                   &B->cmplx_mat[B->imag_part],
                                   real_zero,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
            err_code = RealMatGEMM(op_A_part,
                                   op_B_part,
                                   sign_A_imag*positive_one,
                                   &A->cmplx_mat[A->imag_part],
                                   &B->cmplx_mat[B->real_part],
                                   positive_one,
                                   &T->cmplx_mat[T->imag_part]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGEMM");
#endif
            break;
        default:
            printf("CmplxMatGEMM>> data type of matrix B: %d\n", B->data_type);
            QErrorExit(FILE_AND_LINE, "invalid data type");
        }
        break;
    default:
        printf("CmplxMatGEMM>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    /* sets the data and symmetry types of matrix T */
    T->data_type = A->data_type*B->data_type;
    T->sym_type = QNONSYMMAT;
    /* the scalar beta is zero or C is not assembled, we have C = alpha*op(A)*op(B) */
    if (T==C) {
        /* if the scalar alpha is 1, we do not need to scale */
        if (QAbs(alpha[0]-1)>QZEROTHRSH || QAbs(alpha[1])>QZEROTHRSH) {
            err_code = CmplxMatScale(alpha, C);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatScale");
        }
    }
    /* C = alpha*T+beta*C */
    else {
        /* if the scalar beta is 1, we do not need to scale */
        if (QAbs(beta[0]-1)>QZEROTHRSH || QAbs(beta[1])>QZEROTHRSH) {
            err_code = CmplxMatScale(beta, C);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatScale");
        }
        err_code = CmplxMatAXPY(alpha, T, C);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatAXPY");
        /* cleans */
        err_code = CmplxMatDestroy(T);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDestroy");
        free(T);
        T = NULL;
    }
    return QSUCCESS;
}
