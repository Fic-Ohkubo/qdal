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

   This file tests the function QcMatGEMM().

   2014-03-28, Bin Gao:
   * first version
*/

/* header file of QcMatrix library */
#include "qcmatrix.h"
/* parameters for test suite */
#include "tests/qcmatrix_test_param.h"
/* BLAS routines */
#include "lapack/qcmatrix_c_blas.h"

/*% \brief tests the function QcMatGEMM()
    \author Bin Gao
    \date 2014-03-28
    \param[QcMat:type]{in} A the matrix
    \param[QcMat:type]{in} B the matrix
*/
QVoid test_c_QcMatGEMM(QcMat *A, QcMat *B)
{
    QcMat C;                          /* product matrix */
    QcMat T;                          /* temporary matrix */
    QBool assembled;                  /* indicates if the matrix is assembled or not */
    QInt dim_block;                   /* dimension of blocks */
    QInt dim_mat;                     /* dimension of each block */
    QInt size_values;                 /* number of elements in the matrix */
    QReal *A_real;                    /* values of the real part of the matrix A */
    QReal *A_imag;                    /* values of the imaginary part of the matrix A */
    QReal *B_real;                    /* values of the real part of the matrix B */
    QReal *B_imag;                    /* values of the imaginary part of the matrix B */
    QReal *C_real;                    /* values of the real part of the matrix C */
    QReal *C_imag;                    /* values of the imaginary part of the matrix C */
    QReal *T_real;                    /* temporary values */
    QReal *T_imag;                    /* temporary values */
    QcSymType sym_type[3] = {QSYMMAT,QANTISYMMAT,QNONSYMMAT};  /* all symmetry types */
    QcDataType data_type[3] = {QREALMAT,QIMAGMAT,QCMPLXMAT};   /* all data types */
    QcMatOperation all_mat_operations[4] = {MAT_NO_OPERATION,  /* all matrix operations */
                                            MAT_TRANSPOSE,
                                            MAT_HERM_TRANSPOSE,
                                            MAT_COMPLEX_CONJUGATE};
    QChar trans_A;                    /* operation on A, for BLAS routine */
    QChar trans_B;                    /* operation on B, for BLAS routine */
    QReal alpha[2][2];                /* the scalar number */
    QReal beta[2];                    /* the scalar number */
    QBool is_equal;                   /* indicates if the matrix and array have the same values */
    QReal positive_one = 1.0;         /* positive one */
    QReal negative_one = -1.0;        /* negative one */
    QReal real_zero = 0.0;            /* zero */
    QInt iop, jop, isym, idat, jdat;  /* incremental recorders */
    QInt ival;                        /* incremental recorder over values */
    QErrorCode ierr;                  /* error information */
    /* checks if the matrices A and B are assembled */
    ierr = QcMatIsAssembled(A, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatGEMM>> matrix A is not assembled ...\n");
            printf("test_c_QcMatGEMM>> QcMatGEMM() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatGEMM>> failed to call QcMatIsAssembled(A)\n");
        exit(ierr);
    }
    ierr = QcMatIsAssembled(B, &assembled);
    if (ierr==QSUCCESS) {
        if (assembled!=QTRUE) {
            printf("test_c_QcMatGEMM>> matrix B is not assembled ...\n");
            printf("test_c_QcMatGEMM>> QcMatGEMM() will not be tested ...\n");
            return;
        }
    }
    else {
        printf("test_c_QcMatGEMM>> failed to call QcMatIsAssembled(B)\n");
        exit(ierr);
    }
    /* gets the dimension of blocks */
    ierr = QcMatGetDimBlock(A, &dim_block);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatGEMM>> QcMatGetDimBlock(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatGEMM>> failed to call QcMatGetDimBlock(A)\n");
        exit(ierr);
    }
    /* gets the dimension of each block */
    ierr = QcMatGetDimMat(A, &dim_mat, &dim_mat);
    if (ierr==QSUCCESS) {
        printf("test_c_QcMatGEMM>> QcMatGetDimMat(A) passed ...\n");
    }
    else {
        printf("test_c_QcMatGEMM>> failed to call QcMatGetDimMat(A)\n");
        exit(ierr);
    }
    size_values = dim_block*dim_block*dim_mat*dim_mat;
    /* allocates memory for the elements of the matrices */
    A_real = (QReal *)malloc(sizeof(QReal)*size_values);
    if (A_real==NULL) {
        printf("test_c_QcMatGEMM>> failed to allocate A_real\n");
        exit(QFAILURE);
    }
    A_imag = (QReal *)malloc(sizeof(QReal)*size_values); 
    if (A_imag==NULL) {
        printf("test_c_QcMatGEMM>> failed to allocate A_imag\n");
        exit(QFAILURE);
    }
    B_real = (QReal *)malloc(sizeof(QReal)*size_values); 
    if (B_real==NULL) {
        printf("test_c_QcMatGEMM>> failed to allocate B_real\n");
        exit(QFAILURE);
    }
    B_imag = (QReal *)malloc(sizeof(QReal)*size_values); 
    if (B_imag==NULL) {
        printf("test_c_QcMatGEMM>> failed to allocate B_imag\n");
        exit(QFAILURE);
    }
    C_real = (QReal *)malloc(sizeof(QReal)*size_values); 
    if (C_real==NULL) {
        printf("test_c_QcMatGEMM>> failed to allocate C_real\n");
        exit(QFAILURE);
    }
    C_imag = (QReal *)malloc(sizeof(QReal)*size_values); 
    if (C_imag==NULL) {
        printf("test_c_QcMatGEMM>> failed to allocate C_imag\n");
        exit(QFAILURE);
    }
    T_real = (QReal *)malloc(sizeof(QReal)*size_values); 
    if (T_real==NULL) {
        printf("test_c_QcMatGEMM>> failed to allocate T_real\n");
        exit(QFAILURE);
    }
    /* gets all the values of the matrices */
    ierr = QcMatGetAllValues(A, QFALSE, size_values, A_real, A_imag);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatGEMM>> failed to call QcMatGetAllValues(A)\n");
        exit(ierr);
    }
    ierr = QcMatGetAllValues(B, QFALSE, size_values, B_real, B_imag);
    if (ierr!=QSUCCESS) {
        printf("test_c_QcMatGEMM>> failed to call QcMatGetAllValues(B)\n");
        exit(ierr);
    }
    /* set the alpha as complex number, since other types have been tested in CmplxMatScale() or CmplxMatAXPY() */
    alpha[0][0] = 0.01; alpha[0][1] = 0.01;  /* complex number */
    alpha[1][0] = 0.0; alpha[1][1] = 0.0;    /* zero */
    /* we first test the matrix C which is not assembled */
    beta[0] = 0.0; beta[1] = 0.0;
    for (iop=0; iop<4; iop++) {
        switch (all_mat_operations[iop]) {
        case MAT_NO_OPERATION:
            trans_A = 'N';
            break;
        case MAT_TRANSPOSE:
            trans_A = 'T';
            break;
        case MAT_HERM_TRANSPOSE:
            trans_A = 'T';
            for (ival=0; ival<size_values; ival++) {
                A_imag[ival] = -A_imag[ival];
            }
            break;
        case MAT_COMPLEX_CONJUGATE:
            trans_A = 'N';
            for (ival=0; ival<size_values; ival++) {
                A_imag[ival] = -A_imag[ival];
            }
            break;
        default:
            printf("test_c_QcMatGEMM>> invalid operation on matrix A\n");
            exit(all_mat_operations[iop]);
        }
        for (jop=0; jop<4; jop++) {
            switch (all_mat_operations[jop]) {
            case MAT_NO_OPERATION:
                trans_B = 'N';
                break;
            case MAT_TRANSPOSE:
                trans_B = 'T';
                break;
            case MAT_HERM_TRANSPOSE:
                trans_B = 'T';
                for (ival=0; ival<size_values; ival++) {
                    B_imag[ival] = -B_imag[ival];
                }
                break;
            case MAT_COMPLEX_CONJUGATE:
                trans_B = 'N';
                for (ival=0; ival<size_values; ival++) {
                    B_imag[ival] = -B_imag[ival];
                }
                break;
            default:
                printf("test_c_QcMatGEMM>> invalid operation on matrix B\n");
                exit(all_mat_operations[jop]);
            }
            /* tests different alpha's */
            for (jdat=0; jdat<1; jdat++) {  /* we can not test alpha=0 here, since the matrix C is not assembled */
                ierr = QcMatCreate(&C);
                if (ierr!=QSUCCESS) {
                    printf("test_c_QcMatGEMM>> failed to call QcMatCreate(C)\n");
                    exit(ierr);
#if defined(DEBUG)
                }
                else {
                    printf("test_c_QcMatGEMM>> QcMatCreate(C) passed ...\n");
#endif
                }
                ierr = QcMatGEMM(all_mat_operations[iop],
                                 all_mat_operations[jop],
                                 alpha[jdat],
                                 A,
                                 B,
                                 beta,
                                 &C);
                if (ierr==QSUCCESS) {
                    /* C = (a_{R}+i*a_{I})*(A_{R}+i*A_{I})*(B_{R}+i*B_{I})
                         = (a_{R}+i*a_{I})*(A_{R}*B_{R}-A_{I}*B_{I}+i*(A_{R}*B_{I}+A_{I}*B_{R}))
                         =   a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})
                         -   a_{I}*(A_{R}*B_{I}+A_{I}*B_{R})
                         + i*a_{R}*(A_{R}*B_{I}+A_{I}*B_{R})
                         + i*a_{I}*(A_{R}*B_{R}-A_{I}*B_{I}) */
                    C_BLAS_GEMM(trans_A,            /* A_{R}*B_{R} */
                                trans_B,
                                dim_block*dim_mat,
                                dim_block*dim_mat,
                                dim_block*dim_mat,
                                positive_one,
                                A_real,
                                dim_block*dim_mat,
                                B_real,
                                dim_block*dim_mat,
                                real_zero,
                                C_real,
                                dim_block*dim_mat);
                    C_BLAS_GEMM(trans_A,            /* A_{I}*B_{I} */
                                trans_B,
                                dim_block*dim_mat,
                                dim_block*dim_mat,
                                dim_block*dim_mat,
                                positive_one,
                                A_imag,
                                dim_block*dim_mat,
                                B_imag,
                                dim_block*dim_mat,
                                real_zero,
                                T_real,
                                dim_block*dim_mat);
                    C_BLAS_AXPY(size_values,   /* -A_{I}*B_{I}+A_{R}*B_{R} */
                                negative_one,
                                T_real,
                                1,
                                C_real,
                                1);
                    C_BLAS_GEMM(trans_A,            /* A_{R}*B_{I} */
                                trans_B,
                                dim_block*dim_mat,
                                dim_block*dim_mat,
                                dim_block*dim_mat,
                                positive_one,
                                A_real,
                                dim_block*dim_mat,
                                B_imag,
                                dim_block*dim_mat,
                                real_zero,
                                C_imag,
                                dim_block*dim_mat);
                    C_BLAS_GEMM(trans_A,            /* A_{I}*B_{R} */
                                trans_B,
                                dim_block*dim_mat,
                                dim_block*dim_mat,
                                dim_block*dim_mat,
                                positive_one,
                                A_imag,
                                dim_block*dim_mat,
                                B_real,
                                dim_block*dim_mat,
                                real_zero,
                                T_real,
                                dim_block*dim_mat);
                    C_BLAS_AXPY(size_values,   /* A_{I}*B_{R}+A_{R}*B_{I} */
                                positive_one,
                                T_real,
                                1,
                                C_imag,
                                1);
                    C_BLAS_COPY(size_values,
                                C_real,
                                1,
                                T_real,
                                1);
                    C_BLAS_SCAL(size_values, alpha[jdat][0], C_real, 1);  /* a_{R}*(A_{R}*B_{R}-A_{I}*B_{I}) */
                    C_BLAS_AXPY(size_values,      /* -a_{I}*(A_{R}*B_{I}+A_{I}*B_{R})+a_{R}*(A_{R}*B_{R}-A_{I}*B_{I}) */
                                -alpha[jdat][1],
                                C_imag,
                                1,
                                C_real,
                                1);
                    C_BLAS_SCAL(size_values, alpha[jdat][0], C_imag, 1);  /* a_{R}*(A_{R}*B_{I}+A_{I}*B_{R}) */
                    C_BLAS_AXPY(size_values,     /* a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})+a_{R}*(A_{R}*B_{I}+A_{I}*B_{R}) */
                                alpha[jdat][1],
                                T_real,
                                1,
                                C_imag,
                                1);
                    ierr = QcMatCfArray(&C,
                                        QFALSE,
                                        size_values,
                                        C_real,
                                        C_imag,
                                        &is_equal);
                    if (ierr==QSUCCESS) {
                        if (is_equal!=QTRUE) {
                            /* dumps results to check */
                            printf("test_c_QcMatGEMM>> parameters %d, %d, (%f, %f)\n",
                                   all_mat_operations[iop],
                                   all_mat_operations[jop],
                                   alpha[jdat][0],
                                   alpha[jdat][1]);
#if defined(QCMATRIX_ENABLE_VIEW)
                            ierr = QcMatWrite(A, "QcMatGEMM_A", ASCII_VIEW);
                            if (ierr!=QSUCCESS) {
                                printf("test_c_QcMatGEMM>> failed to call QcMatWrite(A)\n");
                                exit(ierr);
                            }
                            ierr = QcMatWrite(B, "QcMatGEMM_B", ASCII_VIEW);
                            if (ierr!=QSUCCESS) {
                                printf("test_c_QcMatGEMM>> failed to call QcMatWrite(B)\n");
                                exit(ierr);
                            }
                            ierr = QcMatWrite(&C, "QcMatGEMM_C", ASCII_VIEW);
                            if (ierr!=QSUCCESS) {
                                printf("test_c_QcMatGEMM>> failed to call QcMatWrite(C)\n");
                                exit(ierr);
                            }
#endif
                            printf("test_c_QcMatGEMM>> real part of C from BLAS\n");
                            for (ival=0; ival<size_values; ival++) {
                                if (ival%4==3) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                                    printf("%20.12lf\n", C_real[ival]);
#else
                                    printf("%20.12f\n", C_real[ival]);
#endif
                                }
                                else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                                    printf("%20.12lf  ", C_real[ival]);
#else
                                    printf("%20.12f  ", C_real[ival]);
#endif
                                }
                            }
                            printf("test_c_QcMatGEMM>> imaginary part of C from BLAS\n");
                            for (ival=0; ival<size_values; ival++) {
                                if (ival%4==3) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                                    printf("%20.12lf\n", C_imag[ival]);
#else
                                    printf("%20.12f\n", C_imag[ival]);
#endif
                                }
                                else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                                    printf("%20.12lf  ", C_imag[ival]);
#else
                                    printf("%20.12f  ", C_imag[ival]);
#endif
                                }
                            }
                            printf("test_c_QcMatGEMM>> QcMatGEMM(beta=0) failed\n");
                            exit(is_equal);
#if defined(DEBUG)
                        }
                        else {
                            printf("test_c_QcMatGEMM>> QcMatGEMM(beta=0) passed %d, %d, (%f, %f)\n",
                                   all_mat_operations[iop],
                                   all_mat_operations[jop],
                                   alpha[jdat][0]
                                   alpha[jdat][1]);
#endif
                        }
                    }
                    else {
                        printf("test_c_QcMatGEMM>> failed to call QcMatCfArray(C)\n");
                        exit(ierr);
                    }
                }
                else {
                    printf("test_c_QcMatGEMM>> parameters %d, %d, (%f, %f)\n",
                           all_mat_operations[iop],
                           all_mat_operations[jop],
                           alpha[jdat][0],
                           alpha[jdat][1]);
#if defined(QCMATRIX_ENABLE_VIEW)
                    ierr = QcMatWrite(A, "QcMatGEMM_A", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatGEMM>> failed to call QcMatWrite(A)\n");
                        exit(ierr);
                    }
                    ierr = QcMatWrite(B, "QcMatGEMM_B", ASCII_VIEW);
                    if (ierr!=QSUCCESS) {
                        printf("test_c_QcMatGEMM>> failed to call QcMatWrite(B)\n");
                        exit(ierr);
                    }
#endif
                    printf("test_c_QcMatGEMM>> failed to call QcMatGEMM(C)\n");
                    exit(ierr);
                }
                /* cleans the matrix C */
                ierr = QcMatDestroy(&C);
                if (ierr!=QSUCCESS) {
                    printf("test_c_QcMatGEMM>> failed to call QcMatDestroy(C)\n");
                    exit(ierr);
#if defined(DEBUG)
                }
                else {
                    printf("test_c_QcMatGEMM>> QcMatDestroy(C) passed ...\n");
#endif
                }
            }
            if (all_mat_operations[jop]==MAT_HERM_TRANSPOSE ||
                all_mat_operations[jop]==MAT_COMPLEX_CONJUGATE) {
                for (ival=0; ival<size_values; ival++) {
                    B_imag[ival] = -B_imag[ival];
                }
            }
        }
        if (all_mat_operations[iop]==MAT_HERM_TRANSPOSE ||
            all_mat_operations[iop]==MAT_COMPLEX_CONJUGATE) {
            for (ival=0; ival<size_values; ival++) {
                A_imag[ival] = -A_imag[ival];
            }
        }
    }
#if !defined(DEBUG)
    printf("test_c_QcMatGEMM>> QcMatGEMM(beta=0) passed ...\n");
#endif
    /* allocates memory for temporary values */
    T_imag = (QReal *)malloc(sizeof(QReal)*size_values);
    if (T_imag==NULL) {
        printf("test_c_QcMatGEMM>> failed to allocate T_imag\n");
        exit(QFAILURE);
    }
    /* set the beta as complex number, since it is only used in QcMatScale() and/or CmplxMatScale(),
       different beta's have been tested before */
    beta[0] = 0.1; beta[1] = 0.1;
    /* tests different symmetry types (symmetric, anti-symmetric, non-symmetric) for matrix C */
    for (isym=0; isym<3; isym++) {
        /* tests different data types (real, imaginary, complex) for matrix C */
        for (idat=0; idat<3; idat++) {
            /* generates a random matrix C according to its symmetry and data types */
            ierr = QcMatCreate(&C);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatGEMM>> failed to call QcMatCreate(C)\n");
                exit(ierr);
#if defined(DEBUG)
            }
            else {
                printf("test_c_QcMatGEMM>> QcMatCreate(C) passed ...\n");
#endif
            }
            ierr = QcMatSetRandMat(&C,
                                   sym_type[isym],
                                   data_type[idat],
                                   dim_block,
                                   dim_mat,
                                   dim_mat);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatGEMM>> failed to call QcMatSetRandMat(C)\n");
                exit(ierr);
#if defined(DEBUG)
            }
            else {
                printf("test_c_QcMatGEMM>> QcMatSetRandMat(C) passed ...\n");
#endif
            }
            /* temporary matrix as beta*C */
            ierr = QcMatCreate(&T);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatGEMM>> QcMatCreate(T) failed\n");
                exit(ierr);
#if defined(DEBUG)
            }
            else {
                printf("test_c_QcMatGEMM>> QcMatCreate(T) passed ...\n");
#endif
            }
            /* we first test the matrix C which is not assembled */
            for (iop=0; iop<4; iop++) {
                switch (all_mat_operations[iop]) {
                case MAT_NO_OPERATION:
                    trans_A = 'N';
                    break;
                case MAT_TRANSPOSE:
                    trans_A = 'T';
                    break;
                case MAT_HERM_TRANSPOSE:
                    trans_A = 'T';
                    for (ival=0; ival<size_values; ival++) {
                        A_imag[ival] = -A_imag[ival];
                    }
                    break;
                case MAT_COMPLEX_CONJUGATE:
                    trans_A = 'N';
                    for (ival=0; ival<size_values; ival++) {
                        A_imag[ival] = -A_imag[ival];
                    }
                    break;
                default:
                    printf("test_c_QcMatGEMM>> invalid operation on matrix A\n");
                    exit(all_mat_operations[iop]);
                }
                for (jop=0; jop<4; jop++) {
                    switch (all_mat_operations[jop]) {
                    case MAT_NO_OPERATION:
                        trans_B = 'N';
                        break;
                    case MAT_TRANSPOSE:
                        trans_B = 'T';
                        break;
                    case MAT_HERM_TRANSPOSE:
                        trans_B = 'T';
                        for (ival=0; ival<size_values; ival++) {
                            B_imag[ival] = -B_imag[ival];
                        }
                        break;
                    case MAT_COMPLEX_CONJUGATE:
                        trans_B = 'N';
                        for (ival=0; ival<size_values; ival++) {
                            B_imag[ival] = -B_imag[ival];
                        }
                        break;
                    default:
                        printf("test_c_QcMatGEMM>> invalid operation on matrix B\n");
                        exit(all_mat_operations[jop]);
                    }
                    /* tests different alpha's */
                    for (jdat=0; jdat<2; jdat++) {
                        /* T = beta*C */
                        ierr = QcMatDuplicate(&C, COPY_PATTERN_AND_VALUE, &T);
                        if (ierr!=QSUCCESS) {
                            printf("test_c_QcMatGEMM>> QcMatDuplicate(T) failed\n");
                            exit(ierr);
                        }
                        ierr = QcMatScale(beta, &T);
                        if (ierr!=QSUCCESS) {
                            printf("test_c_QcMatGEMM>> QcMatScale(T) failed\n");
                            exit(ierr);
                        }
                        /* calcualtes C = alpha*A*B+beta*C using QcMatGEMM() */
                        ierr = QcMatGEMM(all_mat_operations[iop],
                                         all_mat_operations[jop],
                                         alpha[jdat],
                                         A,
                                         B,
                                         beta,
                                         &C);
                        if (ierr==QSUCCESS) {
                            /* C = (a_{R}+i*a_{I})*(A_{R}+i*A_{I})*(B_{R}+i*B_{I})+beta*C
                                 = (a_{R}+i*a_{I})*(A_{R}*B_{R}-A_{I}*B_{I}+i*(A_{R}*B_{I}+A_{I}*B_{R}))+beta*C
                                 =   a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})
                                 -   a_{I}*(A_{R}*B_{I}+A_{I}*B_{R})
                                 +   (beta*C)_{R}
                                 + i*a_{R}*(A_{R}*B_{I}+A_{I}*B_{R})
                                 + i*a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})
                                 + i*(beta*C)_{I} */
                            C_BLAS_GEMM(trans_A,            /* A_{R}*B_{R} */
                                        trans_B,
                                        dim_block*dim_mat,
                                        dim_block*dim_mat,
                                        dim_block*dim_mat,
                                        positive_one,
                                        A_real,
                                        dim_block*dim_mat,
                                        B_real,
                                        dim_block*dim_mat,
                                        real_zero,
                                        T_real,
                                        dim_block*dim_mat);
                            C_BLAS_GEMM(trans_A,            /* A_{I}*B_{I} */
                                        trans_B,
                                        dim_block*dim_mat,
                                        dim_block*dim_mat,
                                        dim_block*dim_mat,
                                        positive_one,
                                        A_imag,
                                        dim_block*dim_mat,
                                        B_imag,
                                        dim_block*dim_mat,
                                        real_zero,
                                        C_real,
                                        dim_block*dim_mat);
                            C_BLAS_AXPY(size_values,   /* -A_{I}*B_{I}+A_{R}*B_{R} */
                                        negative_one,
                                        C_real,
                                        1,
                                        T_real,
                                        1);
                            C_BLAS_GEMM(trans_A,            /* A_{R}*B_{I} */
                                        trans_B,
                                        dim_block*dim_mat,
                                        dim_block*dim_mat,
                                        dim_block*dim_mat,
                                        positive_one,
                                        A_real,
                                        dim_block*dim_mat,
                                        B_imag,
                                        dim_block*dim_mat,
                                        real_zero,
                                        T_imag,
                                        dim_block*dim_mat);
                            C_BLAS_GEMM(trans_A,            /* A_{I}*B_{R} */
                                        trans_B,
                                        dim_block*dim_mat,
                                        dim_block*dim_mat,
                                        dim_block*dim_mat,
                                        positive_one,
                                        A_imag,
                                        dim_block*dim_mat,
                                        B_real,
                                        dim_block*dim_mat,
                                        real_zero,
                                        C_imag,
                                        dim_block*dim_mat);
                            C_BLAS_AXPY(size_values,   /* A_{I}*B_{R}+A_{R}*B_{I} */
                                        positive_one,
                                        C_imag,
                                        1,
                                        T_imag,
                                        1);
                            /* gets the values of original C matrix scaled by beta */
                            ierr = QcMatGetAllValues(&T,
                                                     QFALSE,
                                                     size_values,
                                                     C_real,
                                                     C_imag);
                            if (ierr!=QSUCCESS) {
                                printf("test_c_QcMatGEMM>> failed to call QcMatGetAllValues(T)\n");
                                exit(ierr);
                            }
                            C_BLAS_AXPY(size_values,     /* a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})+(beta*C)_{R} */
                                        alpha[jdat][0],
                                        T_real,
                                        1,
                                        C_real,
                                        1);
                            C_BLAS_AXPY(size_values,     /* -a_{I}*(A_{R}*B_{I}+A_{I}*B_{R})+a_{R}*(A_{R}*B_{R}-A_{I}*B_{I})+(beta*C)_{R} */
                                        -alpha[jdat][1],
                                        T_imag,
                                        1,
                                        C_real,
                                        1);
                            C_BLAS_AXPY(size_values,     /* a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})+(beta*C)_{I} */
                                        alpha[jdat][1],
                                        T_real,
                                        1,
                                        C_imag,
                                        1);
                            C_BLAS_AXPY(size_values,     /* a_{R}*(A_{R}*B_{I}+A_{I}*B_{R})+a_{I}*(A_{R}*B_{R}-A_{I}*B_{I})+(beta*C)_{I} */
                                        alpha[jdat][0],
                                        T_imag,
                                        1,
                                        C_imag,
                                        1);
                            ierr = QcMatCfArray(&C,
                                                QFALSE,
                                                size_values,
                                                C_real,
                                                C_imag,
                                                &is_equal);
                            if (ierr==QSUCCESS) {
                                if (is_equal!=QTRUE) {
                                    /* dumps results to check */
                                    printf("test_c_QcMatGEMM>> QcMatGEMM(beta=0) passed ...\n");
                                    printf("test_c_QcMatGEMM>> %d, %d, (%f, %f), (%f, %f), %d, %d\n",
                                           all_mat_operations[iop],
                                           all_mat_operations[jop],
                                           alpha[jdat][0],
                                           alpha[jdat][1],
                                           beta[0],
                                           beta[1],
                                           data_type[idat],
                                           sym_type[isym]);
#if defined(QCMATRIX_ENABLE_VIEW)
                                    ierr = QcMatWrite(A, "QcMatGEMM_A", ASCII_VIEW);
                                    if (ierr!=QSUCCESS) {
                                        printf("test_c_QcMatGEMM>> failed to call QcMatWrite(A)\n");
                                        exit(ierr);
                                    }
                                    ierr = QcMatWrite(B, "QcMatGEMM_B", ASCII_VIEW);
                                    if (ierr!=QSUCCESS) {
                                        printf("test_c_QcMatGEMM>> failed to call QcMatWrite(B)\n");
                                        exit(ierr);
                                    }
                                    ierr = QcMatWrite(&C, "QcMatGEMM_C", ASCII_VIEW);
                                    if (ierr!=QSUCCESS) {
                                        printf("test_c_QcMatGEMM>> failed to call QcMatWrite(C)\n");
                                        exit(ierr);
                                    }
                                    ierr = QcMatWrite(&T, "QcMatGEMM_T", ASCII_VIEW);
                                    if (ierr!=QSUCCESS) {
                                        printf("test_c_QcMatGEMM>> failed to call QcMatWrite(T)\n");
                                        exit(ierr);
                                    }
#endif
                                    printf("test_c_QcMatGEMM>> real part of C from BLAS\n");
                                    for (ival=0; ival<size_values; ival++) {
                                        if (ival%4==3) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                                            printf("%20.12lf\n", C_real[ival]);
#else
                                            printf("%20.12f\n", C_real[ival]);
#endif
                                        }
                                        else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                                            printf("%20.12lf  ", C_real[ival]);
#else
                                            printf("%20.12f  ", C_real[ival]);
#endif
                                        }
                                    }
                                    printf("test_c_QcMatGEMM>> imaginary part of C from BLAS\n");
                                    for (ival=0; ival<size_values; ival++) {
                                        if (ival%4==3) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                                            printf("%20.12lf\n", C_imag[ival]);
#else
                                            printf("%20.12f\n", C_imag[ival]);
#endif
                                        }
                                        else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                                            printf("%20.12lf  ", C_imag[ival]);
#else
                                            printf("%20.12f  ", C_imag[ival]);
#endif
                                        }
                                    }
                                    printf("test_c_QcMatGEMM>> QcMatGEMM(beta!=0) failed\n");
                                    exit(is_equal);
#if defined(DEBUG)
                                }
                                else {
                                    printf("test_c_QcMatGEMM>> QcMatGEMM(beta!=0) passed ...\n");
                                    printf("test_c_QcMatGEMM>> %d, %d, (%f, %f), (%f, %f), %d, %d\n",
                                           all_mat_operations[iop],
                                           all_mat_operations[jop],
                                           alpha[jdat][0],
                                           alpha[jdat][1],
                                           beta[0],
                                           beta[1],
                                           data_type[idat],
                                           sym_type[isym]);
#endif
                                }
                            }
                            else {
                                printf("test_c_QcMatGEMM>> failed to call QcMatCfArray(C)\n");
                                exit(ierr);
                            }
                        }
                        else {
                            printf("test_c_QcMatGEMM>> %d, %d, (%f, %f), (%f, %f), %d, %d\n",
                                   all_mat_operations[iop],
                                   all_mat_operations[jop],
                                   alpha[jdat][0],
                                   alpha[jdat][1],
                                   beta[0],
                                   beta[1],
                                   data_type[idat],
                                   sym_type[isym]);
#if defined(QCMATRIX_ENABLE_VIEW)
                            ierr = QcMatWrite(A, "QcMatGEMM_A", ASCII_VIEW);
                            if (ierr!=QSUCCESS) {
                                printf("test_c_QcMatGEMM>> failed to call QcMatWrite(A)\n");
                                exit(ierr);
                            }
                            ierr = QcMatWrite(B, "QcMatGEMM_B", ASCII_VIEW);
                            if (ierr!=QSUCCESS) {
                                printf("test_c_QcMatGEMM>> failed to call QcMatWrite(B)\n");
                                exit(ierr);
                            }
                            ierr = QcMatWrite(&C, "QcMatGEMM_C", ASCII_VIEW);
                            if (ierr!=QSUCCESS) {
                                printf("test_c_QcMatGEMM>> failed to call QcMatWrite(C)\n");
                                exit(ierr);
                            }
#endif
                            printf("test_c_QcMatGEMM>> failed to call QcMatGEMM(C)\n");
                            exit(ierr);
                        }
                    }
                    if (all_mat_operations[jop]==MAT_HERM_TRANSPOSE ||
                        all_mat_operations[jop]==MAT_COMPLEX_CONJUGATE) {
                        for (ival=0; ival<size_values; ival++) {
                            B_imag[ival] = -B_imag[ival];
                        }
                    }
                }
                if (all_mat_operations[iop]==MAT_HERM_TRANSPOSE ||
                    all_mat_operations[iop]==MAT_COMPLEX_CONJUGATE) {
                    for (ival=0; ival<size_values; ival++) {
                        A_imag[ival] = -A_imag[ival];
                    }
                }
            }
            ierr = QcMatDestroy(&T);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatGEMM>> QcMatDestroy(T) failed\n");
                exit(ierr);
#if defined(DEBUG)
            }
            else {
                printf("test_c_QcMatGEMM>> QcMatDestroy(T) passed ...\n");
#endif
            }
            /* cleans the matrix */
            ierr = QcMatDestroy(&C);
            if (ierr!=QSUCCESS) {
                printf("test_c_QcMatGEMM>> failed to call QcMatDestroy(C)\n");
                exit(ierr);
#if defined(DEBUG)
            }
            else {
                printf("test_c_QcMatGEMM>> QcMatDestroy(C) passed ...\n");
#endif
            }
        }
    }
#if !defined(DEBUG)
    printf("test_c_QcMatGEMM>> QcMatGEMM(beta!=0) passed ...\n");
#endif
    /* cleans */
    free(A_real);
    A_real = NULL;
    free(A_imag);
    A_imag = NULL;
    free(B_real);
    B_real = NULL;
    free(B_imag);
    B_imag = NULL;
    free(C_real);
    C_real = NULL;
    free(C_imag);
    C_imag = NULL;
    free(T_real);
    T_real = NULL;
    free(T_imag);
    T_imag = NULL;
}
