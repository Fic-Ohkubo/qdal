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

   This file implements the function CmplxMatAXPY().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatAXPY QcMatAXPY
#define CmplxMatDuplicate QcMatDuplicate
#define CmplxMatScale QcMatScale
#else
#include "impls/cmplx_mat.h"
#endif

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief computes Y = a*X+Y
    \author Bin Gao
    \date 2012-04-04
    \param[QReal:real]{in} multiplier the complex multiplier a with multiplier[0]
        being the real part and multiplier[1] the imaginary part
    \param[CmplxMat:struct]{in} X the first matrix, should be at least assembled
        by CmplxMatAssemble()
    \param[CmplxMat:struct]{inout} Y the second matrix, should be at least created
        by CmplxMatCreate()
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatAXPY(const QReal multiplier[], CmplxMat *X, CmplxMat *Y)
{
    QReal scal_number[2];
    QcDataType num_data_type;  /* data type of the multiplier */
    QErrorCode err_code;  
    /* Y = (a+1)*Y */
    if (X==Y) {
        scal_number[0] = multiplier[0]+1;
        scal_number[1] = multiplier[1];
        err_code = CmplxMatScale(scal_number, Y);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatScale");
    }
    else {
        /* Y is not assembled, we have Y = a*X */
        if (Y->data_type==QNULLMAT) {
            err_code = CmplxMatDuplicate(X, COPY_PATTERN_AND_VALUE, Y);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDuplicate");
            err_code = CmplxMatScale(multiplier, Y);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatScale");
        }
        /* Y is assembled */
        else {
            /* checks the data type of the multiplier */
            if (QAbs(multiplier[0])>QZEROTHRSH) {
                if (QAbs(multiplier[1])>QZEROTHRSH) {
                    num_data_type = QCMPLXMAT;
                }
                else {
                    num_data_type = QREALMAT;
                }
            }
            else if (QAbs(multiplier[1])>QZEROTHRSH) {
                num_data_type = QIMAGMAT;
            }
            /* the multiplier is zero, we have Y = Y */
            else {
                return QSUCCESS;
            }
            switch (X->data_type) {
            case QREALMAT:
                switch (num_data_type) {
                case QREALMAT:
                    /* Y_{R} = a_{R}*X_{R} */
                    if (Y->data_type==QIMAGMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->real_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[0],
                                                &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                        /* sets the data type of Y matrix */
                        Y->data_type = QCMPLXMAT;
                    }
                    /* Y_{R} = a_{R}*X_{R}+Y_{R} */
                    else {
                        err_code = RealMatAXPY(multiplier[0],
                                               &X->cmplx_mat[X->real_part],
                                               &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    break;
                case QIMAGMAT:
                    /* Y_{I} = a_{I}*X_{R} */
                    if (Y->data_type==QREALMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->real_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[1],
                                                &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                        /* sets the data type of Y matrix */
                        Y->data_type = QCMPLXMAT;
                    }
                    /* Y_{I} = a_{I}*X_{R}+Y_{I} */
                    else {
                        err_code = RealMatAXPY(multiplier[1],
                                               &X->cmplx_mat[X->real_part],
                                               &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    break;
                case QCMPLXMAT:
                    /* Y_{R} = a_{R}*X_{R} */
                    if (Y->data_type==QIMAGMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->real_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[0],
                                                &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{R} = a_{R}*X_{R}+Y_{R} */
                    else {
                        err_code = RealMatAXPY(multiplier[0],
                                               &X->cmplx_mat[X->real_part],
                                               &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    /* Y_{I} = a_{I}*X_{R} */
                    if (Y->data_type==QREALMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->real_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[1],
                                                &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{I} = a_{I}*X_{R}+Y_{I} */
                    else {
                        err_code = RealMatAXPY(multiplier[1],
                                               &X->cmplx_mat[X->real_part],
                                               &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    /* sets the data type of Y matrix */
                    Y->data_type = QCMPLXMAT;
                    break;
                default:
                    printf("CmplxMatAXPY>> data type of the number: %d\n",
                           num_data_type);
                    QErrorExit(FILE_AND_LINE, "weird data type");
                }
                break;
            case QIMAGMAT:
                switch (num_data_type) {
                case QREALMAT:
                    /* Y_{I} = a_{R}*X_{I} */
                    if (Y->data_type==QREALMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->imag_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[0],
                                                &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                        /* sets the data type of Y matrix */
                        Y->data_type = QCMPLXMAT;
                    }
                    /* Y_{I} = a_{R}*X_{I}+Y_{I} */
                    else {
                        err_code = RealMatAXPY(multiplier[0],
                                               &X->cmplx_mat[X->imag_part],
                                               &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    break;
                case QIMAGMAT:
                    /* Y_{R} = -a_{I}*X_{I} */
                    if (Y->data_type==QIMAGMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->imag_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(-multiplier[1],
                                                &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                        /* sets the data type of Y matrix */
                        Y->data_type = QCMPLXMAT;
                    }
                    /* Y_{R} = -a_{I}*X_{I}+Y_{R} */
                    else {
                        err_code = RealMatAXPY(-multiplier[1],
                                               &X->cmplx_mat[X->imag_part],
                                               &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    break;
                case QCMPLXMAT:
                    /* Y_{R} = -a_{I}*X_{I} */
                    if (Y->data_type==QIMAGMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->imag_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(-multiplier[1],
                                                &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{R} = -a_{I}*X_{I}+Y_{R} */
                    else {
                        err_code = RealMatAXPY(-multiplier[1],
                                               &X->cmplx_mat[X->imag_part],
                                               &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    /* Y_{I} = a_{R}*X_{I} */
                    if (Y->data_type==QREALMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->imag_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[0],
                                                &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{I} = a_{R}*X_{I}+Y_{I} */
                    else {
                        err_code = RealMatAXPY(multiplier[0],
                                               &X->cmplx_mat[X->imag_part],
                                               &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    /* sets the data type of Y matrix */
                    Y->data_type = QCMPLXMAT;
                    break;
                default:
                    printf("CmplxMatAXPY>> data type of the number: %d\n",
                           num_data_type);
                    QErrorExit(FILE_AND_LINE, "weird data type");
                }
                break;
            case QCMPLXMAT:
                switch (num_data_type) {
                case QREALMAT:
                    /* Y_{R} = a_{R}*X_{R} */
                    if (Y->data_type==QIMAGMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->real_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[0],
                                                &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{R} = a_{R}*X_{R}+Y_{R} */
                    else {
                        err_code = RealMatAXPY(multiplier[0],
                                               &X->cmplx_mat[X->real_part],
                                               &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    /* Y_{I} = a_{R}*X_{I} */
                    if (Y->data_type==QREALMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->imag_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[0],
                                                &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{I} = a_{R}*X_{I}+Y_{I} */
                    else {
                        err_code = RealMatAXPY(multiplier[0],
                                               &X->cmplx_mat[X->imag_part],
                                               &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    break;
                case QIMAGMAT:
                    /* Y_{R} = -a_{I}*X_{I} */
                    if (Y->data_type==QIMAGMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->imag_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(-multiplier[1],
                                                &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{R} = -a_{I}*X_{I}+Y_{R} */
                    else {
                        err_code = RealMatAXPY(-multiplier[1],
                                           &X->cmplx_mat[X->imag_part],
                                           &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    /* Y_{I} = a_{I}*X_{R} */
                    if (Y->data_type==QREALMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->real_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[1],
                                                &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{I} = a_{I}*X_{R}+Y_{I} */
                    else {
                        err_code = RealMatAXPY(multiplier[1],
                                               &X->cmplx_mat[X->real_part],
                                               &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    break;
                case QCMPLXMAT:
                    /* Y_{R} = a_{R}*X_{R} */
                    if (Y->data_type==QIMAGMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->real_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[0],
                                                &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{R} = a_{R}*X_{R}+Y_{R} */
                    else {
                        err_code = RealMatAXPY(multiplier[0],
                                               &X->cmplx_mat[X->real_part],
                                               &Y->cmplx_mat[Y->real_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    /* Y_{R} = -a_{I}*X_{I}+Y_{R} */
                    err_code = RealMatAXPY(-multiplier[1],
                                           &X->cmplx_mat[X->imag_part],
                                           &Y->cmplx_mat[Y->real_part]);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    /* Y_{I} = a_{I}*X_{R} */
                    if (Y->data_type==QREALMAT) {
                        err_code = RealMatDuplicate(&X->cmplx_mat[X->real_part],
                                                    COPY_PATTERN_AND_VALUE,
                                                    &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDuplicate");
                        err_code = RealMatScale(multiplier[1],
                                                &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatScale");
                    }
                    /* Y_{I} = a_{I}*X_{R}+Y_{I} */
                    else {
                        err_code = RealMatAXPY(multiplier[1],
                                               &X->cmplx_mat[X->real_part],
                                               &Y->cmplx_mat[Y->imag_part]);
                        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    }
                    /* Y_{I} = a_{R}*X_{I}+Y_{I} */
                    err_code = RealMatAXPY(multiplier[0],
                                           &X->cmplx_mat[X->imag_part],
                                           &Y->cmplx_mat[Y->imag_part]);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAXPY");
                    break;
                default:
                    printf("CmplxMatAXPY>> data type of the number: %d\n",
                           num_data_type);
                    QErrorExit(FILE_AND_LINE, "weird data type");
                }
                /* sets the data type of Y matrix */
                Y->data_type = QCMPLXMAT;
                break;
            default:
                printf("CmplxMatAXPY>> data type of matrix X: %d\n",
                       X->data_type);
                QErrorExit(FILE_AND_LINE, "invalid data type");
            }
            /* sets the symmetry type of Y matrix */
            if (Y->sym_type!=num_data_type*X->sym_type) Y->sym_type = QNONSYMMAT;
        }
    }
    return QSUCCESS;
}
