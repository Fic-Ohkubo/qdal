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

   This file implements the function CmplxMatGetTrace().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"

/*% \brief gets the traces of the first few diagonal blocks of a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[QcMat:struct]{in} A the matrix, should be at least assembled by QcMatAssemble()
    \param[QInt:int]{in} num_blocks is the number of diagonal blocks
    \param[QReal:real]{out} trace the traces, size is 2*\var{num_blocks}
    \return[QErrorCode:int] error information
*/
QErrorCode QcMatGetTrace(QcMat *A, const QInt num_blocks, QReal *trace)
{
    QReal *part_trace;
    QInt irow, pos_tr;
    QErrorCode err_code;
    if (num_blocks<1) {
        printf("QcMatGetTrace>> input number of diagonal blocks %"QINT_FMT"\n",
               num_blocks);
        QErrorExit(FILE_AND_LINE, "invalid input number of diagonal blocks");
    }
    /* allocates memory for the traces of real or imaginary part */
    part_trace = (QReal *)malloc(num_blocks*sizeof(QReal));
    if (part_trace==NULL) {
        printf("QcMatGetTrace>> input number of diagonal blocks %"QINT_FMT"\n",
               num_blocks);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for part_trace");
    }
    switch (A->data_type) {
    case QREALMAT:
        err_code = RealMatGetTrace(&A->cmplx_mat[A->real_part], num_blocks, part_trace);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetTrace");
        for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
            trace[pos_tr++] = part_trace[irow];  /* real part */
            trace[pos_tr++] = 0;                 /* imaginary part */
        }
        break;
    case QIMAGMAT:
        err_code = RealMatGetTrace(&A->cmplx_mat[A->imag_part], num_blocks, part_trace);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetTrace");
        for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
            trace[pos_tr++] = 0;                 /* real part */
            trace[pos_tr++] = part_trace[irow];  /* imaginary part */
        }
        break;
    case QCMPLXMAT:
        /* real part */
        err_code = RealMatGetTrace(&A->cmplx_mat[A->real_part], num_blocks, part_trace);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetTrace");
        for (irow=0,pos_tr=0; irow<num_blocks; irow++) {
            trace[pos_tr] = part_trace[irow];
            pos_tr += 2;
        }
        /* imaginary part */
        err_code = RealMatGetTrace(&A->cmplx_mat[A->imag_part], num_blocks, part_trace);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetTrace");
        for (irow=0,pos_tr=1; irow<num_blocks; irow++) {
            trace[pos_tr] = part_trace[irow];  
            pos_tr += 2;
        }
        break;
    default:
        free(part_trace);
        part_trace = NULL;
        printf("QcMatGetTrace>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    free(part_trace);
    part_trace = NULL;
    return QSUCCESS;
}
#else
#include "impls/cmplx_mat.h"

/*% \brief gets the trace of a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{in} A the matrix, should be at least assembled by CmplxMatAssemble()
    \param[QReal:real]{out} trace the trace
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatGetTrace(CmplxMat *A, QReal *trace)
{
    QErrorCode err_code;
    switch (A->data_type) {
    case QREALMAT:
        err_code = RealMatGetTrace(&A->cmplx_mat[A->real_part], &trace[0]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetTrace");
        trace[1] = 0;
        break;
    case QIMAGMAT:
        trace[0] = 0;
        err_code = RealMatGetTrace(&A->cmplx_mat[A->imag_part], &trace[1]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetTrace");
        break;
    case QCMPLXMAT:
        err_code = RealMatGetTrace(&A->cmplx_mat[A->real_part], &trace[0]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetTrace");
        err_code = RealMatGetTrace(&A->cmplx_mat[A->imag_part], &trace[1]);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetTrace");
        break;
    default:
        printf("CmplxMatGetTrace>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    return QSUCCESS;
}
#endif  /* defined(ADAPTER_BLOCK_REAL) */
