/*
   QcMatrix: an abstract matrix library
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

   This file implements the function QcMatGetExternalMat().

   2014-06-15, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

#if defined(ADAPTER_BLOCK_CMPLX)
#error "QcMatrix does not support external C square block complex matrix"
#elif defined(ADAPTER_BLOCK_REAL)
QErrorCode QcMatGetExternalMat(QcMat *A, const QcDataType data_type, LANG_C_MATRIX **c_A)
{
    RealMat *A_adapter;
    QErrorCode err_code;
    err_code = QcMatGetAdapterMat(A, data_type, &A_adapter);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetAdapterMat");
    err_code = AdapterMatGetExternalMat(A_adapter, c_A);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling AdapterMatGetExternalMat");
    A_adapter = NULL;
    return QSUCCESS;
}
#elif defined(ADAPTER_CMPLX_MAT)
QErrorCode QcMatGetExternalMat(QcMat *A,
                               const QInt idx_block_row,
                               const QInt idx_block_col,
                               LANG_C_MATRIX **c_A)
{
    CmplxMat *A_adapter;
    QErrorCode err_code;
    err_code = QcMatGetAdapterMat(A, idx_block_row, idx_block_col, &A_adapter);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetAdapterMat");
    err_code = AdapterMatGetExternalMat(A_adapter, c_A);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling AdapterMatGetExternalMat");
    A_adapter = NULL;
    return QSUCCESS;
}
#elif defined(ADAPTER_REAL_MAT)
QErrorCode QcMatGetExternalMat(QcMat *A,
                               const QInt idx_block_row,
                               const QInt idx_block_col,
                               const QcDataType data_type,
                               LANG_C_MATRIX **c_A)
{
    RealMat *A_adapter;
    QErrorCode err_code;
    err_code = QcMatGetAdapterMat(A, idx_block_row, idx_block_col, data_type, &A_adapter);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGetAdapterMat");
    err_code = AdapterMatGetExternalMat(A_adapter, c_A);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling AdapterMatGetExternalMat");
    A_adapter = NULL;
    return QSUCCESS;
}
#endif
