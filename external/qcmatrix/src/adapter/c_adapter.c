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

   This file implements the functions in adapter/mat_adapter.h.

   2014-06-12, Bin Gao:
   * first version
*/

#include "adapter/mat_adapter.h"

QErrorCode AdapterMatCreate(AdapterMat *A)
{
    QErrorCode err_code;
    A->c_mat = (LANG_C_MATRIX *)malloc(sizeof(LANG_C_MATRIX));
    if (A->c_mat==NULL) {
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for A->c_mat");
    }
    err_code = Matrix_Create(A->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Create");
    A->external_mat = QFALSE;
    return QSUCCESS;
}

QErrorCode AdapterMatDestroy(AdapterMat *A)
{
    QErrorCode err_code;
    if (A->external_mat==QFALSE) {
        err_code = Matrix_Destroy(A->c_mat);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Destroy");
        free(A->c_mat);
        A->c_mat = NULL;
    }
    return QSUCCESS;
}

#if defined(ADAPTER_BLOCK_REAL)
QErrorCode AdapterMatBlockCreate(AdapterMat *A, const QInt dim_block)
{
    QErrorCode err_code;
    err_code = Matrix_BlockCreate(A->c_mat, dim_block);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_BlockCreate");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatSetSymType(AdapterMat *A, const QcSymType sym_type)
{
    QErrorCode err_code;
    err_code = Matrix_SetSymType(A->c_mat, sym_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_SetSymType");
    return QSUCCESS;
}

#if defined(ADAPTER_BLOCK_REAL)
QErrorCode AdapterMatSetNonZeroBlocks(AdapterMat *A,
                                      const QInt num_blocks,
                                      const QInt idx_block_row[],
                                      const QInt idx_block_col[])
{
    QErrorCode err_code;
    err_code = Matrix_SetNonZeroBlocks(A->c_mat,
                                       num_blocks,
                                       idx_block_row,
                                       idx_block_col);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_SetNonZeroBlocks");
    return QSUCCESS;
}
#elif defined(ADAPTER_CMPLX_MAT)
QErrorCode AdapterMatSetDataType(AdapterMat *A, const QcDataType data_type)
{
    QErrorCode err_code;
    err_code = Matrix_SetDataType(A->c_mat, data_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_SetDataType");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatSetDimMat(AdapterMat *A, const QInt num_row, const QInt num_col)
{
    QErrorCode err_code;
    err_code = Matrix_SetDimMat(A->c_mat, num_row, num_col);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_SetDimMat");
    return QSUCCESS;
}

#if defined(QCMATRIX_STORAGE_MODE)
QErrorCode AdapterMatSetStorageMode(AdapterMat *A, const QcStorageMode storage_mode)
{
    QErrorCode err_code;
    err_code = Matrix_SetStorageMode(A->c_mat, storage_mode);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_SetStorageMode");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatAssemble(AdapterMat *A)
{
    QErrorCode err_code;
    err_code = Matrix_Assemble(A->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Assemble");
    return QSUCCESS;
}

#if defined(ADAPTER_BLOCK_REAL)
QErrorCode AdapterMatGetDimBlock(AdapterMat *A, QInt *dim_block)
{
    QErrorCode err_code;
    err_code = Matrix_GetDimBlock(A->c_mat, dim_block);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetDimBlock");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatGetSymType(AdapterMat *A, QcSymType *sym_type)
{
    QErrorCode err_code;
    err_code = Matrix_GetSymType(A->c_mat, sym_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetSymType");
    return QSUCCESS;
}

#if defined(ADAPTER_BLOCK_REAL)
QErrorCode AdapterMatGetNonZeroBlocks(AdapterMat *A,
                                      const QInt num_blocks,
                                      const QInt idx_block_row[],
                                      const QInt idx_block_col[],
                                      QBool *is_non_zero)
{
    QErrorCode err_code;
    err_code = Matrix_GetNonZeroBlocks(A->c_mat,
                                       num_blocks,
                                       idx_block_row,
                                       idx_block_col,
                                       is_non_zero);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetNonZeroBlocks");
    return QSUCCESS;
}
#elif defined(ADAPTER_CMPLX_MAT)
QErrorCode AdapterMatGetDataType(AdapterMat *A, QcDataType *data_type)
{
    QErrorCode err_code;
    err_code = Matrix_GetDataType(A->c_mat, data_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetDataType");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatGetDimMat(AdapterMat *A, QInt *num_row, QInt *num_col)
{
    QErrorCode err_code;
    err_code = Matrix_GetDimMat(A->c_mat, num_row, num_col);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetDimMat");
    return QSUCCESS;
}

#if defined(QCMATRIX_STORAGE_MODE)
QErrorCode AdapterMatGetStorageMode(AdapterMat *A, QcStorageMode *storage_mode)
{
    QErrorCode err_code;
    err_code = Matrix_GetStorageMode(A->c_mat, storage_mode);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetStorageMode");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatIsAssembled(AdapterMat *A, QBool *assembled)
{
    QErrorCode err_code;
    err_code = Matrix_IsAssembled(A->c_mat, assembled);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_IsAssembled");
    return QSUCCESS;
}

#if defined(ADAPTER_BLOCK_REAL)
QErrorCode AdapterMatSetValues(AdapterMat *A,
                               const QInt idx_block_row,
                               const QInt idx_block_col,
                               const QInt idx_first_row,
                               const QInt num_row_set,
                               const QInt idx_first_col,
                               const QInt num_col_set,
                               const QReal *values)
{
    QErrorCode err_code;
    err_code = Matrix_SetValues(A->c_mat,
                                idx_block_row,
                                idx_block_col,
                                idx_first_row,
                                num_row_set,
                                idx_first_col,
                                num_col_set,
                                values);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_SetValues");
    return QSUCCESS;
}

QErrorCode AdapterMatAddValues(AdapterMat *A,
                               const QInt idx_block_row,
                               const QInt idx_block_col,
                               const QInt idx_first_row,
                               const QInt num_row_add,
                               const QInt idx_first_col,
                               const QInt num_col_add,
                               const QReal *values)
{
    QErrorCode err_code;
    err_code = Matrix_AddValues(A->c_mat,
                                idx_block_row,
                                idx_block_col,
                                idx_first_row,
                                num_row_add,
                                idx_first_col,
                                num_col_add,
                                values);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_AddValues");
    return QSUCCESS;
}

QErrorCode AdapterMatGetValues(AdapterMat *A,
                               const QInt idx_block_row,
                               const QInt idx_block_col,
                               const QInt idx_first_row,
                               const QInt num_row_get,
                               const QInt idx_first_col,
                               const QInt num_col_get,
                               QReal *values)
{
    QErrorCode err_code;
    err_code = Matrix_GetValues(A->c_mat,
                                idx_block_row,
                                idx_block_col,
                                idx_first_row,
                                num_row_get,
                                idx_first_col,
                                num_col_get,
                                values);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetValues");
    return QSUCCESS;
}
#elif defined(ADAPTER_CMPLX_MAT)
QErrorCode AdapterMatSetValues(AdapterMat *A,
                               const QInt idx_first_row,
                               const QInt num_row_set,
                               const QInt idx_first_col,
                               const QInt num_col_set,
                               const QReal *values_real,
                               const QReal *values_imag)
{
    QErrorCode err_code;
    err_code = Matrix_SetValues(A->c_mat,
                                idx_first_row,
                                num_row_set,
                                idx_first_col,
                                num_col_set,
                                values_real,
                                values_imag);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_SetValues");
    return QSUCCESS;
}

QErrorCode AdapterMatAddValues(AdapterMat *A,
                               const QInt idx_first_row,
                               const QInt num_row_add,
                               const QInt idx_first_col,
                               const QInt num_col_add,
                               const QReal *values_real,
                               const QReal *values_imag)
{
    QErrorCode err_code;
    err_code = Matrix_AddValues(A->c_mat,
                                idx_first_row,
                                num_row_add,
                                idx_first_col,
                                num_col_add,
                                values_real,
                                values_imag);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_AddValues");
    return QSUCCESS;
}

QErrorCode AdapterMatGetValues(AdapterMat *A,
                               const QInt idx_first_row,
                               const QInt num_row_get,
                               const QInt idx_first_col,
                               const QInt num_col_get,
                               QReal *values_real,
                               QReal *values_imag)
{
    QErrorCode err_code;
    err_code = Matrix_GetValues(A->c_mat,
                                idx_first_row,
                                num_row_get,
                                idx_first_col,
                                num_col_get,
                                values_real,
                                values_imag);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetValues");
    return QSUCCESS;
}
#elif defined(ADAPTER_REAL_MAT)
QErrorCode AdapterMatSetValues(AdapterMat *A,
                               const QInt idx_first_row,
                               const QInt num_row_set,
                               const QInt idx_first_col,
                               const QInt num_col_set,
                               const QReal *values)
{
    QErrorCode err_code;
    err_code = Matrix_SetValues(A->c_mat,
                                idx_first_row,
                                num_row_set,
                                idx_first_col,
                                num_col_set,
                                values);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_SetValues");
    return QSUCCESS;
}

QErrorCode AdapterMatAddValues(AdapterMat *A,
                               const QInt idx_first_row,
                               const QInt num_row_add,
                               const QInt idx_first_col,
                               const QInt num_col_add,
                               const QReal *values)
{
    QErrorCode err_code;
    err_code = Matrix_AddValues(A->c_mat,
                                idx_first_row,
                                num_row_add,
                                idx_first_col,
                                num_col_add,
                                values);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_AddValues");
    return QSUCCESS;
}

QErrorCode AdapterMatGetValues(AdapterMat *A,
                               const QInt idx_first_row,
                               const QInt num_row_get,
                               const QInt idx_first_col,
                               const QInt num_col_get,
                               QReal *values)
{
    QErrorCode err_code;
    err_code = Matrix_GetValues(A->c_mat,
                                idx_first_row,
                                num_row_get,
                                idx_first_col,
                                num_col_get,
                                values);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetValues");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatDuplicate(AdapterMat *A,
                               const QcDuplicateOption duplicate_option,
                               AdapterMat *B)
{
    QErrorCode err_code;
    err_code = Matrix_Duplicate(A->c_mat, duplicate_option, B->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Duplicate");
    return QSUCCESS;
}

QErrorCode AdapterMatZeroEntries(AdapterMat *A)
{
    QErrorCode err_code;
    err_code = Matrix_ZeroEntries(A->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_ZeroEntries");
    return QSUCCESS;
}

#if defined(ADAPTER_BLOCK_REAL)
QErrorCode AdapterMatGetTrace(AdapterMat *A, const QInt num_blocks, QReal *trace)
{
    QErrorCode err_code;
    err_code = Matrix_GetTrace(A->c_mat, num_blocks, trace);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetTrace");
    return QSUCCESS;
}

QErrorCode AdapterMatGetMatProdTrace(AdapterMat *A,
                                     AdapterMat *B,
                                     const QcMatOperation op_B,
                                     const QInt num_blocks,
                                     QReal *trace)
{
    QErrorCode err_code;
    err_code = Matrix_GetMatProdTrace(A->c_mat, B->c_mat, op_B, num_blocks, trace);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetMatProdTrace");
    return QSUCCESS;
}
#elif defined(ADAPTER_CMPLX_MAT)
QErrorCode AdapterMatGetTrace(AdapterMat *A, QReal *trace)
{
    QErrorCode err_code;
    err_code = Matrix_GetTrace(A->c_mat, trace);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetTrace");
    return QSUCCESS;
}

QErrorCode AdapterMatGetMatProdTrace(AdapterMat *A,
                                     AdapterMat *B,
                                     const QcMatOperation op_B,
                                     QReal *trace)
{
    QErrorCode err_code;
    err_code = Matrix_GetMatProdTrace(A->c_mat, B->c_mat, op_B, trace);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetMatProdTrace");
    return QSUCCESS;
}
#elif defined(ADAPTER_REAL_MAT)
QErrorCode AdapterMatGetTrace(AdapterMat *A, QReal *trace)
{
    QErrorCode err_code;
    err_code = Matrix_GetTrace(A->c_mat, trace);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetTrace");
    return QSUCCESS;
}

QErrorCode AdapterMatGetMatProdTrace(AdapterMat *A,
                                     AdapterMat *B,
                                     const QcMatOperation op_B,
                                     QReal *trace)
{
    QErrorCode err_code;
    err_code = Matrix_GetMatProdTrace(A->c_mat, B->c_mat, op_B, trace);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GetMatProdTrace");
    return QSUCCESS;
}
#endif

#if defined(QCMATRIX_ENABLE_VIEW)
QErrorCode AdapterMatWrite(AdapterMat *A,
                           const QChar *mat_label,
                           const QcViewOption view_option)
{
    QErrorCode err_code;
    err_code = Matrix_Write(A->c_mat, mat_label, view_option);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Write");
    return QSUCCESS;
}

QErrorCode AdapterMatRead(AdapterMat *A,
                          const QChar *mat_label,
                          const QcViewOption view_option)
{
    QErrorCode err_code;
    err_code = Matrix_Read(A->c_mat, mat_label, view_option);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Read");
    return QSUCCESS;
}
#endif

#if defined(ADAPTER_CMPLX_MAT)
QErrorCode AdapterMatScale(const QReal scal_number[], AdapterMat *A)
{
    QErrorCode err_code;
    err_code = Matrix_Scale(scal_number, A->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Scale");
    return QSUCCESS;
}
#elif defined(ADAPTER_BLOCK_REAL) || defined(ADAPTER_REAL_MAT)
QErrorCode AdapterMatScale(const QReal scal_number, AdapterMat *A)
{
    QErrorCode err_code;
    err_code = Matrix_Scale(scal_number, A->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Scale");
    return QSUCCESS;
}
#endif

#if defined(ADAPTER_CMPLX_MAT)
QErrorCode AdapterMatAXPY(const QReal multiplier[], AdapterMat *X, AdapterMat *Y)
{
    QErrorCode err_code;
    err_code = Matrix_AXPY(multiplier, X->c_mat, Y->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_AXPY");
    return QSUCCESS;
}
#elif defined(ADAPTER_BLOCK_REAL) || defined(ADAPTER_REAL_MAT)
QErrorCode AdapterMatAXPY(const QReal multiplier, AdapterMat *X, AdapterMat *Y)
{
    QErrorCode err_code;
    err_code = Matrix_AXPY(multiplier, X->c_mat, Y->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_AXPY");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatTranspose(const QcMatOperation op_A, AdapterMat *A, AdapterMat *B)
{
    QErrorCode err_code;
    err_code = Matrix_Transpose(op_A, A->c_mat, B->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_Transpose");
    return QSUCCESS;
}

#if defined(ADAPTER_CMPLX_MAT)
QErrorCode AdapterMatGEMM(const QcMatOperation op_A,
                          const QcMatOperation op_B,
                          const QReal alpha[],
                          AdapterMat *A,
                          AdapterMat *B,
                          const QReal beta[],
                          AdapterMat *C)
{
    QErrorCode err_code;
    err_code = Matrix_GEMM(op_A,
                           op_B,
                           alpha,
                           A->c_mat,
                           B->c_mat,
                           beta,
                           C->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GEMM");
    return QSUCCESS;
}
#elif defined(ADAPTER_BLOCK_REAL) || defined(ADAPTER_REAL_MAT)
QErrorCode AdapterMatGEMM(const QcMatOperation op_A,
                          const QcMatOperation op_B,
                          const QReal alpha,
                          AdapterMat *A,
                          AdapterMat *B,
                          const QReal beta,
                          AdapterMat *C)
{
    QErrorCode err_code;
    err_code = Matrix_GEMM(op_A,
                           op_B,
                           alpha,
                           A->c_mat,
                           B->c_mat,
                           beta,
                           C->c_mat);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling Matrix_GEMM");
    return QSUCCESS;
}
#endif

QErrorCode AdapterMatSetExternalMat(AdapterMat *A, LANG_C_MATRIX **c_A)
{
    QErrorCode err_code;
    err_code = AdapterMatDestroy(A);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling AdapterMatDestroy");
    A->c_mat = *c_A;
    A->external_mat = QTRUE;
    return QSUCCESS;
}

QErrorCode AdapterMatGetExternalMat(AdapterMat *A, LANG_C_MATRIX **c_A)
{
    *c_A = A->c_mat;
    return QSUCCESS;
}
