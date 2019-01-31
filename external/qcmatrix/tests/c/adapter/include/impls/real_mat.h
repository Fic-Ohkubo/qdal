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

   This file is the header file of external real matrix library.

   2014-06-17, Bin Gao:
   * first version
*/

#if !defined(REAL_MAT_H)
#define REAL_MAT_H

/* error handling*/
#include "utilities/qcmatrix_error.h"
/* basic data types */
#include "types/qcmatrix_basic_types.h"
/* data types for matrices */
#include "types/qcmatrix_mat_types.h"

/* BLAS routines */
#include "lapack/qcmatrix_c_blas.h"

#if defined(ADAPTER_REAL_MAT)
#define Real_Mat_ Matrix_
#define Real_Mat_Create Matrix_Create
#define Real_Mat_SetSymType Matrix_SetSymType
#define Real_Mat_SetDimMat Matrix_SetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define Real_Mat_SetStorageMode Matrix_SetStorageMode
#endif
#define Real_Mat_Assemble Matrix_Assemble
#define Real_Mat_GetSymType Matrix_GetSymType
#define Real_Mat_GetDimMat Matrix_GetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define Real_Mat_GetStorageMode Matrix_GetStorageMode
#endif
#define Real_Mat_IsAssembled Matrix_IsAssembled
#define Real_Mat_SetValues Matrix_SetValues
#define Real_Mat_AddValues Matrix_AddValues
#define Real_Mat_GetValues Matrix_GetValues
#define Real_Mat_Duplicate Matrix_Duplicate
#define Real_Mat_ZeroEntries Matrix_ZeroEntries
#define Real_Mat_GetTrace Matrix_GetTrace
#define Real_Mat_GetMatProdTrace Matrix_GetMatProdTrace
#define Real_Mat_Destroy Matrix_Destroy
#if defined(QCMATRIX_ENABLE_VIEW)
#define Real_Mat_Write Matrix_Write
#define Real_Mat_Read Matrix_Read
#endif
#define Real_Mat_Scale Matrix_Scale
#define Real_Mat_AXPY Matrix_AXPY
#define Real_Mat_Transpose Matrix_Transpose
#define Real_Mat_GEMM Matrix_GEMM
#endif

#if defined(QCMATRIX_MPI)
/* interface to MPI */
#include <mpi.h>
#endif

/* defines the real matrix */
typedef struct {
#if defined(QCMATRIX_STORAGE_MODE)
    QcStorageMode storage_mode;  /* storage mode, not implemented */
#endif
    QcSymType sym_type;          /* 1 for symmetric, -1 for anti-symmetric */
    QInt num_row;                /* number of rows */
    QInt num_col;                /* number of columns */
    QReal *values;               /* numerical values */
} Real_Mat_;

/* external real matrix library should also implement the following functions */
extern QErrorCode Real_Mat_Create(Real_Mat_*);
extern QErrorCode Real_Mat_Destroy(Real_Mat_*);
extern QErrorCode Real_Mat_SetSymType(Real_Mat_*,const QcSymType);
extern QErrorCode Real_Mat_SetDimMat(Real_Mat_*,const QInt,const QInt);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode Real_Mat_SetStorageMode(Real_Mat_*,const QcStorageMode);
#endif
extern QErrorCode Real_Mat_Assemble(Real_Mat_*);
extern QErrorCode Real_Mat_GetSymType(Real_Mat_*,QcSymType*);
extern QErrorCode Real_Mat_GetDimMat(Real_Mat_*,QInt*,QInt*);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode Real_Mat_GetStorageMode(Real_Mat_*,QcStorageMode*);
#endif
extern QErrorCode Real_Mat_IsAssembled(Real_Mat_*,QBool*);
extern QErrorCode Real_Mat_SetValues(Real_Mat_*,
                                     const QInt,
                                     const QInt,
                                     const QInt,
                                     const QInt,
                                     const QReal*);
extern QErrorCode Real_Mat_AddValues(Real_Mat_*,
                                     const QInt,
                                     const QInt,
                                     const QInt,
                                     const QInt,
                                     const QReal*);
extern QErrorCode Real_Mat_GetValues(Real_Mat_*,
                                     const QInt,
                                     const QInt,
                                     const QInt,
                                     const QInt,
                                     QReal*);
extern QErrorCode Real_Mat_Duplicate(Real_Mat_*,const QcDuplicateOption,Real_Mat_*);
extern QErrorCode Real_Mat_ZeroEntries(Real_Mat_*);
extern QErrorCode Real_Mat_GetTrace(Real_Mat_*,QReal*);
extern QErrorCode Real_Mat_GetMatProdTrace(Real_Mat_*,
                                           Real_Mat_*,
                                           const QcMatOperation,
                                           QReal*);
#if defined(QCMATRIX_ENABLE_VIEW)
extern QErrorCode Real_Mat_Write(Real_Mat_*,const QChar*,const QcViewOption);
extern QErrorCode Real_Mat_Read(Real_Mat_*,const QChar*,const QcViewOption);
#endif
extern QErrorCode Real_Mat_Scale(const QReal,Real_Mat_*);
extern QErrorCode Real_Mat_AXPY(const QReal,Real_Mat_*,Real_Mat_*);
extern QErrorCode Real_Mat_Transpose(const QcMatOperation,Real_Mat_*,Real_Mat_*);
extern QErrorCode Real_Mat_GEMM(const QcMatOperation,
                                const QcMatOperation,
                                const QReal,
                                Real_Mat_*,
                                Real_Mat_*,
                                const QReal,
                                Real_Mat_*);

#endif /* if !defined(REAL_MAT_H) */
