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

   This file is the header file of complex matrix.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(CMPLX_MAT_H)
#define CMPLX_MAT_H

/* error handling*/
#include "utilities/qcmatrix_error.h"
/* basic data types */
#include "types/qcmatrix_basic_types.h"
/* data types for matrices */
#include "types/qcmatrix_mat_types.h"

/* BLAS routines */
#include "lapack/qcmatrix_c_blas.h"

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

#if defined(ADAPTER_CMPLX_MAT)
#define Cmplx_Mat_ Matrix_
#define Cmplx_Mat_Create Matrix_Create
#define Cmplx_Mat_SetSymType Matrix_SetSymType
#define Cmplx_Mat_SetDataType Matrix_SetDataType
#define Cmplx_Mat_SetDimMat Matrix_SetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define Cmplx_Mat_SetStorageMode Matrix_SetStorageMode
#endif
#define Cmplx_Mat_Assemble Matrix_Assemble
#define Cmplx_Mat_GetSymType Matrix_GetSymType
#define Cmplx_Mat_GetDataType Matrix_GetDataType
#define Cmplx_Mat_GetDimMat Matrix_GetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define Cmplx_Mat_GetStorageMode Matrix_GetStorageMode
#endif
#define Cmplx_Mat_IsAssembled Matrix_IsAssembled
#define Cmplx_Mat_SetValues Matrix_SetValues
#define Cmplx_Mat_AddValues Matrix_AddValues
#define Cmplx_Mat_GetValues Matrix_GetValues
#define Cmplx_Mat_Duplicate Matrix_Duplicate
#define Cmplx_Mat_ZeroEntries Matrix_ZeroEntries
#define Cmplx_Mat_GetTrace Matrix_GetTrace
#define Cmplx_Mat_GetMatProdTrace Matrix_GetMatProdTrace
#define Cmplx_Mat_Destroy Matrix_Destroy
#if defined(QCMATRIX_ENABLE_VIEW)
#define Cmplx_Mat_Write Matrix_Write
#define Cmplx_Mat_Read Matrix_Read
#endif
#define Cmplx_Mat_Scale Matrix_Scale
#define Cmplx_Mat_AXPY Matrix_AXPY
#define Cmplx_Mat_Transpose Matrix_Transpose
#define Cmplx_Mat_GEMM Matrix_GEMM
#endif

/* defines the complex matrix */
typedef struct {
    QcSymType sym_type;    /* Hermitian, anti-Hermitian or non-Hermitian */
    QcDataType data_type;  /* real, imaginary or complex */
    QInt real_part;        /* pointer to the real part if allocated, default 0 */
    QInt imag_part;        /* pointer to the imaginary part if allocated, default 1 */
    Real_Mat_ *cmplx_mat;  /* array with size of 2, for the real and imaginary parts */
} Cmplx_Mat_;

/* basic functions of complex matrix */
extern QErrorCode Cmplx_Mat_Create(Cmplx_Mat_*);
extern QErrorCode Cmplx_Mat_SetSymType(Cmplx_Mat_*,const QcSymType);
extern QErrorCode Cmplx_Mat_SetDataType(Cmplx_Mat_*,const QcDataType);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode Cmplx_Mat_SetStorageMode(Cmplx_Mat_*,const QcStorageMode);
#endif
extern QErrorCode Cmplx_Mat_SetDimMat(Cmplx_Mat_*,const QInt,const QInt);
extern QErrorCode Cmplx_Mat_Assemble(Cmplx_Mat_*);
extern QErrorCode Cmplx_Mat_GetSymType(Cmplx_Mat_*,QcSymType*);
extern QErrorCode Cmplx_Mat_GetDataType(Cmplx_Mat_*,QcDataType*);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode Cmplx_Mat_GetStorageMode(Cmplx_Mat_*,QcStorageMode*);
#endif
extern QErrorCode Cmplx_Mat_GetDimMat(Cmplx_Mat_*,QInt*,QInt*);
extern QErrorCode Cmplx_Mat_IsAssembled(Cmplx_Mat_*,QBool*);
extern QErrorCode Cmplx_Mat_SetValues(Cmplx_Mat_*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*,
                                      const QReal*);
extern QErrorCode Cmplx_Mat_AddValues(Cmplx_Mat_*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*,
                                      const QReal*);
extern QErrorCode Cmplx_Mat_GetValues(Cmplx_Mat_*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      QReal*,
                                      QReal*);
extern QErrorCode Cmplx_Mat_Duplicate(Cmplx_Mat_*,const QcDuplicateOption,Cmplx_Mat_*);
extern QErrorCode Cmplx_Mat_ZeroEntries(Cmplx_Mat_*);
extern QErrorCode Cmplx_Mat_GetTrace(Cmplx_Mat_*,QReal*);
extern QErrorCode Cmplx_Mat_GetMatProdTrace(Cmplx_Mat_*,
                                            Cmplx_Mat_*,
                                            const QcMatOperation,
                                            QReal*);
extern QErrorCode Cmplx_Mat_Destroy(Cmplx_Mat_*);
#if defined(QCMATRIX_ENABLE_VIEW)
extern QErrorCode Cmplx_Mat_Write(Cmplx_Mat_*,const QChar*,const QcViewOption);
extern QErrorCode Cmplx_Mat_Read(Cmplx_Mat_*,const QChar*,const QcViewOption);
#endif
/* functions which invoke BLAS routines */
extern QErrorCode Cmplx_Mat_Scale(const QReal[],Cmplx_Mat_*);
extern QErrorCode Cmplx_Mat_AXPY(const QReal[],Cmplx_Mat_*,Cmplx_Mat_*);
extern QErrorCode Cmplx_Mat_Transpose(const QcMatOperation,Cmplx_Mat_*,Cmplx_Mat_*);
/* functions which invoke LAPACK routines */
extern QErrorCode Cmplx_Mat_GEMM(const QcMatOperation,
                                 const QcMatOperation,
                                 const QReal[],
                                 Cmplx_Mat_*,
                                 Cmplx_Mat_*,
                                 const QReal[],
                                 Cmplx_Mat_*);

#endif /* if !defined(CMPLX_MAT_H) */
