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

   This file is the header file of the adapeter between QcMatrix and external
   matrix library.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(MAT_ADAPTER_H)
#define MAT_ADAPTER_H

/* error handling*/
#include "utilities/qcmatrix_error.h"
/* basic data types */
#include "types/qcmatrix_basic_types.h"
/* data types for matrices */
#include "types/qcmatrix_mat_types.h"

/* adapter header file for the C library */
#if defined(ADAPTER_C_LANG)
#include "adapter/c_adapter.h"
/* adapter header file for the Fortran 90 library */
#elif defined(ADAPTER_F90_LANG)
#include "adapter/f90_adapter_c.h"
/* adapter header file for the Fortran 2003 library */
#elif defined(ADAPTER_F03_LANG)
#include "adapter/f03_adapter_c.h"
#endif

/* renames adapter matrix for external square block complex matrix */
#if defined(ADAPTER_BLOCK_CMPLX)
#define AdapterMat QcMat
/* renames adapter matrix for external square block real matrix */
#elif defined(ADAPTER_BLOCK_REAL)
#define AdapterMat RealMat
/* renames adapter matrix for external complex matrix */
#elif defined(ADAPTER_CMPLX_MAT)
#define AdapterMat CmplxMat
/* renames adapter matrix for external real matrix */
#elif defined(ADAPTER_REAL_MAT)
#define AdapterMat RealMat
#else
#error "unknown adapter matrix type"
#endif

/* renames adapter functions for external square block complex matrix */
#if defined(ADAPTER_BLOCK_CMPLX)
#define AdapterMatCreate QcMatCreate
#define AdapterMatDestroy QcMatDestroy
#define AdapterMatBlockCreate QcMatBlockCreate
#define AdapterMatSetSymType QcMatSetSymType
#define AdapterMatSetDataType QcMatSetDataType
#define AdapterMatSetDimMat QcMatSetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define AdapterMatSetStorageMode QcMatSetStorageMode
#endif
#define AdapterMatAssemble QcMatAssemble
#define AdapterMatGetDimBlock QcMatGetDimBlock
#define AdapterMatGetSymType QcMatGetSymType
#define AdapterMatGetDataType QcMatGetDataType
#define AdapterMatGetDimMat QcMatGetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define AdapterMatGetStorageMode QcMatGetStorageMode
#endif
#define AdapterMatIsAssembled QcMatIsAssembled
#define AdapterMatSetValues QcMatSetValues
#define AdapterMatAddValues QcMatAddValues
#define AdapterMatGetValues QcMatGetValues
#define AdapterMatDuplicate QcMatDuplicate
#define AdapterMatZeroEntries QcMatZeroEntries
#define AdapterMatGetTrace QcMatGetTrace
#define AdapterMatGetMatProdTrace QcMatGetMatProdTrace
#if defined(QCMATRIX_ENABLE_VIEW)
#define AdapterMatWrite QcMatWrite
#define AdapterMatRead QcMatRead
#endif
#define AdapterMatScale QcMatScale
#define AdapterMatAXPY QcMatAXPY
#define AdapterMatTranspose QcMatTranspose
#define AdapterMatGEMM QcMatGEMM
/* renames adapter functions for external square block real matrix */
#elif defined(ADAPTER_BLOCK_REAL)
#define AdapterMatCreate RealMatCreate
#define AdapterMatDestroy RealMatDestroy
#define AdapterMatBlockCreate RealMatBlockCreate
#define AdapterMatSetSymType RealMatSetSymType
#define AdapterMatSetNonZeroBlocks RealMatSetNonZeroBlocks
#define AdapterMatSetDimMat RealMatSetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define AdapterMatSetStorageMode RealMatSetStorageMode
#endif
#define AdapterMatAssemble RealMatAssemble
#define AdapterMatGetDimBlock RealMatGetDimBlock
#define AdapterMatGetSymType RealMatGetSymType
#define AdapterMatGetNonZeroBlocks RealMatGetNonZeroBlocks
#define AdapterMatGetDimMat RealMatGetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define AdapterMatGetStorageMode RealMatGetStorageMode
#endif
#define AdapterMatIsAssembled RealMatIsAssembled
#define AdapterMatSetValues RealMatSetValues
#define AdapterMatAddValues RealMatAddValues
#define AdapterMatGetValues RealMatGetValues
#define AdapterMatDuplicate RealMatDuplicate
#define AdapterMatZeroEntries RealMatZeroEntries
#define AdapterMatGetTrace RealMatGetTrace
#define AdapterMatGetMatProdTrace RealMatGetMatProdTrace
#if defined(QCMATRIX_ENABLE_VIEW)
#define AdapterMatWrite RealMatWrite
#define AdapterMatRead RealMatRead
#endif
#define AdapterMatScale RealMatScale
#define AdapterMatAXPY RealMatAXPY
#define AdapterMatTranspose RealMatTranspose
#define AdapterMatGEMM RealMatGEMM
/* renames adapter functions for external complex matrix */
#elif defined(ADAPTER_CMPLX_MAT)
#define AdapterMatCreate CmplxMatCreate
#define AdapterMatDestroy CmplxMatDestroy
#define AdapterMatSetSymType CmplxMatSetSymType
#define AdapterMatSetDataType CmplxMatSetDataType
#define AdapterMatSetDimMat CmplxMatSetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define AdapterMatSetStorageMode CmplxMatSetStorageMode
#endif
#define AdapterMatAssemble CmplxMatAssemble
#define AdapterMatGetSymType CmplxMatGetSymType
#define AdapterMatGetDataType CmplxMatGetDataType
#define AdapterMatGetDimMat CmplxMatGetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define AdapterMatGetStorageMode CmplxMatGetStorageMode
#endif
#define AdapterMatIsAssembled CmplxMatIsAssembled
#define AdapterMatSetValues CmplxMatSetValues
#define AdapterMatAddValues CmplxMatAddValues
#define AdapterMatGetValues CmplxMatGetValues
#define AdapterMatDuplicate CmplxMatDuplicate
#define AdapterMatZeroEntries CmplxMatZeroEntries
#define AdapterMatGetTrace CmplxMatGetTrace
#define AdapterMatGetMatProdTrace CmplxMatGetMatProdTrace
#if defined(QCMATRIX_ENABLE_VIEW)
#define AdapterMatWrite CmplxMatWrite
#define AdapterMatRead CmplxMatRead
#endif
#define AdapterMatScale CmplxMatScale
#define AdapterMatAXPY CmplxMatAXPY
#define AdapterMatTranspose CmplxMatTranspose
#define AdapterMatGEMM CmplxMatGEMM
/* renames adapter functions for external real matrix */
#elif defined(ADAPTER_REAL_MAT)
#define AdapterMatCreate RealMatCreate
#define AdapterMatDestroy RealMatDestroy
#define AdapterMatSetSymType RealMatSetSymType
#define AdapterMatSetDimMat RealMatSetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define AdapterMatSetStorageMode RealMatSetStorageMode
#endif
#define AdapterMatAssemble RealMatAssemble
#define AdapterMatGetSymType RealMatGetSymType
#define AdapterMatGetDimMat RealMatGetDimMat
#if defined(QCMATRIX_STORAGE_MODE)
#define AdapterMatGetStorageMode RealMatGetStorageMode
#endif
#define AdapterMatIsAssembled RealMatIsAssembled
#define AdapterMatSetValues RealMatSetValues
#define AdapterMatAddValues RealMatAddValues
#define AdapterMatGetValues RealMatGetValues
#define AdapterMatDuplicate RealMatDuplicate
#define AdapterMatZeroEntries RealMatZeroEntries
#define AdapterMatGetTrace RealMatGetTrace
#define AdapterMatGetMatProdTrace RealMatGetMatProdTrace
#if defined(QCMATRIX_ENABLE_VIEW)
#define AdapterMatWrite RealMatWrite
#define AdapterMatRead RealMatRead
#endif
#define AdapterMatScale RealMatScale
#define AdapterMatAXPY RealMatAXPY
#define AdapterMatTranspose RealMatTranspose
#define AdapterMatGEMM RealMatGEMM
#else
#error "unknown adapter matrix type"
#endif

/* declaration of adapter functions */
extern QErrorCode AdapterMatCreate(AdapterMat*);
#if defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatBlockCreate(AdapterMat*,const QInt);
#endif
extern QErrorCode AdapterMatSetSymType(AdapterMat*,const QcSymType);
#if defined(ADAPTER_BLOCK_CMPLX)
extern QErrorCode AdapterMatSetDataType(AdapterMat*,
                                        const QInt,
                                        const QInt[],
                                        const QInt[],
                                        const QcDataType[])
#elif defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatSetNonZeroBlocks(AdapterMat*,
                                             const QInt,
                                             const QInt[],
                                             const QInt[]);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode AdapterMatSetDataType(AdapterMat*,const QcDataType);
#endif
extern QErrorCode AdapterMatSetDimMat(AdapterMat*,const QInt,const QInt);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode AdapterMatSetStorageMode(AdapterMat*,const QcStorageMode);
#endif
extern QErrorCode AdapterMatAssemble(AdapterMat*);
#if defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatGetDimBlock(AdapterMat*,QInt*);
#endif
extern QErrorCode AdapterMatGetSymType(AdapterMat*,QcSymType*);
#if defined(ADAPTER_BLOCK_CMPLX)
extern QErrorCode AdapterMatGetDataType(AdapterMat*,
                                        const QInt,
                                        const QInt[],
                                        const QInt[],
                                        QcDataType*)
#elif defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatGetNonZeroBlocks(AdapterMat*,
                                             const QInt,
                                             const QInt[],
                                             const QInt[],
                                             QBool*);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode AdapterMatGetDataType(AdapterMat*,QcDataType*);
#endif
extern QErrorCode AdapterMatGetDimMat(AdapterMat*,QInt*,QInt*);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode AdapterMatGetStorageMode(AdapterMat*,QcStorageMode*);
#endif
extern QErrorCode AdapterMatIsAssembled(AdapterMat*,QBool*);
#if defined(ADAPTER_BLOCK_CMPLX)
extern QErrorCode AdapterMatSetValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*,
                                      const QReal*)
extern QErrorCode AdapterMatAddValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*,
                                      const QReal*)
extern QErrorCode AdapterMatGetValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      QReal*,
                                      QReal*)
#elif defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatSetValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*);
extern QErrorCode AdapterMatAddValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*);
extern QErrorCode AdapterMatGetValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      QReal*);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode AdapterMatSetValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*,
                                      const QReal*);
extern QErrorCode AdapterMatAddValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*,
                                      const QReal*);
extern QErrorCode AdapterMatGetValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      QReal*,
                                      QReal*);
#else
extern QErrorCode AdapterMatSetValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*);
extern QErrorCode AdapterMatAddValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QReal*);
extern QErrorCode AdapterMatGetValues(AdapterMat*,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      const QInt,
                                      QReal*);
#endif
extern QErrorCode AdapterMatDuplicate(AdapterMat*,const QcDuplicateOption,AdapterMat*);
extern QErrorCode AdapterMatZeroEntries(AdapterMat*);
#if defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatGetTrace(AdapterMat*,const QInt,QReal*);
extern QErrorCode AdapterMatGetMatProdTrace(AdapterMat*,
                                            AdapterMat*,
                                            const QcMatOperation,
                                            const QInt,
                                            QReal*);
#else
extern QErrorCode AdapterMatGetTrace(AdapterMat*,QReal*);
extern QErrorCode AdapterMatGetMatProdTrace(AdapterMat*,
                                            AdapterMat*,
                                            const QcMatOperation,
                                            QReal*);
#endif
extern QErrorCode AdapterMatDestroy(AdapterMat*);
#if defined(QCMATRIX_ENABLE_VIEW)
extern QErrorCode AdapterMatWrite(AdapterMat*,const QChar*,const QcViewOption);
extern QErrorCode AdapterMatRead(AdapterMat*,const QChar*,const QcViewOption);
#endif
/* functions which invoke BLAS routines */
#if defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_CMPLX_MAT)
extern QErrorCode AdapterMatScale(const QReal[],AdapterMat*);
extern QErrorCode AdapterMatAXPY(const QReal[],AdapterMat*,AdapterMat*);
#else
extern QErrorCode AdapterMatScale(const QReal,AdapterMat*);
extern QErrorCode AdapterMatAXPY(const QReal,AdapterMat*,AdapterMat*);
#endif
extern QErrorCode AdapterMatTranspose(const QcMatOperation,AdapterMat*,AdapterMat*);
/* functions which invoke LAPACK routines */
#if defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_CMPLX_MAT)
extern QErrorCode AdapterMatGEMM(const QcMatOperation,
                                 const QcMatOperation,
                                 const QReal[],
                                 AdapterMat*,
                                 AdapterMat*,
                                 const QReal[],
                                 AdapterMat*);
#else
extern QErrorCode AdapterMatGEMM(const QcMatOperation,
                                 const QcMatOperation,
                                 const QReal,
                                 AdapterMat*,
                                 AdapterMat*,
                                 const QReal,
                                 AdapterMat*);
#endif

#endif /* if !defined(MAT_ADAPTER_H) */
