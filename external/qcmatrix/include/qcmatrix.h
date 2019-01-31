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

   This file is the header file of square block complex matrix.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(QCMATRIX_H)
#define QCMATRIX_H

/* configure file */
#include "qcmatrix_config.h"

/* external library has implemented square block complex matrix */
#if defined(ADAPTER_BLOCK_CMPLX)
#include "adapter/mat_adapter.h"
#else
/* external library has implemented square block real matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "adapter/mat_adapter.h"
typedef struct {
    QcSymType sym_type;    /* Hermitian, anti-Hermitian or non-Hermitian */
    QcDataType data_type;  /* real, imaginary or complex */
    QInt real_part;        /* pointer to the real part if allocated, default 0 */
    QInt imag_part;        /* pointer to the imaginary part if allocated, default 1 */
    RealMat *cmplx_mat;    /* array with size of 2, for the real and imaginary parts */
} QcMat;
#define CmplxMat QcMat
/* external library has implemented complex matrix */
#elif defined(ADAPTER_CMPLX_MAT)
#include "adapter/mat_adapter.h"
typedef struct {
    QcSymType sym_type;  /* Hermitian, anti-Hermitian or non-Hermitian */
    QInt dim_block;      /* dimension of blocks */
    QBool **assembled;   /* \var(dim_block)x\var(dim_block) array indicating if the blocks are assembled */
    CmplxMat **blocks;   /* \var(dim_block)x\var(dim_block) array for the blocks */
} QcMat;
/* external library has implemented real matrix */
#elif defined(ADAPTER_REAL_MAT)
#include "impls/cmplx_mat.h"
typedef struct {
    QcSymType sym_type;  /* Hermitian, anti-Hermitian or non-Hermitian */
    QInt dim_block;      /* dimension of blocks */
    QBool **assembled;   /* \var(dim_block)x\var(dim_block) array indicating if the blocks are assembled */
    CmplxMat **blocks;   /* \var(dim_block)x\var(dim_block) array for the blocks */
} QcMat;
/* no external matrix library, we use our own simple matrix library */
#else
#include "impls/cmplx_mat.h"
typedef struct {
    QcSymType sym_type;  /* Hermitian, anti-Hermitian or non-Hermitian */
    QInt dim_block;      /* dimension of blocks */
    QBool **assembled;   /* \var(dim_block)x\var(dim_block) array indicating if the blocks are assembled */
    CmplxMat **blocks;   /* \var(dim_block)x\var(dim_block) array for the blocks */
} QcMat;
#endif /* if defined(ADAPTER_BLOCK_REAL) */
/* basic functions of square block complex matrix */
extern QErrorCode QcMatCreate(QcMat*);
extern QErrorCode QcMatBlockCreate(QcMat*,const QInt);
extern QErrorCode QcMatSetSymType(QcMat*,const QcSymType);
extern QErrorCode QcMatSetDataType(QcMat*,
                                   const QInt,
                                   const QInt[],
                                   const QInt[],
                                   const QcDataType[]);
extern QErrorCode QcMatSetDimMat(QcMat*,const QInt,const QInt);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode QcMatSetStorageMode(QcMat*,const QcStorageMode);
#endif
extern QErrorCode QcMatAssemble(QcMat*);
extern QErrorCode QcMatGetDimBlock(QcMat*,QInt*);
extern QErrorCode QcMatGetSymType(QcMat*,QcSymType*);
extern QErrorCode QcMatGetDataType(QcMat*,
                                   const QInt,
                                   const QInt[],
                                   const QInt[],
                                   QcDataType*);
extern QErrorCode QcMatGetDimMat(QcMat*,QInt*,QInt*);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode QcMatGetStorageMode(QcMat*,QcStorageMode*);
#endif
extern QErrorCode QcMatIsAssembled(QcMat*,QBool*);
extern QErrorCode QcMatSetValues(QcMat*,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QReal*,
                                 const QReal*);
extern QErrorCode QcMatAddValues(QcMat*,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QReal*,
                                 const QReal*);
extern QErrorCode QcMatGetValues(QcMat*,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 const QInt,
                                 QReal*,
                                 QReal*);
extern QErrorCode QcMatDuplicate(QcMat*,const QcDuplicateOption,QcMat*);
extern QErrorCode QcMatZeroEntries(QcMat*);
extern QErrorCode QcMatGetTrace(QcMat*,const QInt,QReal*);
extern QErrorCode QcMatGetMatProdTrace(QcMat*,
                                       QcMat*,
                                       const QcMatOperation,
                                       const QInt,
                                       QReal*);
extern QErrorCode QcMatDestroy(QcMat*);
#if defined(QCMATRIX_ENABLE_VIEW)
extern QErrorCode QcMatWrite(QcMat*,const QChar*,const QcViewOption);
extern QErrorCode QcMatRead(QcMat*,const QChar*,const QcViewOption);
#endif
/* functions which invoke BLAS routines */
extern QErrorCode QcMatScale(const QReal[],QcMat*);
extern QErrorCode QcMatAXPY(const QReal[],QcMat*,QcMat*);
extern QErrorCode QcMatTranspose(const QcMatOperation,QcMat*,QcMat*);
/* functions which invoke LAPACK routines */
extern QErrorCode QcMatGEMM(const QcMatOperation,
                            const QcMatOperation,
                            const QReal[],
                            QcMat*,
                            QcMat*,
                            const QReal[],
                            QcMat*);
#endif /* if defined(ADAPTER_BLOCK_CMPLX) */

/* the following functions do not need to implement in external library,
   instead they are provided by QcMatrix */
/* function to retrieve the matrix implemented by the external library */
#if defined(ADAPTER_BLOCK_REAL)
extern QErrorCode QcMatSetAdapterMat(QcMat*,const QcDataType,RealMat**);
extern QErrorCode QcMatGetAdapterMat(QcMat*,const QcDataType,RealMat**);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode QcMatSetAdapterMat(QcMat*,const QInt,const QInt,CmplxMat**);
extern QErrorCode QcMatGetAdapterMat(QcMat*,const QInt,const QInt,CmplxMat**);
#elif defined(ADAPTER_REAL_MAT)
extern QErrorCode QcMatSetAdapterMat(QcMat*,
                                     const QInt,
                                     const QInt,
                                     const QcDataType,
                                     RealMat**);
extern QErrorCode QcMatGetAdapterMat(QcMat*,
                                     const QInt,
                                     const QInt,
                                     const QcDataType,
                                     RealMat**);
#endif
#if defined(ADAPTER_C_LANG)
/* functions to set/get the external C matrix */
#if defined(ADAPTER_BLOCK_REAL)
extern QErrorCode QcMatSetExternalMat(QcMat*,const QcDataType,LANG_C_MATRIX**);
extern QErrorCode QcMatGetExternalMat(QcMat*,const QcDataType,LANG_C_MATRIX**);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode QcMatSetExternalMat(QcMat*,const QInt,const QInt,LANG_C_MATRIX**);
extern QErrorCode QcMatGetExternalMat(QcMat*,const QInt,const QInt,LANG_C_MATRIX**);
#elif defined(ADAPTER_REAL_MAT)
extern QErrorCode QcMatSetExternalMat(QcMat*,
                                      const QInt,
                                      const QInt,
                                      const QcDataType,
                                      LANG_C_MATRIX**);
extern QErrorCode QcMatGetExternalMat(QcMat*,
                                      const QInt,
                                      const QInt,
                                      const QcDataType,
                                      LANG_C_MATRIX**);
#endif
#endif
/* functions relevant to commutator, requires the functions using LAPACK routines */
extern QErrorCode QcMatMatCommutator(QcMat*,QcMat*,QcMat*);
extern QErrorCode QcMatMatSCommutator(QcMat*,QcMat*,QcMat*,QcMat*);
extern QErrorCode QcMatMatHermCommutator(QcMat*,QcMat*,QcMat*);
extern QErrorCode QcMatMatSHermCommutator(QcMat*,QcMat*,QcMat*,QcMat*);
/* functions which may be only used for test suite */
extern QErrorCode QcMatSetRandMat(QcMat*,
                                  const QcSymType,
                                  const QcDataType,
                                  const QInt,
                                  const QInt,
                                  const QInt);
extern QErrorCode QcMatIsEqual(QcMat*,QcMat*,const QBool,QBool*);
extern QErrorCode QcMatCfArray(QcMat*,
                               const QBool,
                               const QInt,
                               const QReal*,
                               const QReal*,
                               QBool*);
extern QErrorCode QcMatGetAllValues(QcMat*,const QBool,const QInt,QReal*,QReal*);
/* functions only for internal use */
#if !defined(ADAPTER_BLOCK_CMPLX) && !defined(ADAPTER_BLOCK_REAL)
extern QErrorCode QcMatIsDiagonal(QcMat*,QBool*);
#endif

#endif /* if !defined(QCMATRIX_H) */
