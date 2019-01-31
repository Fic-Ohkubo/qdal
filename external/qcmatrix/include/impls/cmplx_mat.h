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

/* external library has implemented real matrix */
#if defined(ADAPTER_REAL_MAT)
#include "adapter/mat_adapter.h"
/* no external real matrix library, we use our internal matrix library */
#else
#include "impls/real_mat.h"
#endif

/* defines the complex matrix */
typedef struct {
    QcSymType sym_type;    /* Hermitian, anti-Hermitian or non-Hermitian */
    QcDataType data_type;  /* real, imaginary or complex */
    QInt real_part;        /* pointer to the real part if allocated, default 0 */
    QInt imag_part;        /* pointer to the imaginary part if allocated, default 1 */
    RealMat *cmplx_mat;    /* array with size of 2, for the real and imaginary parts */
} CmplxMat;

/* basic functions of complex matrix */
extern QErrorCode CmplxMatCreate(CmplxMat*);
extern QErrorCode CmplxMatSetSymType(CmplxMat*,const QcSymType);
extern QErrorCode CmplxMatSetDataType(CmplxMat*,const QcDataType);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode CmplxMatSetStorageMode(CmplxMat*,const QcStorageMode);
#endif
extern QErrorCode CmplxMatSetDimMat(CmplxMat*,const QInt,const QInt);
extern QErrorCode CmplxMatAssemble(CmplxMat*);
extern QErrorCode CmplxMatGetSymType(CmplxMat*,QcSymType*);
extern QErrorCode CmplxMatGetDataType(CmplxMat*,QcDataType*);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode CmplxMatGetStorageMode(CmplxMat*,QcStorageMode*);
#endif
extern QErrorCode CmplxMatGetDimMat(CmplxMat*,QInt*,QInt*);
extern QErrorCode CmplxMatIsAssembled(CmplxMat*,QBool*);
extern QErrorCode CmplxMatSetValues(CmplxMat*,
                                    const QInt,
                                    const QInt,
                                    const QInt,
                                    const QInt,
                                    const QReal*,
                                    const QReal*);
extern QErrorCode CmplxMatAddValues(CmplxMat*,
                                    const QInt,
                                    const QInt,
                                    const QInt,
                                    const QInt,
                                    const QReal*,
                                    const QReal*);
extern QErrorCode CmplxMatGetValues(CmplxMat*,
                                    const QInt,
                                    const QInt,
                                    const QInt,
                                    const QInt,
                                    QReal*,
                                    QReal*);
extern QErrorCode CmplxMatDuplicate(CmplxMat*,const QcDuplicateOption,CmplxMat*);
extern QErrorCode CmplxMatZeroEntries(CmplxMat*);
extern QErrorCode CmplxMatGetTrace(CmplxMat*,QReal*);
extern QErrorCode CmplxMatGetMatProdTrace(CmplxMat*,
                                          CmplxMat*,
                                          const QcMatOperation,
                                          QReal*);
extern QErrorCode CmplxMatDestroy(CmplxMat*);
#if defined(QCMATRIX_ENABLE_VIEW)
extern QErrorCode CmplxMatWrite(CmplxMat*,const QChar*,const QcViewOption);
extern QErrorCode CmplxMatRead(CmplxMat*,const QChar*,const QcViewOption);
#endif
/* functions which invoke BLAS routines */
extern QErrorCode CmplxMatScale(const QReal[],CmplxMat*);
extern QErrorCode CmplxMatAXPY(const QReal[],CmplxMat*,CmplxMat*);
extern QErrorCode CmplxMatTranspose(const QcMatOperation,CmplxMat*,CmplxMat*);
/* functions which invoke LAPACK routines */
extern QErrorCode CmplxMatGEMM(const QcMatOperation,
                               const QcMatOperation,
                               const QReal[],
                               CmplxMat*,
                               CmplxMat*,
                               const QReal[],
                               CmplxMat*);
#if defined(ADAPTER_REAL_MAT)
/* function to retrieve the matrix implemented by the external library */
extern QErrorCode CmplxMatSetAdapterMat(CmplxMat*,const QcDataType,RealMat**);
extern QErrorCode CmplxMatGetAdapterMat(CmplxMat*,const QcDataType,RealMat**);
#endif

#endif /* if !defined(CMPLX_MAT_H) */
