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

   This file is the header file of real matrix.

   2012-04-04, Bin Gao:
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
} RealMat;

/* basic functions of real matrix */
extern QErrorCode RealMatCreate(RealMat*);
#if defined(QCMATRIX_MPI)
extern QErrorCode RealMatSetMPIComm(RealMat*,MPI_comm);
#endif
extern QErrorCode RealMatSetSymType(RealMat*,const QcSymType);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode RealMatSetStorageMode(RealMat*,const QcStorageMode);
#endif
extern QErrorCode RealMatSetDimMat(RealMat*,const QInt,const QInt);
extern QErrorCode RealMatAssemble(RealMat*);
extern QErrorCode RealMatGetSymType(RealMat*,QcSymType*);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode RealMatGetStorageMode(RealMat*,QcStorageMode*);
#endif
extern QErrorCode RealMatGetDimMat(RealMat*,QInt*,QInt*);
extern QErrorCode RealMatIsAssembled(RealMat*,QBool*);
extern QErrorCode RealMatSetValues(RealMat*,
                                   const QInt,
                                   const QInt,
                                   const QInt,
                                   const QInt,
                                   const QReal*);
extern QErrorCode RealMatAddValues(RealMat*,
                                   const QInt,
                                   const QInt,
                                   const QInt,
                                   const QInt,
                                   const QReal*);
extern QErrorCode RealMatGetValues(RealMat*,
                                   const QInt,
                                   const QInt,
                                   const QInt,
                                   const QInt,
                                   QReal*);
extern QErrorCode RealMatDuplicate(RealMat*,const QcDuplicateOption,RealMat*);
extern QErrorCode RealMatZeroEntries(RealMat*);
extern QErrorCode RealMatGetTrace(RealMat*,QReal*);
extern QErrorCode RealMatGetMatProdTrace(RealMat*,RealMat*,const QcMatOperation,QReal*);
extern QErrorCode RealMatDestroy(RealMat*);
#if defined(QCMATRIX_ENABLE_VIEW)
extern QErrorCode RealMatWrite(RealMat*,const QChar*,const QcViewOption);
extern QErrorCode RealMatRead(RealMat*,const QChar*,const QcViewOption);
#endif
/* functions which invoke BLAS routines */
extern QErrorCode RealMatScale(const QReal,RealMat*);
extern QErrorCode RealMatAXPY(const QReal,RealMat*,RealMat*);
extern QErrorCode RealMatTranspose(const QcMatOperation,RealMat*,RealMat*);
/* functions which invoke LAPACK routines */
extern QErrorCode RealMatGEMM(const QcMatOperation,
                              const QcMatOperation,
                              const QReal,
                              RealMat*,
                              RealMat*,
                              const QReal,
                              RealMat*);

#endif /* if !defined(REAL_MAT_H) */
