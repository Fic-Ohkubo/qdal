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
   Fortran matrix library.

   2014-06-04, Bin Gao:
   * first version
*/

#if !defined(F03_ADAPTER_C_H)
#define F03_ADAPTER_C_H

/* configure file */
#include "qcmatrix_config.h"

/* The ISO_C_BINDING in Fortran 2003 is used to take care data type conversion and name mangling,
   the pointer f03_mat will be converted to a pointer of matrix type in Fortran 2003 using the
   intrinsic procedure C_F_POINTER during run time */
#if defined(ADAPTER_BLOCK_CMPLX)
typedef struct {
    QVoid *f03_mat;
    QBool external_mat;
} QcMat;
#elif defined(ADAPTER_BLOCK_REAL)
typedef struct {
    QVoid *f03_mat;
    QBool external_mat;
} RealMat;
#elif defined(ADAPTER_CMPLX_MAT)
typedef struct {
    QVoid *f03_mat;
    QBool external_mat;
} CmplxMat;
#elif defined(ADAPTER_REAL_MAT)
typedef struct {
    QVoid *f03_mat;
    QBool external_mat;
} RealMat;
#else
#error "unknown adapter matrix type"
#endif

/* the following two functions set and get the external matrix, and can be useful
   for driver routine and callback functions of application libraries */
#if defined(ADAPTER_BLOCK_CMPLX)
extern QErrorCode AdapterMatSetExternalMat(QcMat*,QVoid**);
extern QErrorCode AdapterMatGetExternalMat(QcMat*,QVoid**);
#elif defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatSetExternalMat(RealMat*,QVoid**);
extern QErrorCode AdapterMatGetExternalMat(RealMat*,QVoid**);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode AdapterMatSetExternalMat(CmplxMat*,QVoid**);
extern QErrorCode AdapterMatGetExternalMat(CmplxMat*,QVoid**);
#elif defined(ADAPTER_REAL_MAT)
extern QErrorCode AdapterMatSetExternalMat(RealMat*,QVoid**);
extern QErrorCode AdapterMatGetExternalMat(RealMat*,QVoid**);
#else
#error "unknown adapter matrix type"
#endif

#endif
