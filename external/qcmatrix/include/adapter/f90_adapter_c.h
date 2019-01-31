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

#if !defined(F90_ADAPTER_C_H)
#define F90_ADAPTER_C_H

/* configure file */
#include "qcmatrix_config.h"

/* SIZEOF_F_TYPE_P is the size in bytes of derived types with a single pointer member,
   could be changed by -DSIZEOF_F_TYPE_P=XX */
#include "adapter/sizeof_f_type_p.h"
/* Fortran 90 library usually defines a matrix using derived type in a module (default
   is lib_matrix, could be set by -DLANG_F_MODULE=XX), which can not be used directly.
   Instead an array f90_imat is used to save the information of the matrix type, which
   will be converted to a Fortran derived type (with a single pointer member pointing
   to the matrix type) during run time by subroutines in src/adapter/f90_adapter_f.F90 */
#if defined(ADAPTER_BLOCK_CMPLX)
typedef struct {
    QInt f90_imat[SIZEOF_F_TYPE_P];
    QBool external_mat;
} QcMat;
#elif defined(ADAPTER_BLOCK_REAL)
typedef struct {
    QInt f90_imat[SIZEOF_F_TYPE_P];
    QBool external_mat;
} RealMat;
#elif defined(ADAPTER_CMPLX_MAT)
typedef struct {
    QInt f90_imat[SIZEOF_F_TYPE_P];
    QBool external_mat;
} CmplxMat;
#elif defined(ADAPTER_REAL_MAT)
typedef struct {
    QInt f90_imat[SIZEOF_F_TYPE_P];
    QBool external_mat;
} RealMat;
#else
#error "unknown adapter matrix type"
#endif

/* the following two functions set and get the external matrix, and can be useful
   for driver routine and callback functions of application libraries */
#if defined(ADAPTER_BLOCK_CMPLX)
extern QErrorCode AdapterMatSetExternalMat(QcMat*,QInt*);
extern QErrorCode AdapterMatGetExternalMat(QcMat*,QInt*);
#elif defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatSetExternalMat(RealMat*,QInt*);
extern QErrorCode AdapterMatGetExternalMat(RealMat*,QInt*);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode AdapterMatSetExternalMat(CmplxMat*,QInt*);
extern QErrorCode AdapterMatGetExternalMat(CmplxMat*,QInt*);
#elif defined(ADAPTER_REAL_MAT)
extern QErrorCode AdapterMatSetExternalMat(RealMat*,QInt*);
extern QErrorCode AdapterMatGetExternalMat(RealMat*,QInt*);
#else
#error "unknown adapter matrix type"
#endif

#endif
