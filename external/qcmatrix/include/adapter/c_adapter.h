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
   C matrix library.

   2014-06-04, Bin Gao:
   * first version
*/

#if !defined(C_ADAPTER_H)
#define C_ADAPTER_H

/* error handling*/
#include "utilities/qcmatrix_error.h"

/* we require a header file from external C library (default is lib_matrix.h),
   which defines the matrix struct (default is matrix_t) and functions named
   as Matrix_...() */
#if defined(LANG_C_HEADER)
#include LANG_C_HEADER
#else
#error "a header file from C library required"
#endif

#if defined(ADAPTER_BLOCK_CMPLX)
/* typedef struct {          */
/*     LANG_C_MATRIX *c_mat; */
/*     QBool external_mat;   */
/* } QcMat;                  */
#error "QcMatrix does not support external C square block complex matrix"
#elif defined(ADAPTER_BLOCK_REAL)
typedef struct {
    LANG_C_MATRIX *c_mat;
    QBool external_mat;
} RealMat;
#elif defined(ADAPTER_CMPLX_MAT)
typedef struct {
    LANG_C_MATRIX *c_mat;
    QBool external_mat;
} CmplxMat;
#elif defined(ADAPTER_REAL_MAT)
typedef struct {
    LANG_C_MATRIX *c_mat;
    QBool external_mat;
} RealMat;
#else
#error "unknown adapter matrix type"
#endif

/* the following two functions set and get the external matrix, and can be useful
   for driver routine and callback functions of application libraries */
#if defined(ADAPTER_BLOCK_CMPLX)
/* extern QErrorCode AdapterMatSetExternalMat(QcMat*,LANG_C_MATRIX**); */
/* extern QErrorCode AdapterMatGetExternalMat(QcMat*,LANG_C_MATRIX**); */
#error "QcMatrix does not support external C square block complex matrix"
#elif defined(ADAPTER_BLOCK_REAL)
extern QErrorCode AdapterMatSetExternalMat(RealMat*,LANG_C_MATRIX**);
extern QErrorCode AdapterMatGetExternalMat(RealMat*,LANG_C_MATRIX**);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode AdapterMatSetExternalMat(CmplxMat*,LANG_C_MATRIX**);
extern QErrorCode AdapterMatGetExternalMat(CmplxMat*,LANG_C_MATRIX**);
#elif defined(ADAPTER_REAL_MAT)
extern QErrorCode AdapterMatSetExternalMat(RealMat*,LANG_C_MATRIX**);
extern QErrorCode AdapterMatGetExternalMat(RealMat*,LANG_C_MATRIX**);
#else
#error "unknown adapter matrix type"
#endif

#endif
