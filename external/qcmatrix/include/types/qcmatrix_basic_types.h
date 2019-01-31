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

   This file defines the basic data types used in QcMatrix.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(QCMATRIX_BASIC_TYPES_H)
#define QCMATRIX_BASIC_TYPES_H

/* definitions of the characteristics of common variable types */
#include <limits.h>

/* configure file */
#include "qcmatrix_config.h"

/* Boolean type */
typedef enum {QFALSE=0, QTRUE=1} QBool;

/* characters */
typedef char QChar;

/* short integer */
typedef short int QShort;
#define QSHRT_MAX SHRT_MAX
#define QSHRT_MIN SHRT_MIN
typedef unsigned short QUShort;
#define QUSHRT_MAX USHRT_MAX
/* basic integer */
#if defined(QCMATRIX_64BIT_INTEGER)
typedef long QInt;
#define QINT_MAX LONG_MAX
#define QINT_MIN LONG_MIN
#define QINT_FMT "ld"
typedef unsigned long QUInt;
#define QUINT_MAX ULONG_MAX
#else
typedef int QInt;
#define QINT_MAX INT_MAX
#define QINT_MIN INT_MIN
#define QINT_FMT "d"
typedef unsigned int QUInt;
#define QUINT_MAX UINT_MAX
#endif
/* long integer */
typedef long QLong;
#define QLONG_MAX LONG_MAX
#define QLONG_MIN LONG_MIN
typedef unsigned long QULong;
#define QULONG_MAX ULONG_MAX

/* real numbers */
#if defined(QCMATRIX_SINGLE_PRECISION)
typedef float QReal;
/* threshold for nearly negligible number */
#define QZEROTHRSH 1.19209290e-07
#if defined(QCMATRIX_ENABLE_HDF5)
/* real numbers of HDF5 library */
#define H5T_NATIVE_REAL H5T_NATIVE_FLOAT
#endif
#else
typedef double QReal;
#define QZEROTHRSH 2.2204460492503131e-16
#if defined(QCMATRIX_ENABLE_HDF5)
#define H5T_NATIVE_REAL H5T_NATIVE_DOUBLE
#endif
#endif

/* size of any object in bytes */
typedef size_t QSizeT;

/* void type */
typedef void QVoid;

#endif
