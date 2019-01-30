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

   This is the header file for error handling.

   2013-11-15, Bin Gao:
   * first version
*/

#if !defined(QCMATRIX_ERROR_H)
#define QCMATRIX_ERROR_H

#include <stdlib.h>
#include <stdio.h>

/* configure file */
#include "qcmatrix_config.h"

/* data type used for return error code from all QcMatrix functions */
typedef int QErrorCode;

/* for returning error code from functions */
#if !defined(QSUCCESS)
#define QSUCCESS 0
#endif
#if !defined(QFAILURE)
#define QFAILURE 1
#endif

/* macros for the location (file and line),
   from http://www.cnblogs.com/-clq/archive/2012/02/16/2354385.html */
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define FILE_AND_LINE __FILE__ "(" TOSTRING(__LINE__) ")"

/* functions for error handling */
#if defined(QCMATRIX_AUTO_ERROR_EXIT)
#define QErrorCheckCode(error_code, location, message) \
    if (error_code!=QSUCCESS) { \
        printf(">> error happened at %s: %s\n", location, message); \
        exit(QFAILURE); \
    }
#define QErrorExit(location, message) \
    printf(">> error happened at %s: %s\n", location, message); \
    exit(QFAILURE);
#else
#define QErrorCheckCode(error_code, location, message) \
    if (error_code!=QSUCCESS) { \
        printf(">> error happened at %s: %s\n", location, message); \
        return error_code; \
    }
#define QErrorExit(location, message) \
    printf(">> error happened at %s: %s\n", location, message); \
    return QFAILURE;
#endif

#endif
