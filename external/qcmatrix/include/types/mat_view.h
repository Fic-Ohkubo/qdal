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

   This file is the header file of functions related to matrix
   writing/reading in QcMatrix.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(MAT_VIEW_H)
#define MAT_VIEW_H

/* configure file */
#include "qcmatrix_config.h"

/* reads and writes a matrix in file using ASCII or binary format,
   should be supported by the external matrix library */
typedef enum {
    BINARY_VIEW=0,
    ASCII_VIEW=1
} QcViewOption;

/* parameters for a matrix on file */
#define TAG_MAT_BLOCK "_block"  /* tag of blocks */
#define TAG_REAL_MAT "_real"    /* tag of the real part of a complex matrix */
#define TAG_IMAG_MAT "_imag"    /* tag of the imaginary part of a complex matrix */
#define TAG_DELIMITER "_"       /* tag of the delimiter */

/* uses strlen() and memcpy() functions */
#include <string.h>
/* uses access() function to check the existence of a file */
#include <unistd.h>

#if defined(QCMATRIX_ENABLE_HDF5)
/* uses HDF5 library to save a matrix */
#include <hdf5.h>
#define QCMATRIX_FILE "qcmatrix.h5"         /* file name for saving the structures of matrices */
#define ATTR_SYM_TYPE "sym_type"          /* attribute name of symmetry type */
#define ATTR_DIM_BLOCK "dim_block"        /* attribute name of the dimension of blocks */
#define ATTR_NUM_ROW "num_row"            /* attribute name of the number of rows of each block */
#define ATTR_NUM_COL "num_col"            /* attribute name of the number of columns of each block */
#if defined(QCMATRIX_STORAGE_MODE)
#define ATTR_STORAGE_MODE "storage_mode"  /* attribute name of the storage mode */
#endif
#endif

#if defined(QCMATRIX_ENABLE_MXML)
/* uses Mini-XML library to save a matrix */
#include <mxml.h>
#endif

#if defined(QCMATRIX_STANDARD_IO)
/* uses standard C functions for reading and writing */
#include <stdio.h>
#endif

#endif
