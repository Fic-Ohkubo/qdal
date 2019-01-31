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

   This file defines different types for the matrix.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(QCMATRIX_MAT_TYPES_H)
#define QCMATRIX_MAT_TYPES_H

/* configure file */
#include "qcmatrix_config.h"

#include "types/mat_symmetry.h"
#include "types/mat_data.h"
#if defined(QCMATRIX_STORAGE_MODE)
#include "types/mat_storage.h"
#endif
#include "types/mat_duplicate.h"
#if defined(QCMATRIX_ENABLE_VIEW)
#include "types/mat_view.h"
#endif
#include "types/mat_operations.h"

#endif
