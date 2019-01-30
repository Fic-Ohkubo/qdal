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

   This file defines the types of operations on the matrix.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(MAT_OPERATIONS_H)
#define MAT_OPERATIONS_H

typedef enum {
    MAT_NO_OPERATION=0,      /* no matrix operation performed */
    MAT_TRANSPOSE=1,         /* transpose */
    MAT_HERM_TRANSPOSE=2,    /* Hermitian transpose */
    MAT_COMPLEX_CONJUGATE=3  /* complex conjugate */
} QcMatOperation;

#endif
