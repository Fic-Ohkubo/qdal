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

   This file defines the matrix storage modes.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(MAT_STORAGE_H)
#define MAT_STORAGE_H

/* matrix storage modes, we here only define an unknown storage mode (for error handling),
   while specific storage modes should be defined and implemented in the external library
 */
typedef int QcStorageMode;
#define UNKNOWN_STORAGE_MODE 0

#endif
