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

   This file defines the duplicate options.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(MAT_DUPLICATE_H)
#define MAT_DUPLICATE_H

/* duplicate option when calling \fn(QcMatDuplicate), in which COPY_PATTERN_AND_VALUE
   indicates copying an entire matrix including its numerical values, while
   COPY_PATTERN_ONLY copies the pattern of a matrix only (previous numerical values
   may be removed depending on the external library) */
typedef enum {
    COPY_PATTERN_ONLY=0,
    COPY_PATTERN_AND_VALUE=1
} QcDuplicateOption;

#endif
