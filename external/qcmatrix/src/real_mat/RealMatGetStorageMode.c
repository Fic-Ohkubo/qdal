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

   This file implements the function RealMatGetStorageMode().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief gets the matrix storage mode
    \param[RealMat:struct]{in} A the matrix, should be at least created by RealMatCreate()
    \param[QcStorageMode:int]{out} storage_mode return matrix storage mode,
        should be defined and implemented in external library
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatGetStorageMode(RealMat *A, QcStorageMode *storage_mode)
{
    *storage_mode = A->storage_mode;
    return QSUCCESS;
}
