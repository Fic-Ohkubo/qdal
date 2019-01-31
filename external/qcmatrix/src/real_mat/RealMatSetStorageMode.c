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

   This file implements the function RealMatSetStorageMode().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief sets the matrix storage mode
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{inout} A the matrix, should be created by RealMatCreate()
    \param[QcStorageMode:int]{in} storage_mode given matrix storage mode,
        should be defined and implemented in external library
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatSetStorageMode(RealMat *A, const QcStorageMode storage_mode)
{
    A->storage_mode = storage_mode;
    return QSUCCESS;
}
