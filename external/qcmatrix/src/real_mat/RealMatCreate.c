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

   This file implements the function RealMatCreate().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*@% \brief creates the context of a matrix, should be invoked at first
     \author Bin Gao
     \date 2014-06-16
     \param[RealMat:struct]{new} A the matrix
     \return[QErrorCode:int] error information
*/
QErrorCode RealMatCreate(RealMat *A)
{
#if defined(QCMATRIX_STORAGE_MODE)
    A->storage_mode = UNKNOWN_STORAGE_MODE;
#endif
    A->sym_type = QNONSYMMAT;
    A->num_row = 0;
    A->num_col = 0;
    A->values = NULL;
    return QSUCCESS;
}
