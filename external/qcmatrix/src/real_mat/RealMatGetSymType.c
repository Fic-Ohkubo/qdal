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

   This file implements the function RealMatGetSymType().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief gets the symmetry type of a matrix
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{in} A the matrix, should be at least created
        by RealMatCreate()
    \param[QcSymType:int]{out} sym_type symmetry type of the matrix, see file
        include/types/mat_symmetry.h
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatGetSymType(RealMat *A, QcSymType *sym_type)
{
    *sym_type = A->sym_type;
    return QSUCCESS;
}
