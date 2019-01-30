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

   This file implements the function QcMatCreate().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief creates the context of a matrix, should be invoked at first
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{new} A the matrix
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatCreate(QcMat *A)
{
    A->sym_type = QNONSYMMAT;
    /* it is important to set A->dim_block as 0, so other functions do not need to
       check if A->assembled or A->blocks is allocated when looping over the blocks */
    A->dim_block = 0;
    A->assembled = NULL;
    A->blocks = NULL;
    return QSUCCESS;
}
