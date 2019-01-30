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

   This file implements the function RealMatAssemble().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/*% \brief assembles a matrix (e.g. allocating memory) so that it could be used
        in further matrix calculations, this function should be invoked after
        RealMatCreate() and RealMatSet...()
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{inout} A the matrix to be assembled
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatAssemble(RealMat *A)
{
    if (A->values!=NULL) {
        free(A->values);
    }
    A->values = (QReal *)malloc(A->num_row*A->num_col*sizeof(QReal));
    if (A->values==NULL) {
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for A->values");
    }
    return QSUCCESS;
}
