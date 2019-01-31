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

   This file implements the function QcMatBlockCreate().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*% \brief sets the dimension of blocks and creates the blocks
    \author Bin Gao
    \date 2012-04-04
    \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate()
    \param[QInt:int]{in} dim_block the dimension of blocks
    \return[QErrorCode:int] error information
*/
QErrorCode QcMatBlockCreate(QcMat *A, const QInt dim_block)
{
    QInt which_part;
    QErrorCode err_code;
    if (dim_block<1) {
        printf("QcMatBlockCreate>> input dimension of blocks %"QINT_FMT"\n", dim_block);
        QErrorExit(FILE_AND_LINE, "invalid dimension of blocks");
    }
    for (which_part=0; which_part<2; which_part++) {
        err_code = RealMatBlockCreate(&A->cmplx_mat[which_part], dim_block);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatBlockCreate");
    }
    return QSUCCESS;
}
