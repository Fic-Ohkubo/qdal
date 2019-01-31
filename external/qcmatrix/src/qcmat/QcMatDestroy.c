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

   This file implements the function QcMatDestroy().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief frees space taken by a matrix
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{delete} A the matrix, should be at least created by QcMatCreate()
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatDestroy(QcMat *A)
{
    QInt irow, icol;
    QErrorCode err_code;
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            err_code = CmplxMatDestroy(&A->blocks[irow][icol]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatDestroy");
        }
    }
    /* prevents repeatedly frees the space */
    if (A->blocks!=NULL) {
        free(A->blocks[0]);
        for (irow=0; irow<A->dim_block; irow++) {
            A->blocks[irow] = NULL;
        }
        free(A->blocks);
        A->blocks = NULL;
    }
    if (A->assembled!=NULL) {
        free(A->assembled[0]);
        for (irow=0; irow<A->dim_block; irow++) {
            A->assembled[irow] = NULL;
        }
        free(A->assembled);
        A->assembled = NULL;
    }
    A->dim_block = 0;
    A->sym_type = QNONSYMMAT;
    return QSUCCESS;
}
