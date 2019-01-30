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

/*@% \brief sets the dimension of blocks and creates the blocks
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate()
     \param[QInt:int]{in} dim_block the dimension of blocks
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatBlockCreate(QcMat *A, const QInt dim_block)
{
    QInt irow, icol;
    QErrorCode err_code;
    /* checks if the blocks are created before */
    if (A->assembled!=NULL || A->blocks!=NULL) {
        QErrorExit(FILE_AND_LINE, "blocks are already created");
    }
    /* checks the dimension of blocks*/
    if (dim_block<1) {
        printf("QcMatBlockCreate>> input dimension of blocks %"QINT_FMT"\n", dim_block);
        QErrorExit(FILE_AND_LINE, "invalid dimension of blocks");
    }
    A->dim_block = dim_block;
    /* creates an array of Boolean pointers */
    A->assembled = (QBool **)malloc(sizeof(QBool *)*A->dim_block);
    if (A->assembled==NULL) {
        printf("QcMatBlockCreate>> input dimension of blocks %"QINT_FMT"\n", dim_block);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for A->assembled");
    }
    /* allocates memory for the two-dimensional Boolean array at one time */
    A->assembled[0] = (QBool *)malloc(sizeof(QBool)*A->dim_block*A->dim_block);
    if (A->assembled[0]==NULL) {
        printf("QcMatBlockCreate>> input dimension of blocks %"QINT_FMT"\n", dim_block);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for A->assembled[0]");
    }
    /* make the pointer of each row point to the correct memory */
    for (irow=1; irow<A->dim_block; irow++) {
        A->assembled[irow] = A->assembled[irow-1]+A->dim_block;
    }
    /* none of the block will be assembled on default */
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            A->assembled[irow][icol] = QFALSE;
        }
    }
    /* creates an array of complex matrix pointers */
    A->blocks = (CmplxMat **)malloc(sizeof(CmplxMat *)*A->dim_block);
    if (A->blocks==NULL) {
        printf("QcMatBlockCreate>> input dimension of blocks %"QINT_FMT"\n", dim_block);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for A->blocks");
    }
    /* allocates memory for the two-dimensional blocks at one time */
    A->blocks[0] = (CmplxMat *)malloc(sizeof(CmplxMat)*A->dim_block*A->dim_block);
    if (A->blocks[0]==NULL) {
        printf("QcMatBlockCreate>> input dimension of blocks %"QINT_FMT"\n", dim_block);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for A->blocks[0]");
    }
    /* make the pointer of each row point to the correct memory */
    for (irow=1; irow<A->dim_block; irow++) {
        A->blocks[irow] = A->blocks[irow-1]+A->dim_block;
    }
    /* creates the context of each block */
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            err_code = CmplxMatCreate(&A->blocks[irow][icol]);
            QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatCreate");
        }
    }
    return QSUCCESS;
}
