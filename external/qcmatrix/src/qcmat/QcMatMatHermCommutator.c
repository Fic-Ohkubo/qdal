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

   This file implements the function QcMatMatHermCommutator().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief calculates the commutator C = A*B-B*A^{\dagger}
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the left matrix, should be at least assembled
         by QcMatAssemble()
     \param[QcMat:struct]{in} B the right matrix, should be at least assembled
         by QcMatAssemble()
     \param[QcMat:struct]{inout} C the result matrix, should be at least created
         by QcMatCreate()
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatMatHermCommutator(QcMat *A, QcMat *B, QcMat *C)
{
    QReal positive_one[2]={1,0};
    QReal negative_one[2]={-1,0};
    QReal real_zero[2]={0,0};
    QErrorCode err_code;
    /* calculates C = A*B */
    err_code = QcMatGEMM(MAT_NO_OPERATION,
                         MAT_NO_OPERATION,
                         positive_one,
                         A,
                         B,
                         real_zero,
                         C);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGEMM");
    /* calculates C = -B*A^{\dagger}+C */
    err_code = QcMatGEMM(MAT_NO_OPERATION,
                         MAT_HERM_TRANSPOSE,
                         negative_one,
                         B,
                         A,
                         positive_one,
                         C);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatGEMM");
    return QSUCCESS;
}
