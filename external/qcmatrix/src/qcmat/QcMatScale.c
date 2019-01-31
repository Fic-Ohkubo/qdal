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

   This file implements the function QcMatScale().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief scales all elements of a matrix by a given (complex) number
     \author Bin Gao
     \date 2012-04-04
     \param[QReal:real]{in} scal_number the scaling number with scal_number[0] being the
         real part and scal_number[1] the imaginary part
     \param[QcMat:struct]{inout} A the matrix to be scaled, should be at least assembled
         by QcMatAssemble()
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatScale(const QReal scal_number[], QcMat *A)
{
    QInt irow, icol;
    QErrorCode err_code;
    for (irow=0; irow<A->dim_block; irow++) {
        for (icol=0; icol<A->dim_block; icol++) {
            if (A->assembled[irow][icol]==QTRUE) {
                err_code = CmplxMatScale(scal_number, &A->blocks[irow][icol]);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatScale");
            }
        }
    }
    return QSUCCESS;
}
