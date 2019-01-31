/* QcMatrix: square block complex matrix for quantum chemistry calculations
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

   This file implements the ECH-TDDFT method for XPS shake-up satellites.

   2014-04-03, Bin Gao:
   * first version
*/

#include "qcmatrix.h"
#include "xray.h"

/*@% \brief calculates the XPS shake-up satellites using the ECH-TDDFT method
     \author Bin Gao
     \date 2014-04-03
     \return[QErrorCode:int] error information
*/
QErrorCode XPS_ECH_TDDFT(QcMat *A)
{
    return QSUCCESS;
}
