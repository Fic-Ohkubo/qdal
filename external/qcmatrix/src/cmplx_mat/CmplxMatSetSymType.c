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

   This file implements the function CmplxMatSetSymType().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatSetSymType QcMatSetSymType
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief sets the symmetry type of a matrix
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{inout} A the matrix, should be created
        by CmplxMatCreate()
    \param[QcSymType:int]{in} sym_type given symmetry type, see file
        include/types/mat_symmetry.h
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatSetSymType(CmplxMat *A, const QcSymType sym_type)
{
    QErrorCode err_code;
    if (sym_type!=QSYMMAT && sym_type!=QANTISYMMAT && sym_type!=QNONSYMMAT) {
        printf("CmplxMatSetSymType>> input symmetry type %d\n", sym_type);
        QErrorExit(FILE_AND_LINE, "invalid symmetry type");
    }
    err_code = RealMatSetSymType(&A->cmplx_mat[A->real_part], sym_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetSymType");
    /* the imaginar part takes the symmetry with opposite sign */
    err_code = RealMatSetSymType(&A->cmplx_mat[A->imag_part], -sym_type);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetSymType");
    A->sym_type = sym_type;
    return QSUCCESS;
}
