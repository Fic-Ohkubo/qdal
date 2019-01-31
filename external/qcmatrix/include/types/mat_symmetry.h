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

   This file defines the symmetry types of matrices.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(MAT_SYMMETRY_H)
#define MAT_SYMMETRY_H

/* symmetry types */
typedef enum {
    QANTISYMMAT=-1,  /* anti-symmetric (anti-Hermitian) matrix */
    QNONSYMMAT=0,    /* non-symmetric (non-Hermitian) matrix */
    QSYMMAT=1        /* symmetric (Hermitian) matrix */
} QcSymType;

/* if this file is changed, the following files need to modify:
   - include/api/qcmatrix_f_mat_data.h90
   - include/api/qcmatrix_f_mat_symmetry.h90
   - include/types/mat_data.h
   - src/adapter/f03_adapter_c.c
   - src/adapter/f90_adapter_c.c
   - src/cmplx_mat/CmplxMatAXPY.c
   - src/cmplx_mat/CmplxMatGEMM.c
   - src/cmplx_mat/CmplxMatRead.c
   - src/cmplx_mat/CmplxMatScale.c
   - src/cmplx_mat/CmplxMatSetSymType.c
   - src/qcmat/f03/f03_api_c.c
   - src/qcmat/f90/f90_api_c.c
   - src/real_mat/RealMatGetMatProdTrace.c
 */

#endif
