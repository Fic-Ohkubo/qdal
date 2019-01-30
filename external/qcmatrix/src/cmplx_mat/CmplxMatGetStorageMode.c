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

   This file implements the function CmplxMatGetStorageMode().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatGetStorageMode QcMatGetStorageMode
#else
#include "impls/cmplx_mat.h"
#endif

/*% \brief gets the matrix storage mode
    \param[CmplxMat:struct]{in} A the matrix, should be at least created by CmplxMatCreate()
    \param[QcStorageMode:int]{out} storage_mode return matrix storage mode,
        should be defined and implemented in external library
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatGetStorageMode(CmplxMat *A, QcStorageMode *storage_mode)
{
    QcStorageMode imag_storage_mode;
    QErrorCode err_code;
    /* the storage modes of the real and imaginary parts should be the same */
    err_code = RealMatGetStorageMode(&A->cmplx_mat[A->real_part], storage_mode);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetStorageMode");
    err_code = RealMatGetStorageMode(&A->cmplx_mat[A->imag_part], &imag_storage_mode);
    QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetStorageMode");
    if (*storage_mode==UNKNOWN_STORAGE_MODE) {
        *storage_mode = imag_storage_mode;
    }
    else if ((*storage_mode!=imag_storage_mode) &&
             (imag_storage_mode!=UNKNOWN_STORAGE_MODE)) {
        printf("dimension of the real part %d\n", *storage_mode);
        printf("dimension of the imaginary part %d\n", imag_storage_mode);
        QErrorExit(FILE_AND_LINE, "invalid storage mode");
    }
    return QSUCCESS;
}
