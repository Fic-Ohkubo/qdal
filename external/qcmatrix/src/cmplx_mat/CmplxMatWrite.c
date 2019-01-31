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

   This file implements the function CmplxMatWrite().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatWrite QcMatWrite
#else
#include "impls/cmplx_mat.h"
#endif

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief writes a matrix to file
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{in} A the matrix, should be at least assembled
        by CmplxMatAssemble()
    \param[QChar:char]{in} mat_label label of the matrix, should be unique
    \param[QcViewOption:int]{in} view_option option of writing, see file
        include/types/mat_view.h
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatWrite(CmplxMat *A,
                         const QChar *mat_label,
                         const QcViewOption view_option)
{
    QChar *part_mat_label;    /* label of the real or imaginary parts */
    QSizeT len_mat_label;     /* length of the matrix label */
    QSizeT len_tag_real_mat;  /* length of the tag of the real part */
    QSizeT len_tag_imag_mat;  /* length of the tag of the imaginary part */
    QErrorCode err_code;
#if defined(QCMATRIX_ENABLE_HDF5)
    /* variables for HDF5 library */
    hid_t file_id;        /* identifier of the QCMATRIX_FILE */
    hid_t dataspace_id ;  /* identifier of the data space */
    hid_t dataset_id;     /* identifier of the dataset */
    /*herr_t err_hdf5;*/      /* error code for the HDF5 */
#endif
#if defined(QCMATRIX_ENABLE_MXML)
    /* variables for Mini-XML library */
#endif
#if defined(QCMATRIX_STANDARD_IO)
    /* standard C functions will be used for reading and writing */
    FILE *fp_mat;  /* file pointer */
#endif
    switch (view_option) {
    case BINARY_VIEW:
#if defined(QCMATRIX_ENABLE_HDF5)
        /* opens the QCMATRIX_FILE */
        if (access(QCMATRIX_FILE, F_OK)==0) {
            file_id = H5Fopen(QCMATRIX_FILE, H5F_ACC_RDWR, H5P_DEFAULT);
        }
        /* creates the QCMATRIX_FILE */
        else {
            file_id = H5Fcreate(QCMATRIX_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        /* creates the data space for the dataset */
        dataspace_id = H5Screate(H5S_SCALAR);
        /* the matrix has been written before, so we open the dataset mat_label in group "/" */
        if (H5Lexists(file_id, mat_label, H5P_DEFAULT)!=QFALSE) {
            dataset_id = H5Dopen(file_id, mat_label, H5P_DEFAULT);
        }
        /* the matrix is written for the first time, we create a dataset mat_label in group "/" */
        else {
            dataset_id = H5Dcreate(file_id,
                                   mat_label,
                                   H5T_NATIVE_INT,
                                   dataspace_id,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT);
        }
        /* writes the data type of the complex matrix */
        H5Dwrite(dataset_id,
                 H5T_NATIVE_INT,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &A->data_type);
        /* closes the data space, dataset and the QCMATRIX_FILE */
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
#else
        fp_mat = fopen(mat_label, "wb");
        if (fp_mat==NULL) {
            printf("CmplxMatWrite>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fwrite(&A->data_type, sizeof(A->data_type), 1, fp_mat);
        fclose(fp_mat);
#endif
        break;
    case ASCII_VIEW:
#if defined(QCMATRIX_ENABLE_MXML)
#else
        fp_mat = fopen(mat_label, "w");
        if (fp_mat==NULL) {
            printf("CmplxMatWrite>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fprintf(fp_mat, "%d\n", A->data_type);
        fclose(fp_mat);
#endif
        break;
    default:
        printf("CmplxMatWrite>> view option: %d\n", view_option);
        QErrorExit(FILE_AND_LINE, "invalid view option");
    }
    /* get the labels for the real and imaginary parts of the matrix */
    len_mat_label = strlen(mat_label);
    len_tag_real_mat = strlen(TAG_REAL_MAT)+1;  /* +1 for the zero-terminator */
    len_tag_imag_mat = strlen(TAG_IMAG_MAT)+1;
    part_mat_label = (QChar *)malloc(len_mat_label+QMax(len_tag_real_mat,len_tag_imag_mat));
    if (part_mat_label==NULL) {
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for part_mat_label");
    }
    memcpy(part_mat_label, mat_label, len_mat_label);
    switch (A->data_type) {
    case QREALMAT:
        memcpy(part_mat_label+len_mat_label, TAG_REAL_MAT, len_tag_real_mat);
        err_code = RealMatWrite(&A->cmplx_mat[A->real_part],
                                part_mat_label,
                                view_option);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatWrite");
        break;
    case QIMAGMAT:
        memcpy(part_mat_label+len_mat_label, TAG_IMAG_MAT, len_tag_imag_mat);
        err_code = RealMatWrite(&A->cmplx_mat[A->imag_part],
                                part_mat_label,
                                view_option);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatWrite");
        break;
    case QCMPLXMAT:
        memcpy(part_mat_label+len_mat_label, TAG_REAL_MAT, len_tag_real_mat);
        err_code = RealMatWrite(&A->cmplx_mat[A->real_part],
                                part_mat_label,
                                view_option);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatWrite");
        memcpy(part_mat_label+len_mat_label, TAG_IMAG_MAT, len_tag_imag_mat);
        err_code = RealMatWrite(&A->cmplx_mat[A->imag_part],
                                part_mat_label,
                                view_option);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatWrite");
        break;
    default:
        printf("CmplxMatWrite>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    free(part_mat_label);
    part_mat_label = NULL;
    return QSUCCESS;
}
