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

   This file implements the function CmplxMatRead().

   2012-04-04, Bin Gao:
   * first version
*/

/* we will implement functions of square block complex matrix if external
   library has implemented real square block matrix */
#if defined(ADAPTER_BLOCK_REAL)
#include "qcmatrix.h"
#define CmplxMatRead QcMatRead
#else
#include "impls/cmplx_mat.h"
#endif

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief reads a matrix from file
    \author Bin Gao
    \date 2012-04-04
    \param[CmplxMat:struct]{inout} A the matrix, should be created by CmplxMatCreate()
    \param[QChar:char]{in} mat_label label of the matrix, should be unique
    \param[QcViewOption:int]{in} view_option option of writing, see file
        include/types/mat_view.h
    \return[QErrorCode:int] error information
*/
QErrorCode CmplxMatRead(CmplxMat *A,
                        const QChar *mat_label,
                        const QcViewOption view_option)
{
#if defined(QCMATRIX_ENABLE_HDF5)
    /* variables for HDF5 library */
    hid_t file_id;        /* identifier of the QCMATRIX_FILE */
    hid_t dataset_id;     /* identifier of the dataset */
    /*herr_t err_hdf5;*/      /* error code for the HDF5 */
#endif
#if defined(QCMATRIX_ENABLE_MXML)
    /* variables for Mini-XML library */
#endif
#if defined(QCMATRIX_STANDARD_IO)
    /* standard C functions will be used for reading and writing */
    FILE *fp_mat;
#endif
    QChar *part_mat_label;    /* label of the real or imaginary parts */
    QSizeT len_mat_label;     /* length of the matrix label */
    QSizeT len_tag_real_mat;  /* length of the tag of the real part */
    QSizeT len_tag_imag_mat;  /* length of the tag of the imaginary part */
    QcSymType imag_sym_type;  /* symmetry of the imaginary part */
    QInt num_row, num_col;    /* dimensions of the real and imaginary parts */
    QErrorCode err_code;
    switch (view_option) {
    case BINARY_VIEW:
#if defined(QCMATRIX_ENABLE_HDF5)
        /* opens the QCMATRIX_FILE */
        file_id = H5Fopen(QCMATRIX_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
        /* opens the dataset mat_label in group "/" */
        dataset_id = H5Dopen(file_id, mat_label, H5P_DEFAULT);
        /* reads the data type of the complex matrix */
        H5Dread(dataset_id,
                H5T_NATIVE_INT,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &A->data_type);
        /* closes the dataset and the QCMATRIX_FILE */
        H5Dclose(dataset_id);
        H5Fclose(file_id);
#else
        fp_mat = fopen(mat_label, "rb");
        if (fp_mat==NULL) {
            printf("CmplxMatRead>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fread(&A->data_type, sizeof(A->data_type), 1, fp_mat);
        fclose(fp_mat);
#endif
        break;
    case ASCII_VIEW:
#if defined(QCMATRIX_ENABLE_MXML)
#else
        fp_mat = fopen(mat_label, "r");
        if (fp_mat==NULL) {
            printf("CmplxMatRead>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        if (fscanf(fp_mat, "%d\n", &A->data_type)!=1) {
            QErrorExit(FILE_AND_LINE, "failed to read A->data_type");
        }
        fclose(fp_mat);
#endif
        break;
    default:
        printf("CmplxMatRead>> view option: %d\n", view_option);
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
    /* reads the real and imaginary parts of the matrix and sets other information */
    switch (A->data_type) {
    case QREALMAT:
        /* reads the real part */
        memcpy(part_mat_label+len_mat_label, TAG_REAL_MAT, len_tag_real_mat);
        err_code = RealMatRead(&A->cmplx_mat[A->real_part],
                               part_mat_label,
                               view_option);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatRead");
        /* gets the symmetry type of the matrix */ 
        err_code = RealMatGetSymType(&A->cmplx_mat[A->real_part], &A->sym_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetSymType");
        /* gets the dimensions of the real part */
        err_code = RealMatGetDimMat(&A->cmplx_mat[A->real_part], &num_row, &num_col);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetDimMat");
        /* sets the dimensions of the imaginary part */
        err_code = RealMatSetDimMat(&A->cmplx_mat[A->imag_part], num_row, num_col);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetDimMat");
        break;
    case QIMAGMAT:
        /* reads the imaginary part */
        memcpy(part_mat_label+len_mat_label, TAG_IMAG_MAT, len_tag_imag_mat);
        err_code = RealMatRead(&A->cmplx_mat[A->imag_part],
                               part_mat_label,
                               view_option);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatRead");
        /* gets the symmetry type of the matrix */
        err_code = RealMatGetSymType(&A->cmplx_mat[A->imag_part], &A->sym_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetSymType");
        A->sym_type = -A->sym_type;
        /* gets the dimensions of the imaginary part */
        err_code = RealMatGetDimMat(&A->cmplx_mat[A->imag_part], &num_row, &num_col);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetDimMat");
        /* sets the dimensions of the real part */
        err_code = RealMatSetDimMat(&A->cmplx_mat[A->real_part], num_row, num_col);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatSetDimMat");
        break;
    case QCMPLXMAT:
        /* reads the real part */
        memcpy(part_mat_label+len_mat_label, TAG_REAL_MAT, len_tag_real_mat);
        err_code = RealMatRead(&A->cmplx_mat[A->real_part],
                               part_mat_label,
                               view_option);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatRead");
        /* gets the symmetry type of the real part */
        err_code = RealMatGetSymType(&A->cmplx_mat[A->real_part], &A->sym_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetSymType");
        /* reads the imaginary part */
        memcpy(part_mat_label+len_mat_label, TAG_IMAG_MAT, len_tag_imag_mat);
        err_code = RealMatRead(&A->cmplx_mat[A->imag_part],
                               part_mat_label,
                               view_option);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatRead");
        /* gets the symmetry type of the imaginary part */
        err_code = RealMatGetSymType(&A->cmplx_mat[A->imag_part], &imag_sym_type);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatGetSymType");
        /* sets the symmetry of the matrix */
        if (A->sym_type!=-imag_sym_type) A->sym_type = QNONSYMMAT;
        break;
    default:
        printf("CmplxMatRead>> data type of matrix A: %d\n", A->data_type);
        QErrorExit(FILE_AND_LINE, "invalid data type");
    }
    free(part_mat_label);
    part_mat_label = NULL;
    return QSUCCESS;
}
