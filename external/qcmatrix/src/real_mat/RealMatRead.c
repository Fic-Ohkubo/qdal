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

   This file implements the function RealMatRead().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief reads a matrix from file
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{inout} A the matrix, should be created by RealMatCreate()
    \param[QChar:char]{in} mat_label label of the matrix, should be unique
    \param[QcViewOption:int]{in} view_option option of writing, see file
        include/types/mat_view.h
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatRead(RealMat *A,
                       const QChar *mat_label,
                       const QcViewOption view_option)
{
#if defined(QCMATRIX_ENABLE_HDF5)
    /* variables for HDF5 library */
    hid_t file_id;        /* identifier of the QCMATRIX_FILE */
    hid_t dataset_id;     /* identifier of the dataset */
    hid_t attr_id;        /* identifier of the attribute */
    /*herr_t err_hdf5;*/      /* error code for the HDF5 */
#endif
#if defined(QCMATRIX_ENABLE_MXML)
    /* variables for Mini-XML library */
#endif
#if defined(QCMATRIX_STANDARD_IO)
    /* standard C functions will be used for reading and writing */
    FILE *fp_mat;   /* file pointer */
    QInt size_val;  /* size of values to write */
    QInt ival;      /* incremental recorder over values */
#endif
    QErrorCode err_code;
    /* destroys previous information */
    if (A->values!=NULL) {
        err_code = RealMatDestroy(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatDestroy");
        err_code = RealMatCreate(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatCreate");
    }
    switch (view_option) {
    case BINARY_VIEW:
#if defined(QCMATRIX_ENABLE_HDF5)
        /* opens the QCMATRIX_FILE */
        file_id = H5Fopen(QCMATRIX_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
        /* opens the dataset mat_label in group "/" */
        dataset_id = H5Dopen(file_id, mat_label, H5P_DEFAULT);
        /* reads the numbers of rows and columns */
        attr_id = H5Aopen(dataset_id, ATTR_NUM_ROW, H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &A->num_row);
        H5Aclose(attr_id);
        attr_id = H5Aopen(dataset_id, ATTR_NUM_COL, H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &A->num_col);
        H5Aclose(attr_id);
#if defined(QCMATRIX_STORAGE_MODE)
        /* reads the storage mode */
        attr_id = H5Aopen(dataset_id, ATTR_STORAGE_MODE, H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &A->storage_mode);
        H5Aclose(attr_id);
#endif
        /* reads the symmetry type */
        attr_id = H5Aopen(dataset_id, ATTR_SYM_TYPE, H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &A->sym_type);
        H5Aclose(attr_id);
        /* assembles the matrix */
        err_code = RealMatAssemble(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
        /* reads the values */
        H5Dread(dataset_id,
                H5T_NATIVE_REAL,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                A->values);
        /* closes the dataset and the QCMATRIX_FILE */
        H5Dclose(dataset_id);
        H5Fclose(file_id);
#else
        fp_mat = fopen(mat_label, "rb");
        if (fp_mat==NULL) {
            printf("RealMatRead>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fread(&A->sym_type, sizeof(A->sym_type), 1, fp_mat);
        fread(&A->num_row, sizeof(A->num_row), 1, fp_mat);
        fread(&A->num_col, sizeof(A->num_col), 1, fp_mat);
#if defined(QCMATRIX_STORAGE_MODE)
        fread(&A->storage_mode, sizeof(A->storage_mode), 1, fp_mat);
#endif
        /* assembles the matrix */
        err_code = RealMatAssemble(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
        fread(A->values, sizeof(A->values[0]), A->num_row*A->num_col, fp_mat);
        fclose(fp_mat);
#endif
        break;
    case ASCII_VIEW:
#if defined(QCMATRIX_ENABLE_MXML)
#else
        fp_mat = fopen(mat_label, "r");
        if (fp_mat==NULL) {
            printf("RealMatRead>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        if (fscanf(fp_mat, "%d\n", &A->sym_type)!=1) {
            QErrorExit(FILE_AND_LINE, "failed to read A->sym_type");
        }
        if (fscanf(fp_mat, "%"QINT_FMT"\n", &A->num_row)!=1) {
            QErrorExit(FILE_AND_LINE, "failed to read A->num_row");
        }
        if (fscanf(fp_mat, "%"QINT_FMT"\n", &A->num_col)!=1) {
            QErrorExit(FILE_AND_LINE, "failed to read A->num_col");
        }
#if defined(QCMATRIX_STORAGE_MODE)
        if (fscanf(fp_mat, "%d\n", &A->storage_mode)!=1) {
            QErrorExit(FILE_AND_LINE, "failed to read A->storage_mode");
        }
#endif
        /* assembles the matrix */
        err_code = RealMatAssemble(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling RealMatAssemble");
        size_val = A->num_row*A->num_col;
        for (ival=0; ival<size_val; ival++) {
            if (ival%6==5) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                if (fscanf(fp_mat, "%f\n", &A->values[ival])!=1) {
                    printf("RealMatRead>> value (%"QINT_FMT")\n", ival);
                    QErrorExit(FILE_AND_LINE, "failed to read A->values");
                }
#else
                if (fscanf(fp_mat, "%lf\n", &A->values[ival])!=1) {
                    printf("RealMatRead>> value (%"QINT_FMT")\n", ival);
                    QErrorExit(FILE_AND_LINE, "failed to read A->values");
                }
#endif
            }
            else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                if (fscanf(fp_mat, "%f  ", &A->values[ival])!=1) {
                    printf("RealMatRead>> value (%"QINT_FMT")\n", ival);
                    QErrorExit(FILE_AND_LINE, "failed to read A->values");
                }
#else
                if (fscanf(fp_mat, "%lf  ", &A->values[ival])!=1) {
                    printf("RealMatRead>> value (%"QINT_FMT")\n", ival);
                    QErrorExit(FILE_AND_LINE, "failed to read A->values");
                }
#endif
            }
        }
        fclose(fp_mat);
#endif
        break;
    default:
        printf("RealMatRead>> view option: %d\n", view_option);
        QErrorExit(FILE_AND_LINE, "invalid view option");
    }
    return QSUCCESS;
}
