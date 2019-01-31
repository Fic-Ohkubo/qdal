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

   This file implements the function RealMatWrite().

   2014-06-16, Bin Gao:
   * first version
*/

#include "impls/real_mat.h"

/* some basic algebraic functions */
#include "utilities/qcmatrix_algebra.h"

/*% \brief writes a matrix to file
    \author Bin Gao
    \date 2014-06-16
    \param[RealMat:struct]{in} A the matrix, should be at least assembled
        by RealMatAssemble()
    \param[QChar:char]{in} mat_label label of the matrix, should be unique
    \param[QcViewOption:int]{in} view_option option of writing, see file
        include/types/mat_view.h
    \return[QErrorCode:int] error information
*/
QErrorCode RealMatWrite(RealMat *A,
                        const QChar *mat_label,
                        const QcViewOption view_option)
{
#if defined(QCMATRIX_ENABLE_HDF5)
    /* variables for HDF5 library */
    QInt num_row, num_col;
#if defined(QCMATRIX_STORAGE_MODE)
    QcStorageMode storage_mode;
#endif
    hid_t file_id;        /* identifier of the QCMATRIX_FILE */
    hid_t dataspace_id ;  /* identifier of the data space */
    hid_t dataset_id;     /* identifier of the dataset */
    hid_t aspace_id;      /* identifier of the data space of the attribute */
    hid_t attr_id;        /* identifier of the attribute */
    hsize_t dims[1];      /* dimensions for the dataset */
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
        dims[0] = A->num_row*A->num_col;
        dataspace_id = H5Screate_simple(1, dims, NULL);
        /* creates the data space for the attribute */
        aspace_id = H5Screate(H5S_SCALAR);
        /* the matrix has been written before, so we open the dataset mat_label in group "/" */
        if (H5Lexists(file_id, mat_label, H5P_DEFAULT)!=QFALSE) {
            dataset_id = H5Dopen(file_id, mat_label, H5P_DEFAULT);
            /* checks the number of rows */
            attr_id = H5Aopen(dataset_id, ATTR_NUM_ROW, H5P_DEFAULT);
            H5Aread(attr_id, H5T_NATIVE_INT, &num_row);
            H5Aclose(attr_id);
            if (num_row!=A->num_row) {
                printf("RealMatWrite>> number of rows %"QINT_FMT"\n", A->num_row);
                printf("RealMatWrite>> number of rows from %s: %"QINT_FMT"\n",
                       QCMATRIX_FILE,
                       num_row);
                QErrorExit(FILE_AND_LINE, "invalid number of rows");
            }
            /* checks the number of columns */
            attr_id = H5Aopen(dataset_id, ATTR_NUM_COL, H5P_DEFAULT);
            H5Aread(attr_id, H5T_NATIVE_INT, &num_col);
            H5Aclose(attr_id);
            if (num_col!=A->num_col) {
                printf("RealMatWrite>> number of columns %"QINT_FMT"\n", A->num_col);
                printf("RealMatWrite>> number of columns from %s: %"QINT_FMT"\n",
                       QCMATRIX_FILE,
                       num_col);
                QErrorExit(FILE_AND_LINE, "invalid number of columns");
            }
#if defined(QCMATRIX_STORAGE_MODE)
            /* checks the storage mode */
            attr_id = H5Aopen(dataset_id, ATTR_STORAGE_MODE, H5P_DEFAULT);
            H5Aread(attr_id, H5T_NATIVE_INT, &storage_mode);
            H5Aclose(attr_id);
            if (storage_mode!=A->storage_mode) {
                printf("RealMatWrite>> storage mode %d\n", A->storage_mode);
                printf("RealMatWrite>> storage mode from %s: %d\n",
                       QCMATRIX_FILE,
                       storage_mode);
                QErrorExit(FILE_AND_LINE, "invalid storage mode");
            }
#endif
            /* opens the attribute of symmetry type */
            attr_id = H5Aopen(dataset_id, ATTR_SYM_TYPE, H5P_DEFAULT);
        }
        /* the matrix is written for the first time, we create a dataset mat_label in group "/"
           and writes the matrix' structure information in the attribute */
        else {
            dataset_id = H5Dcreate(file_id,
                                   mat_label,
                                   H5T_NATIVE_REAL,
                                   dataspace_id,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT);
            /* creates the attribute of number of rows */
            attr_id = H5Acreate(dataset_id,
                                ATTR_NUM_ROW,
                                H5T_NATIVE_INT,
                                aspace_id,
                                H5P_DEFAULT,
                                H5P_DEFAULT);
            H5Awrite(attr_id, H5T_NATIVE_INT, &A->num_row);
            H5Aclose(attr_id);
            /* creates the attribute of number of columns */
            attr_id = H5Acreate(dataset_id,
                                ATTR_NUM_COL,
                                H5T_NATIVE_INT,
                                aspace_id,
                                H5P_DEFAULT,
                                H5P_DEFAULT);
            H5Awrite(attr_id, H5T_NATIVE_INT, &A->num_col);
            H5Aclose(attr_id);
#if defined(QCMATRIX_STORAGE_MODE)
            /* creates the attribute of storage mode */
            attr_id = H5Acreate(dataset_id,
                                ATTR_STORAGE_MODE,
                                H5T_NATIVE_INT,
                                aspace_id,
                                H5P_DEFAULT,
                                H5P_DEFAULT);
            H5Awrite(attr_id, H5T_NATIVE_INT, &A->storage_mode);
            H5Aclose(attr_id);
#endif
            /* creates the attribute of symmetry type */
            attr_id = H5Acreate(dataset_id,
                                ATTR_SYM_TYPE,
                                H5T_NATIVE_INT,
                                aspace_id,
                                H5P_DEFAULT,
                                H5P_DEFAULT);
        }
        /* writes the symmetry type */
        H5Awrite(attr_id, H5T_NATIVE_INT, &A->sym_type);
        H5Aclose(attr_id);
        H5Sclose(aspace_id);
        /* writes the values */
        H5Dwrite(dataset_id,
                 H5T_NATIVE_REAL,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 A->values);
        /* closes the data space, dataset and the QCMATRIX_FILE */
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
#else
        fp_mat = fopen(mat_label, "wb");
        if (fp_mat==NULL) {
            printf("RealMatWrite>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fwrite(&A->sym_type, sizeof(A->sym_type), 1, fp_mat);
        fwrite(&A->num_row, sizeof(A->num_row), 1, fp_mat);
        fwrite(&A->num_col, sizeof(A->num_col), 1, fp_mat);
#if defined(QCMATRIX_STORAGE_MODE)
        fwrite(&A->storage_mode, sizeof(A->storage_mode), 1, fp_mat);
#endif
        fwrite(A->values, sizeof(A->values[0]), A->num_row*A->num_col, fp_mat);
        fclose(fp_mat);
#endif
        break;
    case ASCII_VIEW:
#if defined(QCMATRIX_ENABLE_MXML)
#else
        fp_mat = fopen(mat_label, "w");
        if (fp_mat==NULL) {
            printf("RealMatWrite>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fprintf(fp_mat, "%d\n", A->sym_type);
        fprintf(fp_mat, "%"QINT_FMT"\n", A->num_row);
        fprintf(fp_mat, "%"QINT_FMT"\n", A->num_col);
#if defined(QCMATRIX_STORAGE_MODE)
        fprintf(fp_mat, "%d\n", A->storage_mode);
#endif
        size_val = A->num_row*A->num_col;
        for (ival=0; ival<size_val; ival++) {
            if (ival%6==5) {
#if defined(QCMATRIX_SINGLE_PRECISION)
                fprintf(fp_mat, "%20.12lf\n", A->values[ival]);
#else
                fprintf(fp_mat, "%20.12f\n", A->values[ival]);
#endif
            }
            else {
#if defined(QCMATRIX_SINGLE_PRECISION)
                fprintf(fp_mat, "%20.12lf  ", A->values[ival]);
#else
                fprintf(fp_mat, "%20.12f  ", A->values[ival]);
#endif
            }
        }
        fclose(fp_mat);
#endif
        break;
    default:
        printf("RealMatWrite>> view option: %d\n", view_option);
        QErrorExit(FILE_AND_LINE, "invalid view option");
    }
    return QSUCCESS;
}
