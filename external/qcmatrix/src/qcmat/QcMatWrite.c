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

   This file implements the function QcMatWrite().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief writes a matrix to file
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{in} A the matrix, should be at least assembled by QcMatAssemble()
     \param[QChar:char]{in} mat_label label of the matrix, should be unique
     \param[QcViewOption:int]{in} view_option option of writing, see file
         include/types/mat_view.h
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatWrite(QcMat *A, const QChar *mat_label, const QcViewOption view_option)
{
    QInt dim_block;            /* dimension of blocks of the matrix */
    QInt num_digits;           /* number of digits of the block dimension */
    QSizeT len_mat_label;      /* length of the matrix label */
    QSizeT len_tag_delimiter;  /* length of tag of the delimiter */
    QSizeT len_tag_mat_block;  /* length of tag of blocks */
    QSizeT len_block_label;    /* length of the label of the blocks */
    QSizeT len_indices;        /* length of the label of the row and column indices */
    QChar *block_label;        /* label of blocks as mat_label+TAG_MAT_BLOCK %
                                  +TAG_DELIMITER+"row_index"+TAG_DELIMITER+"colum_index" */
    QSizeT len_row_label;      /* length of the label of the rows */
    QChar *row_label;          /* label of rows as "row_index" */
    QInt irow, icol;
    QErrorCode err_code;
#if defined(QCMATRIX_ENABLE_HDF5)
    /* variables for HDF5 library */
    hid_t file_id;        /* identifier of the QCMATRIX_FILE */
    hid_t dataspace_id ;  /* identifier of the data space */
    hid_t dataset_id;     /* identifier of the dataset */
    hid_t aspace_id;      /* identifier of the data space of the attribute */
    hid_t attr_id;        /* identifier of the attribute */
    hsize_t dims[2];      /* dimensions for the dataset */
    /*herr_t err_hdf5;*/      /* error code for the HDF5 */
#endif
#if defined(QCMATRIX_ENABLE_MXML)
    /* variables for Mini-XML library */
#endif
#if defined(QCMATRIX_STANDARD_IO)
    /* standard C functions will be used for reading and writing */
    FILE *fp_mat;  /* file pointer */
#endif
    /* checks the dimension of blocks */
    if (A->dim_block<1) {
        printf("QcMatWrite>> dimension of blocks %"QINT_FMT"\n", A->dim_block);
        QErrorExit(FILE_AND_LINE, "invalid dimension of blocks");
    }
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
        dims[0] = A->dim_block;
        dims[1] = A->dim_block;
        dataspace_id = H5Screate_simple(2, dims, NULL);
        /* creates the data space for the attribute */
        aspace_id = H5Screate(H5S_SCALAR);
        /* the matrix has been written before, so we open the dataset mat_label in group "/"
           and checks the dimension of blocks */
        if (H5Lexists(file_id, mat_label, H5P_DEFAULT)!=QFALSE) {
            dataset_id = H5Dopen(file_id, mat_label, H5P_DEFAULT);
            attr_id = H5Aopen(dataset_id, ATTR_DIM_BLOCK, H5P_DEFAULT);
            H5Aread(attr_id, H5T_NATIVE_INT, &dim_block);
            H5Aclose(attr_id);
            if (dim_block!=A->dim_block) {
                printf("QcMatWrite>> dimension of blocks %"QINT_FMT"\n", A->dim_block);
                printf("QcMatWrite>> dimension of blocks from %s: %"QINT_FMT"\n",
                       QCMATRIX_FILE,
                       dim_block);
                QErrorExit(FILE_AND_LINE, "invalid dimension of blocks");
            }
            /* opens the attribute of symmetry type */
            attr_id = H5Aopen(dataset_id, ATTR_SYM_TYPE, H5P_DEFAULT);
        }
        /* the matrix is written for the first time, we create a dataset mat_label in group "/"
           and writes the information of dimension, symmetry etc. in the attribute */
        else {
            dataset_id = H5Dcreate(file_id,
                                   mat_label,
                                   H5T_NATIVE_INT,
                                   dataspace_id,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT);
            /* creates the attribute of dimension of blocks */
            attr_id = H5Acreate(dataset_id,
                                ATTR_DIM_BLOCK,
                                H5T_NATIVE_INT,
                                aspace_id,
                                H5P_DEFAULT,
                                H5P_DEFAULT);
            H5Awrite(attr_id, H5T_NATIVE_INT, &A->dim_block);
            H5Aclose(attr_id);
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
        /* writes the structures of the blocks */
        H5Dwrite(dataset_id,
                 H5T_NATIVE_INT,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 A->assembled[0]);
        /* closes the data space, dataset and the QCMATRIX_FILE */
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
#else
        fp_mat = fopen(mat_label, "wb");
        if (fp_mat==NULL) {
            printf("QcMatWrite>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fwrite(&A->sym_type, sizeof(A->sym_type), 1, fp_mat);
        fwrite(&A->dim_block, sizeof(A->dim_block), 1, fp_mat);
        /* A->assembled is a pointer to a pointer, passing its first element (a pointer) to fwrite() */
        fwrite(A->assembled[0],
               sizeof(A->assembled[0][0]),
               A->dim_block*A->dim_block,
               fp_mat);
        fclose(fp_mat);
#endif
        break;
    case ASCII_VIEW:
#if defined(QCMATRIX_ENABLE_MXML)
#else
        fp_mat = fopen(mat_label, "w");
        if (fp_mat==NULL) {
            printf("QcMatWrite>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fprintf(fp_mat, "%d\n", A->sym_type);
        fprintf(fp_mat, "%"QINT_FMT"\n", A->dim_block);
        for (irow=0; irow<A->dim_block; irow++) {
            for (icol=0; icol<A->dim_block-1; icol++) {
                fprintf(fp_mat, "%"QINT_FMT" ", (QInt)A->assembled[irow][icol]);
            }
            fprintf(fp_mat, "%"QINT_FMT"\n", (QInt)A->assembled[irow][A->dim_block-1]);
        }
        fclose(fp_mat);
#endif
        break;
    default:
        printf("QcMatWrite>> view option: %d\n", view_option);
        QErrorExit(FILE_AND_LINE, "invalid view option");
    }
    /* gets the number of digits of the dimension of blocks */
    dim_block = A->dim_block;
    num_digits = 0;
    while (dim_block!=0)
    {
        dim_block /= 10;
        num_digits++;
    }
    /* length of the matrix label */
    len_mat_label = strlen(mat_label);
    /* length of the tag of the delimiter, withouth the zero-terminator */
    len_tag_delimiter = strlen(TAG_DELIMITER);
    /* length of the tag of blocks, withouth the zero-terminator */
    len_tag_mat_block = strlen(TAG_MAT_BLOCK);
    /* length of the label for the blocks, +1 for the zero-terminator */
    len_indices = 2*(num_digits+len_tag_delimiter)+1;
    len_block_label = len_mat_label+len_tag_mat_block+len_indices;
    block_label = (QChar *)malloc(len_block_label);
    if (block_label==NULL) {
        printf("QcMatWrite>> lenght of label for the blocks %"QINT_FMT"\n",
               (QInt)len_block_label);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for block_label");
    }
    /* prepares the label for the blocks of the matrix as mat_label+TAG_MAT_BLOCK+TAG_DELIMITER */
    memcpy(block_label, mat_label, len_mat_label);
    memcpy(block_label+len_mat_label, TAG_MAT_BLOCK, len_tag_mat_block);
    len_block_label = len_mat_label+len_tag_mat_block;
    memcpy(block_label+len_block_label, TAG_DELIMITER, len_tag_delimiter);
    len_block_label += len_tag_delimiter;
    /* allocates memory for the label of the rows, +1 for the zero-terminator */
    row_label = (QChar *)malloc(num_digits+len_tag_delimiter+1);
    if (row_label==NULL) {
        printf("QcMatWrite>> dimension of blocks %"QINT_FMT" (%"QINT_FMT")\n",
               A->dim_block,
               num_digits);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for row_label");
    }
    /* writes non-zero blocks */
    for (irow=0; irow<A->dim_block; irow++) {
        /* generates the label of this row "row_index" */
        snprintf(row_label, num_digits+1, "%"QINT_FMT"", irow);
        len_row_label = strlen(row_label);
        memcpy(row_label+len_row_label, TAG_DELIMITER, len_tag_delimiter);
        len_row_label += len_tag_delimiter;
        /* label for the blocks becomes mat_label+TAG_MAT_BLOCK+TAG_DELIMITER+"row_index"+TAG_DELIMITER */
        memcpy(block_label+len_block_label, row_label, len_row_label);
        len_block_label += len_row_label;
        /* left length for the label of the columns */
        len_mat_label = len_indices-len_row_label;
        for (icol=0; icol<A->dim_block; icol++) {
            if (A->assembled[irow][icol]==QTRUE) {
                /* appends the label of this column */
                snprintf(&block_label[len_block_label],
                         len_mat_label,
                         "%"QINT_FMT"",
                         icol);
                err_code = CmplxMatWrite(&A->blocks[irow][icol],
                                         block_label,
                                         view_option);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatWrite");
            }
        }
        /* changes the length of the label of blocks */
        len_block_label -= len_row_label;
    }
    /* cleans */
    free(row_label);
    row_label = NULL;
    free(block_label);
    block_label = NULL;
    return QSUCCESS;
}
