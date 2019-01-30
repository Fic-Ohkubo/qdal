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

   This file implements the function QcMatRead().

   2012-04-04, Bin Gao:
   * first version
*/

#include "qcmatrix.h"

/*@% \brief reads a matrix from file
     \author Bin Gao
     \date 2012-04-04
     \param[QcMat:struct]{inout} A the matrix, should be created by QcMatCreate()
     \param[QChar:char]{in} mat_label label of the matrix, should be unique
     \param[QcViewOption:int]{in} view_option option of reading, see file
         include/types/mat_view.h
     \return[QErrorCode:int] error information
*/
QErrorCode QcMatRead(QcMat *A, const QChar *mat_label, const QcViewOption view_option)
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
    FILE *fp_mat;    /* file pointer */
    QInt assembled;  /* read from file and assigned to the QBool type A->assembled */
#endif
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
    QInt num_nz_block;         /* number of non-zero blocks */
    QInt num_row, num_col;     /* dimension of each block */
    QInt irow, icol;
    QErrorCode err_code;
    /* destroys previous information */
    if (A->blocks!=NULL) {
        err_code = QcMatDestroy(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatDestroy");
        err_code = QcMatCreate(A);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatCreate");
    }
    switch (view_option) {
    case BINARY_VIEW:
#if defined(QCMATRIX_ENABLE_HDF5)
        /* opens the QCMATRIX_FILE */
        file_id = H5Fopen(QCMATRIX_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
        /* opens the dataset mat_label in group "/" */
        dataset_id = H5Dopen(file_id, mat_label, H5P_DEFAULT);
        /* reads the dimension of blocks */
        attr_id = H5Aopen(dataset_id, ATTR_DIM_BLOCK, H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &A->dim_block);
        H5Aclose(attr_id);
        /* creates the blocks */
        err_code = QcMatBlockCreate(A, A->dim_block);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatBlockCreate");
        /* reads the symmetry type */
        attr_id = H5Aopen(dataset_id, ATTR_SYM_TYPE, H5P_DEFAULT);
        H5Aread(attr_id, H5T_NATIVE_INT, &A->sym_type);
        H5Aclose(attr_id);
        /* reads the structures of the blocks */
        H5Dread(dataset_id,
                H5T_NATIVE_INT,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                A->assembled[0]);
        /* closes the dataset and the QCMATRIX_FILE */
        H5Dclose(dataset_id);
        H5Fclose(file_id);
#else
        fp_mat = fopen(mat_label, "rb");
        if (fp_mat==NULL) {
            printf("QcMatRead>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        fread(&A->sym_type, sizeof(A->sym_type), 1, fp_mat);
        fread(&A->dim_block, sizeof(A->dim_block), 1, fp_mat);
        /* creates the blocks */
        err_code = QcMatBlockCreate(A, A->dim_block);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatBlockCreate");
        /* A->assembled is a pointer to a pointer, passing its first element (a pointer) to fwrite() */
        fread(A->assembled[0],
              sizeof(A->assembled[0][0]),
              A->dim_block*A->dim_block,
              fp_mat);
        fclose(fp_mat);
#endif
        break;
    case ASCII_VIEW:
#if defined(QCMATRIX_ENABLE_MXML)
#else
        fp_mat = fopen(mat_label, "r");
        if (fp_mat==NULL) {
            printf("QcMatRead>> file: %s\n", mat_label);
            QErrorExit(FILE_AND_LINE, "failed to open the file");
        }
        if (fscanf(fp_mat, "%d\n", &A->sym_type)!=1) {
            QErrorExit(FILE_AND_LINE, "failed to read A->sym_type");
        }
        if (fscanf(fp_mat, "%"QINT_FMT"\n", &A->dim_block)!=1) {
            QErrorExit(FILE_AND_LINE, "failed to read A->dim_block");
        }
        /* creates the blocks */
        err_code = QcMatBlockCreate(A, A->dim_block);
        QErrorCheckCode(err_code, FILE_AND_LINE, "calling QcMatBlockCreate");
        for (irow=0; irow<A->dim_block; irow++) {
            for (icol=0; icol<A->dim_block-1; icol++) {
                if (fscanf(fp_mat, "%"QINT_FMT" ", &assembled)==1) {
                    A->assembled[irow][icol] = (QBool)assembled;
                }
                else {
                    printf("QcMatRead>> block (%"QINT_FMT",%"QINT_FMT")\n", irow, icol);
                    QErrorExit(FILE_AND_LINE, "failed to read assembled");
                }
            }
            if (fscanf(fp_mat, "%"QINT_FMT"\n", &assembled)==1) {
                A->assembled[irow][A->dim_block-1] = (QBool)assembled;
            }
            else {
                printf("QcMatRead>> block (%"QINT_FMT",%"QINT_FMT")\n",
                       irow,
                       A->dim_block-1);
                QErrorExit(FILE_AND_LINE, "failed to read assembled");
            }
        }
        fclose(fp_mat);
#endif
        break;
    default:
        printf("QcMatRead>> view option: %d\n", view_option);
        QErrorExit(FILE_AND_LINE, "invalid view option");
    }
    /* gets the number of digits of the block dimension */
    dim_block = A->dim_block;
    num_digits = 0;
    while (dim_block!=0)  /* as read previously, dim_block == A->dim_block */
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
        printf("QcMatRead>> lenght of label for the blocks %"QINT_FMT"\n",
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
        printf("QcMatRead>> dimension of blocks %"QINT_FMT" (%"QINT_FMT")\n",
               A->dim_block,
               num_digits);
        QErrorExit(FILE_AND_LINE, "failed to allocate memory for row_label");
    }
    /* reads non-zero blocks */
    for (irow=0,num_nz_block=0; irow<A->dim_block; irow++) {
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
                num_nz_block++;
                /* appends the label of this column */
                snprintf(&block_label[len_block_label],
                         len_mat_label,
                         "%"QINT_FMT"",
                         icol);
                err_code = CmplxMatRead(&A->blocks[irow][icol],
                                        block_label,
                                        view_option);
                QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatRead");
            }
        }
        /* changes the length of the label of blocks */
        len_block_label -= len_row_label;
    }
    /* sets the dimension of zero blocks */
    if (num_nz_block<A->dim_block*A->dim_block && num_nz_block>0) {
        num_row = 0;
        num_col = 0;
        for (irow=0; irow<A->dim_block; irow++) {
            for (icol=0; icol<A->dim_block; icol++) {
                if (A->assembled[irow][icol]==QTRUE) {
                    err_code = CmplxMatGetDimMat(&A->blocks[irow][icol],
                                                 &num_row,
                                                 &num_col);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatGetDimMat");
                    break;
                }
            }
            if (num_row!=0 && num_col!=0) break;
        }
        for (irow=0; irow<A->dim_block; irow++) {
            for (icol=0; icol<A->dim_block; icol++) {
                if (A->assembled[irow][icol]==QFALSE) {
                    err_code = CmplxMatSetDimMat(&A->blocks[irow][icol],
                                                 num_row,
                                                 num_col);
                    QErrorCheckCode(err_code, FILE_AND_LINE, "calling CmplxMatSetDimMat");
                }
            }
        }
    }
    /* cleans */
    free(row_label);
    row_label = NULL;
    free(block_label);
    block_label = NULL;
    return QSUCCESS;
}
