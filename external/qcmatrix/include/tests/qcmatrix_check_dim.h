#if !defined(QCMATRIX_CHECK_DIM_H)
#define QCMATRIX_CHECK_DIM_H

#include "tests/qcmatrix_test_param.h"
#include "utilities/qcmatrix_error.h"

#define QCheckDimension(dim_block, num_row, num_col, location) \
    if (dim_block>MAX_DIM_BLOCK) { \
        printf("%s>> maximum dimension of blocks %d\n", location, MAX_DIM_BLOCK); \
        printf("%s>> input dimension of blocks %"QINT_FMT"\n", location, dim_block); \
        QErrorExit(location, \
                   "either change tests/qcmatrix_test_param.h or reduce the dimension"); \
    } \
    if (num_row>MAX_DIM_MAT || num_col>MAX_DIM_MAT) { \
        printf("%s>> maximum dimension of each block %d\n", location, MAX_DIM_MAT); \
        printf("%s>> input number of rows %"QINT_FMT"\n", location, num_row); \
        printf("%s>> input number of columns %"QINT_FMT"\n", location, num_col); \
        QErrorExit(location, \
                   "either change tests/qcmatrix_test_param.h or reduce the dimension"); \
    }

#endif
