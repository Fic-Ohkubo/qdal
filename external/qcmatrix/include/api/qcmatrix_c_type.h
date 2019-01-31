#if !defined(QCMATRIX_C_TYPE)
#define QCMATRIX_C_TYPE

#include "qcmatrix_config.h"

#if defined(QCMATRIX_64BIT_INTEGER)
#define C_QINT C_LONG
#else
#define C_QINT C_INT
#endif

#if defined(QCMATRIX_SINGLE_PRECISION)
#define C_QREAL C_FLOAT
#else
#define C_QREAL C_DOUBLE
#endif

#endif
