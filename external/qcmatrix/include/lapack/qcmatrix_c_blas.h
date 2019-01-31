#if !defined(QCMATRIX_C_BLAS_H)
#define QCMATRIX_C_BLAS_H

#include "types/qcmatrix_basic_types.h"

/* declaration of BLAS routines */
extern QVoid C_BLAS_SCAL(const QInt,
                         const QReal,
                         QReal*,
                         const QInt);
extern QVoid C_BLAS_COPY(const QInt,
                         const QReal*,
                         const QInt,
                         QReal*,
                         const QInt);
extern QVoid C_BLAS_AXPY(const QInt,
                         const QReal,
                         const QReal*,
                         const QInt,
                         QReal*,
                         const QInt);
extern QVoid C_BLAS_DOT(const QInt,
                        const QReal*,
                        const QInt,
                        QReal*,
                        const QInt,
                        QReal*);
extern QVoid C_BLAS_GEMM(const QChar,
                         const QChar,
                         const QInt,
                         const QInt,
                         const QInt,
                         const QReal,
                         const QReal*,
                         const QInt,
                         const QReal*,
                         const QInt,
                         const QReal,
                         QReal*,
                         const QInt);

#endif
