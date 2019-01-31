# Source codes implemented by QcMatrix, not depends on the external library
SET(QCMATRIX_SRCS
    ${QCMATRIX_SRCS}
    ${LIB_QCMATRIX_PATH}/src/qcmat/QcMatMatCommutator.c
    ${LIB_QCMATRIX_PATH}/src/qcmat/QcMatMatSCommutator.c
    ${LIB_QCMATRIX_PATH}/src/qcmat/QcMatMatHermCommutator.c
    ${LIB_QCMATRIX_PATH}/src/qcmat/QcMatMatSHermCommutator.c)
# Functions which may be only used for test suite
SET(QCMATRIX_SRCS
    ${QCMATRIX_SRCS}
    ${LIB_QCMATRIX_PATH}/src/qcmat/tests/QcMatSetRandMat.c
    ${LIB_QCMATRIX_PATH}/src/qcmat/tests/QcMatIsEqual.c
    ${LIB_QCMATRIX_PATH}/src/qcmat/tests/QcMatCfArray.c
    ${LIB_QCMATRIX_PATH}/src/qcmat/tests/QcMatGetAllValues.c)
