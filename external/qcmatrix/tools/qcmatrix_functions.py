#!/usr/bin/env python
#
#  QcMatrix: square block complex matrix for quantum chemistry calculations
#  Copyright 2012-2015 Bin Gao
#
#  QcMatrix is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  QcMatrix is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with QcMatrix. If not, see <http://www.gnu.org/licenses/>.
# 
#  This file defines the functions implemented in QcMatrix.
#
#  2014-02-20, Bin Gao:
#  * taken from qcmatrix_fortran.py

#__author__ = "Bin Gao"
#__copyright__ = "Copyright 2012-2015"
#__license__ = "LGPLv3"
#__version__ = "0.1.0"
#__maintainer__ = "Bin Gao"
#__email__ = "bin.gao@uit.no"
#__status__ = "Development"

# functions of the adapter that require memory allocation for the matrix
ADAPTER_MALLOC_FUN = [["MatCreate(Mat *A)", "out"]]
# functions of the adapter that require memory deallocation for the matrix
ADAPTER_FREE_FUN = [["MatDestroy(Mat *A)", "inout"]]
# other functions of the adapter that do not require any memory operation
ADAPTER_NOMEM_FUN = [["MatBlockCreate(Mat *A, QInt dim_block)",                                                               \
                      "in, in",                                                                                               \
                      "defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_BLOCK_REAL)"],                                         \
                     ["MatSetSymType(Mat *A, QSymType sym_type)",                                                             \
                      "in, in"],                                                                                              \
                     ["MatSetDataType(Mat *A, QDataType data_type)",                                                          \
                      "in, in",                                                                                               \
                      "defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_CMPLX_MAT)"],                                          \
                     ["MatSetDimMat(Mat *A, QInt dim_mat)",                                                                   \
                      "in, in"],                                                                                              \
                     ["MatSetStorageMode(Mat *A, QStorageMode storage_mode)",                                                 \
                      "in, in",                                                                                               \
                      "defined(QCMATRIX_STORAGE_MODE)"],                                                                        \
                     ["MatAssemble(Mat *A)",                                                                                  \
                      "in"],                                                                                                  \
                     ["MatGetDimBlock(Mat *A, QInt *dim_block)",                                                              \
                      "in, out",                                                                                              \
                      "defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_BLOCK_REAL)"],                                         \
                     ["MatGetSymType(Mat *A, QSymType *sym_type)",                                                            \
                      "in, out"],                                                                                             \
                     ["MatGetDataType(Mat *A, QDataType *data_type)",                                                         \
                      "in, out",                                                                                              \
                      "defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_CMPLX_MAT)"],                                          \
                     ["MatGetDimMat(Mat *A, QInt *dim_mat)",                                                                  \
                      "in, out"],                                                                                             \
                     ["MatGetStorageMode(Mat *A, QStorageMode *storage_mode)",                                                \
                      "in, out",                                                                                              \
                      "defined(QCMATRIX_STORAGE_MODE)"],                                                                        \
                     ["MatIsAssembled(Mat *A, QInt *assembled)",                                                              \
                      "in, out"],                                                                                             \
                     ["MatDuplicate(Mat *A, QDuplicateOption duplicate_option, Mat *B)",                                      \
                      "in, in, in"],                                                                                          \
                     ["MatZeroEntries(Mat *A)",                                                                               \
                      "in"],                                                                                                  \
                     ["MatGetTrace(Mat *A, QInt num_blocks, QReal *trace)",                                                   \
                      "in, in, out(2*num_blocks)",                                                                            \
                      "defined(ADAPTER_BLOCK_CMPLX)",                                                                         \
                      "MatGetTrace(Mat *A, QInt num_blocks, QReal *trace)",                                                   \
                      "in, in, out(num_blocks)",                                                                              \
                      "defined(ADAPTER_BLOCK_REAL)",                                                                          \
                      "MatGetTrace(Mat *A, QReal *trace)",                                                                    \
                      "in, out(2)",                                                                                           \
                      "defined(ADAPTER_CMPLX_MAT)",                                                                           \
                      "MatGetTrace(Mat *A, QReal *trace)",                                                                    \
                      "in, out"],                                                                                             \
                     ["MatWrite(Mat *A, QChar *mat_label, QViewOption view_option)",                                          \
                      "in, in, in",                                                                                           \
                      "defined(QCMATRIX_ENABLE_VIEW)"],                                                                         \
                     ["MatRead(Mat *A, QChar *mat_label, QViewOption view_option)",                                           \
                      "in, in, in",                                                                                           \
                      "defined(QCMATRIX_ENABLE_VIEW)"],                                                                         \
                     ["MatScale(QReal scal_number[], Mat *A)",                                                                \
                      "in(2), in",                                                                                            \
                      "defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_CMPLX_MAT)",                                           \
                      "MatScale(QReal scal_number, Mat *A)",                                                                  \
                      "in, in"],                                                                                              \
                     ["MatAXPY(QReal multiplier[], Mat *X, Mat *Y)",                                                          \
                      "in(2), in, in",                                                                                        \
                      "defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_CMPLX_MAT)",                                           \
                      "MatAXPY(QReal multiplier, Mat *X, Mat *Y)",                                                            \
                      "in, in, in"],                                                                                          \
                     ["MatTranspose(QMatOperation op_A, Mat *A, Mat *B)",                                                     \
                      "in, in, in"],                                                                                          \
                     ["MatGEMM(QMatOperation op_A, QMatOperation op_B, QReal alpha[], Mat *A, Mat *B, QReal beta[], Mat *C)", \
                      "in, in, in(2), in, in, in(2), in",                                                                     \
                      "defined(ADAPTER_BLOCK_CMPLX) || defined(ADAPTER_CMPLX_MAT)",                                           \
                      "MatGEMM(QMatOperation op_A, QMatOperation op_B, QReal alpha, Mat *A, Mat *B, QReal beta, Mat *C)",     \
                      "in, in, in, in, in, in, in"]                                                                           \
                    ]

# functions of the APIs that require memory allocation for the matrix
API_MALLOC_FUN = [["QcMatCreate(QcMat *A)", "out"]]
# functions of the APIs that require memory deallocation for the matrix
API_FREE_FUN = [["QcMatDestroy(QcMat *A)", "inout"]]
# other functions of the APIs that do not require any memory operation
API_NOMEM_FUN = [["QcMatBlockCreate(QcMat *A, QInt dim_block)",                                                                 \
                  "in, in"],                                                                                                  \
                 ["QcMatSetSymType(QcMat *A, QSymType sym_type)",                                                               \
                  "in, in"],                                                                                                  \
                 ["QcMatSetDataType(QcMat *A, QInt num_blocks, QInt row_ind[], QInt col_ind[], QDataType data_type[])",         \
                  "in, in, in(num_blocks), in(num_blocks), in(num_blocks)"],                                                  \
                 ["QcMatSetDimMat(QcMat *A, QInt dim_mat)",                                                                     \
                  "in, in"],                                                                                                  \
                 ["QcMatSetStorageMode(QcMat *A, QStorageMode storage_mode)",                                                   \
                  "in, in",                                                                                                   \
                  "defined(QCMATRIX_STORAGE_MODE)"],                                                                            \
                 ["QcMatAssemble(QcMat *A)",                                                                                    \
                  "in"],                                                                                                      \
                 ["QcMatGetDimBlock(QcMat *A, QInt *dim_block)",                                                                \
                  "in, out"],                                                                                                 \
                 ["QcMatGetSymType(QcMat *A, QSymType *sym_type)",                                                              \
                  "in, out"],                                                                                                 \
                 ["QcMatGetDataType(QcMat *A, QInt num_blocks, QInt row_ind[], QInt col_ind[], QDataType *data_type)",          \
                  "in, in, in(num_blocks), in(num_blocks), out(num_blocks)"],                                                 \
                 ["QcMatGetDimMat(QcMat *A, QInt *dim_mat)",                                                                    \
                  "in, out"],                                                                                                 \
                 ["QcMatGetStorageMode(QcMat *A, QStorageMode *storage_mode)",                                                  \
                  "in, out",                                                                                                  \
                  "defined(QCMATRIX_STORAGE_MODE)"],                                                                            \
                 ["QcMatIsAssembled(QcMat *A, QInt *assembled)",                                                                \
                  "in, out"],                                                                                                 \
                 ["QcMatDuplicate(QcMat *A, QDuplicateOption duplicate_option, QcMat *B)",                                       \
                  "in, in, in"],                                                                                              \
                 ["QcMatZeroEntries(QcMat *A)",                                                                                 \
                  "in"],                                                                                                      \
                 ["QcMatGetTrace(QcMat *A, QInt num_blocks, QReal *trace)",                                                     \
                  "in, in, out(2*num_blocks)"],                                                                               \
                 ["QcMatWrite(QcMat *A, QChar *mat_label, QViewOption view_option)",                                            \
                  "in, in, in",                                                                                               \
                  "defined(QCMATRIX_ENABLE_VIEW)"],                                                                             \
                 ["QcMatRead(QcMat *A, QChar *mat_label, QViewOption view_option)",                                             \
                  "in, in, in",                                                                                               \
                  "defined(QCMATRIX_ENABLE_VIEW)"],                                                                             \
                 ["QcMatScale(QReal scal_number[], QcMat *A)",                                                                  \
                  "in(2), in"],                                                                                               \
                 ["QcMatAXPY(QReal multiplier[], QcMat *X, QcMat *Y)",                                                           \
                  "in(2), in, in"],                                                                                           \
                 ["QcMatTranspose(QMatOperation op_A, QcMat *A, QcMat *B)",                                                      \
                  "in, in, in"],                                                                                              \
                 ["QcMatGEMM(QMatOperation op_A, QMatOperation op_B, QReal alpha[], QcMat *A, QcMat *B, QReal beta[], QcMat *C)", \
                  "in, in, in(2), in, in, in(2), in"],                                                                        \
                 # QcMatrix functions which are not implemented in the external libraries
                 ["QcMatMatCommutator(QcMat *A, QcMat *B, QcMat *C)",                                                             \
                  "in, in, in"],                                                                                              \
                 ["QcMatMatSCommutator(QcMat *A, QcMat *B, QcMat *S, QcMat *C)",                                                   \
                  "in, in, in, in"],                                                                                          \
                 ["QcMatMatHermCommutator(QcMat *A, QcMat *B, QcMat *C)",                                                         \
                  "in, in, in"],                                                                                              \
                 ["QcMatMatSHermCommutator(QcMat *A, QcMat *B, QcMat *S, QcMat *C)",                                               \
                  "in, in, in, in"]                                                                                           \
                ]

# extracts function names and arguments
# - c_functions: list of C functions
# - c_function_names: C function names
# - c_argument_types: argument types of a Fortran subroutine called from a C function
# - c_argument_names: argument names of a Fortran subroutine called from a C function
def extract_functions(c_functions):
    # loops over functions in \var(c_functions)
    c_function_names = []
    c_argument_types = []
    c_argument_names = []
    for ifun in xrange(len(c_functions)):
        # lists of function names, types of arguments, and names of arguments
        list_function_names = []
        list_argument_types = []
        list_argument_names = []
        for jfun in xrange(0,len(c_functions[ifun]),3):
            # gets the length of the function name
            len_fun_name = c_functions[ifun][jfun].find("(")
            # puts the arguments into a list
            fun_arguments = c_functions[ifun][jfun][len_fun_name+1:-1].split(", ")
            # removes the first three characters "(Q)Mat" in the function name
            if c_functions[ifun][jfun][:4]=="QcMat":
                list_function_names.append(c_functions[ifun][jfun][4:len_fun_name])
            elif c_functions[ifun][jfun][:3]=="Mat":
                list_function_names.append(c_functions[ifun][jfun][3:len_fun_name])
            else:
                raise ValueError(">> function "+c_functions[ifun][jfun]+" is not supported!")
            # tries to find the types and names of arguments in a function
            cfun_arg_types = []
            cfun_arg_names = []
            for iarg in xrange(len(fun_arguments)):
                # gets the length of argument type, which is separated with its name by *one* space
                len_arg_type = fun_arguments[iarg].find(" ")
                argument_type = fun_arguments[iarg][:len_arg_type]
                cfun_arg_types.append(argument_type)
                # character should be transferred by its reference
                if argument_type=="QChar" and fun_arguments[iarg][len_arg_type+1:len_arg_type+2]!="*":
                    raise ValueError(">> "+fun_arguments[iarg][len_arg_type+1:]+" is not supported!")
                # pointer of a pointer is not supported
                if fun_arguments[iarg][len_arg_type+1:len_arg_type+3]=="**":
                    raise ValueError(">> "+fun_arguments[iarg][len_arg_type+1:]+" is not supported!")
                # gets rid of "[]" for an array
                len_arg_name = fun_arguments[iarg][len_arg_type+1:].find("[")+1
                # arrays with more than one dimensional are not supported
                if fun_arguments[iarg][len_arg_type+len_arg_name+1:].find("[")!=-1:
                    raise ValueError(">> "+fun_arguments[iarg][len_arg_type+1:]+" is not supported!")
                if len_arg_name==0:
                    cfun_arg_names.append(fun_arguments[iarg][len_arg_type+1:])
                else:
                    # array of pointers is not supported
                    if fun_arguments[iarg][len_arg_type+1:len_arg_type+2]=="*":
                        raise ValueError(">> "+fun_arguments[iarg][len_arg_type+1:]+" is not supported!")
                    # an array of matrices is not supported
                    if argument_type=="Mat" or argument_type=="QcMat":
                        raise ValueError(">> "+fun_arguments[iarg][len_arg_type+1:]+" is not supported!")
                    cfun_arg_names.append(fun_arguments[iarg][len_arg_type+1:len_arg_type+len_arg_name])
            list_argument_types.append(cfun_arg_types)
            list_argument_names.append(cfun_arg_names)
        # appends the function names, types and names of the arguments
        c_function_names.append(list_function_names)
        c_argument_types.append(list_argument_types)
        c_argument_names.append(list_argument_names)
    return c_function_names,c_argument_types,c_argument_names
