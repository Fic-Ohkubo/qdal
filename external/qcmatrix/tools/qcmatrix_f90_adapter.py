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
#  This file generates the files in Fortran 90 adapter:
#  (1) cmake/f90_adapter.cmake
#  (2) include/private/f90_adapter_c.h
#  (3) src/adapter/f90_adapter_c.c
#  (4) src/adapter/f90_adapter_f.F90
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

# an array used to save the information of the matrix type, and will be converted to
# a Fortran derived type (with a single pointer member pointing to the matrix type)
# during run time
C_F_INT_MATRIX = "f90_imat"
# pointer to the Fortran 90 matrix type
FORTRAN_MAT_PTR = "f90_mat"
# size in bytes of derived types with a single pointer member
SIZEOF_F_TYPE_P = "SIZEOF_F_TYPE_P"
# default size in bytes of the derived types with a single pointer member
DEFAULT_F_TSIZE = "12"

# writes the Fortran 90 subroutines in the adapter to CMake file
# - f_adapter: name of the Fortran adapter
# - fl_cmake: file object of CMake file
# - c_functions: list of C functions
# - c_function_names: names of C functions
def f90_adapter_cmake(f_adapter, fl_cmake, c_functions, c_function_names):
    # code indent in QcMatrix
    from qcmatrix_conventions import DEFAULT_CODE_INDENT
    print ">> writes the Fortran 90 subroutines in the CMake file cmake/"+f_adapter+".cmake ..."
    for ifun in xrange(len(c_functions)):
        len_fun_list = len(c_functions[ifun])
        if len_fun_list>2:
            code_indent = DEFAULT_CODE_INDENT
        else:
            code_indent = ""
        for jfun in xrange(0,len_fun_list,3):
            # we need to process the C preprocessor
            if len_fun_list>2:
                # the first if-else directive
                if jfun==0:
                    fl_cmake.write("IF("+c_functions[ifun][2].replace("defined(", \
                                   "").replace(")","").replace("||","OR")+")\n")
                # there are even numbers of if-else directives
                elif len_fun_list%3==0:
                    fl_cmake.write("ELSEIF("+c_functions[ifun][jfun+2].replace("defined(", \
                                   "").replace(")","").replace("||","OR")+")\n")
                # there are odd numbers of if-else directive, the last one is "#else"
                else:
                    if jfun<len_fun_list-2:
                        fl_cmake.write("ELSEIF("+c_functions[ifun][jfun+2].replace("defined(", \
                                       "").replace(")","").replace("||","OR")+")\n")
                    else:
                        fl_cmake.write("ELSE()\n")
            kfun = jfun/3
            fl_cmake.write(code_indent+"SET(FC_MANGLING_SUB ${FC_MANGLING_SUB} "+f_adapter+ \
                           ":"+"Mat_Ptr_"+c_function_names[ifun][kfun]+")\n")
        if len_fun_list>2:
            fl_cmake.write("ENDIF()\n")
    return

# writes the declaration of Fortran 90 subroutines in src/adapter/f90_adapter_c.c
# - f_adapter: name of the Fortran adapter
# - fl_adapter_c: file object of src/adapter/f90_adapter_c.c
# - c_functions: list of C functions
# - c_function_names: names of C functions
# - c_argument_types: argument types in the C functions
# - c_argument_names: argument names in the C functions
def f90_adapter_c_declaration(f_adapter, fl_adapter_c, c_functions, c_function_names, \
                              c_argument_types, c_argument_names):
    # types of functions in src/adapter/f90_adapter_c.c
    c_adapter_type = {"Mat": "QInt",              \
                      "QDuplicateOption": "QInt", \
                      "QViewOption": "QInt",      \
                      "QSymType": "QInt",         \
                      "QDataType": "QInt",        \
                      "QStorageMode": "QInt",     \
                      "QChar": "QChar",           \
                      "QInt": "QInt",             \
                      "QReal": "QReal",           \
                      "QMatOperation": "QInt"}
    print ">> writes the declaration of Fortran 90 subroutines in src/adapter/"+f_adapter+"_c.c ..."
    for ifun in xrange(len(c_functions)):
        len_fun_list = len(c_functions[ifun])
        for jfun in xrange(0,len_fun_list,3):
            # we need to write the C preprocessor
            if len_fun_list>2:
                # the first if-else directive
                if jfun==0:
                    fl_adapter_c.write("#if "+c_functions[ifun][2]+"\n")
                # there are even numbers of if-else directives
                elif len_fun_list%3==0:
                    fl_adapter_c.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                # there are odd numbers of if-else directive, the last one is "#else"
                else:
                    if jfun<len_fun_list-2:
                        fl_adapter_c.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                    else:
                        fl_adapter_c.write("#else\n")
            # generates a string containing the arguments
            str_fun_args = ""
            kfun = jfun/3
            for iarg in xrange(len(c_argument_types[ifun][kfun])):
                # gets the start of the argument name
                if c_argument_names[ifun][kfun][iarg][0]=="*":
                    start_arg_name = 1
                else:
                    start_arg_name = 0
                # converts matrix to integer
                if c_argument_types[ifun][kfun][iarg]=="Mat":
                    str_fun_args = str_fun_args                                       \
                                 + c_adapter_type[c_argument_types[ifun][kfun][iarg]] \
                                 + " *i"+c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                # characters need to be transferred with the length
                elif c_argument_types[ifun][kfun][iarg]=="QChar":
                    str_fun_args = str_fun_args+"QInt *len_"                                \
                                 + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", " \
                                 + c_adapter_type[c_argument_types[ifun][kfun][iarg]]       \
                                 + " *"+c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                else:
                    str_fun_args = str_fun_args                                       \
                                 + c_adapter_type[c_argument_types[ifun][kfun][iarg]] \
                                 + " *"+c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
            fl_adapter_c.write("extern QVoid "+f_adapter+"_Mat_Ptr_"+ \
                               c_function_names[ifun][kfun]+"("+str_fun_args[:-2]+");\n")
        if len_fun_list>2:
            fl_adapter_c.write("#endif\n")
    return

# writes the functions in src/adapter/f90_adapter_c.c
# - f_adapter: name of the Fortran adapter
# - fl_adapter_c: file object of src/adapter/f90_adapter_c.c
# - c_functions: list of C functions
# - c_function_names: names of C functions
# - c_argument_types: argument types in the C functions
# - c_argument_names: argument names in the C functions
def f90_adapter_c_functions(f_adapter, fl_adapter_c, c_functions, c_function_names, \
                            c_argument_types, c_argument_names):
    # code indent in QcMatrix
    from qcmatrix_conventions import DEFAULT_CODE_INDENT
    print ">> writes the functions in src/adapter/"+f_adapter+"_c.c ..."
    for ifun in xrange(len(c_functions)):
        len_fun_list = len(c_functions[ifun])
        for jfun in xrange(0,len_fun_list,3):
            print ">> writes function "+c_functions[ifun][jfun]+" ..."
            # we need to write the C preprocessor
            if len_fun_list>2:
                # the first if-else directive
                if jfun==0:
                    fl_adapter_c.write("#if "+c_functions[ifun][2]+"\n")
                # there are even numbers of if-else directives
                elif len_fun_list%3==0:
                    fl_adapter_c.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                # there are odd numbers of if-else directive, the last one is "#else"
                else:
                    if jfun<len_fun_list-2:
                        fl_adapter_c.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                    else:
                        fl_adapter_c.write("#else\n")
            fl_adapter_c.write("QErrorCode "+c_functions[ifun][jfun]+"\n")
            fl_adapter_c.write("{\n")
            # a string containing the arguments
            str_fun_args = ""
            # a string containing the character arguments
            str_char_args = ""
            kfun = jfun/3
            for iarg in xrange(len(c_argument_types[ifun][kfun])):
                # gets the start of the argument name
                if c_argument_names[ifun][kfun][iarg][0]=="*":
                    start_arg_name = 1
                    oper_struct_access = "->"
                    oper_reference = ""
                else:
                    start_arg_name = 0
                    oper_struct_access = "."
                    oper_reference = "&"
                # converts matrix to integer
                if c_argument_types[ifun][kfun][iarg]=="Mat":
                    str_fun_args = str_fun_args+"&"                                    \
                                 + c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                                 + oper_struct_access+C_F_INT_MATRIX+"[0], "
                # characters need to be transferred with the length
                elif c_argument_types[ifun][kfun][iarg]=="QChar":
                    str_fun_args = str_fun_args+"&len_"                        \
                                 + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", " \
                                 + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                    str_char_args = str_char_args+DEFAULT_CODE_INDENT+"QInt len_"       \
                                  + c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                                  + " = strlen("                                        \
                                  + c_argument_names[ifun][kfun][iarg][start_arg_name:]+");\n"
                else:
                    str_fun_args = str_fun_args+oper_reference \
                                 + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
            if len(str_char_args)!=0:
                fl_adapter_c.write(str_char_args)
            fl_adapter_c.write(DEFAULT_CODE_INDENT+f_adapter+"_Mat_Ptr_"+ \
                               c_function_names[ifun][kfun]+"("+str_fun_args[:-2]+");\n")
            fl_adapter_c.write(DEFAULT_CODE_INDENT+"return QSUCCESS;\n")
            fl_adapter_c.write("}\n")
        if len_fun_list>2:
            fl_adapter_c.write("#endif\n")
        fl_adapter_c.write("\n")
    return

# writes the subroutines in src/adapter/f90_adapter_f.F90
# - f_adapter: name of the Fortran adapter
# - fl_adapter_f: file object of src/adapter/f90_adapter_f.F90
# - c_functions: list of C functions
# - c_function_names: names of C functions
# - c_argument_types: argument types in the C functions
# - c_argument_names: argument names in the C functions
def f90_adapter_f_subroutines(f_adapter, fl_adapter_f, c_functions, c_function_names, \
                              c_argument_types, c_argument_names):
    # name of external Fortran 90 module and code indent in QcMatrix
    from qcmatrix_conventions import LANG_F_MODULE,DEFAULT_CODE_INDENT
    print ">> writes the subroutines in src/adapter/"+f_adapter+"_f.F90 ..."
    # conversion list for types between C and Fortran
    fc_conversion = {"Mat": "integer",              \
                     "QDuplicateOption": "integer", \
                     "QViewOption": "integer",      \
                     "QSymType": "integer",         \
                     "QDataType": "integer",        \
                     "QStorageMode": "integer",     \
                     "QChar": "character*(*)",      \
                     "QInt": "integer",             \
                     "QReal": "real(QREAL)",        \
                     "QMatOperation": "integer"}
    sub_code_indent = 2*DEFAULT_CODE_INDENT
    for ifun in xrange(len(c_functions)):
        len_fun_list = len(c_functions[ifun])
        for jfun in xrange(0,len_fun_list,3):
            # we need to write the C preprocessor
            if len_fun_list>2:
                # the first if-else directive
                if jfun==0:
                    fl_adapter_f.write("#if "+c_functions[ifun][2]+"\n")
                # there are even numbers of if-else directives
                elif len_fun_list%3==0:
                    fl_adapter_f.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                # there are odd numbers of if-else directive, the last one is "#else"
                else:
                    if jfun<len_fun_list-2:
                        fl_adapter_f.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                    else:
                        fl_adapter_f.write("#else\n")
            # generates a string containing the arguments
            str_adapter_args = ""
            str_eli_f_args = ""
            kfun = jfun/3
            print ">> writes subroutine Mat_Ptr_"+c_function_names[ifun][kfun]+" ..."
            # intent of arguments
            arg_intent_size = c_functions[ifun][jfun+1].split(", ")
            str_arg_intent = ""
            # list of matrix arguments
            mat_argument = []
            mat_alloc = []
            mat_dealloc = []
            # characters arguments
            char_argument = []
            for iarg in xrange(len(c_argument_types[ifun][kfun])):
                len_arg_intent = arg_intent_size[iarg].find("(")
                if len_arg_intent==-1:
                    argument_intent = arg_intent_size[iarg]
                    argument_size = ""
                else:
                    argument_intent = arg_intent_size[iarg][:len_arg_intent]
                    argument_size = arg_intent_size[iarg][len_arg_intent+1:-1]
                # gets the start of the argument name
                if c_argument_names[ifun][kfun][iarg][0]=="*":
                    start_arg_name = 1
                else:
                    start_arg_name = 0
                # converts matrix to integer
                if c_argument_types[ifun][kfun][iarg]=="Mat":
                    if argument_size!="":
                        argument_size = "*"+argument_size
                    str_adapter_args = str_adapter_args+"i" \
                                     + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                    str_eli_f_args = str_eli_f_args                                      \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                                   + "%"+FORTRAN_MAT_PTR+", "
                    str_arg_intent = str_arg_intent+sub_code_indent                      \
                                   + fc_conversion[c_argument_types[ifun][kfun][iarg]]   \
                                   + ", intent("+argument_intent+") :: i"                \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                                   + "("+SIZEOF_F_TYPE_P+argument_size+")\n"
                    mat_argument.append(c_argument_names[ifun][kfun][iarg][start_arg_name:])
                    mat_alloc.append(argument_intent=="out")
                    mat_dealloc.append(argument_intent=="inout")
                # characters need to be transferred with the length
                elif c_argument_types[ifun][kfun][iarg]=="QChar":
                    if argument_size!="":
                        argument_size = "("+argument_size+")"
                    str_adapter_args = str_adapter_args+"len_"                                  \
                                     + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", " \
                                     + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                    str_eli_f_args = str_eli_f_args \
                                   + "f_"+c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                    str_arg_intent = str_arg_intent+sub_code_indent                                \
                                   + fc_conversion["QInt"]+", intent("+argument_intent+") :: len_" \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:]+"\n"      \
                                   + sub_code_indent                                               \
                                   + fc_conversion[c_argument_types[ifun][kfun][iarg]]             \
                                   + ", intent("+argument_intent+") :: "                           \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:]           \
                                   + argument_size+"\n"
                    char_argument.append(c_argument_names[ifun][kfun][iarg][start_arg_name:])
                else:
                    if argument_size!="":
                        argument_size = "("+argument_size+")"
                    str_adapter_args = str_adapter_args \
                                     + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                    str_eli_f_args = str_eli_f_args \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                    str_arg_intent = str_arg_intent+sub_code_indent                      \
                                   + fc_conversion[c_argument_types[ifun][kfun][iarg]]   \
                                   + ", intent("+argument_intent+") :: "                 \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                                   + argument_size+"\n"
            fl_adapter_f.write(DEFAULT_CODE_INDENT+"subroutine Mat_Ptr_"+ \
                               c_function_names[ifun][kfun]+"("+str_adapter_args[:-2]+")\n")
            fl_adapter_f.write(sub_code_indent+"use "+LANG_F_MODULE+", only : Matrix_"+ \
                               c_function_names[ifun][kfun]+"\n")
            fl_adapter_f.write(sub_code_indent+"implicit none\n")
            fl_adapter_f.write(str_arg_intent)
            # writes local variables for matrices, characters and error information
            for imat in xrange(len(mat_argument)):
                fl_adapter_f.write(sub_code_indent+"type(matrix_ptr_t) "+mat_argument[imat]+"\n")
            for ichar in xrange(len(char_argument)):
                fl_adapter_f.write(sub_code_indent+"character(len=len_"+char_argument[ichar]+ \
                                   ") f_"+char_argument[ichar]+"\n")
            if True in mat_alloc:
                fl_adapter_f.write(sub_code_indent+"integer ierr\n")
            # writes the allocation of memory for matrices
            for imat in xrange(len(mat_argument)):
                if mat_alloc[imat]:
                    fl_adapter_f.write(sub_code_indent+"allocate("+mat_argument[imat]+ \
                                       "%"+FORTRAN_MAT_PTR+")\n")
                    #fl_adapter_f.write(sub_code_indent+"if (ierr/=0) then\n")
                    #fl_adapter_f.write(3*DEFAULT_CODE_INDENT+"write(STDOUT,*) \"Mat_Ptr_"+       \
                    #                   c_function_names[ifun][kfun]+">> error when allocating "+ \
                    #                   mat_argument[imat]+"%"+FORTRAN_MAT_PTR+"\", ierr\n")
                    #fl_adapter_f.write(3*DEFAULT_CODE_INDENT+"stop\n")
                    #fl_adapter_f.write(sub_code_indent+"end if\n")
                else:
                    fl_adapter_f.write(sub_code_indent+mat_argument[imat]+" = transfer(i"+ \
                                       mat_argument[imat]+", "+mat_argument[imat]+")\n")
            # saves the C characters into local Fortran characters
            for ichar in xrange(len(char_argument)):
                fl_adapter_f.write(sub_code_indent+"f_"+char_argument[ichar]+" = "+ \
                                   char_argument[ichar]+"(1:len_"+char_argument[ichar]+")\n")
            # calls the subroutine in the external Fortran 90 library
            fl_adapter_f.write(sub_code_indent+"call Matrix_"+c_function_names[ifun][kfun]+ \
                               "("+str_eli_f_args[:-2]+")\n")
            # converts Fortran 90 type matrix to an array of integers, or deallocates memory
            for imat in xrange(len(mat_argument)):
                if mat_alloc[imat]:
                    fl_adapter_f.write(sub_code_indent+"i"+mat_argument[imat]+" = transfer("+ \
                                       mat_argument[imat]+", i"+mat_argument[imat]+")\n")
                elif mat_dealloc[imat]:
                    fl_adapter_f.write(sub_code_indent+"deallocate("+mat_argument[imat]+ \
                                       "%"+FORTRAN_MAT_PTR+")\n")
                    fl_adapter_f.write(sub_code_indent+"nullify("+mat_argument[imat]+ \
                                       "%"+FORTRAN_MAT_PTR+")\n")
                    fl_adapter_f.write(sub_code_indent+"i"+mat_argument[imat]+" = 0\n")
            fl_adapter_f.write(sub_code_indent+"return\n")
            fl_adapter_f.write(DEFAULT_CODE_INDENT+"end subroutine Mat_Ptr_"+ \
                               c_function_names[ifun][kfun]+"\n")
        if len_fun_list>2:
            fl_adapter_f.write("#endif\n")
        fl_adapter_f.write("\n")
    return

# function to generate the files in Fortran 90 adapter:
#  (1) cmake/f90_adapter.cmake
#  (2) include/private/f90_adapter_c.h
#  (3) src/adapter/f90_adapter_c.c
#  (4) src/adapter/f90_adapter_f.F90
#
# we use the following way to name functions:
#
# src/adapter/f90_adapter_c.c   src/adapter/f90_adapter_f.F90   external Fortran library
#     MatSetSymType      -->      Mat_Ptr_SetSymType     -->    Matrix_SetSymType
def generate_f90_adapter(f_adapter):
    # prints creating date
    from datetime import date
    # conventions in QcMatrix
    from qcmatrix_conventions import LANG_F_MODULE,DEFAULT_F_MOD, \
                                    LANG_F_MATRIX,DEFAULT_F_MAT, \
                                    DEFAULT_CODE_INDENT
    # QcMatrix license
    from qcmatrix_license import write_license
    # QcMatrix functions
    from qcmatrix_functions import ADAPTER_MALLOC_FUN, \
                                  ADAPTER_FREE_FUN,   \
                                  ADAPTER_NOMEM_FUN,  \
                                  extract_functions
    print ">> extracts function names and arguments ..."
    malloc_fun_names,malloc_arg_types,malloc_arg_names = extract_functions(ADAPTER_MALLOC_FUN)
    mfree_fun_names,mfree_arg_types,mfree_arg_names = extract_functions(ADAPTER_FREE_FUN)
    nomem_fun_names,nomem_arg_types,nomem_arg_names = extract_functions(ADAPTER_NOMEM_FUN)
    print ">> number of functions", \
          len(ADAPTER_MALLOC_FUN)+len(ADAPTER_FREE_FUN)+len(ADAPTER_NOMEM_FUN)
    # creates file objects
    print ">> open files in the Fortran 90 adapter ..."
    fl_cmake = open(f_adapter+".cmake", "w")
    fl_adapter_h = open(f_adapter+"_c.h", "w")
    fl_adapter_c = open(f_adapter+"_c.c", "w")
    fl_adapter_f = open(f_adapter+"_f.F90", "w")
    # writes license information
    file_description = "   This file implements the functions in private/mat_adapter.h by calling\n" \
                     + "   src/adapter/"+f_adapter+"_f.F90.\n"
    write_license(fl_adapter_c, "C", file_description)
    file_description = "!!  This file provides an adapter between QcMatrix and external Fortran 90\n" \
                     + "!!  matrix library.\n"
    write_license(fl_adapter_f, "Fortran", file_description)
    # CMake generated header file with auto-detected mangling 
    mangling_header = "f90_mangling.h"
    fl_cmake.write("# CMake file for Fortran 90 adapter\n")
    fl_cmake.write("# "+date.today().isoformat()+", Bin Gao:\n")
    fl_cmake.write("# * generated by tools/qcmatrix_fortran.py\n")
    fl_cmake.write("# Sets Fortran compiler flags\n")
    fl_cmake.write("SET(LANG_F_MODULE \""+DEFAULT_F_MOD+"\" CACHE STRING \""+ \
                   "Name of external Fortran 90 matrix module\")\n")
    fl_cmake.write("SET(LANG_F_MATRIX \""+DEFAULT_F_MAT+"\" CACHE STRING \""+ \
                   "Name of external Fortran 90 matrix type\")\n")
    fl_cmake.write("SET(SIZEOF_F_TYPE_P \""+DEFAULT_F_TSIZE+"\" CACHE STRING \"Size (in bytes) of "+ \
                   "Fortran 90 derived types with a single pointer member\")\n")
    fl_cmake.write("SET(ADAPTER_F_FLAGS \"-D"+LANG_F_MODULE+"=${LANG_F_MODULE}"+ \
                   " -D"+LANG_F_MATRIX+"=${LANG_F_MATRIX}"+                      \
                   " -D"+SIZEOF_F_TYPE_P+"=${SIZEOF_F_TYPE_P}\")\n")
    fl_cmake.write("SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS} ${ADAPTER_F_FLAGS}\")\n")
    fl_cmake.write("# Source codes of the adapter for the external Fortran 90 library\n")
    fl_cmake.write("SET(ADAPTER_SRCS\n")
    fl_cmake.write("    ${LIB_QCMATRIX_PATH}/src/adapter/"+f_adapter+"_f.F90\n")
    fl_cmake.write("    ${LIB_QCMATRIX_PATH}/src/adapter/"+f_adapter+"_c.c)\n")
    fl_cmake.write("# Fortran 90 subroutines of the adapter\n")
    # writes include/private/f90_adapter_c.h
    fl_adapter_h.write("/* "+SIZEOF_F_TYPE_P+ \
                       " is the size in bytes of derived types with a single pointer member,\n")
    fl_adapter_h.write("   default is "+DEFAULT_F_TSIZE+ \
                       " and could be changed by -D"+SIZEOF_F_TYPE_P+"=XX */\n")
    fl_adapter_h.write("#if !defined("+SIZEOF_F_TYPE_P+")\n")
    fl_adapter_h.write("#define "+SIZEOF_F_TYPE_P+" "+DEFAULT_F_TSIZE+"\n")
    fl_adapter_h.write("#endif\n")
    fl_adapter_h.write("/* Fortran 90 library usually defines a matrix using derived type in a module (default\n")
    fl_adapter_h.write("   is "+DEFAULT_F_MOD+", could be set by -D"+LANG_F_MODULE+ \
                       "=XX), which can not be used directly.\n")
    fl_adapter_h.write("   Instead an array "+C_F_INT_MATRIX+ \
                       " is used to save the information of the matrix type, which\n")
    fl_adapter_h.write("   will be converted to a Fortran derived type (with a single pointer member pointing\n")
    fl_adapter_h.write("   to the matrix type) during run time by subroutines in src/adapter/"+ \
                       f_adapter+"_f.F90 */\n")
    fl_adapter_h.write("typedef struct {\n")
    fl_adapter_h.write(DEFAULT_CODE_INDENT+"QInt "+C_F_INT_MATRIX+"["+SIZEOF_F_TYPE_P+"];\n")
    fl_adapter_h.write("} Mat;\n")
    fl_adapter_h.close()
    # includes header files in src/adapter/f90_adapter_c.c
    fl_adapter_c.write("\n#include \"private/mat_adapter.h\"\n")
    fl_adapter_c.write("/* uses CMake generated header file with auto-detected mangling */\n")
    fl_adapter_c.write("#include \""+mangling_header+"\"\n\n")
    fl_adapter_c.write("/* declaration of Fortran 90 subroutines */\n")
    # writes the module, derivated type and parameters in src/adapter/f90_adapter_f.F90
    fl_adapter_f.write("\n! macros related to Fortran module and compilers\n")
    fl_adapter_f.write("#if !defined("+LANG_F_MODULE+")\n")
    fl_adapter_f.write("#define "+LANG_F_MODULE+" "+DEFAULT_F_MOD+"\n")
    fl_adapter_f.write("#endif\n")
    fl_adapter_f.write("#if !defined("+LANG_F_MATRIX+")\n")
    fl_adapter_f.write("#define "+LANG_F_MATRIX+" "+DEFAULT_F_MAT+"\n")
    fl_adapter_f.write("#endif\n")
    fl_adapter_f.write("#if !defined("+SIZEOF_F_TYPE_P+")\n")
    fl_adapter_f.write("#define "+SIZEOF_F_TYPE_P+" "+DEFAULT_F_TSIZE+"\n")
    fl_adapter_f.write("#endif\n\n")
    fl_adapter_f.write("module "+f_adapter+"\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"! includes matrix type and relevant subroutines from "+ \
                       "external Fortran 90 module\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"use "+LANG_F_MODULE+", only : "+LANG_F_MATRIX+"\n\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"implicit none\n\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+ \
                       "! derived type in which the single pointer member points to the matrix in\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"! external Fortran 90 module; the idea is from "+ \
                       "Comput. Sci. Eng. 10, 86 (2008)\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"type matrix_ptr_t\n")
    fl_adapter_f.write(2*DEFAULT_CODE_INDENT+"private\n")
    fl_adapter_f.write(2*DEFAULT_CODE_INDENT+"type("+LANG_F_MATRIX+"), pointer :: "+ \
                       FORTRAN_MAT_PTR+"\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"end type matrix_ptr_t\n\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"! data type for real numbers\n")
    fl_adapter_f.write("#if defined(ENABLE_SINGLE_PRECISION)\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"integer, parameter :: QREAL = kind(1.0)\n")
    fl_adapter_f.write("#else\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"integer, parameter :: QREAL = kind(1.0D0)\n")
    fl_adapter_f.write("#endif\n\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"! unit number of standard output\n")
    #fl_adapter_f.write(DEFAULT_CODE_INDENT+"integer, parameter :: STDOUT = 6\n\n")
    fl_adapter_f.write(DEFAULT_CODE_INDENT+"contains\n\n")
    # writes the Fortran 90 subroutines of the adapter to CMake file cmake/f90_adapter.cmake
    f90_adapter_cmake(f_adapter, fl_cmake, ADAPTER_MALLOC_FUN, malloc_fun_names)
    f90_adapter_cmake(f_adapter, fl_cmake, ADAPTER_FREE_FUN, mfree_fun_names)
    f90_adapter_cmake(f_adapter, fl_cmake, ADAPTER_NOMEM_FUN, nomem_fun_names)
    # writes the declaration of Fortran 90 subroutines in src/adapter/f90_adapter_c.c
    f90_adapter_c_declaration(f_adapter, fl_adapter_c, ADAPTER_MALLOC_FUN, malloc_fun_names, \
                              malloc_arg_types, malloc_arg_names)
    f90_adapter_c_declaration(f_adapter, fl_adapter_c, ADAPTER_FREE_FUN, mfree_fun_names, \
                              mfree_arg_types, mfree_arg_names)
    f90_adapter_c_declaration(f_adapter, fl_adapter_c, ADAPTER_NOMEM_FUN, nomem_fun_names, \
                              nomem_arg_types, nomem_arg_names)
    fl_adapter_c.write("\n")
    # writes the functions in src/adapter/f90_adapter_c.c
    f90_adapter_c_functions(f_adapter, fl_adapter_c, ADAPTER_MALLOC_FUN, malloc_fun_names, \
                            malloc_arg_types, malloc_arg_names)
    f90_adapter_c_functions(f_adapter, fl_adapter_c, ADAPTER_FREE_FUN, mfree_fun_names, \
                            mfree_arg_types, mfree_arg_names)
    f90_adapter_c_functions(f_adapter, fl_adapter_c, ADAPTER_NOMEM_FUN, nomem_fun_names, \
                            nomem_arg_types, nomem_arg_names)
    # writes the subroutines in src/adapter/f90_adapter_f.F90
    f90_adapter_f_subroutines(f_adapter, fl_adapter_f, ADAPTER_MALLOC_FUN, malloc_fun_names, \
                              malloc_arg_types, malloc_arg_names)
    f90_adapter_f_subroutines(f_adapter, fl_adapter_f, ADAPTER_FREE_FUN, mfree_fun_names, \
                              mfree_arg_types, mfree_arg_names)
    f90_adapter_f_subroutines(f_adapter, fl_adapter_f, ADAPTER_NOMEM_FUN, nomem_fun_names, \
                              nomem_arg_types, nomem_arg_names)
    # writes left parts of the adapter
    fl_adapter_f.write("end module "+f_adapter+"\n")
    print ">> closes files ..."
    fl_cmake.close()
    fl_adapter_c.close()
    fl_adapter_f.close()
    print
    print ">> please replace the files of the Fortran 90 adapter as:"
    print "   mv "+f_adapter+".cmake ${PATH_QCMATRIX}/cmake/"+f_adapter+".cmake"
    print "   mv "+f_adapter+"_c.h ${PATH_QCMATRIX}/include/private/"+f_adapter+"_c.h"
    print "   mv "+f_adapter+"_c.c ${PATH_QCMATRIX}/src/adapter/"+f_adapter+"_c.c"
    print "   mv "+f_adapter+"_f.F90 ${PATH_QCMATRIX}/src/adapter/"+f_adapter+"_f.F90"
    print
    return
