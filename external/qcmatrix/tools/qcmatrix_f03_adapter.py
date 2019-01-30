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
#  This file generates the files in Fortran 2003 adapter:
#  (1) cmake/f03_adapter.cmake
#  (2) include/private/f03_adapter_c.h
#  (3) src/adapter/f03_adapter_c.c
#  (4) src/adapter/f03_adapter_f.F90
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

# a pointer used to save of the information of the matrix type, and will be converted to
# a pointer of matrix type defined in Fortran 2003 library during run time
C_F_PTR_MATRIX = "f03_mat"

# writes the declaration of Fortran 2003 subroutines in src/adapter/f03_adapter_c.c
# - f_adapter: name of the Fortran adapter
# - fl_adapter_c: file object of src/adapter/f03_adapter_c.c
# - c_functions: list of C functions
# - c_function_names: names of C functions
# - c_argument_types: argument types in the C functions
# - c_argument_names: argument names in the C functions
def f03_adapter_c_declaration(f_adapter, fl_adapter_c, c_functions, c_function_names, \
                              c_argument_types, c_argument_names):
    # types of functions in src/adapter/f03_adapter_c.c
    c_adapter_type = {"Mat": "QVoid",             \
                      "QDuplicateOption": "QInt", \
                      "QViewOption": "QInt",      \
                      "QSymType": "QInt",         \
                      "QDataType": "QInt",        \
                      "QStorageMode": "QInt",     \
                      "QChar": "QChar",           \
                      "QInt": "QInt",             \
                      "QReal": "QReal",           \
                      "QMatOperation": "QInt"}
    print ">> writes the declaration of Fortran 2003 subroutines in src/adapter/"+f_adapter+"_c.c ..."
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
                # converts matrix to C_PTR
                if c_argument_types[ifun][kfun][iarg]=="Mat":
                    str_fun_args = str_fun_args                                       \
                                 + c_adapter_type[c_argument_types[ifun][kfun][iarg]] \
                                 + " **"+c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                # characters need to be transferred with the length
                elif c_argument_types[ifun][kfun][iarg]=="QChar":
                    str_fun_args = str_fun_args                                             \
                                 + "QInt *len_"                                             \
                                 + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", " \
                                 + c_adapter_type[c_argument_types[ifun][kfun][iarg]]       \
                                 + " *"+c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                                 + ", "
                else:
                    str_fun_args = str_fun_args                                       \
                                 + c_adapter_type[c_argument_types[ifun][kfun][iarg]] \
                                 + " *"+c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
            fl_adapter_c.write("extern QVoid Mat_Ptr_"+ \
                               c_function_names[ifun][kfun]+"("+str_fun_args[:-2]+");\n")
        if len_fun_list>2:
            fl_adapter_c.write("#endif\n")
    return

# writes the functions in src/adapter/f03_adapter_c.c
# - f_adapter: name of the Fortran adapter
# - fl_adapter_c: file object of src/adapter/f03_adapter_c.c
# - c_functions: list of C functions
# - c_function_names: names of C functions
# - c_argument_types: argument types in the C functions
# - c_argument_names: argument names in the C functions
def f03_adapter_c_functions(f_adapter, fl_adapter_c, c_functions, c_function_names, \
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
                                 + oper_struct_access+C_F_PTR_MATRIX+", "
                # characters need to be transferred with the length
                elif c_argument_types[ifun][kfun][iarg]=="QChar":
                    str_fun_args = str_fun_args+"&len_"                                     \
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
            fl_adapter_c.write(DEFAULT_CODE_INDENT+"Mat_Ptr_"+c_function_names[ifun][kfun]+ \
                               "("+str_fun_args[:-2]+");\n")
            fl_adapter_c.write(DEFAULT_CODE_INDENT+"return QSUCCESS;\n")
            fl_adapter_c.write("}\n")
        if len_fun_list>2:
            fl_adapter_c.write("#endif\n")
        fl_adapter_c.write("\n")
    return

# writes the subroutines in src/adapter/f03_adapter_f.F90
# - f_adapter: name of the Fortran adapter
# - fl_adapter_f: file object of src/adapter/f03_adapter_f.F90
# - c_functions: list of C functions
# - c_function_names: names of C functions
# - c_argument_types: argument types in the C functions
# - c_argument_names: argument names in the C functions
def f03_adapter_f_subroutines(f_adapter, fl_adapter_f, c_functions, c_function_names, \
                              c_argument_types, c_argument_names):
    # name of external Fortran 2003 module and code indent in QcMatrix
    from qcmatrix_conventions import LANG_F_MODULE,LANG_F_MATRIX,DEFAULT_CODE_INDENT
    print ">> writes the subroutines in src/adapter/"+f_adapter+"_f.F90 ..."
    # conversion list for types between C and Fortran
    fc_conversion = {"Mat": "type(C_PTR)",                 \
                     "QDuplicateOption": "integer(C_INT)", \
                     "QViewOption": "integer(C_INT)",      \
                     "QSymType": "integer(C_INT)",         \
                     "QDataType": "integer(C_INT)",        \
                     "QStorageMode": "integer(C_INT)",     \
                     "QChar": "type(C_PTR)",               \
                     "QInt": "integer(C_INT)",             \
                     "QReal": "real(C_QREAL)",             \
                     "QMatOperation": "integer(C_INT)"}
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
                # converts matrix to C_PTR
                if c_argument_types[ifun][kfun][iarg]=="Mat":
                    if argument_size!="":
                        argument_size = "*"+argument_size
                    str_adapter_args = str_adapter_args \
                                     + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                    str_eli_f_args = str_eli_f_args \
                                   + "f_"+c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                    str_arg_intent = str_arg_intent+DEFAULT_CODE_INDENT                \
                                   + fc_conversion[c_argument_types[ifun][kfun][iarg]] \
                                   + ", intent("+argument_intent+") :: "               \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:]+"\n"
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
                    str_arg_intent = str_arg_intent+DEFAULT_CODE_INDENT                            \
                                   + fc_conversion["QInt"]+", intent("+argument_intent+") :: len_" \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:]+"\n"      \
                                   + DEFAULT_CODE_INDENT                                           \
                                   + fc_conversion[c_argument_types[ifun][kfun][iarg]]             \
                                   + ", intent("+argument_intent+"), value :: "                    \
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
                    str_arg_intent = str_arg_intent+DEFAULT_CODE_INDENT                  \
                                   + fc_conversion[c_argument_types[ifun][kfun][iarg]]   \
                                   + ", intent("+argument_intent+") :: "                 \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                                   + argument_size+"\n"
            fl_adapter_f.write("subroutine Mat_Ptr_"+c_function_names[ifun][kfun]+ \
                               "("+str_adapter_args[:-2]+") bind(C, name=\"Mat_Ptr_"+ \
                               c_function_names[ifun][kfun]+"\")\n")
            fl_adapter_f.write(DEFAULT_CODE_INDENT+"use "+LANG_F_MODULE+", only : "+ \
                               LANG_F_MATRIX+",Matrix_"+c_function_names[ifun][kfun]+"\n")
            fl_adapter_f.write(DEFAULT_CODE_INDENT+"use, intrinsic :: iso_c_binding\n")
            fl_adapter_f.write(DEFAULT_CODE_INDENT+"implicit none\n")
            fl_adapter_f.write(str_arg_intent)
            # writes local variables for matrices, characters and error information
            for imat in xrange(len(mat_argument)):
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"type("+LANG_F_MATRIX+ \
                                   "), pointer :: f_"+mat_argument[imat]+"\n")
            for ichar in xrange(len(char_argument)):
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"character(kind=C_CHAR), pointer :: ptr_"+ \
                                   char_argument[ichar]+"(:)\n")
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"character(len=len_"+ \
                                   char_argument[ichar]+") f_"+char_argument[ichar]+"\n")
            if True in mat_alloc:
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"integer ierr\n")
            if len(char_argument)!=0:
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"integer ic\n")
            # writes the allocation of memory for matrices
            for imat in xrange(len(mat_argument)):
                if mat_alloc[imat]:
                    fl_adapter_f.write(DEFAULT_CODE_INDENT+"allocate(f_"+mat_argument[imat]+")\n")
                else:
                    fl_adapter_f.write(DEFAULT_CODE_INDENT+"call c_f_pointer("+ \
                                       mat_argument[imat]+", f_"+mat_argument[imat]+")\n")
            # saves the C characters into local Fortran characters
            for ichar in xrange(len(char_argument)):
                # it may be better to check the status of C_PTR characters by using
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"if (.not.c_associated("+ \
                                   char_argument[ichar]+")) then\n")
                fl_adapter_f.write(2*DEFAULT_CODE_INDENT+"write(STDOUT,*) \"Mat_Ptr_"+      \
                                   c_function_names[ifun][kfun]+">> "+char_argument[ichar]+ \
                                   " is a C NULL pointer!\"\n")
                fl_adapter_f.write(2*DEFAULT_CODE_INDENT+"stop\n")
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"end if\n")
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"call c_f_pointer("+char_argument[ichar]+ \
                                   ", ptr_"+char_argument[ichar]+", [len_"+char_argument[ichar]+"])\n")
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"do ic = 1, len_"+char_argument[ichar]+"\n")
                fl_adapter_f.write(2*DEFAULT_CODE_INDENT+"f_"+char_argument[ichar]+ \
                                   "(ic:ic) = ptr_"+char_argument[ichar]+"(ic)(1:1)\n")
                fl_adapter_f.write(DEFAULT_CODE_INDENT+"end do\n")
            # calls the subroutine in the external Fortran 2003 library
            fl_adapter_f.write(DEFAULT_CODE_INDENT+"call Matrix_"+c_function_names[ifun][kfun]+ \
                               "("+str_eli_f_args[:-2]+")\n")
            # converts Fortran 2003 type matrix to an array of integers, or deallocates memory
            for imat in xrange(len(mat_argument)):
                if mat_alloc[imat]:
                    fl_adapter_f.write(DEFAULT_CODE_INDENT+mat_argument[imat]+ \
                                       " = c_loc(f_"+mat_argument[imat]+")\n")
                elif mat_dealloc[imat]:
                    fl_adapter_f.write(DEFAULT_CODE_INDENT+"deallocate(f_"+mat_argument[imat]+")\n")
                    fl_adapter_f.write(DEFAULT_CODE_INDENT+"nullify(f_"+mat_argument[imat]+")\n")
                    fl_adapter_f.write(DEFAULT_CODE_INDENT+mat_argument[imat]+" = C_NULL_PTR\n")
            fl_adapter_f.write(DEFAULT_CODE_INDENT+"return\n")
            fl_adapter_f.write("end subroutine Mat_Ptr_"+c_function_names[ifun][kfun]+"\n")
        if len_fun_list>2:
            fl_adapter_f.write("#endif\n")
        fl_adapter_f.write("\n")
    return

# function to generate the files in Fortran 2003 adapter:
# (1) cmake/f03_adapter.cmake
# (2) include/private/f03_adapter_c.h
# (3) src/adapter/f03_adapter_c.c
# (4) src/adapter/f03_adapter_f.F90
#
# we use the following way to name functions:
#
# src/adapter/f03_adapter_c.c   src/adapter/f03_adapter_f.F90   external Fortran library
#     MatSetSymType      -->      Mat_Ptr_SetSymType     -->    Matrix_SetSymType
def generate_f03_adapter(f_adapter):
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
    print ">> open files in the Fortran 2003 adapter ..."
    fl_cmake = open(f_adapter+".cmake", "w")
    fl_adapter_h = open(f_adapter+"_c.h", "w")
    fl_adapter_c = open(f_adapter+"_c.c", "w")
    fl_adapter_f = open(f_adapter+"_f.F90", "w")
    # writes license information
    file_description = "   This file implements the functions in private/mat_adapter.h by calling\n" \
                     + "   src/adapter/"+f_adapter+"_f.F90.\n"
    write_license(fl_adapter_c, "C", file_description)
    file_description = "!!  This file provides an adapter between QcMatrix and external Fortran 2003\n" \
                     + "!!  matrix library.\n"
    write_license(fl_adapter_f, "Fortran", file_description)
    # writes the CMake file cmake/f03_adapter.cmake
    fl_cmake.write("# CMake file for Fortran 2003 adapter\n")
    fl_cmake.write("# "+date.today().isoformat()+", Bin Gao:\n")
    fl_cmake.write("# * generated by tools/qcmatrix_fortran.py\n")
    fl_cmake.write("# Sets Fortran compiler flags\n")
    fl_cmake.write("SET(LANG_F_MODULE \""+DEFAULT_F_MOD+"\" CACHE STRING \""+ \
                   "Name of external Fortran 2003 matrix module\")\n")
    fl_cmake.write("SET(LANG_F_MATRIX \""+DEFAULT_F_MAT+"\" CACHE STRING \""+ \
                   "Name of external Fortran 2003 matrix type\")\n")
    fl_cmake.write("SET(ADAPTER_F_FLAGS \"-D"+LANG_F_MODULE+"=${LANG_F_MODULE}"+ \
                   " -D"+LANG_F_MATRIX+"=${LANG_F_MATRIX}\")\n")
    fl_cmake.write("SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS} ${ADAPTER_F_FLAGS}\")\n")
    fl_cmake.write("# Source codes of the adapter for the external Fortran 2003 library\n")
    fl_cmake.write("SET(ADAPTER_SRCS\n")
    fl_cmake.write("    ${LIB_QCMATRIX_PATH}/src/adapter/"+f_adapter+"_f.F90\n")
    fl_cmake.write("    ${LIB_QCMATRIX_PATH}/src/adapter/"+f_adapter+"_c.c)\n")
    fl_cmake.close()
    # writes include/private/f03_adapter_c.h
    fl_adapter_h.write("/* The ISO_C_BINDING in Fortran 2003 is used to take care data type conversion and name mangling,\n")
    fl_adapter_h.write("   the pointer "+C_F_PTR_MATRIX+ \
                       " will be converted to a pointer of matrix type in Fortran 2003 using the\n")
    fl_adapter_h.write("   intrinsic procedure C_F_POINTER during run time */\n")
    fl_adapter_h.write("typedef struct {\n")
    fl_adapter_h.write(DEFAULT_CODE_INDENT+"QVoid *"+C_F_PTR_MATRIX+";\n")
    fl_adapter_h.write("} Mat;\n")
    fl_adapter_h.close()
    # includes header files in src/adapter/f03_adapter_c.c
    fl_adapter_c.write("\n#include \"private/mat_adapter.h\"\n\n")
    fl_adapter_c.write("/* declaration of Fortran 2003 subroutines */\n")
    # writes macros of data type for real numbers and unit number of standard output
    fl_adapter_f.write("\n! macros related to Fortran module and matrix type\n")
    fl_adapter_f.write("#if !defined("+LANG_F_MODULE+")\n")
    fl_adapter_f.write("#define "+LANG_F_MODULE+" "+DEFAULT_F_MOD+"\n")
    fl_adapter_f.write("#endif\n")
    fl_adapter_f.write("#if !defined("+LANG_F_MATRIX+")\n")
    fl_adapter_f.write("#define "+LANG_F_MATRIX+" "+DEFAULT_F_MAT+"\n")
    fl_adapter_f.write("#endif\n\n")
    fl_adapter_f.write("! data type for real numbers\n")
    fl_adapter_f.write("#if defined(ENABLE_SINGLE_PRECISION)\n")
    fl_adapter_f.write("#define C_QREAL C_FLOAT\n")
    fl_adapter_f.write("#else\n")
    fl_adapter_f.write("#define C_QREAL C_DOUBLE\n")
    fl_adapter_f.write("#endif\n\n")
    fl_adapter_f.write("! unit number of standard output\n")
    fl_adapter_f.write("#if !defined(STDOUT)\n")
    fl_adapter_f.write("#define STDOUT 6\n")
    fl_adapter_f.write("#endif\n\n")
    # writes the declaration of Fortran 2003 subroutines in src/adapter/f03_adapter_c.c
    f03_adapter_c_declaration(f_adapter, fl_adapter_c, ADAPTER_MALLOC_FUN, malloc_fun_names, \
                              malloc_arg_types, malloc_arg_names)
    f03_adapter_c_declaration(f_adapter, fl_adapter_c, ADAPTER_FREE_FUN, mfree_fun_names, \
                              mfree_arg_types, mfree_arg_names)
    f03_adapter_c_declaration(f_adapter, fl_adapter_c, ADAPTER_NOMEM_FUN, nomem_fun_names, \
                              nomem_arg_types, nomem_arg_names)
    fl_adapter_c.write("\n")
    # writes the functions in src/adapter/f03_adapter_c.c
    f03_adapter_c_functions(f_adapter, fl_adapter_c, ADAPTER_MALLOC_FUN, malloc_fun_names, \
                            malloc_arg_types, malloc_arg_names)
    f03_adapter_c_functions(f_adapter, fl_adapter_c, ADAPTER_FREE_FUN, mfree_fun_names, \
                            mfree_arg_types, mfree_arg_names)
    f03_adapter_c_functions(f_adapter, fl_adapter_c, ADAPTER_NOMEM_FUN, nomem_fun_names, \
                            nomem_arg_types, nomem_arg_names)
    # writes the subroutines in src/adapter/f03_adapter_f.F90
    f03_adapter_f_subroutines(f_adapter, fl_adapter_f, ADAPTER_MALLOC_FUN, malloc_fun_names, \
                              malloc_arg_types, malloc_arg_names)
    f03_adapter_f_subroutines(f_adapter, fl_adapter_f, ADAPTER_FREE_FUN, mfree_fun_names, \
                              mfree_arg_types, mfree_arg_names)
    f03_adapter_f_subroutines(f_adapter, fl_adapter_f, ADAPTER_NOMEM_FUN, nomem_fun_names, \
                              nomem_arg_types, nomem_arg_names)
    print ">> closes files ..."
    fl_adapter_c.close()
    fl_adapter_f.close()
    print
    print ">> please replace the files of the Fortran 2003 adapter as:"
    print "   mv "+f_adapter+".cmake ${PATH_QCMATRIX}/cmake/"+f_adapter+".cmake"
    print "   mv "+f_adapter+"_c.h ${PATH_QCMATRIX}/include/private/"+f_adapter+"_c.h"
    print "   mv "+f_adapter+"_c.c ${PATH_QCMATRIX}/src/adapter/"+f_adapter+"_c.c"
    print "   mv "+f_adapter+"_f.F90 ${PATH_QCMATRIX}/src/adapter/"+f_adapter+"_f.F90"
    print
    return
