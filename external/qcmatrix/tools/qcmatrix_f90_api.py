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
#  This file generates the files in Fortran 90 APIs:
#  (1) cmake/f90_api.cmake
#  (2) src/qcmat/f90/f90_api_c.c
#  (3) src/qcmat/f90/f90_api_f.F90
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

# size in bytes of the Fortran 90 integer (for the offset)
SIZEOF_VOID_P = "SIZEOF_VOID_P"
# default size in bytes of the Fortran 90 integer
DEFAULT_SIZEOF_VOID_P = "8"

# writes the declaration of Fortran 90 functions in src/qcmat/f90/f90_api_f.F90
# - f_api: name of the Fortran APIs
# - fl_api_f: file object of src/qcmat/f90/f90_api_f.F90
# - c_functions: list of C functions
# - c_function_names: names of C functions
def f90_api_f_declaration(f_api, fl_api_f, c_functions, c_function_names):
    # code indent in QcMatrix
    from qcmatrix_conventions import DEFAULT_CODE_INDENT
    print ">> writes the declaration of Fortran 90 functions in src/qcmat/f90/"+f_api+"_f.F90 ..."
    for ifun in xrange(len(c_functions)):
        len_fun_list = len(c_functions[ifun])
        for jfun in xrange(0,len_fun_list,3):
            kfun = jfun/3
            # we need to write the C preprocessor
            if len_fun_list>2:
                # the first if-else directive
                if jfun==0:
                    fl_api_f.write("#if "+c_functions[ifun][2]+"\n")
                # there are even numbers of if-else directives
                elif len_fun_list%3==0:
                    fl_api_f.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                # there are odd numbers of if-else directive, the last one is "#else"
                else:
                    if jfun<len_fun_list-2:
                        fl_api_f.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                    else:
                        fl_api_f.write("#else\n")
            fl_api_f.write(DEFAULT_CODE_INDENT+"public :: QcMat"+c_function_names[ifun][kfun]+"\n")
        if len_fun_list>2:
            fl_api_f.write("#endif\n")
    return

# writes the functions in src/qcmat/f90/f90_api_c.c to CMake file
# - f_api: name of the Fortran APIs
# - fl_cmake: file object of CMake file
# - c_functions: list of C functions
# - c_function_names: names of C functions
def f90_api_cmake(f_api, fl_cmake, c_functions, c_function_names):
    # code indent in QcMatrix
    from qcmatrix_conventions import DEFAULT_CODE_INDENT
    print ">> writes the functions in src/qcmat/f90/"+f_api+ \
          "_c.c to the CMake file cmake/"+f_api+".cmake ..."
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
            fl_cmake.write(code_indent+"SET(FC_MANGLING_SUB ${FC_MANGLING_SUB} "+f_api+ \
                           "_"+"QcMat"+c_function_names[ifun][kfun]+")\n")
        if len_fun_list>2:
            fl_cmake.write("ENDIF()\n")
    return

# writes the functions in src/qcmat/f90/f90_api_c.c
# - f_api: name of the Fortran APIs
# - fl_api_c: file object of src/qcmat/f90/f90_api_c.c
# - c_functions: list of C functions
# - c_function_names: names of C functions
# - c_argument_types: argument types in the C functions
# - c_argument_names: argument names in the C functions
def f90_api_c_functions(f_api, fl_api_c, c_functions, c_function_names, \
                        c_argument_types, c_argument_names):
    # code indent in QcMatrix
    from qcmatrix_conventions import DEFAULT_CODE_INDENT
    print ">> writes the functions in src/qcmat/f90/"+f_api+"_c.c ..."
    for ifun in xrange(len(c_functions)):
        len_fun_list = len(c_functions[ifun])
        for jfun in xrange(0,len_fun_list,3):
            print ">> writes function "+c_functions[ifun][jfun]+" ..."
            # intent of arguments
            arg_intent_size = c_functions[ifun][jfun+1].split(", ")
            # we need to write the C preprocessor
            if len_fun_list>2:
                # the first if-else directive
                if jfun==0:
                    fl_api_c.write("#if "+c_functions[ifun][2]+"\n")
                # there are even numbers of if-else directives
                elif len_fun_list%3==0:
                    fl_api_c.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                # there are odd numbers of if-else directive, the last one is "#else"
                else:
                    if jfun<len_fun_list-2:
                        fl_api_c.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                    else:
                        fl_api_c.write("#else\n")
            # string containing the arguments in the function of APIs
            str_api_args = ""
            # string containing the arguments in the function of QcMatrix
            str_fun_args = ""
            kfun = jfun/3
            # list of matrix arguments
            mat_argument = []
            mat_alloc = []
            mat_dealloc = []
            for iarg in xrange(len(c_argument_types[ifun][kfun])):
                # gets the size of argument
                len_arg_intent = arg_intent_size[iarg].find("(")
                if len_arg_intent==-1:
                    argument_intent = arg_intent_size[iarg]
                    argument_size = ""
                else:
                    argument_intent = arg_intent_size[iarg][:len_arg_intent]
                    argument_size = arg_intent_size[iarg][len_arg_intent+1:-1]
                # FIXME: quick solution for QcMat
                if c_argument_types[ifun][kfun][iarg]=="QcMat":
                    str_api_args = str_api_args+c_argument_types[ifun][kfun][iarg]+"_ptr " \
                                 + c_argument_names[ifun][kfun][iarg]+", "
                    str_fun_args = str_fun_args+c_argument_names[ifun][kfun][iarg][1:]+"->f90_mat, "
                    mat_argument.append(c_argument_names[ifun][kfun][iarg][1:]+"->f90_mat")
                    mat_alloc.append(argument_intent=="out")
                    mat_dealloc.append(argument_intent=="inout")
                # argument is a pointer
                elif c_argument_names[ifun][kfun][iarg][0]=="*":
                    str_api_args = str_api_args+c_argument_types[ifun][kfun][iarg]+" " \
                                 + c_argument_names[ifun][kfun][iarg]+", "
                    str_fun_args = str_fun_args+c_argument_names[ifun][kfun][iarg][1:]+", "
                else:
                    # Fortran always passes the address of an argument
                    str_api_args = str_api_args+c_argument_types[ifun][kfun][iarg]+" *" \
                                 + c_argument_names[ifun][kfun][iarg]+", "
                    # argument is an array
                    if len_arg_intent!=-1:
                        str_fun_args = str_fun_args+c_argument_names[ifun][kfun][iarg]+", "
                    else:
                        str_fun_args = str_fun_args+"*"+c_argument_names[ifun][kfun][iarg]+", "
            fl_api_c.write("QVoid "+f_api+"_QcMat"+c_function_names[ifun][kfun]+ \
                           "("+str_api_args+"QErrorCode *ierr)\n")
            fl_api_c.write("{\n")
            # writes the allocation of memory for matrices
            for imat in xrange(len(mat_argument)):
                if mat_alloc[imat]:
                    fl_api_c.write(DEFAULT_CODE_INDENT+mat_argument[imat]+" = (QcMat *)malloc(sizeof(QcMat));\n")
                    fl_api_c.write(DEFAULT_CODE_INDENT+"if ("+mat_argument[imat]+"==NULL) {\n")
                    fl_api_c.write(2*DEFAULT_CODE_INDENT+                                        \
                                   "QErrorExit(FILE_AND_LINE, \"failed to allocate memory for "+ \
                                   mat_argument[imat]+"\");\n")
                    fl_api_c.write(DEFAULT_CODE_INDENT+"}\n")
            fl_api_c.write(DEFAULT_CODE_INDENT+"*ierr = QcMat"+c_function_names[ifun][kfun]+ \
                           "("+str_fun_args[:-2]+");\n")
            # writes the deallocation of memory for matrices
            for imat in xrange(len(mat_argument)):
                if mat_dealloc[imat]:
                    fl_api_c.write(DEFAULT_CODE_INDENT+"free("+mat_argument[imat]+");\n")
                    fl_api_c.write(DEFAULT_CODE_INDENT+mat_argument[imat]+" = NULL;\n")
            fl_api_c.write("}\n")
        if len_fun_list>2:
            fl_api_c.write("#endif\n")
        fl_api_c.write("\n")
    return

# writes the functions in src/qcmat/f90/f90_api_f.F90
# - f_api: name of the Fortran APIs
# - fl_api_f: file object of src/qcmat/f90/f90_api_f.F90
# - c_functions: list of C functions
# - c_function_names: names of C functions
# - c_argument_types: argument types in the C functions
# - c_argument_names: argument names in the C functions
def f90_api_f_functions(f_api, fl_api_f, c_functions, c_function_names, \
                        c_argument_types, c_argument_names):
    # code indent in QcMatrix
    from qcmatrix_conventions import DEFAULT_CODE_INDENT
    print ">> writes the functions in src/qcmat/f90/"+f_api+"_f.F90 ..."
    # conversion list for types between C and Fortran
    fc_conversion = {"QcMat": "type(QcMat)",          \
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
                    fl_api_f.write("#if "+c_functions[ifun][2]+"\n")
                # there are even numbers of if-else directives
                elif len_fun_list%3==0:
                    fl_api_f.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                # there are odd numbers of if-else directive, the last one is "#else"
                else:
                    if jfun<len_fun_list-2:
                        fl_api_f.write("#elif "+c_functions[ifun][jfun+2]+"\n")
                    else:
                        fl_api_f.write("#else\n")
            # generates a string containing the arguments
            str_api_args = ""
            kfun = jfun/3
            print ">> writes function QcMat"+c_function_names[ifun][kfun]+" ..."
            # intent of arguments
            arg_intent_size = c_functions[ifun][jfun+1].split(", ")
            str_arg_intent = ""
            # string containing the arguments for the function in src/qcmat/f90/f90_api_c.c
            str_api_c_args = ""
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
                if argument_size!="":
                    argument_size = "("+argument_size+")"
                str_api_args = str_api_args \
                             + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
                str_arg_intent = str_arg_intent+sub_code_indent                      \
                               + fc_conversion[c_argument_types[ifun][kfun][iarg]]   \
                               + ", intent("+argument_intent+") :: "                 \
                               + c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                               + argument_size+"\n"
                # characters need to be terminated by char(0)
                if c_argument_types[ifun][kfun][iarg]=="QChar":
                    if argument_intent!="in":
                        raise ValueError(">> intent "+argument_intent+" is not supported!")
                    str_api_c_args = str_api_c_args \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:] \
                                   + "//char(0), "
                else:
                    str_api_c_args = str_api_c_args \
                                   + c_argument_names[ifun][kfun][iarg][start_arg_name:]+", "
            fl_api_f.write(DEFAULT_CODE_INDENT+"function QcMat"+c_function_names[ifun][kfun]+ \
                           "("+str_api_args[:-2]+") result(ierr)\n")
            fl_api_f.write(sub_code_indent+"integer :: ierr\n")
            fl_api_f.write(str_arg_intent)
            # calls the function in src/qcmat/f90/f90_api_c.c
            fl_api_f.write(sub_code_indent+"call "+f_api+"_QcMat"+c_function_names[ifun][kfun]+ \
                           "("+str_api_c_args[:-2]+", ierr)\n")

            fl_api_f.write(DEFAULT_CODE_INDENT+"end function QcMat"+c_function_names[ifun][kfun]+"\n")
        if len_fun_list>2:
            fl_api_f.write("#endif\n")
        fl_api_f.write("\n")
    return

# function to generate the files in Fortran 90 APIs:
#  (1) cmake/f90_api.cmake
#  (2) src/qcmat/f90/f90_api_c.c
#  (3) src/qcmat/f90/f90_api_f.F90
def generate_f90_api(f_api):
    # prints creating date
    from datetime import date
    # conventions in QcMatrix
    from qcmatrix_conventions import LANG_F_MODULE,DEFAULT_F_MOD, \
                                    LANG_F_MATRIX,DEFAULT_F_MAT, \
                                    DEFAULT_CODE_INDENT
    # QcMatrix license
    from qcmatrix_license import write_license
    # QcMatrix functions
    from qcmatrix_functions import API_MALLOC_FUN, \
                                  API_FREE_FUN,   \
                                  API_NOMEM_FUN,  \
                                  extract_functions
    print ">> extracts function names and arguments ..."
    malloc_fun_names,malloc_arg_types,malloc_arg_names = extract_functions(API_MALLOC_FUN)
    mfree_fun_names,mfree_arg_types,mfree_arg_names = extract_functions(API_FREE_FUN)
    nomem_fun_names,nomem_arg_types,nomem_arg_names = extract_functions(API_NOMEM_FUN)
    print ">> number of functions", len(API_MALLOC_FUN)+len(API_FREE_FUN)+len(API_NOMEM_FUN)
    # creates file objects
    print ">> open files in the Fortran 90 APIs ..."
    fl_cmake = open(f_api+".cmake", "w")
    fl_api_c = open(f_api+"_c.c", "w")
    fl_api_f = open(f_api+"_f.F90", "w")
    # writes license information
    file_description = "   This file implements the functions of Fortran 90 APIs.\n"
    write_license(fl_api_c, "C", file_description)
    file_description = "!!  This file defines the QcMat type for the Fortran 90 APIs.\n"
    write_license(fl_api_f, "Fortran", file_description)
    # CMake generated header file with auto-detected mangling 
    mangling_header = "f90_mangling.h"
    fl_cmake.write("# CMake file for Fortran 90 APIs\n")
    fl_cmake.write("# "+date.today().isoformat()+", Bin Gao:\n")
    fl_cmake.write("# * generated by tools/qcmatrix_fortran.py\n")
    fl_cmake.write("# Sets Fortran compiler flags\n")
    fl_cmake.write("SET(API_F_FLAGS \"-D"+SIZEOF_VOID_P+"=${CMAKE_SIZEOF_VOID_P}\")\n")
    fl_cmake.write("SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS} ${API_F_FLAGS}\")\n")
    fl_cmake.write("# Source codes of the APIs\n")
    fl_cmake.write("SET(API_SRCS\n")
    fl_cmake.write("    ${LIB_QCMATRIX_PATH}/src/qcmat/f90/"+f_api+"_f.F90\n")
    fl_cmake.write("    ${LIB_QCMATRIX_PATH}/src/qcmat/f90/"+f_api+"_c.c)\n")
    fl_cmake.write("# Name mangling for the functions in src/qcmat/f90/"+f_api+"_c.c\n")
    # includes header files in src/qcmat/f90/f90_api_c.c
    fl_api_c.write("\n#include \"qcmatrix.h\"\n")
    fl_api_c.write("/* uses CMake generated header file with auto-detected mangling */\n")
    fl_api_c.write("#include \""+mangling_header+"\"\n\n")
    fl_api_c.write("/* To be consistent with QcMat type defined in f90_api_f.F90, the struct QcMat_ptr\n")
    fl_api_c.write("   only has one member as the QcMat pointer, so that type(QcMat) can be converted\n")
    fl_api_c.write("   to struct QcMat_prt directly, inspired by the PETSc library */\n")
    fl_api_c.write("typedef struct {\n")
    fl_api_c.write("    QcMat *f90_mat;\n")
    fl_api_c.write("} QcMat_ptr;\n\n")
    # writes the module in src/qcmat/f90/f90_api_f.F90 for the Fortran 90 APIs
    fl_api_f.write("\n! size in bytes of the Fortran 90 integer (for the offset)\n")
    fl_api_f.write("#if !defined("+SIZEOF_VOID_P+")\n")
    fl_api_f.write("#define "+SIZEOF_VOID_P+" "+DEFAULT_SIZEOF_VOID_P+"\n")
    fl_api_f.write("#endif\n\n")
    fl_api_f.write("module qcmatrix_f90\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"implicit none\n\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+ \
                   "! QcMat type (used by the user) whose address will be passed to QcMatrix functions,\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"! in which "+SIZEOF_VOID_P+ \
                   " is the size in bytes of the Fortran 90 integer (for\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"! the offset), inspired by the PETSc library\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"type, public :: QcMat\n")
    fl_api_f.write(2*DEFAULT_CODE_INDENT+"private\n")
    fl_api_f.write(2*DEFAULT_CODE_INDENT+"integer(kind="+SIZEOF_VOID_P+") f90_int\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"end type QcMat\n\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"! data type for real numbers\n")
    fl_api_f.write("#if defined(ENABLE_SINGLE_PRECISION)\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"integer, public, parameter :: QREAL = kind(1.0)\n")
    fl_api_f.write("#else\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"integer, public, parameter :: QREAL = kind(1.0D0)\n")
    fl_api_f.write("#endif\n\n")
    fl_api_f.write("! parameters defined in QcMatrix library that will be used in Fortran APIs\n")
    fl_api_f.write("#include \"private/f_api.h90\"\n\n")
    # writes the declaration of Fortran 90 functions in src/qcmat/f90/f90_api_f.F90
    fl_api_f.write(DEFAULT_CODE_INDENT+"! functions provided by the Fortran 90 APIs\n")
    f90_api_f_declaration(f_api, fl_api_f, API_MALLOC_FUN, malloc_fun_names)
    f90_api_f_declaration(f_api, fl_api_f, API_FREE_FUN, mfree_fun_names)
    f90_api_f_declaration(f_api, fl_api_f, API_NOMEM_FUN, nomem_fun_names)
    fl_api_f.write("\n")
    fl_api_f.write(DEFAULT_CODE_INDENT+"contains\n\n")
    # writes the Fortran 90 subroutines of the APIs to CMake file cmake/f90_api.cmake
    f90_api_cmake(f_api, fl_cmake, API_MALLOC_FUN, malloc_fun_names)
    f90_api_cmake(f_api, fl_cmake, API_FREE_FUN, mfree_fun_names)
    f90_api_cmake(f_api, fl_cmake, API_NOMEM_FUN, nomem_fun_names)
    # writes the functions in src/qcmat/f90/f90_api_c.c
    f90_api_c_functions(f_api, fl_api_c, API_MALLOC_FUN, malloc_fun_names, \
                        malloc_arg_types, malloc_arg_names)
    f90_api_c_functions(f_api, fl_api_c, API_FREE_FUN, mfree_fun_names, \
                        mfree_arg_types, mfree_arg_names)
    f90_api_c_functions(f_api, fl_api_c, API_NOMEM_FUN, nomem_fun_names, \
                        nomem_arg_types, nomem_arg_names)
    # writes the functions in src/qcmat/f90/f90_api_f.F90
    f90_api_f_functions(f_api, fl_api_f, API_MALLOC_FUN, malloc_fun_names, \
                        malloc_arg_types, malloc_arg_names)
    f90_api_f_functions(f_api, fl_api_f, API_FREE_FUN, mfree_fun_names, \
                        mfree_arg_types, mfree_arg_names)
    f90_api_f_functions(f_api, fl_api_f, API_NOMEM_FUN, nomem_fun_names, \
                        nomem_arg_types, nomem_arg_names)
    # writes the left part of Fortran 90 APIs
    fl_api_f.write("end module qcmatrix_f90\n")
    print ">> closes files ..."
    fl_cmake.close()
    fl_api_c.close()
    fl_api_f.close()
    print
    print ">> please replace the files of the Fortran 90 APIs as:"
    print "   mv "+f_api+".cmake ${PATH_QCMATRIX}/cmake/"+f_api+".cmake"
    print "   mv "+f_api+"_c.c ${PATH_QCMATRIX}/src/qcmat/f90/"+f_api+"_c.c"
    print "   mv "+f_api+"_f.F90 ${PATH_QCMATRIX}/src/qcmat/f90/"+f_api+"_f.F90"
    print
    return
