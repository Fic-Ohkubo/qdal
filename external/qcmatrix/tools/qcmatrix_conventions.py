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
#  This file defines the name and coding conventions in QcMatrix.
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

# name of external Fortran 90/2003 module
LANG_F_MODULE = "LANG_F_MODULE"
# default name of the external Fortran 90/2003 module
DEFAULT_F_MOD = "lib_matrix"
# name of matrix type in the external Fortran 90/2003 module
LANG_F_MATRIX = "LANG_F_MATRIX"
# default name of the matrix type
DEFAULT_F_MAT = "matrix_t"
# code indent
DEFAULT_CODE_INDENT = "    "
