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
#  This file generates the files in:
#  (1) Fortran 90 adapter
#  (2) Fortran 03 adapter
#  (3) Fortran 90 APIs
#  (4) Fortran 03 APIs
#
#  2014-02-20, Bin Gao:
#  * moves functions for the adapters and APIs into python modules
#
#  2014-02-01, Bin Gao:
#  * supports Fortran 2003
# 
#  2013-11-02, Bin Gao:
#  * first version

#__author__ = "Bin Gao"
#__copyright__ = "Copyright 2012-2015"
#__license__ = "LGPLv3"
#__version__ = "0.1.0"
#__maintainer__ = "Bin Gao"
#__email__ = "bin.gao@uit.no"
#__status__ = "Development"

if __name__ == '__main__':
    # menu
    while True:
        print "(1) generate Fortran 90 adapter"
        print "(2) generate Fortran 2003 adapter"
        print "(3) generate Fortran 90 APIs"
        print "(4) generate Fortran 2003 APIs"
        print "(5) exit"
        try:
            user_choice = int(raw_input(">> enter your choice: "))
            if (user_choice>=1 and user_choice<=5):
                break
            else:
                continue
        except:
            continue
    if user_choice==1:
        from qcmatrix_f90_adapter import generate_f90_adapter
        print ">> your choice: (1) generate Fortran 90 adapater ..."
        generate_f90_adapter("f90_adapter")
    elif user_choice==2:
        from qcmatrix_f03_adapter import generate_f03_adapter
        print ">> your choice: (2) generate Fortran 2003 adapter ..."
        generate_f03_adapter("f03_adapter")
    elif user_choice==3:
        from qcmatrix_f90_api import generate_f90_api
        print ">> your choice: (3) generate Fortran 90 APIs ..."
        generate_f90_api("f90_api")
    elif user_choice==4:
        from qcmatrix_f03_api import generate_f03_api
        print ">> your choice: (4) generate Fortran 2003 APIs ..."
        generate_f03_api("f03_api")
