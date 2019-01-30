!!  QcMatrix: an abstract matrix library
!!  Copyright 2012-2015 Bin Gao
!!
!!  QcMatrix is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  QcMatrix is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with QcMatrix. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file contains subroutines for error handling.
!!
!!  2014-09-19, Bin Gao:
!!  * first version

    !% \brief prints error message and stops the program
    !  \author Bin Gao
    !  \date 2014-09-19
    !  \param[integer]{in} io_log IO of the logfile
    !  \param[integer]{in} line line number of the error occurred
    !% \param[character]{in} source_code name of source code where the error occurred
    subroutine error_exit(io_log, line, source_code)
        implicit none
        integer(kind=4), intent(in) :: io_log
        integer(kind=4), intent(in) :: line
        character*(*), intent(in) :: source_code
        write(io_log,100) line, source_code
        stop
100     format("error_exit>> error occurred at line ", I4," of ",A)
    end subroutine error_exit
