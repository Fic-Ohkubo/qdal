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
!!  This file returns the time used by calling \fn(cpu_time).
!!
!!  2014-04-07, Bin Gao:
!!  * first version

!% \brief sets current CPU time
!  \author Bin Gao
!  \date 2014-04-07
!% \param[QREAL]{out} curr_time is current time from \fn(cpu_time)
subroutine TimerGet(curr_time)
    implicit none
#include "api/qcmatrix_f_basic.h90"
    real(kind=QREAL), intent(out) :: curr_time
    ! gets the CPU elapsed time
    call cpu_time(curr_time)
    return
end subroutine TimerGet

!% \brief prints the CPU elapsed time
!  \author Bin Gao
!  \date 2014-04-07
!  \param[QREAL]{in} prev_time is the previous time
!  \param[character]{in} msg_timer contains message to print
!% \param[integer]{in} io_log is the logical unit number of the logfile
subroutine TimerView(prev_time, msg_timer, io_log)
    implicit none
#include "api/qcmatrix_f_basic.h90"
    real(kind=QREAL), intent(in) :: prev_time
    character*(*), intent(in) :: msg_timer
    integer(kind=4), intent(in) :: io_log
    real(kind=QREAL) curr_time  !current CPU elapsed time
    ! gets the CPU elapsed time
    call cpu_time(curr_time)
    write(io_log,100) msg_timer, curr_time-prev_time
    return
100 format("TimerView>> CPU TIME of ",A,1X,F16.8," seconds")
end subroutine TimerView
