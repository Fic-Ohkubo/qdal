/*
    Copyright 2013 Jaime Axel Rosal Sandberg

    This file is part of the EFS library.

    The EFS library is free software:  you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The EFS library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the EFS library.  If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef __CHRONO__
#define __CHRONO__

#include <ostream>
#include <sys/time.h>

class Chronometer {
    friend std::ostream & operator<<(std::ostream & os, const Chronometer & chrono);
  private:
    timeval timer[2];
    int total_secs;
    int total_usecs;
    bool running;

  public:
    Chronometer();
    void Start();
    void Stop();
    void Restart();
    double GetTotalTime() const;
};

std::ostream & operator<<(std::ostream & os, const Chronometer & chrono);

#endif
