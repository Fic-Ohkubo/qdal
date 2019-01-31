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

#ifndef __MPIWRAP__
#define __MPIWRAP__

class tensor2;
void Tensor2Broadcast(tensor2 & T);
void Tensor2Reduce   (tensor2 & T);

//variables needed to use this as a library called from fortran code

class Process {
  public:
    char processor_name[128 - 3*4];
    int namelen;
    int numprocs;
    int rank;

    void Init();
    void End();
    bool IsMaster() const;
    void Set1NodeOMPI () const;
    bool CheckMod(int n) const;
};

extern Process ThisNode;

#endif

