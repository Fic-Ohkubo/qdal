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


#include "MPIwrap.hpp"

#ifdef _OPENMP
    #include <omp.h>

#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1

#endif


//when compiling with fortran code, fotran wrappers are needed
#if defined(VAR_MPI) || defined(C_MPI)
    #include "mpi.h"
#endif


Process ThisNode;

void Process::Init() {
  #if defined (C_MPI)
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);
  #elif defined(VAR_MPI)
    //do not initialize, since the code is being called from somewhere else
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);
  #else
    numprocs = 1;
    rank     = 0;
  #endif
}


void Process::End() {
  #ifdef C_MPI
    MPI_Finalize();
  #endif
}

bool Process::IsMaster() const {
  #if defined(VAR_MPI) || defined(C_MPI)
    return (rank==0);
  #else
    return true;
  #endif
}

void Process::Set1NodeOMPI () const {
  #ifdef C_MPI
    int omp_mt = omp_get_max_threads();
    int mpi_mt = numprocs;
    int ths = omp_mt/mpi_mt;
    int xth = omp_mt%mpi_mt;

    if (rank < xth) omp_set_num_threads(ths + 1);
    else            omp_set_num_threads(ths);
  #endif
}

//check whether an integer is congruent to the rank of the process modulo the total number of MPI processes
bool Process::CheckMod(int n) const {
  #if defined(VAR_MPI) || defined(C_MPI)
    return (n%numprocs == rank);
  #else
    return true;
  #endif
}



#include "../math/tensors.hpp"

void Tensor2Broadcast(tensor2 & T) {

  #if defined(VAR_MPI) || defined(C_MPI)
    int ierr;

    MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Bcast (T.c2, T.n*T.m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  #endif
}

void Tensor2Reduce(tensor2 & T) {

  #if defined(VAR_MPI)
    int ierr;

    tensor2 S;
    if (ThisNode.IsMaster()) S.setsize(T.n, T.m);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(T.c2, S.c2, T.n*T.m, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (ThisNode.IsMaster()) T = S;

  #elif defined(C_MPI)
    int ierr;

    tensor2 S;
    S.setsize(T.n, T.m);

    MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Reduce(T.c2, S.c2, T.n*T.m, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Bcast (S.c2, S.n*S.m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    T = S;
  #endif
}


