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


    eigen.cpp:     Just a couple of function wrappers for dsyev_

*/

#include <iostream>
#include "../math/eigen.hpp"
#include "../math/tensors.hpp"
using namespace std;

extern "C" void dsyev_(const char * JOBZ,                    const char * UPLO, const int * N, double * A, const int * LDA, double *W, double *WORK, const int *LWORK, int *INFO);


void DiagonalizeV(tensor2 & A, double * w)  {
    int info, lwork;
    double wkopt;
    double * work;

    int n = A.n;

    lwork = -1;
    dsyev_( "V", "L", &n, A.c2, &n, w, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = new double[lwork];

    dsyev_( "V", "L", &n, A.c2, &n, w,   work, &lwork, &info );

    if (info>0) throw 1412; //meaningless error code

    delete[] work;
}

void DiagonalizeE(const tensor2 & A, double * w)  {
    int info, lwork;
    double wkopt;
    double * work;

    int n = A.n;

    lwork = -1;
    dsyev_( "N", "L", &n, A.c2, &n, w, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = new double[lwork];

    dsyev_( "N", "L", &n, A.c2, &n, w,   work, &lwork, &info );

    if (info>0) throw 1413; //another meaningless error code

    delete[] work;
}




