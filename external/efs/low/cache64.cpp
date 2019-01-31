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


#include "cache64.hpp"
using namespace std;

ostream & operator<<(ostream & os, const cacheline64 & rhs) {
    double a[DOUBLES_PER_CACHE_LINE];
    store(a, rhs);

    for (int i=0; i<DOUBLES_PER_CACHE_LINE; ++i) os << a[i] << " ";
    return os;
}

void PackArrays     (const double * __restrict__  arrays, UI32 ArraySize, cacheline64 * __restrict__  array8) {
    double * darr8 = (double*) array8;

    for (UI32 k=0; k<ArraySize; ++k) {
        for (UI32 n=0; n<DOUBLES_PER_CACHE_LINE; ++n) {
            darr8[DOUBLES_PER_CACHE_LINE*k+n] = arrays[n*ArraySize + k];
        }
    }
}

void UnPackArrays   (const cacheline64 * __restrict__ array8,  UI32 ArraySize, double * arrays) {
    const double * darr8 = (double*) array8;

    for (UI32 n=0; n<DOUBLES_PER_CACHE_LINE; ++n) {
        for (UI32 s=0; s<ArraySize; ++s) {
            arrays[n * ArraySize + s] = darr8[DOUBLES_PER_CACHE_LINE*s + n]; //0; //
        }
    }
}

