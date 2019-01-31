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


#include "cache32.hpp"
using namespace std;

ostream & operator<<(ostream & os, const cacheline32 & rhs) {
    float a[FLOATS_PER_CACHE_LINE];
    store(a, rhs);

    for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) os << a[i] << " ";
    return os;
}

void PackArrays     (const float * __restrict__  arrays, UI32 ArraySize, cacheline32 * __restrict__  array16) {
    float * darr16 = (float*) array16;

    for (UI32 k=0; k<ArraySize; ++k) {
        for (UI32 n=0; n<FLOATS_PER_CACHE_LINE; ++n) {
            darr16[FLOATS_PER_CACHE_LINE*k+n] = arrays[n*ArraySize + k];
        }
    }
}

void UnPackArrays   (const cacheline32 * __restrict__ array16,  UI32 ArraySize, float * arrays) {
    const float * darr16 = (float*) array16;

    for (UI32 n=0; n<FLOATS_PER_CACHE_LINE; ++n) {
        for (UI32 s=0; s<ArraySize; ++s) {
            arrays[n * ArraySize + s] = darr16[FLOATS_PER_CACHE_LINE*s + n]; //0; //
        }
    }
}
