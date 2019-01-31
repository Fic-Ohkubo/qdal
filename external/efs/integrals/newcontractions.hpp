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


#ifndef __NEW_CONTRACTIONS__
#define __NEW_CONTRACTIONS__

#include "../low/cache32.hpp"
#include "../low/cache64.hpp"
#include "../defs.hpp"

class RotationMatrix;

struct RotationMatrices64 {
    cacheline64 PP [3*3];
    cacheline64 DD [5*5];
    cacheline64 FF [7*7];
    cacheline64 GG [9*9];
    //cacheline64 HH[11*11];

    void From (const RotationMatrix * RM, int L);
} __attribute__((aligned(CACHE_LINE_SIZE)));

struct RotationMatrices32 {
    cacheline32 PP [3*3];
    cacheline32 DD [5*5];
    cacheline32 FF [7*7];
    cacheline32 GG [9*9];
    //cacheline32 HH[11*11];

    void From (const RotationMatrix     * RM, int L);
    void From (const RotationMatrices64 * RM, int L);
} __attribute__((aligned(CACHE_LINE_SIZE)));

//dummy class


typedef void (*JXD64) (cacheline64 * F, const cacheline64 * T, const cacheline64 * D);
typedef void (*JXD32) (cacheline32 * F, const cacheline32 * T, const cacheline32 * D);

typedef void (*JXrot64) (cacheline64 * MR, const RotationMatrices64 * __restrict__ RM);
typedef void (*JXrot32) (cacheline32 * MR, const RotationMatrices32 * __restrict__ RM);


struct DMD {

    JXrot64 R64  [LMAX][LMAX];
    JXrot64 RT64 [LMAX][LMAX];

    JXrot32 R32  [LMAX][LMAX];
    JXrot32 RT32 [LMAX][LMAX];


    //coulomb
    //*******
    JXD64 J64abcd[LMAX][LMAX][LMAX][LMAX];
    JXD64 J64cdab[LMAX][LMAX][LMAX][LMAX];

    JXD64 J64abcc[LMAX][LMAX][LMAX];
    JXD64 J64ccab[LMAX][LMAX][LMAX];

    JXD64 J64aacd[LMAX][LMAX][LMAX];
    JXD64 J64cdaa[LMAX][LMAX][LMAX];

    JXD64 J64aacc[LMAX][LMAX];
    JXD64 J64ccaa[LMAX][LMAX];

    JXD64 J64abab[LMAX][LMAX]; //ABAB


    //exchange
    //********
    JXD64 X64acbd[LMAX][LMAX][LMAX][LMAX];
    JXD64 X64adbc[LMAX][LMAX][LMAX][LMAX];
    JXD64 X64bdac[LMAX][LMAX][LMAX][LMAX];
    JXD64 X64bcad[LMAX][LMAX][LMAX][LMAX];

    JXD64 X64acad[LMAX][LMAX][LMAX];
    JXD64 X64adac[LMAX][LMAX][LMAX];

    JXD64 X64acbc[LMAX][LMAX][LMAX];
    JXD64 X64bcac[LMAX][LMAX][LMAX];

    JXD64 X64acac[LMAX][LMAX];

    JXD64 X64aabb[LMAX][LMAX];
    JXD64 X64bbaa[LMAX][LMAX];
    JXD64 X64abba[LMAX][LMAX];



    //coulomb
    //*******
    JXD32 J32abcd[LMAX][LMAX][LMAX][LMAX];
    JXD32 J32cdab[LMAX][LMAX][LMAX][LMAX];

    JXD32 J32abcc[LMAX][LMAX][LMAX];
    JXD32 J32ccab[LMAX][LMAX][LMAX];

    JXD32 J32aacd[LMAX][LMAX][LMAX];
    JXD32 J32cdaa[LMAX][LMAX][LMAX];

    JXD32 J32aacc[LMAX][LMAX];
    JXD32 J32ccaa[LMAX][LMAX];

    JXD32 J32abab[LMAX][LMAX]; //ABAB


    //exchange
    //********
    JXD32 X32acbd[LMAX][LMAX][LMAX][LMAX];
    JXD32 X32adbc[LMAX][LMAX][LMAX][LMAX];
    JXD32 X32bdac[LMAX][LMAX][LMAX][LMAX];
    JXD32 X32bcad[LMAX][LMAX][LMAX][LMAX];

    JXD32 X32acad[LMAX][LMAX][LMAX];
    JXD32 X32adac[LMAX][LMAX][LMAX];

    JXD32 X32acbc[LMAX][LMAX][LMAX];
    JXD32 X32bcac[LMAX][LMAX][LMAX];

    JXD32 X32acac[LMAX][LMAX];

    JXD32 X32aabb[LMAX][LMAX];
    JXD32 X32bbaa[LMAX][LMAX];
    JXD32 X32abba[LMAX][LMAX];




    DMD();
};

#endif
