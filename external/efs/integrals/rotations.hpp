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


#ifndef __ROTATIONS__
#define __ROTATIONS__

#include "../defs.hpp"
#include "../math/affine.hpp"
#include "../low/cache64.hpp"
#include "../low/cache32.hpp"

class AtomProd;

class ERIgeometry {
    public:
    double ABz;
    double CDy;
    double CDz;
    vector3 ACv;

    double c;
    double s;

    void MakeRotation4(const AtomProd & AP12, const AtomProd & AP34);
    void MakeRotation3(const AtomProd & AP12, const AtomProd & AP34);
    void MakeRotation2S(const AtomProd & AP12);
    void MakeRotation2(const point & A, const point & C);
    void MakeRotation1(const point & A);

    void Adjust(bool invAB, bool invCD);
};

struct ERIgeometries64 {
    cacheline64 ABz;
    cacheline64 CDy;
    cacheline64 CDz;
    cacheline64 ACx;
    cacheline64 ACy;
    cacheline64 ACz;
} __attribute__((aligned(CACHE_LINE_SIZE)));

struct ERIgeometries32 {
    cacheline32 ABz;
    cacheline32 CDy;
    cacheline32 CDz;
    cacheline32 ACx;
    cacheline32 ACy;
    cacheline32 ACz;
} __attribute__((aligned(CACHE_LINE_SIZE)));

void PackGeometries (const ERIgeometry     * vars,  ERIgeometries64  & vars8);
void PackGeometries (const ERIgeometry     * vars,  ERIgeometries32 & vars16);
void PackGeometries (const ERIgeometries64 * vars8, ERIgeometries32 & vars16);
void PackGeometries (const ERIgeometries64 & vars8, ERIgeometries32 & vars16);


//dummy struct
struct ERIgeometries {
    cacheline ABz;
    cacheline CDy;
    cacheline CDz;
    cacheline ACx;
    cacheline ACy;
    cacheline ACz;
} __attribute__((aligned(CACHE_LINE_SIZE)));



class RotationMatrix {
  public:
    double PP[3][3];
    double DD[5][5];
    double FF[7][7];
    double GG[9][9];
    //double HH[11][11];

    void From(const RotationMatrix & rhs, const ERIgeometry & rot, int L);
    void From(const ERIgeometry & rot, int L);
    void From(const point & Ap, const point & Cp, int L);

    
    RotationMatrix() {
        for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) PP[i][j] = 0;
        for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) DD[i][j] = 0;
        for (int i=0; i<7; ++i) for (int j=0; j<7; ++j) FF[i][j] = 0;
        for (int i=0; i<9; ++i) for (int j=0; j<9; ++j) GG[i][j] = 0;
    }

};


#endif
