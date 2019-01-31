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


#include "../integrals/rotations.hpp"
#include "../integrals/newcontractions.hpp"
#include "../integrals/atomprod.hpp"
#include "../math/affine.hpp"
#include "../math/angular.hpp"
using namespace LibAngular;


//optimize SphHarRot!!!
static const double MINS  = 1.e-30;
static const double MINS2 = 1.e-15;

inline void ConstructRotHH (vector3 & vx, vector3 & vy, vector3 & vz,  const vector3 & u, const vector3 & v) {
    vx(1,0,0);
    vy(0,1,0);
    vz(0,0,1);

    double m = u.x*u.x + u.y*u.y;

    if (m > MINS) {
        m += u.z*u.z;
        m = sqrt(m);

        //vector3 w = u;
        double wzz;

        double H;

        if (u.z>0) {
            wzz = u.z + m;
            H = m*wzz;
        }
        else {
            wzz = u.z - m;
            H = -m*wzz;
        }

        double iH = 1/H;

        vx.x -= u.x*u.x*iH;
        vx.y -= u.x*u.y*iH;
        vx.z -= u.x*wzz*iH;

        vy.x -= u.x*u.y*iH;
        vy.y -= u.y*u.y*iH;
        vy.z -= u.y*wzz*iH;

        vz.x -= u.x*wzz*iH;
        vz.y -= u.y*wzz*iH;
        vz.z -= wzz*wzz*iH;

        vz.x = -vz.x;
        vz.y = -vz.y;
        vz.z = -vz.z;
    }

    //return;



    double sx, sy;
    sx = vx.x*v.x + vx.y*v.y + vx.z*v.z;
    sy = vy.x*v.x + vy.y*v.y + vy.z*v.z;

    if (fabs(sx) > MINS2) {
        double m = sqrt(sx*sx + sy*sy);

        double H;

        if (sy>0) {
            sy += m;
            H = m*sy;
        }
        else {
            sy -= m;
            H = -m*sy;
        }

        double iH = 1/H;

        double xx = 1 - sx*sx*iH;
        double xy =   - sx*sy*iH;
        double yx =   - sy*sx*iH;
        double yy = 1 - sy*sy*iH;

        double a,b;

        a = xx*vx.x + xy*vy.x;
        b = yx*vx.x + yy*vy.x;
        vx.x = a; vy.x = b;

        a = xx*vx.y + xy*vy.y;
        b = yx*vx.y + yy*vy.y;
        vx.y = a; vy.y = b;

        a = xx*vx.z + xy*vy.z;
        b = yx*vx.z + yy*vy.z;
        vx.z = a; vy.z = b;


        vy.x = -vy.x;
        vy.y = -vy.y;
        vy.z = -vy.z;
    }

}

inline void ConstructRotHH (vector3 & vx, vector3 & vy, vector3 & vz,  const vector3 & u) {
    vx(1,0,0);
    vy(0,1,0);
    vz(0,0,1);

    double m = u.x*u.x + u.y*u.y;

    if (m > MINS) {
        m += u.z*u.z;
        m = sqrt(m);

        //vector3 w = u;
        double wzz;

        double H;

        if (u.z>0) {
            wzz = u.z + m;
            H = m*wzz;
        }
        else {
            wzz = u.z - m;
            H = -m*wzz;
        }

        double iH = 1/H;

        vx.x -= u.x*u.x*iH;
        vx.y -= u.x*u.y*iH;
        vx.z -= u.x*wzz*iH;

        vy.x -= u.x*u.y*iH;
        vy.y -= u.y*u.y*iH;
        vy.z -= u.y*wzz*iH;

        vz.x -= u.x*wzz*iH;
        vz.y -= u.y*wzz*iH;
        vz.z -= wzz*wzz*iH;

        vz.x = -vz.x;
        vz.y = -vz.y;
        vz.z = -vz.z;
    }
}

inline void ConstructRotHH (double & cc, double & ss, const vector3 & vx, const vector3 & vy, const vector3 & vz,  const vector3 & v) {
    double sx, sy;
    sx = vx.x*v.x + vx.y*v.y + vx.z*v.z;
    sy = vy.x*v.x + vy.y*v.y + vy.z*v.z;

    if (fabs(sx) > MINS2) {
        double m = sqrt(sx*sx + sy*sy);

        double H;

        if (sy>0) {
            sy += m;
            H = m*sy;
        }
        else {
            sy -= m;
            H = -m*sy;
        }

        double iH = 1/H;

        cc = 1 - sx*sx*iH;
        ss =   - sx*sy*iH;
    }
    else {
        cc = 1;
        ss = 0;
    }
}


//finds the basis rotation for spherical harmonics
//very inneficient, yet it is only called in an O(N) loop
//pre: vx,vy,vz is an orthonormal basis in R3
template<int L> inline void SphHarRot(double (&M)[2*L+1][2*L+1], const vector3 & vx, const vector3 & vy, const vector3 & vz) {
    double P [2*L+1][2*L+1];
    double Q [2*L+1][2*L+1];
    double S [2*L+1][2*L+1];


    double ca, sa;
    double cb, sb;
    double cc, sc;

    {
        double s = sqrt(vx.z*vx.z + vy.z*vy.z);
        cb =  vz.z;
        sb = -s;

        if (s>MINS2) {
        //if (1) {
            double is = 1/s;
            ca =   vx.z * is;
            sa = - vy.z * is;
            cc = - vz.x * is;
            sc = - vz.y * is;
        }
        else {
            ca = 1;
            sa = 0;
            cc = 1;
            sc = 0;
        }
    }

    double cn, sn;
    double cm, sm;

    {
        for (int j=0; j<=2*L; ++j) {
            P[j][L] = SHrotI[L] [j][L];
        }

        cn = 1; sn = 0;
        for (int l=1; l<=L; ++l) {
            cm = cn*ca - sn*sa;
            sm = cn*sa + sn*ca;

            for (int j=0; j<=2*L; ++j) {
                P[j][L+l] = cm * SHrotI[L] [j][L+l] + sm * SHrotI[L] [j][L-l];
                P[j][L-l] = cm * SHrotI[L] [j][L-l] - sm * SHrotI[L] [j][L+l];
            }

            cn = cm;
            sn = sm;
        }
    }

    //x-y (z-x) in-plane rotation
    {
        for (int j=0; j<=2*L; ++j)
            Q[L][j] = P[L][j];

        cn = 1; sn = 0;
        for (int l=1; l<=L; ++l) {
            cm = cn*cb - sn*sb;
            sm = cn*sb + sn*cb;

            for (int j=0; j<=2*L; ++j) {
                Q[L+l][j] = cm * P[L+l][j] - sm * P[L-l][j];
                Q[L-l][j] = cm * P[L-l][j] + sm * P[L+l][j];
            }

            cn = cm;
            sn = sm;
        }
    }


    //index permutation ^-1
    for (int i=0; i<=2*L; ++i) {
        for (int j=0; j<=2*L; ++j) {
            double sum = 0;
            for (int k=0; k<=2*L; ++k) {
                sum += SHrot[L][i][k] * Q[k][j];
            }
            S[i][j] = sum;
        }
    }


    //x-y in-plane rotation
    {
        for (int j=0; j<=2*L; ++j)
            M[L][j] = S[L][j];

        cn = 1; sn = 0;
        for (int l=1; l<=L; ++l) {
            cm = cn*cc - sn*sc;
            sm = cn*sc + sn*cc;

            for (int j=0; j<=2*L; ++j) {
                M[L+l][j] = cm * S[L+l][j] - sm * S[L-l][j];
                M[L-l][j] = cm * S[L-l][j] + sm * S[L+l][j];
            }

            cn = cm;
            sn = sm;
        }
    }
}

//rotates a spherical harmonic rotation matrix around the z axis
template<int L> inline void SphHarRot(double (&R)[2*L+1][2*L+1], const double (&Rz)[2*L+1][2*L+1], double cc, double ss) {
    double cn = 1;
    double sn = 0;
    double cm, sm;

    for (int j=0; j<=2*L; ++j) {
        R[j][L] =  Rz[j][L];
    }

    for (int l=1; l<=L; ++l) {
        cm = cn*cc - sn*ss;
        sm = cn*ss + sn*cc;

        for (int j=0; j<=2*L; ++j) {
            R[j][L-l] =  cm*Rz[j][L-l] - sm*Rz[j][L+l];
            R[j][L+l] =  sm*Rz[j][L-l] + cm*Rz[j][L+l];
        }

        cn = cm;
        sn = sm;
    }
}



void ERIgeometry::MakeRotation4(const AtomProd & AP12, const AtomProd & AP34) {
    vector3 CD = AP34.B - AP34.A;
    vector3 AC = AP34.A - AP12.A;

    ConstructRotHH(c,s, AP12.vx, AP12.vy, AP12.vz, CD);

    vector3 vx, vy;

    vx =    c*AP12.vx + s*AP12.vy;
    vy = (-s)*AP12.vx + c*AP12.vy;
    const vector3 & vz = AP12.vz;


    ABz = AP12.ABv.z;
    CDy = CD * vy;
    CDz = CD * vz;

    ACv = FrameRot ( vx, vy, vz, AC );
}

void ERIgeometry::MakeRotation3(const AtomProd & AP12, const AtomProd & AP34) {
    vector3 AC = AP34.A - AP12.A;

    ConstructRotHH(c,s, AP34.vx, AP34.vy, AP34.vz, AC);

    vector3 vx, vy;

    vx =    c*AP34.vx + s*AP34.vy;
    vy = (-s)*AP34.vx + c*AP34.vy;
    const vector3 & vz = AP34.vz;


    ABz = 0;
    CDy = 0;
    CDz = AP34.ABv.z;

    ACv = FrameRot ( vx, vy, vz, AC );
}

void ERIgeometry::MakeRotation2S(const AtomProd & AP12) {
    c = 1;
    s = 0;

    ABz = AP12.ABv.z;
    CDy = 0;
    CDz = AP12.ABv.z;

    ACv(0,0,0);
}

void ERIgeometry::MakeRotation2(const point & A, const point & C) {
    vector3 AC = C - A;

    vector3 vx, vy, vz;

    ConstructRotHH(vx,vy,vz, AC);

    //those values are 0 and wont be used, so they serve as storage for the original vector
    ABz = 0; AC.x;
    CDy = 0; AC.y;
    CDz = 0; AC.z;

    ACv = FrameRot ( vx, vy, vz, AC );
}

void ERIgeometry::MakeRotation1(const point & A) {
    ABz = 0;
    CDy = 0;
    CDz = 0;

    ACv(0,0,0);
}

void ERIgeometry::Adjust(bool invAB, bool invCD) {
    if (invAB) {
        ACv.z -=  ABz; //BC
        ABz  = -ABz;   //BA
    }
    if (invCD) {
        ACv.y +=  CDy; //AD
        ACv.z +=  CDz; //AD
        CDy  = -CDy;   //DC
        CDz  = -CDz;   //DC
    }
}

//shouldn't belong here, but whatevs.
void AtomProd::MakeRotation() {
    vector3 AB = B - A;

    ConstructRotHH(vx,vy,vz, AB);
    ABv = FrameRot ( vx, vy, vz, AB );

    if (LMAX>0) SphHarRot<1>(RM.PP, vx, vy, vz);
    if (LMAX>1) SphHarRot<2>(RM.DD, vx, vy, vz);
    if (LMAX>2) SphHarRot<3>(RM.FF, vx, vy, vz);
    if (LMAX>3) SphHarRot<4>(RM.GG, vx, vy, vz);
    //if (LMAX>4) SphHarRot<5>(RM.HH, vx, vy, vz);
}


void RotationMatrix::From(const RotationMatrix & rhs, const ERIgeometry & rot, int L) {
    if (L>0) SphHarRot<1>(PP, rhs.PP, rot.c, rot.s);
    if (L>1) SphHarRot<2>(DD, rhs.DD, rot.c, rot.s);
    if (L>2) SphHarRot<3>(FF, rhs.FF, rot.c, rot.s);
    if (L>3) SphHarRot<4>(GG, rhs.GG, rot.c, rot.s);
    //if (L>4) SphHarRot<5>(HH, rhs.HH, rot.c, rot.s);
}

void RotationMatrix::From(const ERIgeometry & rot, int L) {
    vector3 AC;
    AC.x =  rot.ABz;
    AC.y =  rot.CDy;
    AC.z =  rot.CDz;

    vector3 vx, vy, vz;
    ConstructRotHH(vx,vy,vz, AC);

    if (L>0) SphHarRot<1>(PP, vx, vy, vz);
    if (L>1) SphHarRot<2>(DD, vx, vy, vz);
    if (L>2) SphHarRot<3>(FF, vx, vy, vz);
    if (L>3) SphHarRot<4>(GG, vx, vy, vz);
    //if (L>4) SphHarRot<5>(HH, vx, vy, vz);
}

void RotationMatrix::From(const point & Ap, const point & Cp, int L) {
    vector3 AC;

    AC = Cp - Ap;

    vector3 vx, vy, vz;
    ConstructRotHH(vx,vy,vz, AC);

    if (L>0) SphHarRot<1>(PP, vx, vy, vz);
    if (L>1) SphHarRot<2>(DD, vx, vy, vz);
    if (L>2) SphHarRot<3>(FF, vx, vy, vz);
    if (L>3) SphHarRot<4>(GG, vx, vy, vz);
    //if (L>4) SphHarRot<5>(HH, vx, vy, vz);
}

void PackGeometries (const ERIgeometry * vars, ERIgeometries64 & vars8) {

    vars8.ABz.set ( vars[0].ABz, vars[1].ABz, vars[2].ABz, vars[3].ABz,
                    vars[4].ABz, vars[5].ABz, vars[6].ABz, vars[7].ABz);

    vars8.CDy.set ( vars[0].CDy, vars[1].CDy, vars[2].CDy, vars[3].CDy,
                    vars[4].CDy, vars[5].CDy, vars[6].CDy, vars[7].CDy);

    vars8.CDz.set ( vars[0].CDz, vars[1].CDz, vars[2].CDz, vars[3].CDz,
                    vars[4].CDz, vars[5].CDz, vars[6].CDz, vars[7].CDz);


    vars8.ACx.set ( vars[0].ACv.x, vars[1].ACv.x, vars[2].ACv.x, vars[3].ACv.x,
                    vars[4].ACv.x, vars[5].ACv.x, vars[6].ACv.x, vars[7].ACv.x);

    vars8.ACy.set ( vars[0].ACv.y, vars[1].ACv.y, vars[2].ACv.y, vars[3].ACv.y,
                    vars[4].ACv.y, vars[5].ACv.y, vars[6].ACv.y, vars[7].ACv.y);

    vars8.ACz.set ( vars[0].ACv.z, vars[1].ACv.z, vars[2].ACv.z, vars[3].ACv.z,
                    vars[4].ACv.z, vars[5].ACv.z, vars[6].ACv.z, vars[7].ACv.z);
}

void PackGeometries (const ERIgeometry * vars, ERIgeometries32 & vars16) {

    float ABz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float CDy[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float CDz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACx[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACy[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));

    for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) ABz[i] = vars[i].ABz;
    for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) CDy[i] = vars[i].CDy;
    for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) CDz[i] = vars[i].CDz;
    for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) ACx[i] = vars[i].ACv.x;
    for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) ACy[i] = vars[i].ACv.y;
    for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) ACz[i] = vars[i].ACv.z;

    vars16.ABz.set (ABz);
    vars16.CDy.set (CDy);
    vars16.CDz.set (CDz);
    vars16.ACx.set (ACx);
    vars16.ACy.set (ACy);
    vars16.ACz.set (ACz);
}

const int DPC = DOUBLES_PER_CACHE_LINE;
const int FPC = FLOATS_PER_CACHE_LINE;

void PackGeometries (const ERIgeometries64 * vars, ERIgeometries32 & vars16) {

    float ABz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float CDy[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float CDz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACx[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACy[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));

    for (int i=0; i<DPC; ++i) ABz[i] = vars[0].ABz(i);
    for (int i=0; i<DPC; ++i) CDy[i] = vars[0].CDy(i);
    for (int i=0; i<DPC; ++i) CDz[i] = vars[0].CDz(i);
    for (int i=0; i<DPC; ++i) ACx[i] = vars[0].ACx(i);
    for (int i=0; i<DPC; ++i) ACy[i] = vars[0].ACy(i);
    for (int i=0; i<DPC; ++i) ACz[i] = vars[0].ACz(i);

    for (int i=0; i<DPC; ++i) ABz[DPC+i] = vars[1].ABz(i);
    for (int i=0; i<DPC; ++i) CDy[DPC+i] = vars[1].CDy(i);
    for (int i=0; i<DPC; ++i) CDz[DPC+i] = vars[1].CDz(i);
    for (int i=0; i<DPC; ++i) ACx[DPC+i] = vars[1].ACx(i);
    for (int i=0; i<DPC; ++i) ACy[DPC+i] = vars[1].ACy(i);
    for (int i=0; i<DPC; ++i) ACz[DPC+i] = vars[1].ACz(i);


    vars16.ABz.set (ABz);
    vars16.CDy.set (CDy);
    vars16.CDz.set (CDz);
    vars16.ACx.set (ACx);
    vars16.ACy.set (ACy);
    vars16.ACz.set (ACz);
}

void PackGeometries (const ERIgeometries64 & vars, ERIgeometries32 & vars16) {

    float ABz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float CDy[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float CDz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACx[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACy[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));
    float ACz[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));

    for (int i=0; i<DPC; ++i) ABz[i] = vars.ABz(i);
    for (int i=0; i<DPC; ++i) CDy[i] = vars.CDy(i);
    for (int i=0; i<DPC; ++i) CDz[i] = vars.CDz(i);
    for (int i=0; i<DPC; ++i) ACx[i] = vars.ACx(i);
    for (int i=0; i<DPC; ++i) ACy[i] = vars.ACy(i);
    for (int i=0; i<DPC; ++i) ACz[i] = vars.ACz(i);

    for (int i=0; i<DPC; ++i) ABz[DPC+i] = 0;
    for (int i=0; i<DPC; ++i) CDy[DPC+i] = 0;
    for (int i=0; i<DPC; ++i) CDz[DPC+i] = 0;
    for (int i=0; i<DPC; ++i) ACx[DPC+i] = 0;
    for (int i=0; i<DPC; ++i) ACy[DPC+i] = 0;
    for (int i=0; i<DPC; ++i) ACz[DPC+i] = 0;


    vars16.ABz.set (ABz);
    vars16.CDy.set (CDy);
    vars16.CDz.set (CDz);
    vars16.ACx.set (ACx);
    vars16.ACy.set (ACy);
    vars16.ACz.set (ACz);
}


#include "../integrals/newcontractions.hpp"

void RotationMatrices64::From(const RotationMatrix * RM, int L) {
    if (L>0) {
        for (int i=0; i< 3; ++i)
            for (int j=0; j< 3; ++j) {
                int ij = 3*i+j;

                for (int k=0; k<DPC; ++k)
                    PP[ij](k) = RM[k].PP[i][j];
            }
    }
    if (L>1) {
        for (int i=0; i<5; ++i)
            for (int j=0; j<5; ++j) {
                int ij = 5*i+j;
                for (int k=0; k<DPC; ++k)
                    DD[ij](k) = RM[k].DD[i][j];
            }
    }
    if (L>2) {
        for (int i=0; i<7; ++i)
            for (int j=0; j<7; ++j)  {
                int ij = 7*i+j;
                for (int k=0; k<DPC; ++k)
                    FF[ij](k) = RM[k].FF[i][j];
            }
    }
    if (L>3) {
        for (int i=0; i<9; ++i)
            for (int j=0; j<9; ++j) {
                int ij = 9*i+j;
                for (int k=0; k<DPC; ++k)
                    GG[ij](k) = RM[k].GG[i][j];
            }
    }
    /*
    if (L>4) {
        for (int i=0; i<11; ++i)
            for (int j=0; j<11; ++j) {
                int ij = 11*i+j;
                for (int k=0; k<DPC; ++k)
                    HH[ij](k) = RM[k].HH[i][j];
            }
    }
    */
}


void RotationMatrices32::From(const RotationMatrix * RM, int L) {
    if (L>0) {
        for (int i=0; i< 3; ++i)
            for (int j=0; j< 3; ++j) {
                int ij = 3*i+j;

                for (int k=0; k<FPC; ++k)
                    PP[ij](k) = RM[k].PP[i][j];
            }
    }
    if (L>1) {
        for (int i=0; i<5; ++i)
            for (int j=0; j<5; ++j) {
                int ij = 5*i+j;
                for (int k=0; k<FPC; ++k)
                    DD[ij](k) = RM[k].DD[i][j];
            }
    }
    if (L>2) {
        for (int i=0; i<7; ++i)
            for (int j=0; j<7; ++j)  {
                int ij = 7*i+j;
                for (int k=0; k<FPC; ++k)
                    FF[ij](k) = RM[k].FF[i][j];
            }
    }
    if (L>3) {
        for (int i=0; i<9; ++i)
            for (int j=0; j<9; ++j) {
                int ij = 9*i+j;
                for (int k=0; k<FPC; ++k)
                    GG[ij](k) = RM[k].GG[i][j];
            }
    }
    /*
    if (L>4) {
        for (int i=0; i<11; ++i)
            for (int j=0; j<11; ++j) {
                int ij = 11*i+j;
                for (int k=0; k<FPC; ++k)
                    HH[ij](k) = RM[k].HH[i][j];
            }
    }
    */
}

void RotationMatrices32::From(const RotationMatrices64 * RM, int L) {
    if (L>0) {
        for (int i=0; i< 3; ++i)
            for (int j=0; j< 3; ++j) {
                int ij = 3*i+j;
                PP[ij].set(RM[0].PP[ij], RM[1].PP[ij]);
            }
    }
    if (L>1) {
        for (int i=0; i<5; ++i)
            for (int j=0; j<5; ++j) {
                int ij = 5*i+j;
                DD[ij].set(RM[0].DD[ij], RM[1].DD[ij]);
            }
    }
    if (L>2) {
        for (int i=0; i<7; ++i)
            for (int j=0; j<7; ++j)  {
                int ij = 7*i+j;
                FF[ij].set(RM[0].FF[ij], RM[1].FF[ij]);
            }
    }
    if (L>3) {
        for (int i=0; i<9; ++i)
            for (int j=0; j<9; ++j) {
                int ij = 9*i+j;
                GG[ij].set(RM[0].GG[ij], RM[1].GG[ij]);
            }
    }

}



