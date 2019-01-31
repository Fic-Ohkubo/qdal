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


        Implements the Spherical Kernel Sieve algorithm
        based on the recurrence relations by Fortunelli and Salvetti
*/

#include "integrals.hpp"
#include <set>
using namespace std;

const SI8 mask = 127; // = 0xEF = 0111 1111

inline SI8 abs(SI8 rhs) {
    return rhs&mask;
}

inline void absms(integralSph & I) {
    I.ma &= mask;
    I.mb &= mask;
    I.mc &= mask;
    I.md &= mask;
}

inline bool IsZero(const integralSph & I) {
    if (I.la<0 || I.lb<0 || I.lc<0 || I.ld<0) return true;
    if ((abs(I.ma)>I.la) || (abs(I.mb)>I.lb) || (abs(I.mc)>I.lc) || (abs(I.md)>I.ld)) return true;
    return false;
}

// ABx, bool ABy, bool ABz, bool CDx, bool CDy, bool CDz, bool ACx, bool ACy, bool ACz) {

void SKSa(integralSph I, set<integralSph> & iset, set<integralSph> & fset, bool vbool[9]) {

    //first block
    const integralSph IR1  (-1,  0,  0, 0, 0, 0, 0, 0,  1, 0, 0, 0,  0); // ABz
    const integralSph IR2  (-1, -1,  0, 0, 0, 0, 0, 0,  1, 0, 0, 0,  0); // NOT USED
    const integralSph IR3  (-1,  0, -1, 0, 0, 0, 0, 0,  0, 1, 0, 0,  0);
    const integralSph IR4  (-1, -1, -1, 1, 0, 0, 0, 0,  0, 1, 0, 0,  0);

    //AC
    const integralSph IR5  (-1,  0,  0, 0, 0, 0, 0, 0,  0, 1, 0, 0,  1);
    const integralSph IR6  (-1, -1,  0, 0, 0, 0, 0, 0,  0, 1, 0, 0,  1);
    const integralSph IR7  (-1,  0, -1, 0, 0, 0, 0, 0,  0, 2, 0, 0,  1);
    const integralSph IR8  (-1, -1, -1, 1, 0, 0, 0, 0,  0, 2, 0, 0,  1);

    //AB
    const integralSph IR9  (-1,  0, 0, 0,  0, 0, 0, 0,  1, 1, 0, 0,  1);
    const integralSph IR10 (-1, -1, 0, 0,  0, 0, 0, 0,  1, 1, 0, 0,  1); // NOT USED
    const integralSph IR11 (-1,  0, 0, 0, -1, 0, 0, 0,  0, 1, 0, 1,  1);
    const integralSph IR12 (-1, -1, 0, 0, -1, 1, 0, 0,  0, 1, 0, 1,  1);

    //CD
    const integralSph IR13 (-1,  0, 0, 0, 0, 0,  0, 0,  0, 1, 1, 0,  1);
    const integralSph IR14 (-1, -1, 0, 0, 0, 0,  0, 0,  0, 1, 1, 0,  1);
    const integralSph IR15 (-1,  0, 0, 0, 0, 0, -1, 0,  0, 1, 0, 1,  1);
    const integralSph IR16 (-1, -1, 0, 0, 0, 0, -1, 1,  0, 1, 0, 1,  1);

    absms(I);

    if (IsZero(I)) return; // out of boundaries (counts as 0)



    if (iset.count(I)) return; // don't perform redundant work
    iset.insert(I); // register the integral

    // iff end of recursion, just insert  the  integral in the set
    if (I.la==0) {
        fset.insert(I);
        return;
    }

    bool ABxy = vbool[0] || vbool[1];
    bool ABz  = vbool[2];

    bool CDxy = vbool[3] || vbool[4];
    bool CDz  = vbool[5];

    bool ACxy = vbool[6] || vbool[7];
    bool ACz  = vbool[8];


    if (ABz)  SKSa(I+IR1, iset, fset, vbool);
    if (ABxy) SKSa(I+IR2, iset, fset, vbool);
    SKSa(I+IR3, iset, fset, vbool);
    SKSa(I+IR4, iset, fset, vbool);

    if (ACz)  SKSa(I+IR5, iset, fset, vbool);
    if (ACxy) SKSa(I+IR6, iset, fset, vbool);
    SKSa(I+IR7, iset, fset, vbool);
    SKSa(I+IR8, iset, fset, vbool);

    if (ABz)  SKSa(I+IR9 , iset, fset, vbool);
    if (ABxy) SKSa(I+IR10, iset, fset, vbool);
    SKSa(I+IR11, iset, fset, vbool);
    SKSa(I+IR12, iset, fset, vbool);

    if (CDz)  SKSa(I+IR13, iset, fset, vbool);
    if (CDxy) SKSa(I+IR14, iset, fset, vbool);
    SKSa(I+IR15, iset, fset, vbool);
    SKSa(I+IR16, iset, fset, vbool);
}

void SKSb(integralSph I, set<integralSph> & iset, set<integralSph> & fset, bool vbool[9]) {

    //first block
    const integralSph IR1   (-1,  0,  0, 0, 0, 0, 0, 0,  1, 0, 0, 0,  0, true); // ABz
    const integralSph IR2   (-1, -1,  0, 0, 0, 0, 0, 0,  1, 0, 0, 0,  0, true); // NOT USED
    const integralSph IR3   (-1,  0, -1, 0, 0, 0, 0, 0,  0, 1, 0, 0,  0, true);
    const integralSph IR4   (-1, -1, -1, 1, 0, 0, 0, 0,  0, 1, 0, 0,  0, true);

    //AC
    const integralSph IR5   (-1,  0,  0, 0, 0, 0, 0, 0,  0, 1, 0, 0,  1, true);
    const integralSph IR6   (-1, -1,  0, 0, 0, 0, 0, 0,  0, 1, 0, 0,  1, true);
    const integralSph IR7   (-1,  0, -1, 0, 0, 0, 0, 0,  0, 2, 0, 0,  1, true);
    const integralSph IR8   (-1, -1, -1, 1, 0, 0, 0, 0,  0, 2, 0, 0,  1, true);

    //AB
    const integralSph IR9   (-1,  0, 0, 0,  0, 0, 0, 0,  1, 1, 0, 0,  1, true);
    const integralSph IR10  (-1, -1, 0, 0,  0, 0, 0, 0,  1, 1, 0, 0,  1, true); // NOT USED
    const integralSph IR11  (-1,  0, 0, 0, -1, 0, 0, 0,  0, 1, 0, 1,  1, true);
    const integralSph IR12  (-1, -1, 0, 0, -1, 1, 0, 0,  0, 1, 0, 1,  1, true);

    //CD
    const integralSph IR13  (-1,  0, 0, 0, 0, 0,  0, 0,  0, 1, 1, 0,  1, true);
    const integralSph IR14  (-1, -1, 0, 0, 0, 0,  0, 0,  0, 1, 1, 0,  1, true);
    const integralSph IR15  (-1,  0, 0, 0, 0, 0, -1, 0,  0, 1, 0, 1,  1, true);
    const integralSph IR16  (-1, -1, 0, 0, 0, 0, -1, 1,  0, 1, 0, 1,  1, true);

    // needed specially for B
    const integralSph IR1b  (-1,  0,  0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, true); // ABz
    const integralSph IR2b  (-1, -1,  0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, true); // NOT USED
    const integralSph IR9b  (-1,  0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,  1, true);
    const integralSph IR10b (-1, -1, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,  1, true); // NOT USED

    absms(I);

    if (IsZero(I)) return; // out of boundaries (counts as 0)



    if (iset.count(I)) return; // don't perform redundant work
    iset.insert(I); // register the integral

    // iff end of recursion, just insert  the  integral in the set
    if (I.lb==0) {
        fset.insert(I);
        return;
    }

    bool ABxy = vbool[0] || vbool[1];
    bool ABz  = vbool[2];

    bool CDxy = vbool[3] || vbool[4];
    bool CDz  = vbool[5];

    bool ACxy = vbool[6] || vbool[7];
    bool ACz  = vbool[8];


    if (ABz)  SKSb(I+IR1,  iset, fset, vbool);
    if (ABz)  SKSb(I+IR1b, iset, fset, vbool);
    if (ABxy) SKSb(I+IR2,  iset, fset, vbool);
    if (ABxy) SKSb(I+IR2b, iset, fset, vbool);
    SKSb(I+IR3, iset, fset, vbool);
    SKSb(I+IR4, iset, fset, vbool);

    if (ACz)  SKSb(I+IR5, iset, fset, vbool);
    if (ACxy) SKSb(I+IR6, iset, fset, vbool);
    SKSb(I+IR7, iset, fset, vbool);
    SKSb(I+IR8, iset, fset, vbool);

    if (ABz)  SKSb(I+IR9  , iset, fset, vbool);
    if (ABz)  SKSb(I+IR9b , iset, fset, vbool);
    if (ABxy) SKSb(I+IR10 , iset, fset, vbool);
    if (ABxy) SKSb(I+IR10b, iset, fset, vbool);
    SKSb(I+IR11, iset, fset, vbool);
    SKSb(I+IR12, iset, fset, vbool);

    if (CDz)  SKSb(I+IR13, iset, fset, vbool);
    if (CDxy) SKSb(I+IR14, iset, fset, vbool);
    SKSb(I+IR15, iset, fset, vbool);
    SKSb(I+IR16, iset, fset, vbool);
}

void SKSc(integralSph I, set<integralSph> & iset, set<integralSph> & fset, bool vbool[9]) {

    //first block
    const integralSph IR1  (-1,  0,  0, 0, 0, 0, 0, 0,  1, 0, 0, 0,  0, false, true); // ABz
    const integralSph IR2  (-1, -1,  0, 0, 0, 0, 0, 0,  1, 0, 0, 0,  0, false, true); // NOT USED
    const integralSph IR3  (-1,  0, -1, 0, 0, 0, 0, 0,  0, 1, 0, 0,  0, false, true);
    const integralSph IR4  (-1, -1, -1, 1, 0, 0, 0, 0,  0, 1, 0, 0,  0, false, true);

    //AC
    const integralSph IR5  (-1,  0,  0, 0, 0, 0, 0, 0,  0, 1, 0, 0,  1, false, true);
    const integralSph IR6  (-1, -1,  0, 0, 0, 0, 0, 0,  0, 1, 0, 0,  1, false, true);
    const integralSph IR7  (-1,  0, -1, 0, 0, 0, 0, 0,  0, 2, 0, 0,  1, false, true);
    const integralSph IR8  (-1, -1, -1, 1, 0, 0, 0, 0,  0, 2, 0, 0,  1, false, true);

    //AB
    const integralSph IR9  (-1,  0, 0, 0,  0, 0, 0, 0,  1, 1, 0, 0,  1, false, true);
    const integralSph IR10 (-1, -1, 0, 0,  0, 0, 0, 0,  1, 1, 0, 0,  1, false, true); // NOT USED
    const integralSph IR11 (-1,  0, 0, 0, -1, 0, 0, 0,  0, 1, 0, 1,  1, false, true);
    const integralSph IR12 (-1, -1, 0, 0, -1, 1, 0, 0,  0, 1, 0, 1,  1, false, true);

    //CD
    const integralSph IR13 (-1,  0, 0, 0, 0, 0,  0, 0,  0, 1, 1, 0,  1, false, true);
    const integralSph IR14 (-1, -1, 0, 0, 0, 0,  0, 0,  0, 1, 1, 0,  1, false, true);
    const integralSph IR15 (-1,  0, 0, 0, 0, 0, -1, 0,  0, 1, 0, 1,  1, false, true);
    const integralSph IR16 (-1, -1, 0, 0, 0, 0, -1, 1,  0, 1, 0, 1,  1, false, true);

    absms(I);

    if (IsZero(I)) return; // out of boundaries (counts as 0)



    if (iset.count(I)) return; // don't perform redundant work
    iset.insert(I); // register the integral

    // iff end of recursion, just insert  the  integral in the set
    if (I.lc==0) {
        fset.insert(I);
        return;
    }

    bool ABxy = vbool[0] || vbool[1];
    bool ABz  = vbool[2];

    bool CDxy = vbool[3] || vbool[4];
    bool CDz  = vbool[5];

    bool ACxy = vbool[6] || vbool[7];
    bool ACz  = vbool[8];


    if (CDz)  SKSc(I+IR1, iset, fset, vbool);
    if (CDxy) SKSc(I+IR2, iset, fset, vbool);
    SKSc(I+IR3, iset, fset, vbool);
    SKSc(I+IR4, iset, fset, vbool);

    if (ACz)  SKSc(I+IR5, iset, fset, vbool);
    if (ACxy) SKSc(I+IR6, iset, fset, vbool);
    SKSc(I+IR7, iset, fset, vbool);
    SKSc(I+IR8, iset, fset, vbool);

    if (CDz)  SKSc(I+IR9 , iset, fset, vbool);
    if (CDxy) SKSc(I+IR10, iset, fset, vbool);
    SKSc(I+IR11, iset, fset, vbool);
    SKSc(I+IR12, iset, fset, vbool);

    if (ABz)  SKSc(I+IR13, iset, fset, vbool);
    if (ABxy) SKSc(I+IR14, iset, fset, vbool);
    SKSc(I+IR15, iset, fset, vbool);
    SKSc(I+IR16, iset, fset, vbool);
}

void SKSd(integralSph I, set<integralSph> & iset, set<integralSph> & fset, bool vbool[9]) {

    //first block
    const integralSph IR1  (-1, 0,  0, 0,  0, 0,  0, 0,   1, 0, 0, 0,  0, true, true); // ABz
    const integralSph IR2  (-1,-1,  0, 0,  0, 0,  0, 0,   1, 0, 0, 0,  0, true, true); // NOT USED
    const integralSph IR3  (-1, 0, -1, 0,  0, 0,  0, 0,   0, 1, 0, 0,  0, true, true);
    const integralSph IR4  (-1,-1, -1, 1,  0, 0,  0, 0,   0, 1, 0, 0,  0, true, true);

    //AC
    const integralSph IR5  (-1, 0,  0, 0,  0, 0,  0, 0,   0, 1, 0, 0,  1, true, true);
    const integralSph IR6  (-1,-1,  0, 0,  0, 0,  0, 0,   0, 1, 0, 0,  1, true, true);
    const integralSph IR7  (-1, 0, -1, 0,  0, 0,  0, 0,   0, 2, 0, 0,  1, true, true);
    const integralSph IR8  (-1,-1, -1, 1,  0, 0,  0, 0,   0, 2, 0, 0,  1, true, true);

    //AB
    const integralSph IR9  (-1, 0,  0, 0,  0, 0,  0, 0,   1, 1, 0, 0,  1, true, true);
    const integralSph IR10 (-1,-1,  0, 0,  0, 0,  0, 0,   1, 1, 0, 0,  1, true, true); // NOT USED
    const integralSph IR11 (-1, 0,  0, 0, -1, 0,  0, 0,   0, 1, 0, 1,  1, true, true);
    const integralSph IR12 (-1,-1,  0, 0, -1, 1,  0, 0,   0, 1, 0, 1,  1, true, true);

    //CD
    const integralSph IR13 (-1, 0,  0, 0,  0, 0,  0, 0,   0, 1, 1, 0,  1, true, true);
    const integralSph IR14 (-1,-1,  0, 0,  0, 0,  0, 0,   0, 1, 1, 0,  1, true, true);
    const integralSph IR15 (-1, 0,  0, 0,  0, 0, -1, 0,   0, 1, 0, 1,  1, true, true);
    const integralSph IR16 (-1,-1,  0, 0,  0, 0, -1, 1,   0, 1, 0, 1,  1, true, true);



    // needed specially for D
    const integralSph IR1b  (-1, 0,  0, 0,  0, 0,  0, 0,   0, 0, 0, 0,  0, true, true); // ABz
    const integralSph IR2b  (-1,-1,  0, 0,  0, 0,  0, 0,   0, 0, 0, 0,  0, true, true); // NOT USED
    const integralSph IR9b  (-1, 0,  0, 0,  0, 0,  0, 0,   0, 1, 0, 0,  1, true, true);
    const integralSph IR10b (-1,-1,  0, 0,  0, 0,  0, 0,   0, 1, 0, 0,  1, true, true); // NOT USED



    absms(I);

    if (IsZero(I)) return; // out of boundaries (counts as 0)



    if (iset.count(I)) return; // don't perform redundant work
    iset.insert(I); // register the integral

    // iff end of recursion, just insert  the  integral in the set
    if (I.ld==0) {
        fset.insert(I);
        return;
    }

    bool ABxy = vbool[0] || vbool[1];
    bool ABz  = vbool[2];

    bool CDxy = vbool[3] || vbool[4];
    bool CDz  = vbool[5];

    bool ACxy = vbool[6] || vbool[7];
    bool ACz  = vbool[8];


    if (CDz)  SKSd(I+IR1,  iset, fset, vbool);
    if (CDz)  SKSd(I+IR1b, iset, fset, vbool);
    if (CDxy) SKSd(I+IR2,  iset, fset, vbool);
    if (CDxy) SKSd(I+IR2b, iset, fset, vbool);
    SKSd(I+IR3, iset, fset, vbool);
    SKSd(I+IR4, iset, fset, vbool);

    if (ACz)  SKSd(I+IR5, iset, fset, vbool);
    if (ACxy) SKSd(I+IR6, iset, fset, vbool);
    SKSd(I+IR7, iset, fset, vbool);
    SKSd(I+IR8, iset, fset, vbool);

    if (CDz)  SKSd(I+IR9  , iset, fset, vbool);
    if (CDz)  SKSd(I+IR9b , iset, fset, vbool);
    if (CDxy) SKSd(I+IR10 , iset, fset, vbool);
    if (CDxy) SKSd(I+IR10b, iset, fset, vbool);
    SKSd(I+IR11, iset, fset, vbool);
    SKSd(I+IR12, iset, fset, vbool);

    if (ABz)  SKSd(I+IR13, iset, fset, vbool);
    if (ABxy) SKSd(I+IR14, iset, fset, vbool);
    SKSd(I+IR15, iset, fset, vbool);
    SKSd(I+IR16, iset, fset, vbool);
}


void SetVbool (GEOM geometry, bool vbool[9]) {

    bool ABx, ABy, ABz;
    bool CDx, CDy, CDz;
    bool ACx, ACy, ACz;

    switch(geometry) {
        case ABCD:
            ABx = ABy = CDx = false;
            ABz = CDy = CDz = true;
            ACx = ACy = ACz = true;
        break;

        case ABAD:
        /*
            ABx = ABy = CDx = false;
            ABz = CDy = CDz = true;
            ACx = ACy = ACz = false;
            */
            ABx = ABy = CDx = false;
            ABz = CDy = CDz = true;
            ACx = ACy = ACz = true;
        break;

        case AACD:
            ABx = ABy = ABz = CDx = CDy = false;
            CDz = true;
            ACx = false;
            ACy = ACz = true;
        break;

        case AACC:
            ABx = ABy = ABz = CDx = CDy = CDz = false;
            ACx = ACy = false;
            ACz = true;
        break;

        case ABAB:
            ABx = ABy = CDx = CDy = false;
            ABz = CDz = true;
            ACx = ACy = false;
            ACz = true;
        break;

        case AAAD:
            ABx = ABy = ABz = CDx = CDy = false;
            CDz = true;
            ACx = ACy = ACz = false;
        break;

        case AAAA:
            ABx = ABy = ABz = CDx = CDy = CDz = false;
            ACx = ACy = ACz = false;
        break;


        default:
            ABx = ABy = ABz = CDx = CDy = CDz = false;
            ACx = ACy = ACz = false;
        break;
    }


    vbool[0] = ABx;
    vbool[1] = ABy;
    vbool[2] = ABz;
    vbool[3] = CDx;
    vbool[4] = CDy;
    vbool[5] = CDz;
    vbool[6] = ACx;
    vbool[7] = ACy;
    vbool[8] = ACz;
}

void TestSKS (set<ikernel> & SphKernelSet, GEOM geometry, int La, int Lb, int Lc, int Ld) {
    set<integralSph> SphEriSet;
    set<integralSph> IES; // to store intermediate ERIs temporally

    set<integralSph>::iterator it;


    for (int ma=0; ma<=La; ++ma) {
        for (int mb=0; mb<=Lb; ++mb) {
            for (int mc=0; mc<=Lc; ++mc) {
                for (int md=0; md<=Ld; ++md) {
                    integralSph IS(La,ma, Lb,mb, Lc,mc, Ld,md);

                    SphEriSet.insert(IS);
                }
            }
        }
    }


    bool vbool[9];

    SetVbool(geometry, vbool);

    /*
     = {    false, false, true,
                        false,  true,  true,
                        true,   true,   true};
    */

    set<integralSph> SphEriSet1, SphEriSet2, SphEriSet3, SphEriSet4;

    //original set; (la+1)^4
    //cout << SphEriSet.size() << " ";

    // exhaust RR over SphHar on A
    for (it=SphEriSet.begin(); it!=SphEriSet.end(); ++it) SKSa(*it, IES, SphEriSet1, vbool);

    //cout << SphEriSet1.size() << " ";

    IES.clear(); // not that it matters, really

    for (it=SphEriSet1.begin(); it!=SphEriSet1.end(); ++it) SKSb(*it, IES, SphEriSet2, vbool);

    //cout << SphEriSet2.size() << " ";

    IES.clear();

    for (it=SphEriSet2.begin(); it!=SphEriSet2.end(); ++it) SKSc(*it, IES, SphEriSet3, vbool);

    //cout << SphEriSet3.size() << " ";

    IES.clear();

    for (it=SphEriSet3.begin(); it!=SphEriSet3.end(); ++it) SKSd(*it, IES, SphEriSet4, vbool);

    //should contain only kernels
    //cout << SphEriSet4.size() << endl;


    // transform everything to kernels
    for (it=SphEriSet4.begin(); it!=SphEriSet4.end(); ++it) {
        int u = it->b;
        int v = it->p;
        int s = it->d;
        int t = it->q;
        int m = it->m;

        ikernel K(u,v,s,t, m);

        SphKernelSet.insert(K);
    }
}

