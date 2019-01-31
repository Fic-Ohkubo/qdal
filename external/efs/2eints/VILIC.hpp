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



    VILIC: Variable Instruction Length Interpreted Code
*/

#ifndef __VILIC__
#define __VILIC__

#include "../defs.hpp"

/*

AERR4/MOVE      F0e[s2->dest] = F0e[s2->op1] * iKpq;


K4E,K4F,K3E ... dest, op1, aux

AERR3/AERR2/AERR1/AERR0        store(dest, load(op1) + load(op2)

CDR4          F0f[s2->dest] = (R2 * F0f[s2->op1] + F0e[s2->ope]) * im2[s2->aux];
CDR3/CDR2         store(dest, (AQ2 * load(op1) - AQABb * load(op2) + AB2b2 * load(op3) + load(ope)) * im2[aux]);
CDR1/CDR0         store(dest, (AC2 * load(op1) - ACAB * load(op2) + ACCDd * load(op3) +
                         AB2 * load(op4) - ABCDd * load(op5) + CD2d2 * load(op6) +
                        load(ope)) * im2[aux]);

*/


class Op_K4;


struct VILsimple {
    UI16 dest;
    UI16 op1;
};

struct VILaerr {
    UI16 dest;
    UI16 op1;
    UI16 op2;
};

struct VILcdr4 {
    UI16 dest;
    UI16 op1;
    UI16 ope;
    UI16 aux;
};

struct VILcdr23 {
    UI16 dest;
    UI16 op1;
    UI16 op2;
    UI16 op3;
    UI16 ope;
    UI16 aux;
};

struct VILcdr01 {
    UI16 dest;
    UI16 op1;
    UI16 op2;
    UI16 op3;
    UI16 op4;
    UI16 op5;
    UI16 op6;
    UI16 ope;
    UI16 aux;
};

struct VILfma {
    UI16 dest;
    UI16 op1;
    UI16 aux;
};


struct VILcode {

    UI16 * code;
    int    size;
    int    ninstr;

    UI32 nbi[32];


//BOYS,
//ADRR4, AERR4, CDR4, K4D, K4E, K4F,
    VILsimple * pAERR4;
    VILcdr4   * pCDR4;

    VILfma    * pK4E;
    VILfma    * pK4F;

//ADRR3, AERR3, CDR3, K3D, K3E, K3F,
    VILaerr   * pAERR3;
    VILcdr23  * pCDR3;

    VILfma    * pK3E;
    VILfma    * pK3F;

//ADRR2, AERR2, CDR2, K2D, K2E, K2F,
    VILaerr   * pAERR2;
    VILcdr23  * pCDR2;

    VILfma    * pK2E;
    VILfma    * pK2F;

//ADRR1, AERR1, CDR1, K1D, K1E, K1F,

    VILaerr   * pAERR1;
    VILcdr01  * pCDR1;

    VILfma    * pK1E;
    VILfma    * pK1F;

//ADRR0, AERR0, CDR0,

    VILaerr   * pAERR0;
    VILcdr01  * pCDR0;

    //VILsimple * pMOVE;

//NADA

    VILcode() {
        code = 0;
        size = 0;
        ninstr = 0;

        for (int i=0; i<32; ++i) nbi[i] = 0;
    }

    ~VILcode() {
        delete[] code;
    }

    void Init (const UI32 ninstrK4[32], const Op_K4 * const nseqK4[32]);
};




#endif
