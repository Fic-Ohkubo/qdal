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


#include "VILIC.hpp"
#include "IICinit.hpp"


/*
enum RRTYPE2 {BOYS,
                ADRR4, AERR4, CDR4, K4D, K4E, K4F,
                ADRR3, AERR3, CDR3, K3D, K3E, K3F,
                ADRR2, AERR2, CDR2, K2D, K2E, K2F,
                ADRR1, AERR1, CDR1, K1D, K1E, K1F,
                ADRR0, AERR0, CDR0,
                NADA};  // 29 < 32
*/


/*
    KERNELS  dest op1
    MMDZ     dest op1 op2 op3 op4 aux
    MMDY           x
    MMDX           x   x
    CTEBZ    dest op1 op2 op3 aux
    CTEKZ
    CTEKY
    CTEKX              x
    CTEBY              x
    CTEBX              x
    HRR      dest op1 op2
    HRR0     dest op1
    SPH      dest op1 op2 op3 op4 op5 op6 aux
    REORDER  dest op1
*/


void VILcode::Init (const UI32 ninstrK4[32], const Op_K4 * const nseqK4[32]) {
    size   = 0;
    ninstr = 0;

    for (UI8 iblock=BOYS; iblock<=NADA; ++iblock) {
        //offsets[iblock] = size;

        int ss = 0;

        if (iblock==AERR4)
            ss = 2;
        else if (iblock==AERR3  || iblock==AERR2  || iblock==AERR1  || iblock==AERR0)
            ss = 3;
        else if (iblock==CDR4)
            ss = 4;
        else if (iblock==CDR2 || iblock==CDR3)
            ss = 6;
        else if (iblock==CDR0 || iblock==CDR1)
            ss = 9;
        else if (iblock==K4E || iblock==K4F || iblock==K3E || iblock==K3F || iblock==K2E || iblock==K2F || iblock==K1E || iblock==K1F)
            ss = 3;

        size   += ninstrK4[iblock] * ss;
        ninstr += ninstrK4[iblock];

        nbi[iblock] = ninstrK4[iblock];
    }


    code = new UI16[size];

    int p = 0;

    for (UI8 iblock=BOYS; iblock<=NADA; ++iblock) {

        if (iblock==AERR4) {

            VILsimple * pKK = (VILsimple*)(code+p);

            pAERR4 = pKK;

            for (int i=0; i<ninstrK4[iblock]; ++i) {
                pKK[i].dest = nseqK4[iblock][i].dest;
                pKK[i].op1  = nseqK4[iblock][i].op1;
                p += 2;
            }
        }

        else if (iblock==AERR3  || iblock==AERR2  || iblock==AERR1  || iblock==AERR0) {

            VILaerr * pKK = (VILaerr*)(code+p);

            if      (iblock==AERR3) pAERR3 = pKK;
            else if (iblock==AERR2) pAERR2 = pKK;
            else if (iblock==AERR1) pAERR1 = pKK;
            else if (iblock==AERR0) pAERR0 = pKK;

            for (int i=0; i<ninstrK4[iblock]; ++i) {
                pKK[i].dest = nseqK4[iblock][i].dest;
                pKK[i].op1  = nseqK4[iblock][i].op1;
                pKK[i].op2  = nseqK4[iblock][i].op2;
                p += 3;
            }

        }

        else if (iblock==CDR4) {

            VILcdr4 * pKK = (VILcdr4*)(code+p);

            pCDR4 = pKK;

            for (int i=0; i<ninstrK4[iblock]; ++i) {
                pKK[i].dest = nseqK4[iblock][i].dest;
                pKK[i].op1  = nseqK4[iblock][i].op1;
                pKK[i].ope  = nseqK4[iblock][i].ope;
                pKK[i].aux  = nseqK4[iblock][i].aux;
                p += 4;
            }
        }

        else if (iblock==CDR2 || iblock==CDR3) {

            VILcdr23 * pKK = (VILcdr23*)(code+p);

            if      (iblock==CDR3) pCDR3 = pKK;
            else if (iblock==CDR2) pCDR2 = pKK;

            for (int i=0; i<ninstrK4[iblock]; ++i) {
                pKK[i].dest = nseqK4[iblock][i].dest;
                pKK[i].op1  = nseqK4[iblock][i].op1;
                pKK[i].op2  = nseqK4[iblock][i].op2;
                pKK[i].op3  = nseqK4[iblock][i].op3;
                pKK[i].ope  = nseqK4[iblock][i].ope;
                pKK[i].aux  = nseqK4[iblock][i].aux;
                p += 6;
            }
        }

        else if (iblock==CDR0 || iblock==CDR1) {

            VILcdr01 * pKK = (VILcdr01*)(code+p);

            if      (iblock==CDR1) pCDR1 = pKK;
            else if (iblock==CDR0) pCDR0 = pKK;

            for (int i=0; i<ninstrK4[iblock]; ++i) {
                pKK[i].dest = nseqK4[iblock][i].dest;
                pKK[i].op1  = nseqK4[iblock][i].op1;
                pKK[i].op2  = nseqK4[iblock][i].op2;
                pKK[i].op3  = nseqK4[iblock][i].op3;
                pKK[i].op4  = nseqK4[iblock][i].op4;
                pKK[i].op5  = nseqK4[iblock][i].op5;
                pKK[i].op6  = nseqK4[iblock][i].op6;
                pKK[i].ope  = nseqK4[iblock][i].ope;
                pKK[i].aux  = nseqK4[iblock][i].aux;
                p += 9;
            }
        }

        else if (iblock==K4E || iblock==K4F || iblock==K3E || iblock==K3F || iblock==K2E || iblock==K2F || iblock==K1E || iblock==K1F) {

            VILfma * pKK = (VILfma*)(code+p); //new VILfma[ninstrK4[iblock]]; //

            if      (iblock==K4E) pK4E = pKK;
            else if (iblock==K3E) pK3E = pKK;
            else if (iblock==K2E) pK2E = pKK;
            else if (iblock==K1E) pK1E = pKK;

            if      (iblock==K4F) pK4F = pKK;
            else if (iblock==K3F) pK3F = pKK;
            else if (iblock==K2F) pK2F = pKK;
            else if (iblock==K1F) pK1F = pKK;

            for (int i=0; i<ninstrK4[iblock]; ++i) {
                pKK[i].dest = nseqK4[iblock][i].dest;
                pKK[i].op1  = nseqK4[iblock][i].op1;
                pKK[i].aux  = nseqK4[iblock][i].aux;
                p += 3;
            }
        }

    }

}

