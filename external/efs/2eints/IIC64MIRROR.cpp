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



    This is an implementation of the K4+MIRROR algorithm to compute 2-electron integrals.
    If you use it in your research, please cite the original  papers:

    "An algorithm for the efficient evaluation of two-electron repulsion integrals over contracted Gaussian-type basis functions",
    Jaime Axel Rosal Sandberg, Zilvinas Rinkevicius, J. Chem. Phys. 137, 234105 (2012); http://dx.doi.org/10.1063/1.4769730

    "New recurrence relations for analytic evaluation of two-electron repulsion integrals over highly contracted gaussian-type orbitals",
    Jaime Axel Rosal Sandberg, Zilvinas Rinkevicius, In preparation
*/

#include <string.h>
#include "../defs.hpp"
#include "../libquimera/libquimera.hpp"
#include "../basis/SPprototype.hpp"
#include "../basis/shellpair.hpp"
#include "../integrals/rotations.hpp"
#include "../integrals/atomprod.hpp"
#include "../2eints/IICinit.hpp"
#include "../2eints/quimera.hpp"

static const UI8 PF = 8; //prefetch distance

static const int DPC = DOUBLES_PER_CACHE_LINE;
static const int FPC = FLOATS_PER_CACHE_LINE;


void ERIroutine::TransformABCD(const cacheline64 * __restrict__  mem, const ERIgeometries64 & vars8, cacheline64 * __restrict__ ERI8, p_ERIbuffer & buffer) const {
    const cacheline64 & ABz = vars8.ABz;
    const cacheline64 & CDy = vars8.CDy;
    const cacheline64 & CDz = vars8.CDz;

    const cacheline64 & ACx = vars8.ACx;
    const cacheline64 & ACy = vars8.ACy;
    const cacheline64 & ACz = vars8.ACz;


    register Op_MIRROR * s2 = eseq;

    cacheline64 * buffer64 = (cacheline64*)buffer.bufferMIRROR;

    memset((void*)buffer64, 0, sizeof(cacheline64)*MaxMem);


    //copy kernels to new memory
    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        //__builtin_prefetch(buffer64 + s2[PF].op1);
        //__builtin_prefetch(mem            + s2[PF].dest, 1, 1);

        buffer64[i] = mem[i];
    }


    for (; s2<nseq[MMDZ]; ++s2) {

        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int zzz = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  ACz[i] * op3[i] - ABz[i] * op1[i] + CDz[i] * op2[i] - z * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp = ACz * load(op3) - ABz * load(op1) + CDz * load(op2);

        if (zzz>0) {
            double z = zzz;
            tmp -= load(op4) * z;
        }

        buffer64[s2->dest] = tmp;
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double z = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = z * op1[i] + ABz[i] * op2[i] + op3[i];
        */
        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + ABz * load(op2) + load(op3));
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double z = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -z * op1[i] + CDz[i] * op2[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + CDz * load(op2) - load(op3));
    }

    for (; s2<nseq[MMDY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int yyy = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACy[i] * op3[i] + CDy[i] * op2[i] - y * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp = ACy * load(op3) + CDy * load(op2);

        if (yyy>0) {
            double y = yyy;
            tmp -= load(op4) * y;
        }

        buffer64[s2->dest] = tmp;
    }

    for (; s2<nseq[MMDX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int xxx = s2->aux;
        //double x = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACx[i] * op3[i] - x * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp = ACx * load(op3);

        if (xxx>0) {
            double x = xxx;
            tmp -= load(op4) * x;
        }

        buffer64[s2->dest] = tmp;
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double y = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = CDy[i] * op2[i] - y * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + CDy * load(op2) - load(op3));
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double x = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -x * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x - load(op3));
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double y = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  y * op1[i] + op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + load(op3));
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double x = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = x * op1[i] + op3[i];
        */
        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x + load(op3));
    }

    for (; s2<nseq[HRRBZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - ABz[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - ABz * load(op2));
    }

    for (; s2<nseq[HRRBY]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBX]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHA]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (la>1) {
            cacheline64 sum;
            if (Ca[m][0]!=0) sum  = buffer64[s2->op1] * Ca[m][0];
            if (Ca[m][1]!=0) sum += buffer64[s2->op2] * Ca[m][1];
            if (Ca[m][2]!=0) sum += buffer64[s2->op3] * Ca[m][2];
            if (Ca[m][3]!=0) sum += buffer64[s2->op4] * Ca[m][3];
            if (Ca[m][4]!=0) sum += buffer64[s2->op5] * Ca[m][4];
            if (Ca[m][5]!=0) sum += buffer64[s2->op6] * Ca[m][5];
            buffer64[s2->dest] = sum;
        }
        else
            buffer64[s2->dest] = buffer64[s2->op1];
    }

    for (; s2<nseq[SPHB]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (lb>1) {
            cacheline64 sum;
            if (Cb[m][0]!=0) sum  = buffer64[s2->op1] * Cb[m][0];
            if (Cb[m][1]!=0) sum += buffer64[s2->op2] * Cb[m][1];
            if (Cb[m][2]!=0) sum += buffer64[s2->op3] * Cb[m][2];
            if (Cb[m][3]!=0) sum += buffer64[s2->op4] * Cb[m][3];
            if (Cb[m][4]!=0) sum += buffer64[s2->op5] * Cb[m][4];
            if (Cb[m][5]!=0) sum += buffer64[s2->op6] * Cb[m][5];
            buffer64[s2->dest] = sum;
        }
        else
            buffer64[s2->dest] = buffer64[s2->op1];
    }

    for (; s2<nseq[HRRKZ]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDz[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDz * load(op2));
    }

    for (; s2<nseq[HRRKY]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDy[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDy * load(op2));
    }

    for (; s2<nseq[HRRKX]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHC]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (lc>1) {
            cacheline64 sum;
            if (Cc[m][0]!=0) sum  = buffer64[s2->op1] * Cc[m][0];
            if (Cc[m][1]!=0) sum += buffer64[s2->op2] * Cc[m][1];
            if (Cc[m][2]!=0) sum += buffer64[s2->op3] * Cc[m][2];
            if (Cc[m][3]!=0) sum += buffer64[s2->op4] * Cc[m][3];
            if (Cc[m][4]!=0) sum += buffer64[s2->op5] * Cc[m][4];
            if (Cc[m][5]!=0) sum += buffer64[s2->op6] * Cc[m][5];
            buffer64[s2->dest] = sum;
        }
        else
            buffer64[s2->dest] = buffer64[s2->op1];
    }

    for (; s2<nseq[SPHD]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (ld>1) {
            cacheline64 sum;
            if (Cd[m][0]!=0) sum  = buffer64[s2->op1] * Cd[m][0];
            if (Cd[m][1]!=0) sum += buffer64[s2->op2] * Cd[m][1];
            if (Cd[m][2]!=0) sum += buffer64[s2->op3] * Cd[m][2];
            if (Cd[m][3]!=0) sum += buffer64[s2->op4] * Cd[m][3];
            if (Cd[m][4]!=0) sum += buffer64[s2->op5] * Cd[m][4];
            if (Cd[m][5]!=0) sum += buffer64[s2->op6] * Cd[m][5];
            buffer64[s2->dest] = sum;
        }
        else
            buffer64[s2->dest] = buffer64[s2->op1];
    }

    for (; s2<nseq[REORDER]; ++s2) {
        double * dest = (double*)(ERI8 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(ERI8 + s2[PF].dest, 1, 1);

        ERI8[s2->dest] = buffer64[s2->op1];
    }
}

void ERIroutine::TransformAACD(const cacheline64 * __restrict__  mem, const ERIgeometries64 & vars8, cacheline64 * __restrict__ ERI8, p_ERIbuffer & buffer) const {

    const cacheline64 & CDz = vars8.CDz;

    const cacheline64 & ACy = vars8.ACy;
    const cacheline64 & ACz = vars8.ACz;


    register Op_MIRROR * s2 = eseq;

    cacheline64 * buffer64 = (cacheline64*)buffer.bufferMIRROR;

    memset((void*)buffer64, 0, sizeof(cacheline64)*MaxMem);

    //copy kernels to new memory
    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        //__builtin_prefetch(buffer64 + s2[PF].op1);
        //__builtin_prefetch(mem            + s2[PF].dest, 1, 1);

        buffer64[i] = mem[i];
    }


    for (; s2<nseq[MMDZ]; ++s2) {

        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int zzz = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  ACz[i] * op3[i] - ABz[i] * op1[i] + CDz[i] * op2[i] - z * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp = ACz * load(op3) + CDz * load(op2);

        if (zzz>0) {
            double z = zzz;
            tmp -= load(op4) * z;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double z = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = z * op1[i] + ABz[i] * op2[i] + op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + load(op3));
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double z = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -z * op1[i] + CDz[i] * op2[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + CDz * load(op2) - load(op3));
    }

    for (; s2<nseq[MMDY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int yyy = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACy[i] * op3[i] + CDy[i] * op2[i] - y * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp = ACy * load(op3);

        if (yyy>0) {
            double y = yyy;
            tmp -= load(op4) * y;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[MMDX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int xxx = s2->aux;
        //double x = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACx[i] * op3[i] - x * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp;
        tmp.set(0);

        if (xxx>0) {
            double x = xxx;
            tmp -= load(op4) * x;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double y = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = CDy[i] * op2[i] - y * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y - load(op3));
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double x = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -x * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x - load(op3));
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double y = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  y * op1[i] + op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + load(op3));
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double x = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = x * op1[i] + op3[i];
        */
        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x + load(op3));
    }

    for (; s2<nseq[HRRBZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - ABz[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBY]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBX]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHA]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 sum = buffer64[s2->op1];

        if (la>1) {
            sum *= Ca[m][0];

            if (Ca[m][1]!=0) sum += buffer64[s2->op2] * Ca[m][1];
            if (Ca[m][2]!=0) sum += buffer64[s2->op3] * Ca[m][2];
            if (Ca[m][3]!=0) sum += buffer64[s2->op4] * Ca[m][3];
            if (Ca[m][4]!=0) sum += buffer64[s2->op5] * Ca[m][4];
            if (Ca[m][5]!=0) sum += buffer64[s2->op6] * Ca[m][5];
        }

        buffer64[s2->dest] = sum;
    }

    for (; s2<nseq[SPHB]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 sum = buffer64[s2->op1];

        if (lb>1) {
            sum *= Cb[m][0];

            if (Cb[m][1]!=0) sum += buffer64[s2->op2] * Cb[m][1];
            if (Cb[m][2]!=0) sum += buffer64[s2->op3] * Cb[m][2];
            if (Cb[m][3]!=0) sum += buffer64[s2->op4] * Cb[m][3];
            if (Cb[m][4]!=0) sum += buffer64[s2->op5] * Cb[m][4];
            if (Cb[m][5]!=0) sum += buffer64[s2->op6] * Cb[m][5];
        }

        buffer64[s2->dest] = sum;
    }

    for (; s2<nseq[HRRKZ]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDz[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDz * load(op2));
    }

    for (; s2<nseq[HRRKY]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDy[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRKX]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHC]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 sum = buffer64[s2->op1];

        if (lc>1) {
            sum *= Cc[m][0];

            if (Cc[m][1]!=0) sum += buffer64[s2->op2] * Cc[m][1];
            if (Cc[m][2]!=0) sum += buffer64[s2->op3] * Cc[m][2];
            if (Cc[m][3]!=0) sum += buffer64[s2->op4] * Cc[m][3];
            if (Cc[m][4]!=0) sum += buffer64[s2->op5] * Cc[m][4];
            if (Cc[m][5]!=0) sum += buffer64[s2->op6] * Cc[m][5];
        }

        buffer64[s2->dest] = sum;
    }

    for (; s2<nseq[SPHD]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 sum = buffer64[s2->op1];

        if (ld>1) {
            sum *= Cd[m][0];

            if (Cd[m][1]!=0) sum += buffer64[s2->op2] * Cd[m][1];
            if (Cd[m][2]!=0) sum += buffer64[s2->op3] * Cd[m][2];
            if (Cd[m][3]!=0) sum += buffer64[s2->op4] * Cd[m][3];
            if (Cd[m][4]!=0) sum += buffer64[s2->op5] * Cd[m][4];
            if (Cd[m][5]!=0) sum += buffer64[s2->op6] * Cd[m][5];
        }

        buffer64[s2->dest] = sum;
    }

    for (; s2<nseq[REORDER]; ++s2) {
        double * dest = (double*)(ERI8 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(ERI8 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

}

void ERIroutine::TransformAACC(const cacheline64 * __restrict__  mem, const ERIgeometries64 & vars8, cacheline64 * __restrict__ ERI8, p_ERIbuffer & buffer) const {


    register Op_MIRROR * s2 = eseq;

    cacheline64 * buffer64 = (cacheline64*)buffer.bufferMIRROR;

    memset((void*)buffer64, 0, sizeof(cacheline64)*MaxMem);

    //copy kernels to new memory
    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        //__builtin_prefetch(buffer64 + s2[PF].op1);
        //__builtin_prefetch(mem            + s2[PF].dest, 1, 1);

        buffer64[i] = mem[i];
    }


    for (; s2<nseq[MMDZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int zzz = s2->aux;
        //double x = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp;
        tmp = buffer64[s2->op3] * vars8.ACz;

        if (zzz>0) {
            double z = zzz;
            tmp -= load(op4) * z;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double z = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = z * op1[i] + ABz[i] * op2[i] + op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + load(op3));
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double z = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -z * op1[i] + CDz[i] * op2[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z  - load(op3));
    }

    for (; s2<nseq[MMDY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int yyy = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp;
        tmp.set(0);

        if (yyy>0) {
            double y = yyy;
            tmp -= load(op4) * y;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[MMDX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int xxx = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp;
        tmp.set(0);

        if (xxx>0) {
            double x = xxx;
            tmp -= load(op4) * x;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double y = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = CDy[i] * op2[i] - y * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y - load(op3));
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double x = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -x * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x - load(op3));
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double y = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  y * op1[i] + op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + load(op3));
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double x = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = x * op1[i] + op3[i];
        */
        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x + load(op3));
    }

    for (; s2<nseq[HRRBZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - ABz[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBY]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBX]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHA]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 sum = buffer64[s2->op1];

        if (la>1) {
            sum *= Ca[m][0];

            if (Ca[m][1]!=0) sum += buffer64[s2->op2] * Ca[m][1];
            if (Ca[m][2]!=0) sum += buffer64[s2->op3] * Ca[m][2];
            if (Ca[m][3]!=0) sum += buffer64[s2->op4] * Ca[m][3];
            if (Ca[m][4]!=0) sum += buffer64[s2->op5] * Ca[m][4];
            if (Ca[m][5]!=0) sum += buffer64[s2->op6] * Ca[m][5];
        }

        buffer64[s2->dest] = sum;
    }

    for (; s2<nseq[SPHB]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 sum = buffer64[s2->op1];

        if (lb>1) {
            sum *= Cb[m][0];

            if (Cb[m][1]!=0) sum += buffer64[s2->op2] * Cb[m][1];
            if (Cb[m][2]!=0) sum += buffer64[s2->op3] * Cb[m][2];
            if (Cb[m][3]!=0) sum += buffer64[s2->op4] * Cb[m][3];
            if (Cb[m][4]!=0) sum += buffer64[s2->op5] * Cb[m][4];
            if (Cb[m][5]!=0) sum += buffer64[s2->op6] * Cb[m][5];
        }

        buffer64[s2->dest] = sum;
    }

    for (; s2<nseq[HRRKZ]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDz[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRKY]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDy[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRKX]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHC]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 sum = buffer64[s2->op1];

        if (lc>1) {
            sum *= Cc[m][0];

            if (Cc[m][1]!=0) sum += buffer64[s2->op2] * Cc[m][1];
            if (Cc[m][2]!=0) sum += buffer64[s2->op3] * Cc[m][2];
            if (Cc[m][3]!=0) sum += buffer64[s2->op4] * Cc[m][3];
            if (Cc[m][4]!=0) sum += buffer64[s2->op5] * Cc[m][4];
            if (Cc[m][5]!=0) sum += buffer64[s2->op6] * Cc[m][5];
        }

        buffer64[s2->dest] = sum;
    }

    for (; s2<nseq[SPHD]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 sum = buffer64[s2->op1];

        if (ld>1) {
            sum *= Cd[m][0];

            if (Cd[m][1]!=0) sum += buffer64[s2->op2] * Cd[m][1];
            if (Cd[m][2]!=0) sum += buffer64[s2->op3] * Cd[m][2];
            if (Cd[m][3]!=0) sum += buffer64[s2->op4] * Cd[m][3];
            if (Cd[m][4]!=0) sum += buffer64[s2->op5] * Cd[m][4];
            if (Cd[m][5]!=0) sum += buffer64[s2->op6] * Cd[m][5];
        }

        buffer64[s2->dest] = sum;
    }

    for (; s2<nseq[REORDER]; ++s2) {
        double * dest = (double*)(ERI8 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(ERI8 + s2[PF].dest, 1, 1);

        ERI8[s2->dest] = buffer64[s2->op1];
    }

}

void ERIroutine::TransformAAAA(const double * __restrict__ uv_m_st, double * W, p_ERIbuffer & buffer)  const {

    //fill memory needed for the transformations with 0s, since some variables referenced with no RR associated are 0
    for (int i=0; i<MaxMem; ++i) buffer.buffer[i] = 0;


    register Op_MIRROR * s2 = eseq;

    //copy kernels to new memory
    //for (int i=0; i<ninstr[KERNELS]; ++i)      buffer.buffer[i] = uv_m_st[i];

    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        buffer.buffer[i] = uv_m_st[i];
    }


    for (; s2<nseq[MMDZ]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op4  = s2->op4;
        double z = s2->aux;
        if (s2->aux!=0) buffer.buffer[dest] =  - z * buffer.buffer[op4];
        else            buffer.buffer[dest] = 0;
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op3  = s2->op3;
        double z = s2->aux;

        buffer.buffer[dest] = z * buffer.buffer[op1] + buffer.buffer[op3];
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op3  = s2->op3;
        double z = s2->aux;

        buffer.buffer[dest] = -z * buffer.buffer[op1] - buffer.buffer[op3];
    }

    for (; s2<nseq[MMDY]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op4  = s2->op4;
        double y = s2->aux;

        if (s2->aux!=0) buffer.buffer[dest] = - y * buffer.buffer[op4];
        else            buffer.buffer[dest] = 0;
    }

    for (; s2<nseq[MMDX]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op4  = s2->op4;
        double x = s2->aux;

        if (s2->aux!=0) buffer.buffer[dest] = - x * buffer.buffer[op4];
        else            buffer.buffer[dest] = 0;
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op3  = s2->op3;
        double y = s2->aux;

        buffer.buffer[dest] = - y * buffer.buffer[op1] - buffer.buffer[op3];
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op3  = s2->op3;
        double x = s2->aux;

        buffer.buffer[dest] = -x * buffer.buffer[op1] - buffer.buffer[op3];
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op3  = s2->op3;
        double y = s2->aux;

        buffer.buffer[dest] =  y * buffer.buffer[op1] + buffer.buffer[op3];
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op3  = s2->op3;
        double x = s2->aux;

        buffer.buffer[dest] = x * buffer.buffer[op1] + buffer.buffer[op3];
    }

    //this three should be NULL if simplified
    for (; s2<nseq[HRRBZ]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        buffer.buffer[dest] = buffer.buffer[op1];
    }
    for (; s2<nseq[HRRBY]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        buffer.buffer[dest] = buffer.buffer[op1];
    }
    for (; s2<nseq[HRRBX]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        buffer.buffer[dest] = buffer.buffer[op1];
    }



    for (; s2<nseq[SPHA]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op2  = s2->op2;
        UI32 op3  = s2->op3;
        UI32 op4  = s2->op4;
        UI32 op5  = s2->op5;
        UI32 op6  = s2->op6;
        UI16 m    = s2->aux;

        if (la>1) {
            double           sum  = buffer.buffer[op1] * Ca[m][0];
            if (Ca[m][1]!=0) sum += buffer.buffer[op2] * Ca[m][1];
            if (Ca[m][2]!=0) sum += buffer.buffer[op3] * Ca[m][2];
            if (Ca[m][3]!=0) sum += buffer.buffer[op4] * Ca[m][3];
            if (Ca[m][4]!=0) sum += buffer.buffer[op5] * Ca[m][4];
            if (Ca[m][5]!=0) sum += buffer.buffer[op6] * Ca[m][5];
            buffer.buffer[dest] = sum;
        }
        else {
            buffer.buffer[dest] = buffer.buffer[op1];
        }
    }

    for (; s2<nseq[SPHB]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op2  = s2->op2;
        UI32 op3  = s2->op3;
        UI32 op4  = s2->op4;
        UI32 op5  = s2->op5;
        UI32 op6  = s2->op6;
        UI16 m    = s2->aux;

        if (lb>1) {
            double           sum  = buffer.buffer[op1] * Cb[m][0];
            if (Cb[m][1]!=0) sum += buffer.buffer[op2] * Cb[m][1];
            if (Cb[m][2]!=0) sum += buffer.buffer[op3] * Cb[m][2];
            if (Cb[m][3]!=0) sum += buffer.buffer[op4] * Cb[m][3];
            if (Cb[m][4]!=0) sum += buffer.buffer[op5] * Cb[m][4];
            if (Cb[m][5]!=0) sum += buffer.buffer[op6] * Cb[m][5];
            buffer.buffer[dest] = sum;
        }
        else {
            buffer.buffer[dest] = buffer.buffer[op1];
        }
    }

    //and again, this other three should be NULL if simplified
    for (; s2<nseq[HRRKZ]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        buffer.buffer[dest] = buffer.buffer[op1];
    }
    for (; s2<nseq[HRRKY]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        buffer.buffer[dest] = buffer.buffer[op1];
    }
    for (; s2<nseq[HRRKX]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        buffer.buffer[dest] = buffer.buffer[op1];
    }


    for (; s2<nseq[SPHC]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op2  = s2->op2;
        UI32 op3  = s2->op3;
        UI32 op4  = s2->op4;
        UI32 op5  = s2->op5;
        UI32 op6  = s2->op6;
        UI16 m    = s2->aux;

        if (lc>1) {
            double           sum  = buffer.buffer[op1] * Cc[m][0];
            if (Cc[m][1]!=0) sum += buffer.buffer[op2] * Cc[m][1];
            if (Cc[m][2]!=0) sum += buffer.buffer[op3] * Cc[m][2];
            if (Cc[m][3]!=0) sum += buffer.buffer[op4] * Cc[m][3];
            if (Cc[m][4]!=0) sum += buffer.buffer[op5] * Cc[m][4];
            if (Cc[m][5]!=0) sum += buffer.buffer[op6] * Cc[m][5];
            buffer.buffer[dest] = sum;
        }
        else {
            buffer.buffer[dest] = buffer.buffer[op1];
        }
    }

    for (; s2<nseq[SPHD]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;
        UI32 op2  = s2->op2;
        UI32 op3  = s2->op3;
        UI32 op4  = s2->op4;
        UI32 op5  = s2->op5;
        UI32 op6  = s2->op6;
        UI16 m    = s2->aux;

        if (ld>1) {
            double           sum  = buffer.buffer[op1] * Cd[m][0];
            if (Cd[m][1]!=0) sum += buffer.buffer[op2] * Cd[m][1];
            if (Cd[m][2]!=0) sum += buffer.buffer[op3] * Cd[m][2];
            if (Cd[m][3]!=0) sum += buffer.buffer[op4] * Cd[m][3];
            if (Cd[m][4]!=0) sum += buffer.buffer[op5] * Cd[m][4];
            if (Cd[m][5]!=0) sum += buffer.buffer[op6] * Cd[m][5];
            buffer.buffer[dest] = sum;
        }
        else {
            buffer.buffer[dest] = buffer.buffer[op1];
        }
    }


    for (; s2<nseq[REORDER]; ++s2) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        W[dest] = buffer.buffer[op1];
    }

}


void ERIroutine::TransformABAB(const cacheline64 * __restrict__  mem, const ERIgeometries64 & vars8, cacheline64 * __restrict__ ERI8, p_ERIbuffer & buffer) const {
    const cacheline64 & ABz = vars8.ABz;
    const cacheline64 & CDy = vars8.CDy;
    const cacheline64 & CDz = vars8.CDz;

    const cacheline64 & ACx = vars8.ACx;
    const cacheline64 & ACy = vars8.ACy;
    const cacheline64 & ACz = vars8.ACz;


    register Op_MIRROR * s2 = eseq;

    cacheline64 * buffer64 = (cacheline64*)buffer.bufferMIRROR;

    memset((void*)buffer64, 0, sizeof(cacheline64)*MaxMem);


    //copy kernels to new memory
    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        //__builtin_prefetch(buffer64 + s2[PF].op1);
        //__builtin_prefetch(mem            + s2[PF].dest, 1, 1);

        buffer64[i] = mem[i];
    }


    for (; s2<nseq[MMDZ]; ++s2) {

        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int zzz = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  ACz[i] * op3[i] - ABz[i] * op1[i] + CDz[i] * op2[i] - z * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        cacheline64 tmp = ACz * load(op3) - ABz * load(op1) + CDz * load(op2);

        if (zzz>0) {
            double z = zzz;
            tmp -= load(op4) * z;
        }

        buffer64[s2->dest] = tmp;
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double z = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = z * op1[i] + ABz[i] * op2[i] + op3[i];
        */
        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + ABz * load(op2) + load(op3));
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double z = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -z * op1[i] + CDz[i] * op2[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + CDz * load(op2) - load(op3));
    }

    for (; s2<nseq[MMDY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int yyy = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACy[i] * op3[i] + CDy[i] * op2[i] - y * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (yyy>0) {
            double y = yyy;
            buffer64[s2->dest] = load(op4)*(-y);
        }
        else buffer64[s2->dest].set(0);
    }

    for (; s2<nseq[MMDX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op4  = (double*)(buffer64 + s2->op4);
        unsigned short int xxx = s2->aux;
        //double x = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACx[i] * op3[i] - x * op4[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op4);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (xxx>0) {
            double x = xxx;
            buffer64[s2->dest] = load(op4) * (-x);
        }
        else buffer64[s2->dest].set(0);
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double y = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = CDy[i] * op2[i] - y * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y - load(op3));
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double x = -double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -x * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x - load(op3));
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double y = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  y * op1[i] + op3[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + load(op3));
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op3  = (double*)(buffer64 + s2->op3);
        double x = double(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = x * op1[i] + op3[i];
        */
        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op3);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x + load(op3));
    }

    for (; s2<nseq[HRRBZ]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - ABz[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - ABz * load(op2));
    }

    for (; s2<nseq[HRRBY]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBX]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHA]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (la>1) {
            cacheline64 sum;
            if (Ca[m][0]!=0) sum  = buffer64[s2->op1] * Ca[m][0];
            if (Ca[m][1]!=0) sum += buffer64[s2->op2] * Ca[m][1];
            if (Ca[m][2]!=0) sum += buffer64[s2->op3] * Ca[m][2];
            if (Ca[m][3]!=0) sum += buffer64[s2->op4] * Ca[m][3];
            if (Ca[m][4]!=0) sum += buffer64[s2->op5] * Ca[m][4];
            if (Ca[m][5]!=0) sum += buffer64[s2->op6] * Ca[m][5];
            buffer64[s2->dest] = sum;
        }
        else
            buffer64[s2->dest] = buffer64[s2->op1];
    }

    for (; s2<nseq[SPHB]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (lb>1) {
            cacheline64 sum;
            if (Cb[m][0]!=0) sum  = buffer64[s2->op1] * Cb[m][0];
            if (Cb[m][1]!=0) sum += buffer64[s2->op2] * Cb[m][1];
            if (Cb[m][2]!=0) sum += buffer64[s2->op3] * Cb[m][2];
            if (Cb[m][3]!=0) sum += buffer64[s2->op4] * Cb[m][3];
            if (Cb[m][4]!=0) sum += buffer64[s2->op5] * Cb[m][4];
            if (Cb[m][5]!=0) sum += buffer64[s2->op6] * Cb[m][5];
            buffer64[s2->dest] = sum;
        }
        else
            buffer64[s2->dest] = buffer64[s2->op1];
    }

    for (; s2<nseq[HRRKZ]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDz[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].op2);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDz * load(op2));
    }

    for (; s2<nseq[HRRKY]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDy[i] * op2[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRKX]; ++s2) {
        const double * op1  = (double*)(buffer64 + s2->op1);
        double * dest = (double*)(buffer64 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHC]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (lc>1) {
            cacheline64 sum;
            if (Cc[m][0]!=0) sum  = buffer64[s2->op1] * Cc[m][0];
            if (Cc[m][1]!=0) sum += buffer64[s2->op2] * Cc[m][1];
            if (Cc[m][2]!=0) sum += buffer64[s2->op3] * Cc[m][2];
            if (Cc[m][3]!=0) sum += buffer64[s2->op4] * Cc[m][3];
            if (Cc[m][4]!=0) sum += buffer64[s2->op5] * Cc[m][4];
            if (Cc[m][5]!=0) sum += buffer64[s2->op6] * Cc[m][5];
            buffer64[s2->dest] = sum;
        }
        else
            buffer64[s2->dest] = buffer64[s2->op1];
    }

    for (; s2<nseq[SPHD]; ++s2) {
        double * dest = (double*)(buffer64 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);
        const double * op2  = (double*)(buffer64 + s2->op2);
        const double * op3  = (double*)(buffer64 + s2->op3);
        const double * op4  = (double*)(buffer64 + s2->op4);
        const double * op5  = (double*)(buffer64 + s2->op5);
        const double * op6  = (double*)(buffer64 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(buffer64 + s2[PF].dest, 1, 1);

        if (ld>1) {
            cacheline64 sum;
            if (Cd[m][0]!=0) sum  = buffer64[s2->op1] * Cd[m][0];
            if (Cd[m][1]!=0) sum += buffer64[s2->op2] * Cd[m][1];
            if (Cd[m][2]!=0) sum += buffer64[s2->op3] * Cd[m][2];
            if (Cd[m][3]!=0) sum += buffer64[s2->op4] * Cd[m][3];
            if (Cd[m][4]!=0) sum += buffer64[s2->op5] * Cd[m][4];
            if (Cd[m][5]!=0) sum += buffer64[s2->op6] * Cd[m][5];
            buffer64[s2->dest] = sum;
        }
        else
            buffer64[s2->dest] = buffer64[s2->op1];
    }

    for (; s2<nseq[REORDER]; ++s2) {
        double * dest = (double*)(ERI8 + s2->dest);
        const double * op1  = (double*)(buffer64 + s2->op1);

        __builtin_prefetch(buffer64 + s2[PF].op1);
        __builtin_prefetch(ERI8 + s2[PF].dest, 1, 1);

        ERI8[s2->dest] = buffer64[s2->op1];
    }
}
