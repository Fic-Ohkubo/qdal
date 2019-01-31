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


#include <pmmintrin.h>
#include <xmmintrin.h>
#include <iostream>
#include <set>
#include <string.h>


#include "../defs.hpp"
#include "../libquimera/libquimera.hpp"
#include "../basis/SPprototype.hpp"
#include "../basis/shellpair.hpp"
#include "../integrals/rotations.hpp"
#include "../integrals/atomprod.hpp"
#include "../2eints/IICinit.hpp"
#include "../2eints/quimera.hpp"
#include "../math/angular.hpp"
#include "../math/gamma.hpp"
using namespace LibIGamma;


#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

using namespace std;
using namespace LibQuimera;
using namespace LibAngular;



static const UI8 PF = 8; //prefetch distance

static const int DPC = DOUBLES_PER_CACHE_LINE;
static const int FPC = FLOATS_PER_CACHE_LINE;




// VECTOR FLOAT ROUTINES
// *********************

void ERIroutine::TransformABCD(const cacheline32 * __restrict__  mem, const ERIgeometries32 & vars16, cacheline32 * __restrict__ ERI16, p_ERIbuffer & buffer) const {
    const cacheline32 & ABz = vars16.ABz;
    const cacheline32 & CDy = vars16.CDy;
    const cacheline32 & CDz = vars16.CDz;

    const cacheline32 & ACx = vars16.ACx;
    const cacheline32 & ACy = vars16.ACy;
    const cacheline32 & ACz = vars16.ACz;


    register Op_MIRROR * s2 = eseq;

    cacheline32 * buffer32 = (cacheline32*)buffer.bufferMIRROR;


    //copy kernels to new memory
    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        //__builtin_prefetch(buffer32 + s2[PF].op1);
        //__builtin_prefetch(mem            + s2[PF].dest, 1, 1);

        buffer32[i] = mem[i];
    }

    for (; s2<nseq[MMDZ]; ++s2) {

        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int zzz = s2->aux;

        //for (int i=0; i<8; ++i)
        //dest[i] =  ACz[i] * op3[i] - ABz[i] * op1[i] + CDz[i] * op2[i] - z * op4[i];

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp = ACz * load(op3) - ABz * load(op1) + CDz * load(op2);

        if (zzz>0) {
            float z = zzz;
            tmp -= load(op4) * z;
        }

        buffer32[s2->dest] = tmp;
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float z = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = z * op1[i] + ABz[i] * op2[i] + op3[i];
        */
        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + ABz * load(op2) + load(op3));
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float z = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -z * op1[i] + CDz[i] * op2[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + CDz * load(op2) - load(op3));
    }

    for (; s2<nseq[MMDY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int yyy = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACy[i] * op3[i] + CDy[i] * op2[i] - y * op4[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp = ACy * load(op3) + CDy * load(op2);

        if (yyy>0) {
            float y = yyy;
            tmp -= load(op4) * y;
        }

        buffer32[s2->dest] = tmp;
    }

    for (; s2<nseq[MMDX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int xxx = s2->aux;
        //float x = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACx[i] * op3[i] - x * op4[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp = ACx * load(op3);

        if (xxx>0) {
            float x = xxx;
            tmp -= load(op4) * x;
        }

        buffer32[s2->dest] = tmp;
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float y = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = CDy[i] * op2[i] - y * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + CDy * load(op2) - load(op3));
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float x = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -x * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x - load(op3));
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float y = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  y * op1[i] + op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + load(op3));
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float x = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = x * op1[i] + op3[i];
        */
        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x + load(op3));
    }

    for (; s2<nseq[HRRBZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - ABz[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - ABz * load(op2));
    }

    for (; s2<nseq[HRRBY]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBX]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHA]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        if (la>1) {
            cacheline32 sum;
            if (fCa[m][0]!=0) sum  = buffer32[s2->op1] * fCa[m][0];
            if (fCa[m][1]!=0) sum += buffer32[s2->op2] * fCa[m][1];
            if (fCa[m][2]!=0) sum += buffer32[s2->op3] * fCa[m][2];
            if (fCa[m][3]!=0) sum += buffer32[s2->op4] * fCa[m][3];
            if (fCa[m][4]!=0) sum += buffer32[s2->op5] * fCa[m][4];
            if (fCa[m][5]!=0) sum += buffer32[s2->op6] * fCa[m][5];
            buffer32[s2->dest] = sum;
        }
        else
            buffer32[s2->dest] = buffer32[s2->op1];
    }

    for (; s2<nseq[SPHB]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        if (lb>1) {
            cacheline32 sum;
            if (fCb[m][0]!=0) sum  = buffer32[s2->op1] * fCb[m][0];
            if (fCb[m][1]!=0) sum += buffer32[s2->op2] * fCb[m][1];
            if (fCb[m][2]!=0) sum += buffer32[s2->op3] * fCb[m][2];
            if (fCb[m][3]!=0) sum += buffer32[s2->op4] * fCb[m][3];
            if (fCb[m][4]!=0) sum += buffer32[s2->op5] * fCb[m][4];
            if (fCb[m][5]!=0) sum += buffer32[s2->op6] * fCb[m][5];
            buffer32[s2->dest] = sum;
        }
        else
            buffer32[s2->dest] = buffer32[s2->op1];
    }

    for (; s2<nseq[HRRKZ]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDz[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDz * load(op2));
    }

    for (; s2<nseq[HRRKY]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDy[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDy * load(op2));
    }

    for (; s2<nseq[HRRKX]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHC]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        if (lc>1) {
            cacheline32 sum;
            if (fCc[m][0]!=0) sum  = buffer32[s2->op1] * fCc[m][0];
            if (fCc[m][1]!=0) sum += buffer32[s2->op2] * fCc[m][1];
            if (fCc[m][2]!=0) sum += buffer32[s2->op3] * fCc[m][2];
            if (fCc[m][3]!=0) sum += buffer32[s2->op4] * fCc[m][3];
            if (fCc[m][4]!=0) sum += buffer32[s2->op5] * fCc[m][4];
            if (fCc[m][5]!=0) sum += buffer32[s2->op6] * fCc[m][5];
            buffer32[s2->dest] = sum;
        }
        else
            buffer32[s2->dest] = buffer32[s2->op1];
    }

    for (; s2<nseq[SPHD]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        if (ld>1) {
            cacheline32 sum;
            if (fCd[m][0]!=0) sum  = buffer32[s2->op1] * fCd[m][0];
            if (fCd[m][1]!=0) sum += buffer32[s2->op2] * fCd[m][1];
            if (fCd[m][2]!=0) sum += buffer32[s2->op3] * fCd[m][2];
            if (fCd[m][3]!=0) sum += buffer32[s2->op4] * fCd[m][3];
            if (fCd[m][4]!=0) sum += buffer32[s2->op5] * fCd[m][4];
            if (fCd[m][5]!=0) sum += buffer32[s2->op6] * fCd[m][5];
            buffer32[s2->dest] = sum;
        }
        else
            buffer32[s2->dest] = buffer32[s2->op1];
    }

    for (; s2<nseq[REORDER]; ++s2) {
        float * dest = (float*)(ERI16 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(ERI16 + s2[PF].dest, 1, 1);

        ERI16[s2->dest] = buffer32[s2->op1];
    }
}

void ERIroutine::TransformAACD(const cacheline32 * __restrict__  mem, const ERIgeometries32 & vars16, cacheline32 * __restrict__ ERI16, p_ERIbuffer & buffer) const {

    const cacheline32 & CDz = vars16.CDz;

    const cacheline32 & ACy = vars16.ACy;
    const cacheline32 & ACz = vars16.ACz;


    register Op_MIRROR * s2 = eseq;

    cacheline32 * buffer32 = (cacheline32*)buffer.bufferMIRROR;

    //copy kernels to new memory
    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        //__builtin_prefetch(buffer32 + s2[PF].op1);
        //__builtin_prefetch(mem            + s2[PF].dest, 1, 1);

        buffer32[i] = mem[i];
    }


    for (; s2<nseq[MMDZ]; ++s2) {

        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int zzz = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  ACz[i] * op3[i] - ABz[i] * op1[i] + CDz[i] * op2[i] - z * op4[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp = ACz * load(op3) + CDz * load(op2);

        if (zzz>0) {
            float z = zzz;
            tmp -= load(op4) * z;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float z = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = z * op1[i] + ABz[i] * op2[i] + op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + load(op3));
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float z = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -z * op1[i] + CDz[i] * op2[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + CDz * load(op2) - load(op3));
    }

    for (; s2<nseq[MMDY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int yyy = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACy[i] * op3[i] + CDy[i] * op2[i] - y * op4[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp = ACy * load(op3);

        if (yyy>0) {
            float y = yyy;
            tmp -= load(op4) * y;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[MMDX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int xxx = s2->aux;
        //float x = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACx[i] * op3[i] - x * op4[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp;
        tmp.set(0.f);

        if (xxx>0) {
            float x = xxx;
            tmp -= load(op4) * x;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float y = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = CDy[i] * op2[i] - y * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y - load(op3));
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float x = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -x * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x - load(op3));
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float y = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  y * op1[i] + op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + load(op3));
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float x = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = x * op1[i] + op3[i];
        */
        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x + load(op3));
    }

    for (; s2<nseq[HRRBZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - ABz[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBY]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBX]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHA]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 sum = buffer32[s2->op1];

        if (la>1) {
            sum *= fCa[m][0];

            if (fCa[m][1]!=0) sum += buffer32[s2->op2] * fCa[m][1];
            if (fCa[m][2]!=0) sum += buffer32[s2->op3] * fCa[m][2];
            if (fCa[m][3]!=0) sum += buffer32[s2->op4] * fCa[m][3];
            if (fCa[m][4]!=0) sum += buffer32[s2->op5] * fCa[m][4];
            if (fCa[m][5]!=0) sum += buffer32[s2->op6] * fCa[m][5];
        }

        buffer32[s2->dest] = sum;
    }

    for (; s2<nseq[SPHB]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 sum = buffer32[s2->op1];

        if (lb>1) {
            sum *= fCb[m][0];

            if (fCb[m][1]!=0) sum += buffer32[s2->op2] * fCb[m][1];
            if (fCb[m][2]!=0) sum += buffer32[s2->op3] * fCb[m][2];
            if (fCb[m][3]!=0) sum += buffer32[s2->op4] * fCb[m][3];
            if (fCb[m][4]!=0) sum += buffer32[s2->op5] * fCb[m][4];
            if (fCb[m][5]!=0) sum += buffer32[s2->op6] * fCb[m][5];
        }

        buffer32[s2->dest] = sum;
    }

    for (; s2<nseq[HRRKZ]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDz[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDz * load(op2));
    }

    for (; s2<nseq[HRRKY]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDy[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRKX]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHC]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 sum = buffer32[s2->op1];

        if (lc>1) {
            sum *= fCc[m][0];

            if (fCc[m][1]!=0) sum += buffer32[s2->op2] * fCc[m][1];
            if (fCc[m][2]!=0) sum += buffer32[s2->op3] * fCc[m][2];
            if (fCc[m][3]!=0) sum += buffer32[s2->op4] * fCc[m][3];
            if (fCc[m][4]!=0) sum += buffer32[s2->op5] * fCc[m][4];
            if (fCc[m][5]!=0) sum += buffer32[s2->op6] * fCc[m][5];
        }

        buffer32[s2->dest] = sum;
    }

    for (; s2<nseq[SPHD]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 sum = buffer32[s2->op1];

        if (ld>1) {
            sum *= fCd[m][0];

            if (fCd[m][1]!=0) sum += buffer32[s2->op2] * fCd[m][1];
            if (fCd[m][2]!=0) sum += buffer32[s2->op3] * fCd[m][2];
            if (fCd[m][3]!=0) sum += buffer32[s2->op4] * fCd[m][3];
            if (fCd[m][4]!=0) sum += buffer32[s2->op5] * fCd[m][4];
            if (fCd[m][5]!=0) sum += buffer32[s2->op6] * fCd[m][5];
        }

        buffer32[s2->dest] = sum;
    }

    for (; s2<nseq[REORDER]; ++s2) {
        float * dest = (float*)(ERI16 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(ERI16 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

}

void ERIroutine::TransformAACC(const cacheline32 * __restrict__  mem, const ERIgeometries32 & vars16, cacheline32 * __restrict__ ERI16, p_ERIbuffer & buffer) const {


    register Op_MIRROR * s2 = eseq;

    cacheline32 * buffer32 = (cacheline32*)buffer.bufferMIRROR;

    //copy kernels to new memory
    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        //__builtin_prefetch(buffer32 + s2[PF].op1);
        //__builtin_prefetch(mem            + s2[PF].dest, 1, 1);

        buffer32[i] = mem[i];
    }


    for (; s2<nseq[MMDZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int zzz = s2->aux;
        //float x = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp;
        tmp = buffer32[s2->op3] * vars16.ACz;

        if (zzz>0) {
            float z = zzz;
            tmp -= load(op4) * z;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float z = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = z * op1[i] + ABz[i] * op2[i] + op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + load(op3));
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float z = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -z * op1[i] + CDz[i] * op2[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z  - load(op3));
    }

    for (; s2<nseq[MMDY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int yyy = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp;
        tmp.set(0.f);

        if (yyy>0) {
            float y = yyy;
            tmp -= load(op4) * y;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[MMDX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int xxx = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp;
        tmp.set(0.f);

        if (xxx>0) {
            float x = xxx;
            tmp -= load(op4) * x;
        }

        store(dest, tmp);
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float y = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = CDy[i] * op2[i] - y * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y - load(op3));
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float x = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -x * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x - load(op3));
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float y = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  y * op1[i] + op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + load(op3));
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float x = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = x * op1[i] + op3[i];
        */
        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x + load(op3));
    }

    for (; s2<nseq[HRRBZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - ABz[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBY]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBX]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHA]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 sum = buffer32[s2->op1];

        if (la>1) {
            sum *= fCa[m][0];

            if (fCa[m][1]!=0) sum += buffer32[s2->op2] * fCa[m][1];
            if (fCa[m][2]!=0) sum += buffer32[s2->op3] * fCa[m][2];
            if (fCa[m][3]!=0) sum += buffer32[s2->op4] * fCa[m][3];
            if (fCa[m][4]!=0) sum += buffer32[s2->op5] * fCa[m][4];
            if (fCa[m][5]!=0) sum += buffer32[s2->op6] * fCa[m][5];
        }

        buffer32[s2->dest] = sum;
    }

    for (; s2<nseq[SPHB]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 sum = buffer32[s2->op1];

        if (lb>1) {
            sum *= fCb[m][0];

            if (fCb[m][1]!=0) sum += buffer32[s2->op2] * fCb[m][1];
            if (fCb[m][2]!=0) sum += buffer32[s2->op3] * fCb[m][2];
            if (fCb[m][3]!=0) sum += buffer32[s2->op4] * fCb[m][3];
            if (fCb[m][4]!=0) sum += buffer32[s2->op5] * fCb[m][4];
            if (fCb[m][5]!=0) sum += buffer32[s2->op6] * fCb[m][5];
        }

        buffer32[s2->dest] = sum;
    }

    for (; s2<nseq[HRRKZ]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDz[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRKY]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDy[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRKX]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHC]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 sum = buffer32[s2->op1];

        if (lc>1) {
            sum *= fCc[m][0];

            if (fCc[m][1]!=0) sum += buffer32[s2->op2] * fCc[m][1];
            if (fCc[m][2]!=0) sum += buffer32[s2->op3] * fCc[m][2];
            if (fCc[m][3]!=0) sum += buffer32[s2->op4] * fCc[m][3];
            if (fCc[m][4]!=0) sum += buffer32[s2->op5] * fCc[m][4];
            if (fCc[m][5]!=0) sum += buffer32[s2->op6] * fCc[m][5];
        }

        buffer32[s2->dest] = sum;
    }

    for (; s2<nseq[SPHD]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 sum = buffer32[s2->op1];

        if (ld>1) {
            sum *= fCd[m][0];

            if (fCd[m][1]!=0) sum += buffer32[s2->op2] * fCd[m][1];
            if (fCd[m][2]!=0) sum += buffer32[s2->op3] * fCd[m][2];
            if (fCd[m][3]!=0) sum += buffer32[s2->op4] * fCd[m][3];
            if (fCd[m][4]!=0) sum += buffer32[s2->op5] * fCd[m][4];
            if (fCd[m][5]!=0) sum += buffer32[s2->op6] * fCd[m][5];
        }

        buffer32[s2->dest] = sum;
    }

    for (; s2<nseq[REORDER]; ++s2) {
        float * dest = (float*)(ERI16 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(ERI16 + s2[PF].dest, 1, 1);

        ERI16[s2->dest] = buffer32[s2->op1];
    }

}

void ERIroutine::TransformABAB(const cacheline32 * __restrict__  mem, const ERIgeometries32 & vars16, cacheline32 * __restrict__ ERI16, p_ERIbuffer & buffer) const {
    const cacheline32 & ABz = vars16.ABz;
    const cacheline32 & CDy = vars16.CDy;
    const cacheline32 & CDz = vars16.CDz;

    const cacheline32 & ACx = vars16.ACx;
    const cacheline32 & ACy = vars16.ACy;
    const cacheline32 & ACz = vars16.ACz;


    register Op_MIRROR * s2 = eseq;

    cacheline32 * buffer32 = (cacheline32*)buffer.bufferMIRROR;


    //copy kernels to new memory
    for (int i=0; s2<nseq[KERNELS]; ++s2,++i) {
        UI32 dest = s2->dest;
        UI32 op1  = s2->op1;

        //__builtin_prefetch(buffer32 + s2[PF].op1);
        //__builtin_prefetch(mem            + s2[PF].dest, 1, 1);

        buffer32[i] = mem[i];
    }


    for (; s2<nseq[MMDZ]; ++s2) {

        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int zzz = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  ACz[i] * op3[i] - ABz[i] * op1[i] + CDz[i] * op2[i] - z * op4[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp = ACz * load(op3) - ABz * load(op1) + CDz * load(op2);

        if (zzz>0) {
            float z = zzz;
            tmp -= load(op4) * z;
        }

        buffer32[s2->dest] = tmp;
    }

    for (; s2<nseq[CTEBZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float z = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = z * op1[i] + ABz[i] * op2[i] + op3[i];
        */
        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + ABz * load(op2) + load(op3));
    }

    for (; s2<nseq[CTEKZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float z = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -z * op1[i] + CDz[i] * op2[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * z + CDz * load(op2) - load(op3));
    }

    for (; s2<nseq[MMDY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int yyy = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACy[i] * op3[i] + CDy[i] * op2[i] - y * op4[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp = ACy * load(op3) + CDy * load(op2);

        if (yyy>0) {
            float y = yyy;
            tmp -= load(op4) * y;
        }

        buffer32[s2->dest] = tmp;
    }

    for (; s2<nseq[MMDX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        unsigned short int xxx = s2->aux;
        //float x = s2->aux;

        /*
        for (int i=0; i<8; ++i)
        dest[i] = ACx[i] * op3[i] - x * op4[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op4);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        cacheline32 tmp = ACx * load(op3);

        if (xxx>0) {
            float x = xxx;
            tmp -= load(op4) * x;
        }

        buffer32[s2->dest] = tmp;
    }

    for (; s2<nseq[CTEKY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float y = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = CDy[i] * op2[i] - y * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + CDy * load(op2) - load(op3));
    }

    for (; s2<nseq[CTEKX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float x = -float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = -x * op1[i] - op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x - load(op3));
    }

    for (; s2<nseq[CTEBY]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float y = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] =  y * op1[i] + op3[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * y + load(op3));
    }

    for (; s2<nseq[CTEBX]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op3  = (float*)(buffer32 + s2->op3);
        float x = float(s2->aux);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = x * op1[i] + op3[i];
        */
        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op3);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) * x + load(op3));
    }

    for (; s2<nseq[HRRBZ]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - ABz[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - ABz * load(op2));
    }

    for (; s2<nseq[HRRBY]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[HRRBX]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHA]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        if (la>1) {
            cacheline32 sum;
            if (fCa[m][0]!=0) sum  = buffer32[s2->op1] * fCa[m][0];
            if (fCa[m][1]!=0) sum += buffer32[s2->op2] * fCa[m][1];
            if (fCa[m][2]!=0) sum += buffer32[s2->op3] * fCa[m][2];
            if (fCa[m][3]!=0) sum += buffer32[s2->op4] * fCa[m][3];
            if (fCa[m][4]!=0) sum += buffer32[s2->op5] * fCa[m][4];
            if (fCa[m][5]!=0) sum += buffer32[s2->op6] * fCa[m][5];
            buffer32[s2->dest] = sum;
        }
        else
            buffer32[s2->dest] = buffer32[s2->op1];
    }

    for (; s2<nseq[SPHB]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        if (lb>1) {
            cacheline32 sum;
            if (fCb[m][0]!=0) sum  = buffer32[s2->op1] * fCb[m][0];
            if (fCb[m][1]!=0) sum += buffer32[s2->op2] * fCb[m][1];
            if (fCb[m][2]!=0) sum += buffer32[s2->op3] * fCb[m][2];
            if (fCb[m][3]!=0) sum += buffer32[s2->op4] * fCb[m][3];
            if (fCb[m][4]!=0) sum += buffer32[s2->op5] * fCb[m][4];
            if (fCb[m][5]!=0) sum += buffer32[s2->op6] * fCb[m][5];
            buffer32[s2->dest] = sum;
        }
        else
            buffer32[s2->dest] = buffer32[s2->op1];
    }

    for (; s2<nseq[HRRKZ]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDz[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDz * load(op2));
    }

    for (; s2<nseq[HRRKY]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i] - CDy[i] * op2[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].op2);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1) - CDy * load(op2));
    }

    for (; s2<nseq[HRRKX]; ++s2) {
        const float * op1  = (float*)(buffer32 + s2->op1);
        float * dest = (float*)(buffer32 + s2->dest);

        /*
        for (int i=0; i<8; ++i)
        dest[i] = op1[i];
        */

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        store(dest, load(op1));
    }

    for (; s2<nseq[SPHC]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        if (lc>1) {
            cacheline32 sum;
            if (fCc[m][0]!=0) sum  = buffer32[s2->op1] * fCc[m][0];
            if (fCc[m][1]!=0) sum += buffer32[s2->op2] * fCc[m][1];
            if (fCc[m][2]!=0) sum += buffer32[s2->op3] * fCc[m][2];
            if (fCc[m][3]!=0) sum += buffer32[s2->op4] * fCc[m][3];
            if (fCc[m][4]!=0) sum += buffer32[s2->op5] * fCc[m][4];
            if (fCc[m][5]!=0) sum += buffer32[s2->op6] * fCc[m][5];
            buffer32[s2->dest] = sum;
        }
        else
            buffer32[s2->dest] = buffer32[s2->op1];
    }

    for (; s2<nseq[SPHD]; ++s2) {
        float * dest = (float*)(buffer32 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);
        const float * op2  = (float*)(buffer32 + s2->op2);
        const float * op3  = (float*)(buffer32 + s2->op3);
        const float * op4  = (float*)(buffer32 + s2->op4);
        const float * op5  = (float*)(buffer32 + s2->op5);
        const float * op6  = (float*)(buffer32 + s2->op6);
        UI16 m    = s2->aux;

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(buffer32 + s2[PF].dest, 1, 1);

        if (ld>1) {
            cacheline32 sum;
            if (fCd[m][0]!=0) sum  = buffer32[s2->op1] * fCd[m][0];
            if (fCd[m][1]!=0) sum += buffer32[s2->op2] * fCd[m][1];
            if (fCd[m][2]!=0) sum += buffer32[s2->op3] * fCd[m][2];
            if (fCd[m][3]!=0) sum += buffer32[s2->op4] * fCd[m][3];
            if (fCd[m][4]!=0) sum += buffer32[s2->op5] * fCd[m][4];
            if (fCd[m][5]!=0) sum += buffer32[s2->op6] * fCd[m][5];
            buffer32[s2->dest] = sum;
        }
        else
            buffer32[s2->dest] = buffer32[s2->op1];
    }

    for (; s2<nseq[REORDER]; ++s2) {
        float * dest = (float*)(ERI16 + s2->dest);
        const float * op1  = (float*)(buffer32 + s2->op1);

        __builtin_prefetch(buffer32 + s2[PF].op1);
        __builtin_prefetch(ERI16 + s2[PF].dest, 1, 1);

        ERI16[s2->dest] = buffer32[s2->op1];
    }
}


static const float iii[] =
{0.5, 0.75, 1.875, 6.5625, 29.53125, 162.421875, 1055.7421875, 7918.06640625, 67303.564453125,
639383.8623046875, 6.7135305541992188e6, 7.7205601373291016e7,
 9.6507001716613770e8, 1.3028445231742859e10,
 1.8891245586027145e11, 2.9281430658342075e12,
 4.8314360586264424e13, 8.4550131025962743e14,
 1.5641774239803107e16, 3.0501459767616059e17,
 6.2527992523612922e18, 1.3443518392576778e20,
 3.0247916383297751e21, 7.1082603500749715e22,
 1.7415237857683680e24, 4.4408856537093384e25,
 1.1768346982329747e27, 3.2362954201406804e28,
 9.2234419474009391e29, 2.7209153744832770e31,
 8.2987918921739949e32, 2.6141194460348084e34};

static const float im2[32] = {
    1.,     1./3.,  1./5.,  1./7.,  1./9.,  1./11., 1./13., 1./15.,
    1./17., 1./19., 1./21., 1./23., 1./25., 1./27., 1./29., 1./31.,
    1./33., 1./35., 1./37., 1./39., 1./41., 1./43., 1./45., 1./47.,
    1./49., 1./51., 1./53., 1./55., 1./57., 1./59., 1./61., 1./63.,
};



template <int labcd> inline void GammaF(float & e, float & m, float KpqR2) {

    double ** gamma = IncompleteGammas[labcd].gamma_table;

    if (KpqR2+fvg_step>fvg_max) {
        float ir2 = 1.f/(KpqR2);
        m = fPI3 * sqrt(ir2);

        if (labcd>0) {
            float ir4  = ir2  * ir2;
            float ir8  = ir4  * ir4;
            float ir16 = ir8  * ir8;
            float ir32 = ir16 * ir16;
            float ir64 = ir32 * ir32;

            if (labcd&1)  m *= ir2;
            if (labcd&2)  m *= ir4;
            if (labcd&4)  m *= ir8;
            if (labcd&8)  m *= ir16;
            if (labcd&16) m *= ir32;
            if (labcd&32) m *= ir64;

            m *= iii[labcd-1];
        }

        e = 0;
    }
    else {
        float p = fivg_step * KpqR2;
        int pos = int(p+0.5f);
        float x0 = fvg_step*float(pos);

        float Ax1 = x0 - KpqR2;
        float Ax2 = 0.5f * Ax1*Ax1;
        float Ax3 = 0.33333333333333333333f * Ax1*Ax2;

        //m = PI54 * (gamma[pos][labcd+1] + gamma[pos][labcd+2] * Ax1 + gamma[pos][labcd+3] * Ax2 + gamma[pos][labcd+4] * Ax3);
        m = fPI54 * (float(gamma[pos][1]) + float(gamma[pos][2]) * Ax1 + float(gamma[pos][3]) * Ax2 + float(gamma[pos][4]) * Ax3);
        e = fPI54 * float(gamma[pos][0]) * (1.f+Ax1+Ax2+Ax3);
    }
}

template <int Labcd> void Gamma(cacheline32 * e, cacheline32 * m, float iKpq, const cacheline32 & R2, const cacheline32 & W) {

    float Kpq = 0.5f/iKpq;

    cacheline32 KpqR2;
    KpqR2 = R2 * Kpq;

    //only non-vector part of the algorithm
    for (int k=0; k<FPC; ++k) {
        float ed, md;

        //double ed, md;
        GammaF<Labcd>(ed, md, KpqR2(k));

        (*e)(k) = ed;
        (*m)(k) = md;
    }

    //this could be precomputed and stored
    float Kpqn = sqrt(Kpq); {
        Kpq += Kpq;

        float Kpq2  = Kpq   * Kpq;
        float Kpq4  = Kpq2  * Kpq2;
        float Kpq8  = Kpq4  * Kpq4;
        float Kpq16 = Kpq8  * Kpq8;
        float Kpq32 = Kpq16 * Kpq16;

        if (Labcd&1)  Kpqn *= Kpq;
        if (Labcd&2)  Kpqn *= Kpq2;
        if (Labcd&4)  Kpqn *= Kpq4;
        if (Labcd&8)  Kpqn *= Kpq8;
        if (Labcd&16) Kpqn *= Kpq16;
        if (Labcd&32) Kpqn *= Kpq32;
    }

    cacheline32 WW;
    WW = W * Kpqn;

    if (Labcd>0) *e *= WW;
    *m *= WW;

    //AER downward recursion
    //for (int i=Labcd-1; i>=0; --i) e[i] = e[i+1] * iKpq;

    //downward recursion
    //for (int i=Labcd-1; i>=0; --i) m[i] = (R2 * m[i+1] + e[i]) * (1./float(2*i+1));
}

void ERIroutine::CalcGammas     (const ERIgeometries32 & vars16, const ERITile32 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, p_ERIbuffer & buffer, bool OnlyJ) const {

    void (*cGamma)(cacheline32 * e, cacheline32 * m, float iKpq, const cacheline32 & R2, const cacheline32 & W) ;

    if (Lt== 0) cGamma = Gamma< 0>;
    if (Lt== 1) cGamma = Gamma< 1>;
    if (Lt== 2) cGamma = Gamma< 2>;
    if (Lt== 3) cGamma = Gamma< 3>;
    if (Lt== 4) cGamma = Gamma< 4>;
    if (Lt== 5) cGamma = Gamma< 5>;
    if (Lt== 6) cGamma = Gamma< 6>;
    if (Lt== 7) cGamma = Gamma< 7>;
    if (Lt== 8) cGamma = Gamma< 8>;
    if (Lt== 9) cGamma = Gamma< 9>;
    if (Lt==10) cGamma = Gamma<10>;
    if (Lt==11) cGamma = Gamma<11>;
    if (Lt==12) cGamma = Gamma<12>;
    if (Lt==13) cGamma = Gamma<13>;
    if (Lt==14) cGamma = Gamma<14>;
    if (Lt==15) cGamma = Gamma<15>;
    if (Lt==16) cGamma = Gamma<16>;
    if (Lt==17) cGamma = Gamma<17>;
    if (Lt==18) cGamma = Gamma<18>;
    if (Lt==19) cGamma = Gamma<19>;
    if (Lt==20) cGamma = Gamma<20>;


    const int JA = ABp.Ja;
    const int JB = ABp.Jb;
    const int JC = CDp.Ja;
    const int JD = CDp.Jb;

    const PrimitiveSet & PSab = ABp.Psets[ET.nKab];
    const PrimitiveSet & PScd = CDp.Psets[ET.nKcd];

    cacheline32 * KAB = (cacheline32*)buffer.KAB;
    cacheline32 * KCD = (cacheline32*)buffer.KCD;

    {
        int addr[FLOATS_PER_CACHE_LINE];
        for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) addr[i] = ET.ap12[i] *ET.nK2ab;

        float ws[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));

        int ab=0;
        for (int b=0; b<PSab.nKb; ++b) {
            for (int a=0; a<PSab.nKa[b]; ++a) {
                int ba = b*ABp.Ka + a;
                for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) ws[i] = ET.wSP12[addr[i]+ba];
                KAB[ab].set(ws);
                ++ab;
            }
        }
    }


    {
        int addr[FLOATS_PER_CACHE_LINE];
        for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) addr[i] = ET.ap34[i] *ET.nK2cd;

        float ws[FLOATS_PER_CACHE_LINE] __attribute__((aligned(CACHE_LINE_SIZE)));

        int cd = 0;
        for (int d=0; d<PScd.nKb; ++d) {
            for (int c=0; c<PScd.nKa[d]; ++c) {
                int dc = d*CDp.Ka + c;

                for (int i=0; i<FLOATS_PER_CACHE_LINE; ++i) ws[i] = ET.wSP34[addr[i]+dc];

                KCD[cd].set(ws);
                ++cd;
            }
        }
    }




    const cacheline32 & ABz = vars16.ABz;
    const cacheline32 & CDy = vars16.CDy;
    const cacheline32 & CDz = vars16.CDz;

    const cacheline32 & ACx = vars16.ACx;
    const cacheline32 & ACy = vars16.ACy;
    const cacheline32 & ACz = vars16.ACz;


    cacheline32 AC2;  AC2  = ACx*ACx;


    cacheline32 * F0 = (cacheline32*)buffer.F0;

    int dc = 0;
    for (int d=0; d<PScd.nKb; ++d) {
        for (int c=0; c<PScd.nKa[d]; ++c) {

            float rcd  = CDp.BCP.r[d][c];

            cacheline32 AQy; AQy = ACy + CDy * rcd;
            cacheline32 AQz; AQz = ACz + CDz * rcd;

            cacheline32 AQ2; AQ2 = AC2 + AQy*AQy;


            int ba = 0;
            for (int b=0; b<PSab.nKb; ++b) {
                for (int a=0; a<PSab.nKa[b]; ++a) {

                    float rab  = ABp.BCP.r[b][a];

                    cacheline32 ww;

                    ww.set(ABp.BCP.K[b][a]*CDp.BCP.K[d][c]);


                    cacheline32 PQz;
                    cacheline32 PQ2;

                    PQz = AQz - ABz*rab;
                    PQ2 = AQ2 + PQz*PQz; //gPQ2 = gAQ2 - gAQAB*rab + gAB2*(rab*rab);
                    ww *=  KAB[ba] * KCD[dc];

                    float iKpq  = (ABp.BCP.k[b][a] + CDp.BCP.k[d][c]);
                    cGamma (F0, F0+1, iKpq, PQ2, ww);

                    F0 += 2;

                    ++ba;
                }
            }

            ++dc;
        }
    }

}

inline void LoopsContraction2(cacheline32 * __restrict__ D, const cacheline32 * __restrict__ S, const float * __restrict__ Vs,
                                const VILfma * __restrict__ IIC, int ninstr) {

    for (int i=0; i<ninstr; ++i) {
        cacheline32       * dest = D + IIC[i].dest;
        const cacheline32 * op1  = S + IIC[i].op1;
        UI16 aux                 =     IIC[i].aux;

        //__builtin_prefetch(D + s2[PF].dest, 1, 1);

        store(dest, load(dest) + load(op1) * Vs[aux]);
    }

}

void ERIroutine::InnerContractCDR_GC(cacheline32 * (&F0), const ERIgeometries32 & vars8, const PrimitiveSet & PSab,
                                      const ShellPairPrototype & ABp, float ikcd, float rcd, p_ERIbuffer & buffer,
                                      cacheline32 & AB2, cacheline32 & X2) const {

    const bool SAB = true;
    const bool SCD = true;

    const cacheline32 & ABz = vars8.ABz;
    const cacheline32 & CDy = vars8.CDy;
    const cacheline32 & CDz = vars8.CDz;

    const cacheline32 & ACx = vars8.ACx;
    const cacheline32 & ACy = vars8.ACy;
    const cacheline32 & ACz = vars8.ACz;

    const int JA = ABp.Ja;
    const int JB = ABp.Jb;

    const int nJ1 = JA;
    const int nJ2 = JA*JB;

    const int mK3f = nJ2*memK3J2f;
    const int mK4f = nJ1*memK4J1f;

    const int mK3e = nJ2*memK3J2e;
    const int mK4e = nJ1*memK4J1e;

    {
        memset(buffer.K3J2e, 0, mK3e*sizeof(cacheline32));
        memset(buffer.K3J2f, 0, mK3f*sizeof(cacheline32));

        cacheline32 AQy;
        cacheline32 AQz;

        cacheline32 X2Y2;
        cacheline32 AQ2;

        cacheline32 AQAB;


        if (SAB && SCD) {
            AQy  = ACy + CDy * rcd;
            AQz  = ACz + CDz * rcd;

            X2Y2 = X2   + AQy*AQy;
            AQ2  = X2Y2 + AQz*AQz;

            AQAB = (AQz*ABz); AQAB += AQAB;
        }
        else if (!SAB && SCD) {
            AQz  = ACz + CDz * rcd;
            AQ2  = X2 + AQz*AQz;
        }

        for (int b=0; b<PSab.nKb; ++b) {
            memset(buffer.K4J1e, 0, mK4e*sizeof(cacheline32));
            memset(buffer.K4J1f, 0, mK4f*sizeof(cacheline32));

            for (int a=0; a<PSab.nKa[b]; ++a) {

                float ikab = ABp.BCP.k[b][a];

                float iKpq  = (ikab + ikcd);

                cacheline32 PQz;
                cacheline32 R2;

                if (SAB && SCD) {
                    float rab  = ABp.BCP.r[b][a];
                    PQz = AQz - ABz*rab;
                    R2 = X2Y2 + PQz*PQz;
                }

                cacheline32 * F0e = (cacheline32*)buffer.F0e;
                cacheline32 * F0f = (cacheline32*)buffer.F0f;

                //copy the gamma of higher L and the contracted exponential
                F0e[0] = F0[0];
                F0f[0] = F0[1];

                //AERR K4 (uncontracted)
                {

                    for (int i=1; i<K4VILcode.nbi[AERR4]; ++i) {
                        cacheline32       * dest = F0e + K4VILcode.pAERR4[i].dest;
                        const cacheline32 * op1  = F0e + K4VILcode.pAERR4[i].op1;

                        *dest = *op1 * iKpq;
                    }

                }

                //CDR K4 (uncontracted downward recursion)
                {

                    for (int i=1; i<K4VILcode.nbi[CDR4]; ++i) {
                        cacheline32       * dest = F0f + K4VILcode.pCDR4[i].dest;
                        const cacheline32 * op1  = F0f + K4VILcode.pCDR4[i].op1;
                        const cacheline32 * ope  = F0e + K4VILcode.pCDR4[i].ope;
                        UI16 aux                 = K4VILcode.pCDR4[i].aux;

                        if      ( SAB &&  SCD) *dest = ( R2 * *op1 + *ope) * im2[aux];
                        else if (!SAB &&  SCD) *dest = (AQ2 * *op1 + *ope) * im2[aux];
                        else if (!SAB && !SCD) *dest = ( X2 * *op1 + *ope) * im2[aux];
                    }

                }

                //K4 E + F
                {
                    const cacheline32 * J0e = (cacheline32*)F0e;
                    cacheline32 * J1e = (cacheline32*)buffer.K4J1e;

                    const cacheline32 * J0f = (cacheline32*)F0f;
                    cacheline32 * J1f = (cacheline32*)buffer.K4J1f;

                    float Vv[maxJ][32] __attribute__((aligned(CACHE_LINE_SIZE)));

                    for (int j=0; j<JA; ++j) {
                        Vv[j][0] = ABp.BCP.Na[a][j];

                        for (int v=1; v<=maxV; ++v)
                            Vv[j][v] = Vv[j][v-1] * ABp.BCP.k[b][a];
                    }

                    for (int j=0; j<min(JA, a+1); ++j) {
                        LoopsContraction2(J1e+j*memK4J1e, J0e, Vv[j], K4VILcode.pK4E, K4VILcode.nbi[K4E]);
                        LoopsContraction2(J1f+j*memK4J1f, J0f, Vv[j], K4VILcode.pK4F, K4VILcode.nbi[K4F]);
                    }
                }

                F0 += 2;
            }

            //K3J1 //AERR K3
            {
                cacheline32 * Jf = (cacheline32*)buffer.K4J1f;
                cacheline32 * Je = (cacheline32*)buffer.K4J1e;

                //__m128 ikcd_v = _mm_load1_ps(&ikcd);
                //u_atom ikcd_v(ikcd);

                for (int nj1=0; nj1<nJ1; ++nj1) {
                    for (int i=0; i<K4VILcode.nbi[AERR3]; ++i) {
                        cacheline32       * dest = Je + K4VILcode.pAERR3[i].dest;
                        const cacheline32 * op1  = Je + K4VILcode.pAERR3[i].op1;
                        const cacheline32 * op2  = Je + K4VILcode.pAERR3[i].op2;

                        store(dest, load(op1) + load(op2) * ikcd);
                    }

                    Je += memK4J1e;
                }
            }

            //CDR3
            {

                cacheline32 AQABb;
                cacheline32 AB2b2;

                if (SAB && SCD) {
                    AQABb = AQAB * ABp.BCP.b1[b];
                    AB2b2 = AB2  * ABp.BCP.b2[b];
                }

                cacheline32 * Jf = (cacheline32*)buffer.K4J1f;
                cacheline32 * Je = (cacheline32*)buffer.K4J1e;

                for (int nj1=0; nj1<nJ1; ++nj1) {
                    for (int i=0; i<K4VILcode.nbi[CDR3]; ++i) {
                        cacheline32       * dest = Jf + K4VILcode.pCDR3[i].dest;
                        const cacheline32 * op1  = Jf + K4VILcode.pCDR3[i].op1;
                        const cacheline32 * op2  = Jf + K4VILcode.pCDR3[i].op2;
                        const cacheline32 * op3  = Jf + K4VILcode.pCDR3[i].op3;
                        const cacheline32 * ope  = Je + K4VILcode.pCDR3[i].ope;
                        UI16 aux  = K4VILcode.pCDR3[i].aux;

                        if (SAB && SCD) {
                            store(dest, (AQ2 * load(op1) - AQABb * load(op2) + AB2b2 * load(op3) + load(ope)) * im2[aux]);
                        }
                        else if (!SAB && SCD) {
                            store(dest, (AQ2 * load(op1)                                         + load(ope)) * im2[aux]);
                        }
                        else if (!SAB && !SCD) {
                            store(dest, (X2 * load(op1)                                         + load(ope)) * im2[aux]);
                        }
                    }


                    Jf += memK4J1f;
                    Je += memK4J1e;
                }
            }

            //K3 E + F
            {
                const cacheline32 * J1e = (cacheline32*)buffer.K4J1e;
                cacheline32 * J2e = (cacheline32*)buffer.K3J2e;

                const cacheline32 * J1f = (cacheline32*)buffer.K4J1f;
                cacheline32 * J2f = (cacheline32*)buffer.K3J2f;

                float Uu[maxJ][32] __attribute__((aligned(CACHE_LINE_SIZE)));

                for (int j=0; j<JB; ++j) {
                    Uu[j][0] = ABp.BCP.Nb[b][j];

                    for (int u=1; u<=maxU; ++u)
                        Uu[j][u] = Uu[j][u-1] * ABp.BCP.b1[b];
                }


                for (int nj1=0; nj1<nJ1; ++nj1) {
                    for (int j=0; j<min(JB, b+1); ++j) {
                        LoopsContraction2(J2e+j*memK3J2e, J1e, Uu[j], K4VILcode.pK3E, K4VILcode.nbi[K3E]);
                        LoopsContraction2(J2f+j*memK3J2f, J1f, Uu[j], K4VILcode.pK3F, K4VILcode.nbi[K3F]);
                    }

                    J2e += JB*memK3J2e;
                    J1e += memK4J1e;

                    J2f += JB*memK3J2f;
                    J1f += memK4J1f;
                }
            }

        }

        //AERR K2
        {
            cacheline32 * J2 = (cacheline32*)buffer.K3J2e;

            //__m128 ikcd_v = _mm_load1_ps(&ikcd);
            //u_atom ikcd_v(ikcd);

            //K2J2
            for (int nj2=0; nj2<nJ2; ++nj2) {
                for (int i=0; i<K4VILcode.nbi[AERR2]; ++i) {
                    cacheline32       * dest = J2 + K4VILcode.pAERR2[i].dest;
                    const cacheline32 * op1  = J2 + K4VILcode.pAERR2[i].op1;
                    const cacheline32 * op2  = J2 + K4VILcode.pAERR2[i].op2;

                    store(dest, load(op1) + load(op2) * ikcd);
                }

                J2 += memK3J2e;
            }
        }

        //CDR K2
        {
            cacheline32 * J2 = (cacheline32*)buffer.K3J2f;
            cacheline32 * Je = (cacheline32*)buffer.K3J2e;

            for (int nj2=0; nj2<nJ2; ++nj2) {
                for (int i=0; i<K4VILcode.nbi[CDR2]; ++i) {
                    cacheline32       * dest = J2 + K4VILcode.pCDR2[i].dest;
                    const cacheline32 * op1  = J2 + K4VILcode.pCDR2[i].op1;
                    const cacheline32 * op2  = J2 + K4VILcode.pCDR2[i].op2;
                    const cacheline32 * op3  = J2 + K4VILcode.pCDR2[i].op3;
                    const cacheline32 * ope  = Je + K4VILcode.pCDR2[i].ope;
                    UI16 aux  = K4VILcode.pCDR2[i].aux;

                    if (SAB && SCD) {
                        store(dest, (AQ2 * load(op1) - AQAB * load(op2) + AB2 * load(op3) + load(ope)) * im2[aux]);
                    }
                    else if (!SAB && SCD) {
                        store(dest, (AQ2 * load(op1)                                      + load(ope)) * im2[aux]);
                    }
                    else if (!SAB && !SCD) {
                        store(dest, (X2 * load(op1)                                      + load(ope)) * im2[aux]);
                    }
                }

                J2 += memK3J2f;
                Je += memK3J2e;
            }
        }

    }

}

void ERIroutine::ContractCDR_GC(const ERIgeometries32 & vars8, const ERITile32 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline32 * __restrict__ uv_m_st8, p_ERIbuffer & buffer) const {

    const bool SAB = true;
    const bool SCD = true;

    const PrimitiveSet & PSab = ABp.Psets[ET.nKab];
    const PrimitiveSet & PScd = CDp.Psets[ET.nKcd];

    const cacheline32 & ABz = vars8.ABz;
    const cacheline32 & CDy = vars8.CDy;
    const cacheline32 & CDz = vars8.CDz;

    const cacheline32 & ACx = vars8.ACx;
    const cacheline32 & ACy = vars8.ACy;
    const cacheline32 & ACz = vars8.ACz;

    const int JA = ABp.Ja;
    const int JB = ABp.Jb;
    const int JC = CDp.Ja;
    const int JD = CDp.Jb;

    const int nJ1 = JA;
    const int nJ2 = JA*JB;
    const int nJ3 = JA*JB*JC;
    const int nJ4 = JA*JB*JC*JD;


    const int mK1f = nJ4*memK1J4f;
    const int mK2f = nJ3*memK2J3f;
    const int mK3f = nJ2*memK3J2f;
    const int mK4f = nJ1*memK4J1f;

    const int mK1e = nJ4*memK1J4e;
    const int mK2e = nJ3*memK2J3e;
    const int mK3e = nJ2*memK3J2e;
    const int mK4e = nJ1*memK4J1e;


    cacheline32 CD2;
    cacheline32 ACCD;
    cacheline32 ABCD;

    cacheline32 AC2;
    cacheline32 ACAB;
    cacheline32 AB2;

    cacheline32 X2;

    cacheline32 X2Y2;


    if (SAB && SCD) {
        CD2  =        CDy*CDy + CDz*CDz;
        ACCD =        ACy*CDy + ACz*CDz; ACCD += ACCD;
        ABCD =                  ABz*CDz; ABCD += ABCD;

        X2   = ACx*ACx;

        AC2  =      X2+ACy*ACy+ACz*ACz;
        ACAB =                 ABz*ACz; ACAB += ACAB;
        AB2  =                 ABz*ABz;
    }
    else if (!SAB && SCD) {
        CD2  =                  CDz*CDz;
        ACCD =                  ACz*CDz; ACCD += ACCD;
        X2   = ACx*ACx;
        X2Y2 = X2 + ACy*ACy;
        AC2  =      X2Y2+ACz*ACz;
    }
    else if (!SAB && !SCD) {
        AC2  =      ACx*ACx+ACy*ACy+ACz*ACz;
    }



    cacheline32 * F0 = (cacheline32*)buffer.F0;


    memset(buffer.K1J4e, 0, mK1e*sizeof(cacheline32));
    memset(buffer.K1J4f, 0, mK1f*sizeof(cacheline32));

    for (int d=0; d<PScd.nKb; ++d) {
        memset(buffer.K2J3e, 0, mK2e*sizeof(cacheline32));
        memset(buffer.K2J3f, 0, mK2f*sizeof(cacheline32));

        for (int c=0; c<PScd.nKa[d]; ++c) {

            float ikcd  = CDp.BCP.k[d][c];
            float rcd   = CDp.BCP.r[d][c];

            /*

            if      (InnerContractionRoutine!=NULL &&  SAB &&  SCD)
                InnerContractionRoutine      ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2);
            else if (InnerContractionRoutine!=NULL && !SAB &&  SCD)
                InnerContractionRoutine      ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2Y2);
            else if (InnerContractionRoutine!=NULL && !SAB &&  !SCD)
                InnerContractionRoutine      ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, AC2);

            else

            if (SAB && SCD)
                InnerContractCDR_GC<SAB,SCD> ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2);
            else if (!SAB && SCD)
                InnerContractCDR_GC<SAB,SCD> ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2Y2);
            else if (!SAB && !SCD)
                InnerContractCDR_GC<SAB,SCD> ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, AC2);
            */

            InnerContractCDR_GC ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2);

            {
                const cacheline32 * J2e = (cacheline32*)buffer.K3J2e;
                cacheline32 * J3e = (cacheline32*)buffer.K2J3e;

                const cacheline32 * J2f = (cacheline32*)buffer.K3J2f;
                cacheline32 * J3f = (cacheline32*)buffer.K2J3f;

                float Tt[maxJ][32] __attribute__((aligned(CACHE_LINE_SIZE)));

                for (int j=0; j<JC; ++j) {
                    Tt[j][0] = CDp.BCP.Na[c][j];

                    for (int t=1; t<=maxT; ++t)
                        Tt[j][t] = Tt[j][t-1] * CDp.BCP.k[d][c];
                }

                //K2J3e+f contraction
                for (int nj2=0; nj2<nJ2; ++nj2) {
                    for (int j=0; j<min(JC, c+1); ++j) {
                        LoopsContraction2(J3e+j*memK2J3e, J2e, Tt[j], K4VILcode.pK2E, K4VILcode.nbi[K2E]);
                        LoopsContraction2(J3f+j*memK2J3f, J2f, Tt[j], K4VILcode.pK2F, K4VILcode.nbi[K2F]);
                    }

                    J3e += JC*memK2J3e;
                    J2e += memK3J2e;

                    J3f += JC*memK2J3f;
                    J2f += memK3J2f;
                }

            }


        }

        //AERR K1
        {
            cacheline32 * J3 = (cacheline32*)buffer.K2J3e;

            //K1J3
            for (int nj3=0; nj3<nJ3; ++nj3) {
                for (int i=0; i<K4VILcode.nbi[AERR1]; ++i) {
                    cacheline32       * dest = J3 + K4VILcode.pAERR1[i].dest;
                    const cacheline32 * op1  = J3 + K4VILcode.pAERR1[i].op1;
                    const cacheline32 * op2  = J3 + K4VILcode.pAERR1[i].op2;

                    store(dest, load(op1) + load(op2));
                }


                J3 += memK2J3e;
            }
        }

        cacheline32 CD2d2;
        cacheline32 ACCDd;
        cacheline32 ABCDd;

        if (SAB && SCD) {
            CD2d2 = CD2  * CDp.BCP.b2[d];
            ACCDd = ACCD * CDp.BCP.b1[d];
            ABCDd = ABCD * CDp.BCP.b1[d];
        }
        else if (!SAB && SCD) {
            CD2d2 = CD2  * CDp.BCP.b2[d];
            ACCDd = ACCD * CDp.BCP.b1[d];
        }

        //CDR K1
        {
            cacheline32 * J3 = (cacheline32*)buffer.K2J3f;
            cacheline32 * Je = (cacheline32*)buffer.K2J3e;

            for (int nj3=0; nj3<nJ3; ++nj3) {
                for (int i=0; i<K4VILcode.nbi[CDR1]; ++i) {
                    cacheline32       * dest = J3 + K4VILcode.pCDR1[i].dest;
                    const cacheline32 * op1  = J3 + K4VILcode.pCDR1[i].op1;
                    const cacheline32 * op2  = J3 + K4VILcode.pCDR1[i].op2;
                    const cacheline32 * op3  = J3 + K4VILcode.pCDR1[i].op3;
                    const cacheline32 * op4  = J3 + K4VILcode.pCDR1[i].op4;
                    const cacheline32 * op5  = J3 + K4VILcode.pCDR1[i].op5;
                    const cacheline32 * op6  = J3 + K4VILcode.pCDR1[i].op6;
                    const cacheline32 * ope  = Je + K4VILcode.pCDR1[i].ope;
                    UI16 aux  = K4VILcode.pCDR1[i].aux;

                    if (SAB && SCD) {
                        store(dest, (AC2 * load(op1) - ACAB * load(op2) + ACCDd * load(op3) +
                         AB2 * load(op4) - ABCDd * load(op5) + CD2d2 * load(op6) +
                        load(ope)) * im2[aux]);
                    }
                    else if (!SAB && SCD) {
                        store(dest, (AC2 * load(op1)                    + ACCDd * load(op3) +
                                                               CD2d2 * load(op6) +
                        load(ope)) * im2[aux]);
                    }
                    else if (!SAB && !SCD) {
                        store(dest, (AC2 * load(op1) +
                        load(ope)) * im2[aux]);
                    }
                }

                J3 += memK2J3f;
                Je += memK2J3e;
            }

        }

        //K1J4 E + F
        {
            const cacheline32 * J3e = (cacheline32*)buffer.K2J3e;
            cacheline32 * J4e = (cacheline32*)buffer.K1J4e;

            const cacheline32 * J3f = (cacheline32*)buffer.K2J3f;
            cacheline32 * J4f = (cacheline32*)buffer.K1J4f;

            float Ss[maxJ][32] __attribute__((aligned(CACHE_LINE_SIZE)));

            for (int j=0; j<JD; ++j) {
                Ss[j][0] = CDp.BCP.Nb[d][j];

                for (int s=1; s<=maxS; ++s)
                    Ss[j][s] = Ss[j][s-1] * CDp.BCP.b1[d];
            }


            for (int nj3=0; nj3<nJ3; ++nj3) {

                for (int j=0; j<min(JD, d+1); ++j) {
                    LoopsContraction2(J4e+j*memK1J4e, J3e, Ss[j], K4VILcode.pK1E, K4VILcode.nbi[K1E]);
                    LoopsContraction2(J4f+j*memK1J4f, J3f, Ss[j], K4VILcode.pK1F, K4VILcode.nbi[K1F]);
                }

                J4e += JD*memK1J4e;
                J3e += memK2J3e;

                J4f += JD*memK1J4f;
                J3f += memK2J3f;
            }
        }
    }


    //AERR K0
    {
        cacheline32 * J4 = (cacheline32*)buffer.K1J4e;

        //K0J4
        for (int nj4=0; nj4<nJ4; ++nj4) {
            for (int i=0; i<K4VILcode.nbi[AERR0]; ++i) {
                cacheline32       * dest = J4 + K4VILcode.pAERR0[i].dest;
                const cacheline32 * op1  = J4 + K4VILcode.pAERR0[i].op1;
                const cacheline32 * op2  = J4 + K4VILcode.pAERR0[i].op2;

                store(dest, load(op1) + load(op2));
            }

            J4 += memK1J4e;
        }
    }

    //CDR K0
    {
        cacheline32 * J4 = (cacheline32*)buffer.K1J4f;
        cacheline32 * Je = (cacheline32*)buffer.K1J4e;

        for (int nj4=0; nj4<nJ4; ++nj4) {
            for (int i=0; i<K4VILcode.nbi[CDR0]; ++i) {
                cacheline32       * dest = J4 + K4VILcode.pCDR0[i].dest;
                const cacheline32 * op1  = J4 + K4VILcode.pCDR0[i].op1;
                const cacheline32 * op2  = J4 + K4VILcode.pCDR0[i].op2;
                const cacheline32 * op3  = J4 + K4VILcode.pCDR0[i].op3;
                const cacheline32 * op4  = J4 + K4VILcode.pCDR0[i].op4;
                const cacheline32 * op5  = J4 + K4VILcode.pCDR0[i].op5;
                const cacheline32 * op6  = J4 + K4VILcode.pCDR0[i].op6;
                const cacheline32 * ope  = Je + K4VILcode.pCDR0[i].ope;
                UI16 aux  = K4VILcode.pCDR0[i].aux;

                if (SAB && SCD) {
                    store(dest, (AC2 * load(op1) - ACAB * load(op2) + ACCD * load(op3) +
                     AB2 * load(op4) - ABCD * load(op5) + CD2 * load(op6) +
                    load(ope)) * im2[aux]);
                }
                else if (!SAB && SCD) {
                    store(dest, (AC2 * load(op1)                    + ACCD * load(op3) +
                                                          CD2 * load(op6) +
                    load(ope)) * im2[aux]);
                }
                else if (!SAB && !SCD) {
                    store(dest, (AC2 * load(op1) +
                    load(ope)) * im2[aux]);
                }
            }


            J4 += memK1J4f;
            Je += memK1J4e;
        }
    }

    //copy to the kernel buffer
    {
        const cacheline32 * J4 = (cacheline32*)buffer.K1J4f;
        cacheline32 * I = uv_m_st8;

        for (int nj4=0; nj4<nJ4; ++nj4) {
            register Op_MIRROR * s2 = eseq;

            for (int i=0; s2<nseq[KERNELS]; ++i,++s2) {
                cacheline32       * dest = I + s2->dest;
                const cacheline32 * op1 = J4 + s2->op1;

                __builtin_prefetch(J4 + s2[PF].op1);
                __builtin_prefetch(I  + s2[PF].dest, 1, 1);

                if (s2->aux == 0) store(dest, load(op1));
                else              memset(dest, 0, sizeof(cacheline32)); // zero the kernel
            }

            J4 += memK1J4f;
            I += nKernels;
        }
    }


}



