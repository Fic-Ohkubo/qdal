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


#ifndef __TRANSFORM__
#define __TRANSFORM__

#include <cstdlib>
#include <string>
#include <fstream>

#include "../defs.hpp"
#include "../low/chrono.hpp"
#include "../libquimera/libquimera.hpp"
#include "VILIC.hpp"

using namespace LibQuimera;


class cacheline64;
class ERIgeometries64;
class ERIgeometry;
class AtomProd;
class ShellPairPrototype;
class ShellPair;
class p_ERIbuffer;
class PrimitiveSet;


enum RRTYPE  {KERNELS, MMDZ, CTEBZ, CTEKZ, MMDY, MMDX, CTEKY, CTEKX, CTEBY, CTEBX, HRRBZ, HRRBY, HRRBX, SPHA, SPHB, HRRKZ, HRRKY, HRRKX, SPHC, SPHD, REORDER, OTHER};

enum RRTYPE2 {BOYS,
                ADRR4, AERR4, CDR4, K4D, K4E, K4F,
                ADRR3, AERR3, CDR3, K3D, K3E, K3F,
                ADRR2, AERR2, CDR2, K2D, K2E, K2F,
                ADRR1, AERR1, CDR1, K1D, K1E, K1F,
                ADRR0, AERR0, CDR0,
                NADA};  // 29 < 32




//this is HUGE!!!
class Op_K4 {
  public:
    UI16 dest;
    UI16 op1;
    UI16 op2;
    UI16 op3;
    UI16 op4;
    UI16 op5;
    UI16 op6;
    UI16 ope;
    UI16 aux;

    //RRTYPE ETYPE;

    bool operator==(const Op_K4 & rhs) const {
        return (dest==rhs.dest && op1==rhs.op1 && op2==rhs.op2 && op3==rhs.op3 && op4==rhs.op4 && op5==rhs.op5 && op6==rhs.op6 && ope==rhs.ope && aux==rhs.aux);
    }

    bool operator!=(const Op_K4 & rhs) const {
        return !(*this==rhs);
    }

} __attribute__((aligned(32)));

class Op_MIRROR16 {
  public:
    UI16 dest;
    UI16 op1;
    UI16 op2;
    UI16 op3;
    UI16 op4;
    UI16 op5;
    UI16 op6;
    UI16 aux; //also ope
    //RRTYPE ETYPE;

    bool operator==(const Op_MIRROR16 & rhs) const {
        return (dest==rhs.dest && op1==rhs.op1 && op2==rhs.op2 && op3==rhs.op3 && op4==rhs.op4 && op5==rhs.op5 && op6==rhs.op6 && aux==rhs.aux);
    }

    bool operator!=(const Op_MIRROR16 & rhs) const {
        return !(*this==rhs);
    }

} __attribute__((aligned(16)));


class Op_MIRROR {
  public:
    UI32 dest;
    UI32 op1;
    UI32 op2;
    UI32 op3;
    UI32 op4;
    UI32 op5;
    UI32 op6;
    UI32 aux; //also ope
    //RRTYPE ETYPE;

    bool operator==(const Op_MIRROR & rhs) const {
        return (dest==rhs.dest && op1==rhs.op1 && op2==rhs.op2 && op3==rhs.op3 && op4==rhs.op4 && op5==rhs.op5 && op6==rhs.op6 && aux==rhs.aux);
    }

    bool operator!=(const Op_MIRROR & rhs) const {
        return !(*this==rhs);
    }

} __attribute__((aligned(32)));



class ERIroutine {
    friend class p_Qalgorithm;
    friend class p_Quimera;

  private:
    //state
    bool IsSet;
    bool IsInitialized;


    //type of ERI batch
    UI8 la;
    UI8 lb;
    UI8 lc;
    UI8 ld;
    UI8 Am;
    GEOM geometry;
    bool useCDR;
    bool useGC;

    //instructions/pointer to function
    UI32 ninstrK4[32];
    UI32 ninstr  [32];

    Op_K4 * nseqK4[32];
    Op_K4 * eseqK4;

    Op_MIRROR * nseq[32];
    Op_MIRROR * eseq;

    VILcode K4VILcode;
    VILcode M0VILcode;


    //copy the coefficients to a local, more compact representation
    double Ca[11][6];
    double Cb[11][6];
    double Cc[11][6];
    double Cd[11][6];

    float fCa[11][6];
    float fCb[11][6];
    float fCc[11][6];
    float fCd[11][6];

    //useful info
    UI32 MaxMem;
    UI32 K4Mem;
    UI32 nKernels;

    UI8 Lt;


    int memF0;
    int memF0e;
    int memF0f;
    int memK4J1e;
    int memK4J1f;
    int memK3J2e;
    int memK3J2f;
    int memK2J3e;
    int memK2J3f;
    int memK1J4e;
    int memK1J4f;



    UI16 maxV;
    UI16 maxU;
    UI16 maxT;
    UI16 maxS;



    //for benchmarking
    //****************
    UI64 NFLOPS;
    UI64 NFLOPSK;
    UI64 NK4;


    //FLOPs in each loop
    int nK4J0;
    int nK4J1;
    int nK3J1;
    int nK3J2;
    int nK2J2;
    int nK2J3;
    int nK1J3;
    int nK1J4;
    int nK0J4;

    int nK4J0e;
    int nK4J1e;
    int nK3J1e;
    int nK3J2e;
    int nK2J2e;
    int nK2J3e;
    int nK1J3e;
    int nK1J4e;
    int nK0J4e;

    int nK4J0f;
    int nK4J1f;
    int nK3J1f;
    int nK3J2f;
    int nK2J2f;
    int nK2J3f;
    int nK1J3f;
    int nK1J4f;
    int nK0J4f;

    void (*InnerContractionRoutine)  (cacheline64 * (&F0), const ERIgeometries64 & vars8, const PrimitiveSet & PSab, const ShellPairPrototype & ABp, double ikcd, double rcd, p_ERIbuffer & buffer, cacheline64 & AB2, cacheline64 & X2);


    //routines not accessed outside the class
    std::string IdString(bool fuse=false) const;

    void ReadICfromFile();
    void WriteIC2File() const;


    template <bool SAB, bool SCD> void InnerContractCDR_GC (const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, double ikcd, double rcd, p_ERIbuffer & buffer) const;
    template <bool SAB, bool SCD> void InnerContractCDR_GC (cacheline64 * (&F0), const ERIgeometries64 & vars8, const PrimitiveSet & PSab,
                                                            const ShellPairPrototype & ABp, double ikcd, double rcd, p_ERIbuffer & buffer,
                                                            cacheline64 & AB2, cacheline64 & X2) const;

    void                               InnerContractCDR_GC (cacheline32 * (&F0), const ERIgeometries32 & vars8, const PrimitiveSet & PSab,
                                                            const ShellPairPrototype & ABp, float ikcd, float rcd, p_ERIbuffer & buffer,
                                                            cacheline32 & AB2, cacheline32 & X2) const;


    void LoadK2C();

    void WriteInnerContraction();

    void WriteLoopsContraction2(std::ofstream & file, const std::string & D, const std::string & S, const Op_K4 * IICbegin, const Op_K4 * IICend, bool firstloop);

    void WriteLoopsContraction(std::ofstream & file, const std::string & D, const std::string & S,
                                int memOffset,
                                const std::string & W, int maxW,
                                const std::string & NN,
                                const std::string & JJ,
                                const Op_K4 * IICbegin,
                                const Op_K4 * IICend,
                                bool firstloop=false
                                );

  public:

    ERIroutine() {
        NFLOPSK = 0;
        NK4 = 0;

        IsSet = false;
        IsInitialized = false;
        useCDR = false;
    }

    ~ERIroutine() {
    }

    void Init();
    void Write();
    void Read();
    void Set(UI8 La, UI8 Lb, UI8 Lc, UI8 Ld, GEOM geom, bool cdr);

    //sets the memory buffers and returns the size used
    size_t SetMemoryBuffers(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * mem, p_ERIbuffer & buffer) const;


    //ONE-CENTER ROUTINES
    //*******************

    void TransformAAAA(const double    *  mem,                             double    * ERI , p_ERIbuffer & buffer) const; //only need one for each atom
    void ContractCDR_GC (const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, double * __restrict__ uv_m_st, p_ERIbuffer & buffer, bool OnlyJ = false) const;

    //PACKED DOUBLE ROUTINES
    //**********************

    //MIRROR transformations
    void TransformABCD(const cacheline64 *  mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, p_ERIbuffer & buffer) const;
    void TransformAACD(const cacheline64 *  mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, p_ERIbuffer & buffer) const;
    void TransformAACC(const cacheline64 *  mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, p_ERIbuffer & buffer) const;
    void TransformABAB(const cacheline64 *  mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, p_ERIbuffer & buffer) const;

    //K4/CDR general ontractions
    template <bool SAB, bool SCD> void ContractCDR_GC (const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * __restrict__ uv_m_st8, p_ERIbuffer & buffer) const;
    template <bool SAB, bool SCD> void CalcGammas     (const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, p_ERIbuffer & buffer, bool OnlyJ = false) const;



    //PACKED FLOAT ROUTINES
    //*********************

    //MIRROR transformations
    void TransformABCD(const cacheline32 *  mem, const ERIgeometries32 & vars16, cacheline32 * ERI16, p_ERIbuffer & buffer) const;
    void TransformAACD(const cacheline32 *  mem, const ERIgeometries32 & vars16, cacheline32 * ERI16, p_ERIbuffer & buffer) const;
    void TransformAACC(const cacheline32 *  mem, const ERIgeometries32 & vars16, cacheline32 * ERI16, p_ERIbuffer & buffer) const;
    void TransformABAB(const cacheline32 *  mem, const ERIgeometries32 & vars16, cacheline32 * ERI16, p_ERIbuffer & buffer) const;

    void TransformCUDA(const cacheline32 *  mem, const ERIgeometries32 * vars16, cacheline32 * ERIs, int NTILES) const;


    //K4/CDR general ontractions
    void ContractCDR_GC (const ERIgeometries32 & vars16, const ERITile32 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline32 * __restrict__ uv_m_st8, p_ERIbuffer & buffer) const;
    void CalcGammas     (const ERIgeometries32 & vars16, const ERITile32 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, p_ERIbuffer & buffer, bool OnlyJ = false) const;
};

#endif
