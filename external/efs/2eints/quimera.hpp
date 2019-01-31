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


#ifndef __QUIMERA__
#define __QUIMERA__

#include <map>
#include <queue>

#include "../defs.hpp"
#include "../libquimera/libquimera.hpp"
#include "../2eints/IICinit.hpp"
using namespace LibQuimera;


class ShellPairPrototype;
class cacheline64;
class cacheline32;
class cacheline;

class ERIgeometries64;
class ShellPair;
class ERIroutine;
class ERItype;
class K4benchmark;

class p_Qalgorithm;
class p_ERIbuffer;
class p_Quimera;

typedef void (*K4routine)  (p_ERIbuffer &, const ERIgeometries64  &, const ERITile64 &, const ShellPairPrototype &, const ShellPairPrototype &, cacheline64 *);
typedef void (*K4routineF) (p_ERIbuffer &, const ERIgeometries32 &, const ERITile32 &, const ShellPairPrototype &, const ShellPairPrototype &, cacheline32 *);


struct ERItype {
    //basic information
    GEOM geometry;
    UI8  La;
    UI8  Lb;
    UI8  Lc;
    UI8  Ld;

    bool isGC;  //general contraction
    bool isLR;  //long range integrals
    bool isCDR; //CDR in K4

    ERItype & operator()(GEOM geo, UI8 la, UI8 lb, UI8 lc, UI8 ld, bool gc, bool lr, bool cdr) {
        geometry = geo;
        La = la;
        Lb = lb;
        Lc = lc;
        Ld = ld;
        isGC  = gc;
        isLR  = lr;
        isCDR = cdr;

        return *this;
    }

    //to enumerate them
    inline bool operator<(const ERItype & rhs) const {
        if (geometry!=rhs.geometry) return geometry<rhs.geometry;
        if (La      !=rhs.La      ) return La      <rhs.La;
        if (Lb      !=rhs.Lb      ) return Lb      <rhs.Lb;
        if (Lc      !=rhs.Lc      ) return Lc      <rhs.Lc;
        if (Ld      !=rhs.Ld      ) return Ld      <rhs.Ld;

        if (isGC  != rhs.isGC)  return isGC;
        if (isCDR != rhs.isCDR) return isCDR;
        if (isLR  != rhs.isLR)  return isLR;

        return false;
    }
};


//some class should have performance data
//also, variable precision
struct K4benchmark {

    ERItype type;


    double deltaK4;
    UI64 ncallsK4;
    UI64 tflopsK4;

    double deltaMIRROR;
    UI64 ncallsMIRROR;
    UI64 tflopsMIRROR;


    //average Ks
    UI64 aK4;

    UI64 aKa;
    UI64 aKb;
    UI64 aKc;
    UI64 aKd;

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

    int nMIRROR;



    K4benchmark();
    void   AddCall(const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp);
    double Statistics() const;
    void   Reset();
};


class p_ERIbuffer {
    friend class ERIroutine;
    friend class LibQuimera::ERIbuffer;

    double * buffer;
    double * bufferk4;
    double * v_m;
    double * uv_m;
    double * uv_m_t;

    /*
    double * f0;
    double * j0;
    double * j1;
    double * j2;
    double * j3;
    double * j4;
    */

    double * f0;
    double * f0e;
    double * f0f;
    double * k4j1e;
    double * k4j1f;
    double * k3j2e;
    double * k3j2f;
    double * k2j3e;
    double * k2j3f;
    double * k1j4e;
    double * k1j4f;


    cacheline * bufferMIRROR;
    cacheline * bufferK4;
    cacheline * v_m8;
    cacheline * uv_m8;
    cacheline * uv_m_t8;

    cacheline * KAB;
    cacheline * KCD;

  public:

    cacheline * F0;
    cacheline * F0e;
    cacheline * F0f;
    cacheline * K4J1e;
    cacheline * K4J1f;
    cacheline * K3J2e;
    cacheline * K3J2f;
    cacheline * K2J3e;
    cacheline * K2J3f;
    cacheline * K1J4e;
    cacheline * K1J4f;

/*
    cacheline * F0;
    cacheline * J0;
    cacheline * J1;
    cacheline * J2;
    cacheline * J3;
    cacheline * J4;
*/
    size_t SetBuffers (cacheline * mem, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp,  const p_Qalgorithm & routine);
};

class p_Qalgorithm {
    friend class p_ERIbuffer;
    friend class p_Quimera;
    friend class LibQuimera::Qalgorithm;

    //private:
    ERItype eritype;

    ERIroutine * IC; //interpreted code
    K4routine  * SC; //actual routine
    bool IsDynamic;

    //K4benchmark * bench;

    //for K4buffer
    UI32 MaxMem;
    UI32 K4Mem;

    UI32 memF0;
    UI32 memF0e;
    UI32 memF0f;
    UI32 memK4J1e;
    UI32 memK4J1f;
    UI32 memK3J2e;
    UI32 memK3J2f;
    UI32 memK2J3e;
    UI32 memK2J3f;
    UI32 memK1J4e;
    UI32 memK1J4f;

/*
    UI32 memF0;
    UI32 memJ1;
    UI32 memJ2;
    UI32 memJ3;
    UI32 memJ4;
*/

    UI32 nKernels;

    //init
    p_Qalgorithm(K4routine & routine);
    void SetStatic (ERItype   & type, K4routine & routine);
    void SetIC     (ERItype   & type);
    void Initialize();
    void Clear();

    //public
    p_Qalgorithm();
    ~p_Qalgorithm();

    UI32 GetNKernels() const;

    size_t MemSize    (const ShellPairPrototype & ABp, const ShellPairPrototype & CDp) const;

    //benchmarked (only double packed)
    void K4bench(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * uv_m_st8, p_ERIbuffer & buffer, bool UseCase = false) const;
    void MIRRORbench(const cacheline64 *   mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, p_ERIbuffer & buffer) const;

    //1 center routines (non benchmarked)
    void K4(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, double * uv_m_st, p_ERIbuffer & buffer, bool UseCase = false) const;
    void MIRROR(const double *  uv_m_st, double * W, p_ERIbuffer & buffer) const;

    //packed double routines
    void K4(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * uv_m_st8, p_ERIbuffer & buffer, bool UseCase = false) const;
    void MIRROR(const cacheline64 *   mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, p_ERIbuffer & buffer) const;

    //packed single routines
    void K4(const ERIgeometries32 & vars16, const ERITile32 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline32 * uv_m_st16, p_ERIbuffer & buffer, bool UseCase = false) const;
    void MIRROR(const cacheline32 *   mem, const ERIgeometries32 & vars16, cacheline32 * ERI16, p_ERIbuffer & buffer) const;

};

class p_Quimera {

  private:

    friend class LibQuimera::Quimera;
    friend class p_Qalgorithm;

    //private:
    std::map             <ERItype, Qalgorithm*>  Solvers;
    std::map             <ERItype, K4benchmark*> Benchmarks;

    std::priority_queue <ERItype>              SetLater;

    void SetStatic      (GEOM geo, UI8 la, UI8 lb, UI8 lc, UI8 ld, K4routine routine, bool overwrite);
    void SetInterpreted (GEOM geo, UI8 la, UI8 lb, UI8 lc, UI8 ld,                    bool overwrite);

    void LinkStatic     (bool overwrite);
    void LinkLibrary    (bool overwrite);

    const Qalgorithm * SelectAlgorithm (GEOM geometry, UI8  La, UI8  Lb, UI8  Lc, UI8  Ld);

    void ListNeeded     (GEOM geometry, UI8  La, UI8  Lb, UI8  Lc, UI8  Ld);
    void GenerateNeeded (bool overwrite);

    void Statistics() const;

  public:

    p_Quimera();
    ~p_Quimera();

} extern p_Q;



#endif


