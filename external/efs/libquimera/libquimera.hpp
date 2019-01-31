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





    QUIMERA library;
    Query-based Unified Interface for Managing Electron-Repulsion Algorithms

    Much like the mythological animal and like the actual biological chimeras, the QUIMERA module
    is composed by potentially very different implementations of ERI algorithms in a single functional
    entity.

    The QUIMERA module is query-based. That is, once initialized, the function 'SelectAlgorithm'
    returns a pointer to an object which implements both the K4() and MIRROR() functions, to generate
    ERI kernels and transform them into batches, respectively.

    The module grants an interface to a working algorithm, but its implementation is hidden from the
    user. The object returned by 'SelectAlgorithm' can work by interpreting ERI Intermediate Code, by
    providing a pointer to a hard-coded static routine, or to a dynamically linked routine, depending
    on whether efficient hard-coded libraries are present.

    In a future it might also encapsulate OTF-generated specialized code and code executed in MIC
    coprocessors or CUDA-able cards. This, however, should be transparent to the user.


    2013.01.20, Jaime Axel Rosal Sandberg
    * initial version, product of major refactoring of earlier code
*/


#ifndef __LIB_QUIMERA__
#define __LIB_QUIMERA__

#include <cstdlib>
#include <string>
#include "../low/plotters.hpp"
#include "../defs.hpp"

class p_ERIbuffer;
class p_Qalgorithm;
class p_Quimera;
class cacheline64;
class cacheline32;
class cacheline;
class ShellPairPrototype;
class ERIgeometries64;
class ERIgeometries32;


namespace LibQuimera __attribute__ ((visibility ("default"))) {
    class Qalgorithm;

    //ERI computational quantum
    struct ERITile64 {
        UI32 ap12[DOUBLES_PER_CACHE_LINE];  // 32 bytes
        UI32 ap34[DOUBLES_PER_CACHE_LINE];  // 32 bytes

        const double * wSP12;  // 8 bytes
        const double * wSP34;  // 8 bytes

        // these are the offsets within the atom pair block for Sparse tensors;
        // they are the same  for all integrals involving the same elemental GTO quadruplet,
        // but storing them makes it possible to decouple the OC and DL parts of the computation
        // later on, so that all transformations/contraction involving same angular momenta can be
        // processed in a batch; also, the structure is already >64 bytes, so it doesnt harm performance
        UI32 offAB;
        UI32 offCD;
        UI32 offAC;
        UI32 offAD;

        UI32 offBC;
        UI32 offBD;

        UI16 nKab;  //number of gaussian products to be evaluated -1
        UI16 nKcd;

        UI16 nK2ab; // total number of gaussian products in the elemental shellpair
        UI16 nK2cd;

        UI8 used;      //number of used elements (0-DOUBLES_PER_CACHE_LINE)
        UI8 distance;  //distance block

        UI8 use[DOUBLES_PER_CACHE_LINE];

        UI8 lowgeom;   // for avoiding 3-center ABAD-type integrals when using CDR

        UI8 align[5];

        //UI8 align[34]; // maybe i'll find something more useful later
    } __attribute__((aligned(128)));

    struct ERITile32 {
        UI32 ap12[FLOATS_PER_CACHE_LINE];  // 64 bytes
        UI32 ap34[FLOATS_PER_CACHE_LINE];  // 64 bytes

        const double * wSP12;  // 8 bytes
        const double * wSP34;  // 8 bytes

        // these are the offsets within the atom pair block for Sparse tensors;
        // they are the same  for all integrals involving the same elemental GTO quadruplet,
        // but storing them makes it possible to decouple the OC and DL parts of the computation
        // later on, so that all transformations/contraction involving same angular momenta can be
        // processed in a batch; also, the structure is already >64 bytes, so it doesnt harm performance
        UI32 offAB;
        UI32 offCD;
        UI32 offAC;
        UI32 offAD;

        UI32 offBC;
        UI32 offBD;

        UI16 nKab;  //number of gaussian products to be evaluated -1
        UI16 nKcd;

        UI16 nK2ab; // total number of gaussian products in the elemental shellpair
        UI16 nK2cd;

        UI8 used;      //number of used elements (0-FLOATS_PER_CACHE_LINE)
        UI8 distance;  //distance block

        UI8 use[FLOATS_PER_CACHE_LINE];

        UI8 align2[62];

        //UI8 align[34]; // maybe i'll find something more useful later
    } __attribute__((aligned(256)));


    class ERIbuffer {
        friend class Qalgorithm;

      private:
        p_ERIbuffer * pERIbuffer; //private implementation

      public:
        ERIbuffer();
        ~ERIbuffer();
        size_t SetBuffers (cacheline * mem, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp,  const Qalgorithm & routine);
    };

    class Qalgorithm {
        friend class ::p_Quimera;
        friend class ERIbuffer;

      private:
        p_Qalgorithm * pQalgorithm; //private implementation

      public:
        Qalgorithm();
        ~Qalgorithm();

        UI32 GetNKernels() const;

        size_t MemSize (const ShellPairPrototype & ABp, const ShellPairPrototype & CDp) const;

        //benchmarked (only double packed)
        void K4bench(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * uv_m_st8, ERIbuffer & buffer, bool UseCase = false) const;
        void MIRRORbench(const cacheline64 *   mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, ERIbuffer & buffer) const;

        //1 center routines (non benchmarked)
        void K4(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, double * uv_m_st, ERIbuffer & buffer, bool UseCase = false) const;
        void MIRROR(const double *  uv_m_st, double * W, ERIbuffer & buffer) const;

        //packed double routines
        void K4(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * uv_m_st8, ERIbuffer & buffer, bool UseCase = false) const;
        void MIRROR(const cacheline64 *   mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, ERIbuffer & buffer) const;

        //packed single routines
        void K4(const ERIgeometries32 & vars16, const ERITile32 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline32 * uv_m_st16, ERIbuffer & buffer, bool UseCase = false) const;
        void MIRROR(const cacheline32 *   mem, const ERIgeometries32 & vars16, cacheline32 * ERI16, ERIbuffer & buffer) const;
    };

    class Quimera {
        friend class Qalgorithm;

      public:
        static MessagePlotter QMessenger;
        static MessagePlotter QBenchmarker;

      private:
        p_Quimera * pQuimera; //private implementation

      public:
        Quimera();
        ~Quimera();

        const Qalgorithm * SelectAlgorithm (GEOM geometry, UI8  La, UI8  Lb, UI8  Lc, UI8  Ld) const ;

        void ListNeeded     (GEOM geometry, UI8  La, UI8  Lb, UI8  Lc, UI8  Ld);
        void GenerateNeeded (bool overwrite);

        void Statistics() const;
    } extern Q;

    int InitLibQuimera();

    extern bool    Initialized;

    extern std::string IC_DIR;

    //not used in DALTON
    extern double  case_w2;
    extern bool    UseCASE;
    extern bool    UseCDR;
    extern bool    UseGC;
}



#endif
