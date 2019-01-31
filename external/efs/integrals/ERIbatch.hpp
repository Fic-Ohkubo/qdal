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


#ifndef __BATCH_EVALUATOR__
#define __BATCH_EVALUATOR__

#include "ERIblocks.hpp"
#include "../low/chrono.hpp"
#include "../libquimera/libquimera.hpp"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif


class BatchEvaluator;
class ERIblock;
class RotationMatrices32;
class RotationMatrices64;

//a memory buffer for the
struct BlockBuffer {
    size_t MaxBlockMem;
    BYTE * BlockMem;


    ERIgeometry   * E0d;
    ERIgeometry   * E0s;

    ERIgeometries64 * E8d;
    ERIgeometries32 * E8s;

    cacheline64     * I8d;
    cacheline32     * I8s;

    cacheline64     * T8d;
    cacheline32     * T8s;

    RotationMatrices64 * RM8d;
    RotationMatrices32 * RM8s;


    cacheline64 * Dab;
    cacheline64 * Dcd;

    cacheline64 * Dac;
    cacheline64 * Dad;
    cacheline64 * Dbc;
    cacheline64 * Dbd;


    cacheline64 * Fab;
    cacheline64 * Fcd;

    cacheline64 * Fac;
    cacheline64 * Fad;
    cacheline64 * Fbc;
    cacheline64 * Fbd;


    void SetBuffers(const ERIBatch & TheBatch, LibQuimera::ERIbuffer & ERIbuff);
};


class ERIbatchBuffer {
  private:
    //LibQuimera::ERITile64  * TheTiles;

    BYTE * mem;

    UI64 TotalMem;

    std::map<UI64, UI64> memFree;

    UI32 count;
    UI64 UsedMem;

    #ifdef _OPENMP
        omp_lock_t Lock;
    #endif

  public:
    ERIbatchBuffer();
    void Init(size_t MEM);
    void * Allocate (UI64 & Size);
    void Free(void * List, UI64 Size);

    void Check();
    ~ERIbatchBuffer();
};


//info about the kind of integral (BatchInfo)
class ERIBatch:public BatchInfo {
    friend class BatchEvaluator;
    friend class ERIblock;
    friend class Fock2e;
    friend class BlockBuffer;

  private:
    LibQuimera::ERITile64  * TileList64;
    LibQuimera::ERITile32  * TileList32;
    UI64       Ntiles64;
    UI64       Ntiles32;
    UI64       memAllocated;
    bool      J;
    bool      X;

    UI64 Nprescreened;


  public:
    void EvaluateOCLD (const Sparse & Ds);

    const ERIBatch & operator=(const ERIblock & rhs);

    ERIBatch() {
        TileList64 = NULL;
        Ntiles64   = 0;
        TileList32 = NULL;
        Ntiles32   = 0;
    }

    void clear(ERIbatchBuffer * EbatchBuffer) {
        if (TileList64!=NULL) EbatchBuffer->Free(TileList64, memAllocated);
        if (TileList32!=NULL) EbatchBuffer->Free(TileList32, memAllocated);
        TileList64 = NULL;
        Ntiles64   = 0;
        TileList32 = NULL;
        Ntiles32   = 0;
    }

};


class BatchEvaluator {

    //memory buffers and pointers for ERI evaluation
  private:
    LibQuimera::ERIbuffer    b4ERIs; //a buffer for the intermediate ERI values
    BlockBuffer  b4block;


    //GEOMETRY AND ROTATION MATRIX GENERATION
    //***************************************
    void Evaluate_LP_ABCD(ERIBatch & TheBatch);
    void Evaluate_LP_AACD(ERIBatch & TheBatch);
    void Evaluate_LP_AACC(ERIBatch & TheBatch);
    void Evaluate_LP_ABAB(ERIBatch & TheBatch);

    void Evaluate_RM_ABCD(ERIBatch & TheBatch);
    void Evaluate_RM_AACD(ERIBatch & TheBatch);
    void Evaluate_RM_AACC(ERIBatch & TheBatch);
    void Evaluate_RM_ABAB(ERIBatch & TheBatch);

    void MoveLists(ERIBatch & TheBatch);
    void KillLists(ERIBatch & TheBatch);


    //COLD PRISM STEPS
    //****************
    void Evaluate_OC(ERIBatch & TheBatch);

    void Evaluate_L(ERIBatch & TheBatch);
    void Evaluate_L_sameF(ERIBatch & TheBatch);
    template<bool J, bool X> void Evaluate_D       (ERIBatch & TheBatch, const Sparse & Ds);
    template<bool J, bool X> void Evaluate_D_sameF (ERIBatch & TheBatch, const Sparse & Ds);


  //functions and variables to be accessed publicly
  public:

    //chronometers
    //************
    Chronometer chronoG;
    Chronometer chronoF;
    Chronometer chronoK;
    Chronometer chronoT;
    Chronometer chronoC;

    Sparse Jcont;
    Sparse Xcont;
    //Sparse JX;

    BatchEvaluator();
    ~BatchEvaluator();




    //FUNCTIONS FOR BLOCK EVALUATION OF ERIs
    //**************************************

    void Init(size_t MEM);
    void SetBuffers(const ERIBatch & TheBatch);
    void Evaluate (ERIBatch & TheBatch, const Sparse & Ds); // evaluates coulomb and contraction separately


    cacheline64 CalcCauchySchwarz(const ERIBatch & TheBatch, const ShellPair & ABs, const ERIgeometries64 & eriG8, const LibQuimera::ERITile64 & ET, UI8 lenk);
    double   CalcCauchySchwarz(const ERIBatch & TheBatch, const ShellPair & ABs, const ERIgeometry & eriG);
    void     ConstructW       (const ERIBatch & TheBatch, const ShellPair & ABs, const ShellPair & CDs, const ERIgeometry & eriG, double * WWW, UI32 wa, bool OnlyJ);

    void ConstructWS(const ERIBatch & TheBatch, double w, const ShellPair & ABs, const ShellPair & CDs, const ERIgeometry & eriG, double * WWW, UI32 wa);
    void ConstructWL(const ERIBatch & TheBatch, double w, const ShellPair & ABs, const ShellPair & CDs, const ERIgeometry & eriG, double * WWW, UI32 wa);

};




#endif
