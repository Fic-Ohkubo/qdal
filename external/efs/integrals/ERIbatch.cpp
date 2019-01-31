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


#include "ERIbatch.hpp"
#include "ERIblocks.hpp"
#include "../math/angular.hpp"
#include "../libechidna/libechidna.hpp"
#include "../basis/shellpair.hpp"
#include "../integrals/newcontractions.hpp"


#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

using namespace LibAngular;


const int DPC = DOUBLES_PER_CACHE_LINE;
const int FPC = FLOATS_PER_CACHE_LINE;


static inline BYTE * align2cacheline64(BYTE * p) {
    PSIZE addr = (PSIZE)(p); //convert to integer
    if (addr % CACHE_LINE_SIZE != 0) addr += CACHE_LINE_SIZE - addr % CACHE_LINE_SIZE; //align to next address multiple of 64
    return (BYTE*)addr; //cast back to pointer
}

static inline size_t align2cacheline64(size_t offset) {
    if (offset % CACHE_LINE_SIZE != 0)
        offset += CACHE_LINE_SIZE - offset % CACHE_LINE_SIZE; //align to next address multiple of 64
    return offset;
}


ERIbatchBuffer::ERIbatchBuffer() {
    //TheTiles = NULL;
    //MaxTiles = 0;
    mem = NULL;
    count   = 0;
    UsedMem  = 0;
}

void ERIbatchBuffer::Init(size_t MEM) {
    TotalMem = MEM;
    UI32 error;

    do {
        //aligning to cache line size avoids different cores fighting for the same cache line
        error = posix_memalign((void**)(&mem), 4096, TotalMem);

        if (error) TotalMem/=2;
    } while (error!=0);

    memFree[0] = TotalMem;
    #ifdef _OPENMP
    omp_init_lock (&Lock);
    #endif
}

ERIbatchBuffer::~ERIbatchBuffer() {
    free(mem);
    #ifdef _OPENMP
    omp_destroy_lock (&Lock);
    #endif
}

void * ERIbatchBuffer::Allocate (UI64 & Size) {
    std::map<UI64, UI64>::iterator it;

    BYTE * ret;

    #ifdef _OPENMP
    omp_set_lock (&Lock);
    #endif
    {
        // best second choice
        std::map<UI64, UI64>::iterator itc = memFree.begin();

        for (it=memFree.begin(); it!=memFree.end(); ++it) {
            if (it->second>=Size) break;                  //found some memory
            else if (it->second > itc->second) itc = it;  //found a better second choice
        }

        //couldn't find the memory requested
        if (it==memFree.end()) {
            //no memory at all, actually
            if (itc==memFree.end())
                ret = NULL;
            else {
                //return as many as found
                Size = itc->second;

                UI64 offset = itc->first;
                memFree.erase(itc);

                ret = mem + offset;
                ++count;
                UsedMem += Size;
            }
        }
        //found enough memory
        else {
            UI64 offset = it->first;
            UI64 tsize    = it->second;
            UI64 left = tsize - Size;
            memFree.erase(it);

            if (left>0) memFree[offset+Size] = left; //memory left

            ret = mem + offset;
            ++count;
            UsedMem += Size;
        }
    }
    #ifdef _OPENMP
    omp_unset_lock (&Lock);
    #endif


    return (void*)ret;
}

void ERIbatchBuffer::Free(void * List, UI64 Size) {
    UI64 offset  = (BYTE*)List - mem;

    #ifdef _OPENMP
    omp_set_lock (&Lock);
    #endif
    {
        std::map<UI64, UI64>::iterator it, itn;

        //insert the freed memory
        memFree[offset] = Size;

        //merge with next if possible
        {
            it = memFree.find(offset);
            itn = it;
            ++itn;
            if ( itn!=memFree.end() &&  offset + Size == itn->first ) {
                //merge the memory
                it->second += itn->second;
                //delete next
                memFree.erase(itn);
            }
        }

         //marge with previous, if possible
        {
            it = memFree.find(offset);

            if ( it!=memFree.begin() ) {
                itn = it;
                --it;
                if ( it->first + it->second == offset ) {
                    //merge the memory
                    it->second += itn->second;
                    //delete next
                    memFree.erase(itn);
                }
            }
        }

        --count;
        UsedMem -= Size;
    }
    #ifdef _OPENMP
    omp_unset_lock (&Lock);
    #endif

}

void ERIbatchBuffer::Check()  {
    std::map<UI64, UI64>::const_iterator it;

    if (count>0 || UsedMem>0) {
        Echidna::EDebugger << "ERI batch buffer was not completely freed after computation: " << std::endl;
        Echidna::EDebugger << "Imbalance: " << count << std::endl;
        Echidna::EDebugger << "Lost mem : " << UsedMem << std::endl;

        for (it=memFree.begin(); it!=memFree.end(); ++it) {
            Echidna::EDebugger << it->first << " " << it->second << std::endl;
        }

        Echidna::EDebugger << "Erasing" << std::endl;
        count = 0;
        UsedMem = 0;
        memFree.clear();
        memFree[0] = TotalMem;
    }
}




void BlockBuffer::SetBuffers(const ERIBatch & TheBatch,  LibQuimera::ERIbuffer & ERIbuff) {

    //compute and initialize the (fixed) amount of buffer memory required for contraction and evaluation
    size_t ERIoffset = ERIbuff.SetBuffers((cacheline*)BlockMem, *TheBatch.ABp, *TheBatch.CDp,  *TheBatch.ERIalgorithm);


    UI32 NTiles64 = TheBatch.Ntiles64;
    UI32 NTiles32 = TheBatch.Ntiles32;

    //UI32 NInts  = NTiles64 * DPC + NTiles32 * FPC;
    //if (NTiles32==0) NTiles32 = NTiles64;


    //distribute memory accordingly, align to cache line
    BYTE * p = BlockMem + ERIoffset;

    E0d = (ERIgeometry*)p; p += NTiles64 * DPC * sizeof(ERIgeometry);
    E0s = (ERIgeometry*)p; p += NTiles32 * FPC * sizeof(ERIgeometry);

    E8d = (ERIgeometries64*)p; p += NTiles64 * sizeof(ERIgeometries64);
    E8s = (ERIgeometries32*)p; p += NTiles32 * sizeof(ERIgeometries32);

    I8d = (cacheline64*)p; p += NTiles64 * TheBatch.msize4 * sizeof(cacheline64);
    I8s = (cacheline32*)p; p += NTiles32 * TheBatch.msize4 * sizeof(cacheline32);

    T8d = (cacheline64*)p; p += NTiles64 * TheBatch.wsize4 * sizeof(cacheline64);
    T8s = (cacheline32*)p; p += NTiles32 * TheBatch.wsize4 * sizeof(cacheline32);


    RM8d = (RotationMatrices64*)p; p += NTiles64 * sizeof(RotationMatrices64);
    RM8s = (RotationMatrices32*)p; p += NTiles32 * sizeof(RotationMatrices32);


    int JA = TheBatch.ABp->Ja;
    int JB = TheBatch.ABp->Jb;
    int JC = TheBatch.CDp->Ja;
    int JD = TheBatch.CDp->Jb;

    Dab  = (cacheline64*)p; p += JA*JB* (2*TheBatch.la+1) * (2*TheBatch.lb+1) * sizeof(cacheline64);
    Dcd  = (cacheline64*)p; p += JC*JD* (2*TheBatch.lc+1) * (2*TheBatch.ld+1) * sizeof(cacheline64);

    Dac  = (cacheline64*)p; p += JA*JC* (2*TheBatch.la+1) * (2*TheBatch.lc+1) * sizeof(cacheline64);
    Dad  = (cacheline64*)p; p += JA*JD* (2*TheBatch.la+1) * (2*TheBatch.ld+1) * sizeof(cacheline64);
    Dbc  = (cacheline64*)p; p += JB*JC* (2*TheBatch.lb+1) * (2*TheBatch.lc+1) * sizeof(cacheline64);
    Dbd  = (cacheline64*)p; p += JB*JD* (2*TheBatch.lb+1) * (2*TheBatch.ld+1) * sizeof(cacheline64);


    Fab  = (cacheline64*)p; p += JA*JB* (2*TheBatch.la+1) * (2*TheBatch.lb+1) * sizeof(cacheline64);
    Fcd  = (cacheline64*)p; p += JC*JD* (2*TheBatch.lc+1) * (2*TheBatch.ld+1) * sizeof(cacheline64);

    Fac  = (cacheline64*)p; p += JA*JC* (2*TheBatch.la+1) * (2*TheBatch.lc+1) * sizeof(cacheline64);
    Fad  = (cacheline64*)p; p += JA*JD* (2*TheBatch.la+1) * (2*TheBatch.ld+1) * sizeof(cacheline64);
    Fbc  = (cacheline64*)p; p += JB*JC* (2*TheBatch.lb+1) * (2*TheBatch.lc+1) * sizeof(cacheline64);
    Fbd  = (cacheline64*)p; p += JB*JD* (2*TheBatch.lb+1) * (2*TheBatch.ld+1) * sizeof(cacheline64);


    if (p > BlockMem + MaxBlockMem) {
        //Echidna::EMessenger << "Error in ERI block memory allocation: final offset is larger than available memory ";
        //Echidna::EMessenger << (p-BlockMem) << "  " << MaxBlockMem << std::endl;
        //Echidna::EMessenger << ERIoffset << std::endl;

        std::cout << "Error in ERI block memory allocation: final offset is larger than available memory ";
        std::cout << (p-BlockMem) << "  " << MaxBlockMem << std::endl;
        std::cout << ERIoffset << std::endl;

        throw(523);
    }
}




//copy everything
const ERIBatch & ERIBatch::operator=(const ERIblock & rhs) {
    ABp = rhs.ABp;
    CDp = rhs.CDp;

    SkipSame  = rhs.SkipSame;
    SameShell = rhs.SameShell;
    SameList  = rhs.SameList;

    geometry  = rhs.geometry;

    la   = rhs.la;
    lb   = rhs.lb;
    lc   = rhs.lc;
    ld   = rhs.ld;
    Lt   = rhs.Lt;
    Lmax = rhs.Lmax;

    //size of the ERI batch
    J4     = rhs.J4;
    fsize  = rhs.fsize;

    msize  = rhs.msize;
    msize4 = rhs.msize4;

    wsize  = rhs.wsize;
    wsize4 = rhs.wsize4;


    //funciones necesarias
    //********************
    ERIalgorithm = rhs.ERIalgorithm;

    //window arguments
    //****************

    AP12list = rhs.AP12list;
    AP34list = rhs.AP34list;

    SP12 = rhs.SP12;
    SP34 = rhs.SP34;

    tmax12 = rhs.tmax12;
    tmax34 = rhs.tmax34;

    State   = rhs.State;
}


BatchEvaluator::BatchEvaluator() {
}

BatchEvaluator::~BatchEvaluator() {
}

void BatchEvaluator::Init(size_t MEM) {
    b4block.MaxBlockMem = MEM;

    UI32 error;

    do {
        //aligning to cache line size avoids different cores fighting for the same cache line
        error = posix_memalign((void**)(&b4block.BlockMem), sizeof(cacheline64), b4block.MaxBlockMem);

        if (error) b4block.MaxBlockMem/=2;
    } while (error!=0);

    if (b4block.MaxBlockMem!=MEM) {
        std::cout << "Couldn't allocate " << MEM << " bytes of memory for ERI buffering; allocating " << b4block.MaxBlockMem << " instead" << std::endl;
    }

};

void BatchEvaluator::SetBuffers(const ERIBatch & TheBatch) {
    b4block.SetBuffers(TheBatch, b4ERIs);
}


void BatchEvaluator::Evaluate (ERIBatch & TheBatch, const Sparse & Ds) {

    GEOM geometry = TheBatch.geometry;

    bool J = TheBatch.J;
    bool X = TheBatch.X;
    bool S = (TheBatch.geometry==ABAB) && TheBatch.SameShell;

    b4block.SetBuffers(TheBatch, b4ERIs);

    //geometries
    chronoG.Start();
    if      (geometry==ABCD) Evaluate_LP_ABCD(TheBatch);
    else if (geometry==AACD) Evaluate_LP_AACD(TheBatch);
    else if (geometry==AACC) Evaluate_LP_AACC(TheBatch);
    else if (geometry==ABAB) Evaluate_LP_ABAB(TheBatch);
    chronoG.Stop();

    //OC steps: gamma function generation and K4 kernel contraction
    Evaluate_OC(TheBatch);

    //L step: MIRROR transformations
    if (!S) Evaluate_L      (TheBatch);
    else    Evaluate_L_sameF(TheBatch);

    //generate vector packed rotation matrices for spherical harmonics
    chronoG.Start();
    if      (geometry==ABCD) Evaluate_RM_ABCD(TheBatch);
    else if (geometry==AACD) Evaluate_RM_AACD(TheBatch);
    else if (geometry==AACC) Evaluate_RM_AACC(TheBatch);
    else if (geometry==ABAB) Evaluate_RM_ABAB(TheBatch);
    chronoG.Stop();

    //D step: ERI Digestion with the DM
    if (!S) {
        if      ( J &&  X) Evaluate_D <true,  true>  (TheBatch, Ds);
        else if ( J && !X) Evaluate_D <true,  false> (TheBatch, Ds);
        else if (!J &&  X) Evaluate_D <false, true>  (TheBatch, Ds);
        else               Evaluate_D <false, false> (TheBatch, Ds);
    }
    else {
        if      ( J &&  X) Evaluate_D_sameF <true,  true>  (TheBatch, Ds);
        else if ( J && !X) Evaluate_D_sameF <true,  false> (TheBatch, Ds);
        else if (!J &&  X) Evaluate_D_sameF <false, true>  (TheBatch, Ds);
        else               Evaluate_D_sameF <false, false> (TheBatch, Ds);
    }
};


void BatchEvaluator::Evaluate_LP_ABCD(ERIBatch & TheBatch) {

    {
        ERIgeometry     * erig  = b4block.E0d;
        ERIgeometries64 * erig8 = b4block.E8d;

        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {
            for (UI32 j=0; j<DPC; ++j) {
                UI32 ap12 = TheBatch.TileList64[i].ap12[j]; const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                UI32 ap34 = TheBatch.TileList64[i].ap34[j]; const AtomProd & AP34 = (*TheBatch.AP34list)[ap34];
                erig[j].MakeRotation4  (AP12, AP34);
                erig[j].Adjust(TheBatch.ABp->inverted, TheBatch.CDp->inverted);
            }
            PackGeometries (erig, *erig8);

            erig8++;
            erig += DPC;
        }
    }

    {
        ERIgeometry     * erig  = b4block.E0s;
        ERIgeometries32 * erig8 = b4block.E8s;

        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {
            for (UI32 j=0; j<FPC; ++j) {
                UI32 ap12 = TheBatch.TileList32[i].ap12[j]; const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                UI32 ap34 = TheBatch.TileList32[i].ap34[j]; const AtomProd & AP34 = (*TheBatch.AP34list)[ap34];
                erig[j].MakeRotation4  (AP12, AP34);
                erig[j].Adjust(TheBatch.ABp->inverted, TheBatch.CDp->inverted);
            }
            PackGeometries (erig, *erig8);

            erig8++;
            erig += FPC;
        }
    }

}

void BatchEvaluator::Evaluate_LP_AACD(ERIBatch & TheBatch) {


    {
        ERIgeometry     * erig  = b4block.E0d;
        ERIgeometries64 * erig8 = b4block.E8d;

        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {
            for (UI32 j=0; j<DPC; ++j) {
                UI32 ap12 = TheBatch.TileList64[i].ap12[j]; const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                UI32 ap34 = TheBatch.TileList64[i].ap34[j]; const AtomProd & AP34 = (*TheBatch.AP34list)[ap34];

                erig[j].MakeRotation3  (AP12, AP34);
                erig[j].Adjust(TheBatch.ABp->inverted, TheBatch.CDp->inverted);
            }

            PackGeometries (erig, *erig8);
            erig8++;
            erig += DPC;
        }
    }

    {
        ERIgeometry     * erig  = b4block.E0s;
        ERIgeometries32 * erig8 = b4block.E8s;

        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {
            for (UI32 j=0; j<FPC; ++j) {
                UI32 ap12 = TheBatch.TileList32[i].ap12[j]; const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                UI32 ap34 = TheBatch.TileList32[i].ap34[j]; const AtomProd & AP34 = (*TheBatch.AP34list)[ap34];

                erig[j].MakeRotation3  (AP12, AP34);
                erig[j].Adjust(TheBatch.ABp->inverted, TheBatch.CDp->inverted);
            }

            PackGeometries (erig, *erig8);
            erig8++;
            erig += FPC;
        }
    }

}

void BatchEvaluator::Evaluate_LP_AACC(ERIBatch & TheBatch) {

    {
        ERIgeometry     * erig  = b4block.E0d;
        ERIgeometries64 * erig8 = b4block.E8d;

        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {
            for (UI32 k=0; k<DPC; ++k) {
                UI32 ap12 = TheBatch.TileList64[i].ap12[k]; const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k]; const AtomProd & AP34 = (*TheBatch.AP34list)[ap34];
                erig[k].MakeRotation2  (AP12.A, AP34.A);
            }
            PackGeometries (erig, *erig8);
            erig8++;
            erig += DPC;
        }
    }

    {
        ERIgeometry     * erig  = b4block.E0s;
        ERIgeometries32 * erig8 = b4block.E8s;

        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {
            for (UI32 k=0; k<FPC; ++k) {
                UI32 ap12 = TheBatch.TileList32[i].ap12[k]; const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k]; const AtomProd & AP34 = (*TheBatch.AP34list)[ap34];
                erig[k].MakeRotation2  (AP12.A, AP34.A);
            }
            PackGeometries (erig, *erig8);
            erig8++;
            erig += FPC;
        }
    }

}

void BatchEvaluator::Evaluate_LP_ABAB(ERIBatch & TheBatch) {

    {
        ERIgeometry     * erig  = b4block.E0d;
        ERIgeometries64 * erig8 = b4block.E8d;

        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {
            for (UI32 k=0; k<DPC; ++k) {
                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                erig[k].MakeRotation2S ( AP12 );
                erig[k].Adjust(TheBatch.ABp->inverted, TheBatch.CDp->inverted);
            }
            PackGeometries (erig, *erig8);
            erig8++;
            erig += DPC;
        }
    }

    {
        ERIgeometry     * erig  = b4block.E0s;
        ERIgeometries32 * erig8 = b4block.E8s;

        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {
            for (UI32 k=0; k<FPC; ++k) {
                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                erig[k].MakeRotation2S ( AP12 );
                erig[k].Adjust(TheBatch.ABp->inverted, TheBatch.CDp->inverted);
            }
            PackGeometries (erig, *erig8);
            erig8++;
            erig += FPC;
        }
    }


}

//decoupled paths
void BatchEvaluator::Evaluate_OC(ERIBatch & TheBatch) {
    chronoK.Start(); {

        bool OnlyJ = TheBatch.J && !TheBatch.X;

        {
            cacheline64     * i8   = b4block.I8d;
            ERIgeometries64 * e8   = b4block.E8d;

            //regular kernel contraction
            for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {
                TheBatch.ERIalgorithm->K4bench (e8[i], TheBatch.TileList64[i], *TheBatch.ABp, *TheBatch.CDp, i8, b4ERIs, OnlyJ);

                i8 += TheBatch.msize4;
            }
        }

        {
            cacheline32     * i8   = b4block.I8s;
            ERIgeometries32 * e8   = b4block.E8s;

            //regular kernel contraction
            for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {
                TheBatch.ERIalgorithm->K4 (e8[i], TheBatch.TileList32[i], *TheBatch.ABp, *TheBatch.CDp, i8, b4ERIs, OnlyJ);

                i8 += TheBatch.msize4;
            }
        }
    }

    chronoK.Stop();
}


void BatchEvaluator::Evaluate_L(ERIBatch & TheBatch) {
    const UI8 JA = TheBatch.ABp->Ja;
    const UI8 JB = TheBatch.ABp->Jb;
    const UI8 JC = TheBatch.CDp->Ja;
    const UI8 JD = TheBatch.CDp->Jb;

    chronoT.Start(); {

        {
            ERIgeometries64 * e8  = b4block.E8d;
            cacheline64     * i8  = b4block.I8d;
            cacheline64     * t8  = b4block.T8d;

            //packed
            for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {

                for (int ja=0; ja<TheBatch.ABp->Ja; ++ja) {
                    for (int jb=0; jb<TheBatch.ABp->Jb; ++jb) {
                        for (int jc=0; jc<TheBatch.CDp->Ja; ++jc) {
                            for (int jd=0; jd<TheBatch.CDp->Jb; ++jd) {
                                TheBatch.ERIalgorithm->MIRRORbench(i8, e8[i], t8, b4ERIs);

                                i8 += TheBatch.msize;
                                t8 += TheBatch.wsize;
                            }
                        }
                    }
                }
            }
        }

        {
            ERIgeometries32 * e8  = b4block.E8s;
            cacheline32     * i8  = b4block.I8s;
            cacheline32     * t8  = b4block.T8s;

            //packed
            for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

                for (int ja=0; ja<TheBatch.ABp->Ja; ++ja) {
                    for (int jb=0; jb<TheBatch.ABp->Jb; ++jb) {
                        for (int jc=0; jc<TheBatch.CDp->Ja; ++jc) {
                            for (int jd=0; jd<TheBatch.CDp->Jb; ++jd) {
                                TheBatch.ERIalgorithm->MIRROR(i8, e8[i], t8, b4ERIs);

                                i8 += TheBatch.msize;
                                t8 += TheBatch.wsize;
                            }
                        }
                    }
                }
            }
        }

    }
    chronoT.Stop();

}

void BatchEvaluator::Evaluate_L_sameF(ERIBatch & TheBatch) {
    const UI8 JA = TheBatch.ABp->Ja;
    const UI8 JB = TheBatch.ABp->Jb;
    const UI8 JC = TheBatch.CDp->Ja;
    const UI8 JD = TheBatch.CDp->Jb;

    //APPLY TRANSFORMATIONS

    chronoT.Start(); {
        {
            ERIgeometries64 * e8 = b4block.E8d;
            cacheline64     * i8 = b4block.I8d;
            cacheline64     * t8 = b4block.T8d;

            //packed
            for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {

                for (int ja=0; ja<TheBatch.ABp->Ja; ++ja) {
                    for (int jb=0; jb<TheBatch.ABp->Jb; ++jb) {
                        for (int jc=0; jc<TheBatch.CDp->Ja; ++jc) {
                            for (int jd=0; jd<TheBatch.CDp->Jb; ++jd) {
                                //skip unnecessary (repeated) transformations
                                if (ja*JB+jb >= jc*JD+jd) {
                                    TheBatch.ERIalgorithm->MIRRORbench(i8, e8[i], t8, b4ERIs);
                                }
                                i8 += TheBatch.msize;
                                t8 += TheBatch.wsize;
                            }
                        }
                    }
                }
            }
        }

        {
            ERIgeometries32 * e8 = b4block.E8s;
            cacheline32     * i8 = b4block.I8s;
            cacheline32     * t8 = b4block.T8s;

            //packed
            for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

                for (int ja=0; ja<TheBatch.ABp->Ja; ++ja) {
                    for (int jb=0; jb<TheBatch.ABp->Jb; ++jb) {
                        for (int jc=0; jc<TheBatch.CDp->Ja; ++jc) {
                            for (int jd=0; jd<TheBatch.CDp->Jb; ++jd) {
                                //skip unnecessary (repeated) transformations
                                if (ja*JB+jb >= jc*JD+jd) {
                                    TheBatch.ERIalgorithm->MIRROR(i8, e8[i], t8, b4ERIs);
                                }
                                i8 += TheBatch.msize;
                                t8 += TheBatch.wsize;
                            }
                        }
                    }
                }
            }
        }

    }
    chronoT.Stop();

}




void BatchEvaluator::Evaluate_RM_ABCD(ERIBatch & TheBatch) {

    {
        ERIgeometry        * erig = b4block.E0d;
        RotationMatrices64 * rm8  = b4block.RM8d;

        RotationMatrix RM[DPC];

        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {

            //pack the rotation matrices in a cacheline64-array
            for (UI32 k=0; k<DPC; ++k) {
                if (TheBatch.TileList64[i].use[k] == 0) continue;
                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k];

                RM[k].From((*TheBatch.AP12list)[ap12].RM,  erig[k], TheBatch.Lmax);
            }

            erig += DPC;
            rm8[i].From(RM, TheBatch.Lmax);
        }
    }

    {
        ERIgeometry        * erig = b4block.E0s;
        RotationMatrices32 * rm8  = b4block.RM8s;

        RotationMatrix RM[FPC];

        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

            //pack the rotation matrices in a cacheline32 array
            for (UI32 k=0; k<FPC; ++k) {
                if (TheBatch.TileList32[i].use[k] == 0) continue;
                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k];

                RM[k].From((*TheBatch.AP12list)[ap12].RM,  erig[k], TheBatch.Lmax);
            }

            erig += FPC;
            rm8[i].From(RM, TheBatch.Lmax);
        }
    }


}

void BatchEvaluator::Evaluate_RM_AACD(ERIBatch & TheBatch) {

    {
        ERIgeometry        * erig = b4block.E0d;
        RotationMatrices64 * rm8  = b4block.RM8d;

        RotationMatrix RM[DPC];

        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {

            //pack the rotation matrices in a cacheline64-array
            for (UI32 k=0; k<DPC; ++k) {
                if (TheBatch.TileList64[i].use[k] == 0) continue;
                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k];

                RM[k].From((*TheBatch.AP34list)[ap34].RM,  erig[k], TheBatch.Lmax);
            }

            erig += DPC;
            rm8[i].From(RM, TheBatch.Lmax);
        }
    }

    {
        ERIgeometry        * erig = b4block.E0s;
        RotationMatrices32 * rm8  = b4block.RM8s;

        RotationMatrix RM[FPC];

        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

            //pack the rotation matrices in a cacheline64-array
            for (UI32 k=0; k<FPC; ++k) {
                if (TheBatch.TileList32[i].use[k] == 0) continue;
                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k];

                RM[k].From((*TheBatch.AP34list)[ap34].RM,  erig[k], TheBatch.Lmax);
            }

            erig += FPC;
            rm8[i].From(RM, TheBatch.Lmax);
        }
    }


}

void BatchEvaluator::Evaluate_RM_AACC(ERIBatch & TheBatch) {

    {
        ERIgeometry        * erig = b4block.E0d;
        RotationMatrices64 * rm8  = b4block.RM8d;

        RotationMatrix RM[DPC];

        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {

            //pack the rotation matrices in a cacheline64-array
            for (UI32 k=0; k<DPC; ++k) {
                if (TheBatch.TileList64[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k];

                const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                const AtomProd & AP34 = (*TheBatch.AP34list)[ap34];

                RM[k].From(AP12.A, AP34.A, TheBatch.Lmax);
            }

            erig += DPC;
            rm8[i].From(RM, TheBatch.Lmax);
        }

    }

    {
        ERIgeometry        * erig = b4block.E0s;
        RotationMatrices32 * rm8  = b4block.RM8s;

        RotationMatrix RM[FPC];

        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

            //pack the rotation matrices in a cacheline64-array
            for (UI32 k=0; k<FPC; ++k) {
                if (TheBatch.TileList32[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k];

                const AtomProd & AP12 = (*TheBatch.AP12list)[ap12];
                const AtomProd & AP34 = (*TheBatch.AP34list)[ap34];

                RM[k].From(AP12.A, AP34.A, TheBatch.Lmax);
            }

            erig += FPC;
            rm8[i].From(RM, TheBatch.Lmax);
        }

    }


}

void BatchEvaluator::Evaluate_RM_ABAB(ERIBatch & TheBatch) {

    {
        ERIgeometry        * erig = b4block.E0d;
        RotationMatrices64 * rm8  = b4block.RM8d;

        RotationMatrix RM[DPC];

        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {

            //pack the rotation matrices in a cacheline64-array
            for (UI32 k=0; k<DPC; ++k) {
                if (TheBatch.TileList64[i].use[k] == 0) continue;
                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k];

                RM[k] = (*TheBatch.AP12list)[ap12].RM;
            }

            erig += DPC;
            rm8[i].From(RM, TheBatch.Lmax);
        }
    }

    {
        ERIgeometry        * erig = b4block.E0s;
        RotationMatrices32 * rm8  = b4block.RM8s;

        RotationMatrix RM[FPC];

        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

            //pack the rotation matrices in a cacheline64-array
            for (UI32 k=0; k<FPC; ++k) {
                if (TheBatch.TileList32[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k];

                RM[k] = (*TheBatch.AP12list)[ap12].RM;
            }

            erig += FPC;
            rm8[i].From(RM, TheBatch.Lmax);
        }
    }

}



void BatchEvaluator::MoveLists(ERIBatch & TheBatch) {

    TheBatch.Ntiles32 = TheBatch.Ntiles64/2;
    TheBatch.TileList32 = new LibQuimera::ERITile32[TheBatch.Ntiles32+1];

    //fuse tiles
    for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

        LibQuimera::ERITile32 & ET32 = TheBatch.TileList32[i];

        LibQuimera::ERITile64 & ET64a = TheBatch.TileList64[2*i];
        LibQuimera::ERITile64 & ET64b = TheBatch.TileList64[2*i+1];

        //fuse tiles
        for (UI32 k=0; k<DPC; ++k) {
            ET32.ap12[k]     = ET64a.ap12[k];
            ET32.ap34[k]     = ET64a.ap34[k];

            ET32.use[k]      = ET64a.use[k];


            ET32.ap12[DPC+k] = ET64b.ap12[k];
            ET32.ap34[DPC+k] = ET64b.ap34[k];

            ET32.use[DPC+k]  = ET64b.use[k];
        }

        ET32.used = ET64a.used + ET64b.used;

        ET32.wSP12 = ET64a.wSP12;
        ET32.wSP34 = ET64a.wSP34;

        ET32.offAB = ET64a.offAB;
        ET32.offCD = ET64a.offCD;
        ET32.offAC = ET64a.offAC;
        ET32.offAD = ET64a.offAD;
        ET32.offBC = ET64a.offBC;
        ET32.offBD = ET64a.offBD;

        ET32.nKab = ET64a.nKab;
        ET32.nKcd = ET64a.nKcd;

        ET32.nK2ab = ET64a.nK2ab;
        ET32.nK2cd = ET64a.nK2cd;
    }

    if (TheBatch.Ntiles64%2) {

        LibQuimera::ERITile32 & ET32 = TheBatch.TileList32[TheBatch.Ntiles32];
        LibQuimera::ERITile64 & ET64 = TheBatch.TileList64[2*TheBatch.Ntiles32];

        //fuse tiles
        for (UI32 k=0; k<DPC; ++k) {
            ET32.ap12[k]     = ET64.ap12[k];
            ET32.ap34[k]     = ET64.ap34[k];

            ET32.use[k]      = ET64.use[k];

            ET32.ap12[DPC+k] = ET64.ap12[k];
            ET32.ap34[DPC+k] = ET64.ap34[k];

            ET32.use[DPC+k]  = 0;
        }

        ET32.used = ET64.used;

        ET32.wSP12 = ET64.wSP12;
        ET32.wSP34 = ET64.wSP34;

        ET32.offAB = ET64.offAB;
        ET32.offCD = ET64.offCD;
        ET32.offAC = ET64.offAC;
        ET32.offAD = ET64.offAD;
        ET32.offBC = ET64.offBC;
        ET32.offBD = ET64.offBD;

        ET32.nKab = ET64.nKab;
        ET32.nKcd = ET64.nKcd;

        ET32.nK2ab = ET64.nK2ab;
        ET32.nK2cd = ET64.nK2cd;
    }

    TheBatch.Ntiles32 += TheBatch.Ntiles64%2;

    /*
    for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {
        PackGeometries(&b4block.E8d[2*i], b4block.E8s[i]);
    }
    if (TheBatch.Ntiles64%2) PackGeometries(b4block.E8d[2*TheBatch.Ntiles32], b4block.E8s[TheBatch.Ntiles32]);

    TheBatch.Ntiles32 += TheBatch.Ntiles64%2;


    ERIgeometry        * ed = b4block.E0d;
    ERIgeometry        * es = b4block.E0s;

    for (UI32 i=0; i<TheBatch.Ntiles64*DPC; ++i) {
        es[i] = ed[i];
    }
    */


    /*
    //fuse cachelines
    for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {
        for (int p=0; p<TheBatch.wsize4; ++p) {
            b4block.T8s[(i) * TheBatch.wsize4 + p].set(
            b4block.T8d[(2*i) * TheBatch.wsize4 + p],
            b4block.T8d[(2*i+1) * TheBatch.wsize4 + p]);
        }
    }
    */

    /*
    for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {
        for (int p=0; p<TheBatch.msize4; ++p) {
            b4block.I8s[(i) * TheBatch.msize4 + p].set(
            b4block.I8d[(2*i) * TheBatch.msize4 + p],
            b4block.I8d[(2*i+1) * TheBatch.msize4 + p]);
        }
    }
    */


    TheBatch.Ntiles64 = 0;
}

void BatchEvaluator::KillLists(ERIBatch & TheBatch) {

    delete[] TheBatch.TileList32;
    TheBatch.TileList32 = NULL;
}


void jxd64NULL(cacheline64 * F, const cacheline64 * T, const cacheline64 * D) {}
void jxd32NULL(cacheline32 * F, const cacheline32 * T, const cacheline32 * D) {}

#include <string.h>

template<bool J, bool X>
void BatchEvaluator::Evaluate_D      (ERIBatch & TheBatch, const Sparse & Ds) {


    const UI8 JA = TheBatch.ABp->Ja;
    const UI8 JB = TheBatch.ABp->Jb;
    const UI8 JC = TheBatch.CDp->Ja;
    const UI8 JD = TheBatch.CDp->Jb;

    const UI8 la = TheBatch.la;
    const UI8 lb = TheBatch.lb;
    const UI8 lc = TheBatch.lc;
    const UI8 ld = TheBatch.ld;

    const UI8 MA = 2*TheBatch.la+1;
    const UI8 MB = 2*TheBatch.lb+1;
    const UI8 MC = 2*TheBatch.lc+1;
    const UI8 MD = 2*TheBatch.ld+1;

    int sAB = JA*JB*MA*MB;
    int sCD = JC*JD*MC*MD;
    int sAC = JA*JC*MA*MC;
    int sAD = JA*JD*MA*MD;
    int sBC = JB*JC*MB*MC;
    int sBD = JB*JD*MB*MD;

    //CONTRACT WITH DENSITY MATRIX

    chronoC.Start();
    {
        JXD64 Jabcd, Jcdab, Xacbd, Xadbc, Xbdac, Xbcad;

        JXrot64 Rab,  Rac,  Rad,  Rbc,  Rbd,  Rcd;
        JXrot64 RTab, RTac, RTad, RTbc, RTbd, RTcd;

        // select functions
        {
            bool AB = !TheBatch.ABp->samef;
            bool CD = !TheBatch.CDp->samef;


            if (AB && CD) {
                Jabcd = Echidna::dmd.J64abcd[la][lb][lc][ld];
                Jcdab = Echidna::dmd.J64cdab[la][lb][lc][ld];
                Xacbd = Echidna::dmd.X64acbd[la][lb][lc][ld];
                Xadbc = Echidna::dmd.X64adbc[la][lb][lc][ld];
                Xbdac = Echidna::dmd.X64bdac[la][lb][lc][ld];
                Xbcad = Echidna::dmd.X64bcad[la][lb][lc][ld];
            }
            else if (AB && !CD) {
                Jabcd = Echidna::dmd.J64abcc[la][lb] [lc];
                Jcdab = Echidna::dmd.J64ccab[la][lb] [lc];
                Xacbd = Echidna::dmd.X64acbc[la][lb] [lc];
                Xadbc = jxd64NULL;
                Xbdac = Echidna::dmd.X64bcac[la][lb] [lc];
                Xbcad = jxd64NULL;
            }
            else if (!AB && CD) {
                Jabcd = Echidna::dmd.J64aacd[la] [lc][ld];
                Jcdab = Echidna::dmd.J64cdaa[la] [lc][ld];
                Xacbd = Echidna::dmd.X64acad[la] [lc][ld];
                Xadbc = Echidna::dmd.X64adac[la] [lc][ld];
                Xbdac = jxd64NULL;
                Xbcad = jxd64NULL;
            }
            else if (!AB && !CD) {
                Jabcd = Echidna::dmd.J64aacc[la] [lc];
                Jcdab = Echidna::dmd.J64ccaa[la] [lc];
                Xacbd = Echidna::dmd.X64acac[la] [lc];
                Xadbc = jxd64NULL;
                Xbdac = jxd64NULL;
                Xbcad = jxd64NULL;
            }

            Rab = Echidna::dmd.R64[la][lb];
            Rac = Echidna::dmd.R64[la][lc];
            Rad = Echidna::dmd.R64[la][ld];
            Rbc = Echidna::dmd.R64[lb][lc];
            Rbd = Echidna::dmd.R64[lb][ld];
            Rcd = Echidna::dmd.R64[lc][ld];

            RTab = Echidna::dmd.RT64[la][lb];
            RTac = Echidna::dmd.RT64[la][lc];
            RTad = Echidna::dmd.RT64[la][ld];
            RTbc = Echidna::dmd.RT64[lb][lc];
            RTbd = Echidna::dmd.RT64[lb][ld];
            RTcd = Echidna::dmd.RT64[lc][ld];
        }


        cacheline64        * t8  = b4block.T8d;
        RotationMatrices64 * rm8 = b4block.RM8d;

        cacheline64 *Dab, *Dac, *Dad, *Dbc, *Dbd, *Dcd;
        cacheline64 *Fab, *Fac, *Fad, *Fbc, *Fbd, *Fcd;

        Dab = b4block.Dab;
        Dcd = b4block.Dcd;
        Dac = b4block.Dac;
        Dad = b4block.Dad;
        Dbc = b4block.Dbc;
        Dbd = b4block.Dbd;

        Fab = b4block.Fab;
        Fcd = b4block.Fcd;
        Fac = b4block.Fac;
        Fad = b4block.Fad;
        Fbc = b4block.Fbc;
        Fbd = b4block.Fbd;


        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {

            RotationMatrices64 & RMV = rm8[i];

            //compute offsets, digest
            UI32 offAB = TheBatch.TileList64[i].offAB;
            UI32 offCD = TheBatch.TileList64[i].offCD;

            UI32 offAC = TheBatch.TileList64[i].offAC;
            UI32 offAD = TheBatch.TileList64[i].offAD;
            UI32 offBC = TheBatch.TileList64[i].offBC;
            UI32 offBD = TheBatch.TileList64[i].offBD;

            if (J) {
                memset(Fab, 0, sAB*sizeof(cacheline64));
                memset(Fcd, 0, sCD*sizeof(cacheline64));
            }

            if (X) {
                memset(Fac, 0, sAC*sizeof(cacheline64));
                memset(Fad, 0, sAD*sizeof(cacheline64));
                memset(Fbc, 0, sBC*sizeof(cacheline64));
                memset(Fbd, 0, sBD*sizeof(cacheline64));
            }

            // set density matrix blocks
            for (UI32 k=0; k<DPC; ++k) {
                if (TheBatch.TileList64[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k];

                UI32 ta = TheBatch.SP12[ap12].ata;
                UI32 tb = TheBatch.SP12[ap12].atb;
                UI32 tc = TheBatch.SP34[ap34].ata;
                UI32 td = TheBatch.SP34[ap34].atb;

                if (J) {
                    const double * dDab = Ds(ta,tb) + offAB;
                    const double * dDcd = Ds(tc,td) + offCD;
                    for (int mm=0; mm<sAB; ++mm) Dab[mm](k) = dDab[mm];
                    for (int mm=0; mm<sCD; ++mm) Dcd[mm](k) = dDcd[mm];
                }

                if (X) {
                    const double * dDac = Ds(ta,tc) + offAC;
                    const double * dDad = Ds(ta,td) + offAD;
                    const double * dDbc = Ds(tb,tc) + offBC;
                    const double * dDbd = Ds(tb,td) + offBD;
                    for (int mm=0; mm<sAC; ++mm) Dac[mm](k) = dDac[mm];
                    for (int mm=0; mm<sAD; ++mm) Dad[mm](k) = dDad[mm];
                    for (int mm=0; mm<sBC; ++mm) Dbc[mm](k) = dDbc[mm];
                    for (int mm=0; mm<sBD; ++mm) Dbd[mm](k) = dDbd[mm];
                }

            }

            // rotate density matrices
            {
                if (J) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jb=0; jb<JB; ++jb) {
                            UI32 oAB = (ja*JB + jb)*(MA*MB);

                            Rab(Dab + oAB, &RMV);
                        }
                    }

                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            Rcd(Dcd + oCD, &RMV);
                        }
                    }
                }

                if (X) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oAC = (ja*JC + jc)*(MA*MC);

                            Rac(Dac + oAC, &RMV);
                        }
                    }

                    for (int ja=0; ja<JA; ++ja) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oAD = (ja*JD + jd)*(MA*MD);

                            Rad(Dad + oAD, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oBC = (jb*JC + jc)*(MB*MC);

                            Rbc(Dbc + oBC, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            Rbd(Dbd + oBD, &RMV);
                        }
                    }
                }

            }

            // digestion (contraction)
            // ***********************

            for (int ja=0; ja<JA; ++ja) {
                for (int jb=0; jb<JB; ++jb) {
                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {

                            UI32 oAB = (ja*JB + jb)*(MA*MB);
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            UI32 oAC = (ja*JC + jc)*(MA*MC);
                            UI32 oAD = (ja*JD + jd)*(MA*MD);
                            UI32 oBC = (jb*JC + jc)*(MB*MC);
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            if (J) {
                                Jabcd(Fab + oAB, t8, Dcd + oCD);
                                Jcdab(Fcd + oCD, t8, Dab + oAB);
                            }

                            if (X) {
                                Xacbd(Fac + oAC, t8, Dbd + oBD);
                                Xadbc(Fad + oAD, t8, Dbc + oBC);
                                Xbdac(Fbd + oBD, t8, Dac + oAC);
                                Xbcad(Fbc + oBC, t8, Dad + oAD);
                            }

                            t8 += TheBatch.wsize;
                        }
                    }
                }
            }

            // rotate fock matrices
            {
                if (J) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jb=0; jb<JB; ++jb) {
                            UI32 oAB = (ja*JB + jb)*(MA*MB);

                            RTab(Fab + oAB, &RMV);
                        }
                    }

                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            RTcd(Fcd + oCD, &RMV);
                        }
                    }
                }

                if (X) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oAC = (ja*JC + jc)*(MA*MC);

                            RTac(Fac + oAC, &RMV);
                        }
                    }

                    for (int ja=0; ja<JA; ++ja) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oAD = (ja*JD + jd)*(MA*MD);

                            RTad(Fad + oAD, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oBC = (jb*JC + jc)*(MB*MC);

                            RTbc(Fbc + oBC, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            RTbd(Fbd + oBD, &RMV);
                        }
                    }
                }

            }

            // update fock matrix
            for (UI32 k=0; k<DPC; ++k) {
                if (TheBatch.TileList64[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k];

                UI32 ta = TheBatch.SP12[ap12].ata;
                UI32 tb = TheBatch.SP12[ap12].atb;
                UI32 tc = TheBatch.SP34[ap34].ata;
                UI32 td = TheBatch.SP34[ap34].atb;

                //coulomb
                if (J) {
                    double * dFab = Jcont(ta,tb) + offAB;
                    double * dFcd = Jcont(tc,td) + offCD;
                    for (int mm=0; mm<sAB; ++mm) dFab[mm] += Fab[mm](k);
                    for (int mm=0; mm<sCD; ++mm) dFcd[mm] += Fcd[mm](k);
                }

                //exchange
                if (X) {
                    double * dFac = Xcont(ta,tc) + offAC;
                    double * dFad = Xcont(ta,td) + offAD;
                    double * dFbc = Xcont(tb,tc) + offBC;
                    double * dFbd = Xcont(tb,td) + offBD;
                    for (int mm=0; mm<sAC; ++mm) dFac[mm] -= Fac[mm](k);
                    for (int mm=0; mm<sAD; ++mm) dFad[mm] -= Fad[mm](k);
                    for (int mm=0; mm<sBC; ++mm) dFbc[mm] -= Fbc[mm](k);
                    for (int mm=0; mm<sBD; ++mm) dFbd[mm] -= Fbd[mm](k);
                }

            }

        }
    }

    {
        JXD32 Jabcd, Jcdab, Xacbd, Xadbc, Xbdac, Xbcad;

        JXrot32 Rab,  Rac,  Rad,  Rbc,  Rbd,  Rcd;
        JXrot32 RTab, RTac, RTad, RTbc, RTbd, RTcd;

        // select functions
        {
            bool AB = !TheBatch.ABp->samef;
            bool CD = !TheBatch.CDp->samef;


            if (AB && CD) {
                Jabcd = Echidna::dmd.J32abcd[la][lb][lc][ld];
                Jcdab = Echidna::dmd.J32cdab[la][lb][lc][ld];
                Xacbd = Echidna::dmd.X32acbd[la][lb][lc][ld];
                Xadbc = Echidna::dmd.X32adbc[la][lb][lc][ld];
                Xbdac = Echidna::dmd.X32bdac[la][lb][lc][ld];
                Xbcad = Echidna::dmd.X32bcad[la][lb][lc][ld];
            }
            else if (AB && !CD) {
                Jabcd = Echidna::dmd.J32abcc[la][lb] [lc];
                Jcdab = Echidna::dmd.J32ccab[la][lb] [lc];
                Xacbd = Echidna::dmd.X32acbc[la][lb] [lc];
                Xadbc = jxd32NULL;
                Xbdac = Echidna::dmd.X32bcac[la][lb] [lc];
                Xbcad = jxd32NULL;
            }
            else if (!AB && CD) {
                Jabcd = Echidna::dmd.J32aacd[la] [lc][ld];
                Jcdab = Echidna::dmd.J32cdaa[la] [lc][ld];
                Xacbd = Echidna::dmd.X32acad[la] [lc][ld];
                Xadbc = Echidna::dmd.X32adac[la] [lc][ld];
                Xbdac = jxd32NULL;
                Xbcad = jxd32NULL;
            }
            else if (!AB && !CD) {
                Jabcd = Echidna::dmd.J32aacc[la] [lc];
                Jcdab = Echidna::dmd.J32ccaa[la] [lc];
                Xacbd = Echidna::dmd.X32acac[la] [lc];
                Xadbc = jxd32NULL;
                Xbdac = jxd32NULL;
                Xbcad = jxd32NULL;
            }

            Rab = Echidna::dmd.R32[la][lb];
            Rac = Echidna::dmd.R32[la][lc];
            Rad = Echidna::dmd.R32[la][ld];
            Rbc = Echidna::dmd.R32[lb][lc];
            Rbd = Echidna::dmd.R32[lb][ld];
            Rcd = Echidna::dmd.R32[lc][ld];

            RTab = Echidna::dmd.RT32[la][lb];
            RTac = Echidna::dmd.RT32[la][lc];
            RTad = Echidna::dmd.RT32[la][ld];
            RTbc = Echidna::dmd.RT32[lb][lc];
            RTbd = Echidna::dmd.RT32[lb][ld];
            RTcd = Echidna::dmd.RT32[lc][ld];
        }


        cacheline32        * t8  = b4block.T8s;
        RotationMatrices32 * rm8 = b4block.RM8s;

        cacheline32 *Dab, *Dac, *Dad, *Dbc, *Dbd, *Dcd;
        cacheline32 *Fab, *Fac, *Fad, *Fbc, *Fbd, *Fcd;

        Dab = (cacheline32*)b4block.Dab;
        Dcd = (cacheline32*)b4block.Dcd;
        Dac = (cacheline32*)b4block.Dac;
        Dad = (cacheline32*)b4block.Dad;
        Dbc = (cacheline32*)b4block.Dbc;
        Dbd = (cacheline32*)b4block.Dbd;

        Fab = (cacheline32*)b4block.Fab;
        Fcd = (cacheline32*)b4block.Fcd;
        Fac = (cacheline32*)b4block.Fac;
        Fad = (cacheline32*)b4block.Fad;
        Fbc = (cacheline32*)b4block.Fbc;
        Fbd = (cacheline32*)b4block.Fbd;


        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

            RotationMatrices32 & RMV = rm8[i];

            //compute offsets, digest
            UI32 offAB = TheBatch.TileList32[i].offAB;
            UI32 offCD = TheBatch.TileList32[i].offCD;

            UI32 offAC = TheBatch.TileList32[i].offAC;
            UI32 offAD = TheBatch.TileList32[i].offAD;
            UI32 offBC = TheBatch.TileList32[i].offBC;
            UI32 offBD = TheBatch.TileList32[i].offBD;

            if (J) {
                memset(Fab, 0, sAB*sizeof(cacheline32));
                memset(Fcd, 0, sCD*sizeof(cacheline32));
            }

            if (X) {
                memset(Fac, 0, sAC*sizeof(cacheline32));
                memset(Fad, 0, sAD*sizeof(cacheline32));
                memset(Fbc, 0, sBC*sizeof(cacheline32));
                memset(Fbd, 0, sBD*sizeof(cacheline32));
            }

            // set density matrix blocks
            for (UI32 k=0; k<FPC; ++k) {
                if (TheBatch.TileList32[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k];

                UI32 ta = TheBatch.SP12[ap12].ata;
                UI32 tb = TheBatch.SP12[ap12].atb;
                UI32 tc = TheBatch.SP34[ap34].ata;
                UI32 td = TheBatch.SP34[ap34].atb;

                if (J) {
                    const double * dDab = Ds(ta,tb) + offAB;
                    const double * dDcd = Ds(tc,td) + offCD;
                    for (int mm=0; mm<sAB; ++mm) Dab[mm](k) = dDab[mm];
                    for (int mm=0; mm<sCD; ++mm) Dcd[mm](k) = dDcd[mm];
                }

                if (X) {
                    const double * dDac = Ds(ta,tc) + offAC;
                    const double * dDad = Ds(ta,td) + offAD;
                    const double * dDbc = Ds(tb,tc) + offBC;
                    const double * dDbd = Ds(tb,td) + offBD;
                    for (int mm=0; mm<sAC; ++mm) Dac[mm](k) = dDac[mm];
                    for (int mm=0; mm<sAD; ++mm) Dad[mm](k) = dDad[mm];
                    for (int mm=0; mm<sBC; ++mm) Dbc[mm](k) = dDbc[mm];
                    for (int mm=0; mm<sBD; ++mm) Dbd[mm](k) = dDbd[mm];
                }

            }

            // rotate density matrices
            {
                if (J) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jb=0; jb<JB; ++jb) {
                            UI32 oAB = (ja*JB + jb)*(MA*MB);

                            Rab(Dab + oAB, &RMV);
                        }
                    }

                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            Rcd(Dcd + oCD, &RMV);
                        }
                    }
                }

                if (X) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oAC = (ja*JC + jc)*(MA*MC);

                            Rac(Dac + oAC, &RMV);
                        }
                    }

                    for (int ja=0; ja<JA; ++ja) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oAD = (ja*JD + jd)*(MA*MD);

                            Rad(Dad + oAD, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oBC = (jb*JC + jc)*(MB*MC);

                            Rbc(Dbc + oBC, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            Rbd(Dbd + oBD, &RMV);
                        }
                    }
                }

            }

            // digestion (contraction)
            // ***********************

            for (int ja=0; ja<JA; ++ja) {
                for (int jb=0; jb<JB; ++jb) {
                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {

                            UI32 oAB = (ja*JB + jb)*(MA*MB);
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            UI32 oAC = (ja*JC + jc)*(MA*MC);
                            UI32 oAD = (ja*JD + jd)*(MA*MD);
                            UI32 oBC = (jb*JC + jc)*(MB*MC);
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            if (J) {
                                Jabcd(Fab + oAB, t8, Dcd + oCD);
                                Jcdab(Fcd + oCD, t8, Dab + oAB);
                            }

                            if (X) {
                                Xacbd(Fac + oAC, t8, Dbd + oBD);
                                Xadbc(Fad + oAD, t8, Dbc + oBC);
                                Xbdac(Fbd + oBD, t8, Dac + oAC);
                                Xbcad(Fbc + oBC, t8, Dad + oAD);
                            }

                            t8 += TheBatch.wsize;
                        }
                    }
                }
            }

            // rotate fock matrices
            {
                if (J) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jb=0; jb<JB; ++jb) {
                            UI32 oAB = (ja*JB + jb)*(MA*MB);

                            RTab(Fab + oAB, &RMV);
                        }
                    }

                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            RTcd(Fcd + oCD, &RMV);
                        }
                    }
                }

                if (X) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oAC = (ja*JC + jc)*(MA*MC);

                            RTac(Fac + oAC, &RMV);
                        }
                    }

                    for (int ja=0; ja<JA; ++ja) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oAD = (ja*JD + jd)*(MA*MD);

                            RTad(Fad + oAD, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oBC = (jb*JC + jc)*(MB*MC);

                            RTbc(Fbc + oBC, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            RTbd(Fbd + oBD, &RMV);
                        }
                    }
                }

            }

            // update fock matrix
            for (UI32 k=0; k<FPC; ++k) {
                if (TheBatch.TileList32[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k];

                UI32 ta = TheBatch.SP12[ap12].ata;
                UI32 tb = TheBatch.SP12[ap12].atb;
                UI32 tc = TheBatch.SP34[ap34].ata;
                UI32 td = TheBatch.SP34[ap34].atb;

                //coulomb
                if (J) {
                    double * dFab = Jcont(ta,tb) + offAB;
                    double * dFcd = Jcont(tc,td) + offCD;
                    for (int mm=0; mm<sAB; ++mm) dFab[mm] += Fab[mm](k);
                    for (int mm=0; mm<sCD; ++mm) dFcd[mm] += Fcd[mm](k);
                }

                //exchange
                if (X) {
                    double * dFac = Xcont(ta,tc) + offAC;
                    double * dFad = Xcont(ta,td) + offAD;
                    double * dFbc = Xcont(tb,tc) + offBC;
                    double * dFbd = Xcont(tb,td) + offBD;
                    for (int mm=0; mm<sAC; ++mm) dFac[mm] -= Fac[mm](k);
                    for (int mm=0; mm<sAD; ++mm) dFad[mm] -= Fad[mm](k);
                    for (int mm=0; mm<sBC; ++mm) dFbc[mm] -= Fbc[mm](k);
                    for (int mm=0; mm<sBD; ++mm) dFbd[mm] -= Fbd[mm](k);
                }

            }

        }
    }



    chronoC.Stop();

}


template<bool J, bool X>
void BatchEvaluator::Evaluate_D_sameF(ERIBatch & TheBatch, const Sparse & Ds) {

    const UI8 JA = TheBatch.ABp->Ja;
    const UI8 JB = TheBatch.ABp->Jb;
    const UI8 JC = TheBatch.CDp->Ja;
    const UI8 JD = TheBatch.CDp->Jb;

    const UI8 la = TheBatch.la;
    const UI8 lb = TheBatch.lb;
    const UI8 lc = TheBatch.lc;
    const UI8 ld = TheBatch.ld;

    const UI8 MA = 2*TheBatch.la+1;
    const UI8 MB = 2*TheBatch.lb+1;
    const UI8 MC = 2*TheBatch.lc+1;
    const UI8 MD = 2*TheBatch.ld+1;

    int sAB = JA*JB*MA*MB;
    int sCD = JC*JD*MC*MD;
    int sAC = JA*JC*MA*MC;
    int sAD = JA*JD*MA*MD;
    int sBC = JB*JC*MB*MC;
    int sBD = JB*JD*MB*MD;


    //CONTRACT WITH DENSITY MATRIX
    chronoC.Start();
    {
        JXD64 Jabcd, Jcdab, Xacbd, Xadbc, Xbdac, Xbcad;
        JXD64 Jabab, Xaabb, Xabba, Xbbaa;

        JXrot64 Rab,  Rac,  Rad,  Rbc,  Rbd,  Rcd;
        JXrot64 RTab, RTac, RTad, RTbc, RTbd, RTcd;

        // select functions
        {
            {
                Jabab = Echidna::dmd.J64abab[la][lb];
                Xaabb = Echidna::dmd.X64aabb[la][lb];
                Xabba = Echidna::dmd.X64abba[la][lb];
                Xbbaa = Echidna::dmd.X64bbaa[la][lb];
            }

            {
                Jabcd = Echidna::dmd.J64abcd[la][lb][lc][ld];
                Jcdab = Echidna::dmd.J64cdab[la][lb][lc][ld];
                Xacbd = Echidna::dmd.X64acbd[la][lb][lc][ld];
                Xadbc = Echidna::dmd.X64adbc[la][lb][lc][ld];
                Xbdac = Echidna::dmd.X64bdac[la][lb][lc][ld];
                Xbcad = Echidna::dmd.X64bcad[la][lb][lc][ld];
            }

            Rab = Echidna::dmd.R64[la][lb];
            Rac = Echidna::dmd.R64[la][lc];
            Rad = Echidna::dmd.R64[la][ld];
            Rbc = Echidna::dmd.R64[lb][lc];
            Rbd = Echidna::dmd.R64[lb][ld];
            Rcd = Echidna::dmd.R64[lc][ld];

            RTab = Echidna::dmd.RT64[la][lb];
            RTac = Echidna::dmd.RT64[la][lc];
            RTad = Echidna::dmd.RT64[la][ld];
            RTbc = Echidna::dmd.RT64[lb][lc];
            RTbd = Echidna::dmd.RT64[lb][ld];
            RTcd = Echidna::dmd.RT64[lc][ld];
        }


        cacheline64        * t8  = b4block.T8d;
        RotationMatrices64 * rm8 = b4block.RM8d;

        cacheline64 *Dab, *Dac, *Dad, *Dbc, *Dbd, *Dcd;
        cacheline64 *Fab, *Fac, *Fad, *Fbc, *Fbd, *Fcd;

        Dab = b4block.Dab;
        Dcd = b4block.Dcd;
        Dac = b4block.Dac;
        Dad = b4block.Dad;
        Dbc = b4block.Dbc;
        Dbd = b4block.Dbd;

        Fab = b4block.Fab;
        Fcd = b4block.Fcd;
        Fac = b4block.Fac;
        Fad = b4block.Fad;
        Fbc = b4block.Fbc;
        Fbd = b4block.Fbd;


        for (UI32 i=0; i<TheBatch.Ntiles64; ++i) {

            RotationMatrices64 & RMV = rm8[i];

            //compute offsets, digest
            UI32 offAB = TheBatch.TileList64[i].offAB;
            UI32 offCD = TheBatch.TileList64[i].offCD;

            UI32 offAC = TheBatch.TileList64[i].offAC;
            UI32 offAD = TheBatch.TileList64[i].offAD;
            UI32 offBC = TheBatch.TileList64[i].offBC;
            UI32 offBD = TheBatch.TileList64[i].offBD;

            if (J) {
                memset(Fab, 0, sAB*sizeof(cacheline64));
                memset(Fcd, 0, sCD*sizeof(cacheline64));
            }

            if (X) {
                memset(Fac, 0, sAC*sizeof(cacheline64));
                memset(Fad, 0, sAD*sizeof(cacheline64));
                memset(Fbc, 0, sBC*sizeof(cacheline64));
                memset(Fbd, 0, sBD*sizeof(cacheline64));
            }

            // set density matrix blocks
            for (UI32 k=0; k<DPC; ++k) {
                if (TheBatch.TileList64[i].use[k] == 0) continue;
                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k];

                UI32 ta = TheBatch.SP12[ap12].ata;
                UI32 tb = TheBatch.SP12[ap12].atb;
                UI32 tc = TheBatch.SP34[ap34].ata;
                UI32 td = TheBatch.SP34[ap34].atb;

                if (J) {
                    const double * dDab = Ds(ta,tb) + offAB;
                    const double * dDcd = Ds(tc,td) + offCD;
                    for (int mm=0; mm<sAB; ++mm) Dab[mm](k) = dDab[mm];
                    for (int mm=0; mm<sCD; ++mm) Dcd[mm](k) = dDcd[mm];
                }

                if (X) {
                    const double * dDac = Ds(ta,tc) + offAC;
                    const double * dDad = Ds(ta,td) + offAD;
                    const double * dDbc = Ds(tb,tc) + offBC;
                    const double * dDbd = Ds(tb,td) + offBD;
                    for (int mm=0; mm<sAC; ++mm) Dac[mm](k) = dDac[mm];
                    for (int mm=0; mm<sAD; ++mm) Dad[mm](k) = dDad[mm];
                    for (int mm=0; mm<sBC; ++mm) Dbc[mm](k) = dDbc[mm];
                    for (int mm=0; mm<sBD; ++mm) Dbd[mm](k) = dDbd[mm];
                }

            }

            // rotate density matrices
            {
                if (J) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jb=0; jb<JB; ++jb) {
                            UI32 oAB = (ja*JB + jb)*(MA*MB);

                            Rab(Dab + oAB, &RMV);
                        }
                    }

                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            Rcd(Dcd + oCD, &RMV);
                        }
                    }
                }

                if (X) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oAC = (ja*JC + jc)*(MA*MC);

                            Rac(Dac + oAC, &RMV);
                        }
                    }

                    for (int ja=0; ja<JA; ++ja) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oAD = (ja*JD + jd)*(MA*MD);

                            Rad(Dad + oAD, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oBC = (jb*JC + jc)*(MB*MC);

                            Rbc(Dbc + oBC, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            Rbd(Dbd + oBD, &RMV);
                        }
                    }
                }

            }

            // digestion (contraction)
            // ***********************

            for (int ja=0; ja<JA; ++ja) {
                for (int jb=0; jb<JB; ++jb) {
                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {

                            UI32 oAB = (ja*JB + jb)*(MA*MB);
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            UI32 oAC = (ja*JC + jc)*(MA*MC);
                            UI32 oAD = (ja*JD + jd)*(MA*MD);
                            UI32 oBC = (jb*JC + jc)*(MB*MC);
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            if (ja*JB+jb > jc*JD+jd) {

                                if (J) {
                                    Jabcd(Fab + oAB, t8, Dcd + oCD);
                                    Jcdab(Fcd + oCD, t8, Dab + oAB);
                                }
                                if (X) {
                                    Xacbd(Fac + oAC, t8, Dbd + oBD);
                                    Xadbc(Fad + oAD, t8, Dbc + oBC);
                                    Xbdac(Fbd + oBD, t8, Dac + oAC);
                                    Xbcad(Fbc + oBC, t8, Dad + oAD);
                                }
                            }

                            else if (ja==jc && jb==jd) {

                                if (J) {
                                    Jabab(Fab + oAB, t8, Dcd + oCD);
                                }

                                if (X) {
                                    Xaabb(Fac + oAC, t8, Dbd + oBD);
                                    Xbbaa(Fbd + oBD, t8, Dac + oAC);

                                    Xabba(Fad + oAD, t8, Dad + oAD);
                                }
                            }

                            t8 += TheBatch.wsize;
                        }
                    }
                }
            }

            // rotate fock matrices
            {
                if (J) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jb=0; jb<JB; ++jb) {
                            UI32 oAB = (ja*JB + jb)*(MA*MB);

                            RTab(Fab + oAB, &RMV);
                        }
                    }

                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            RTcd(Fcd + oCD, &RMV);
                        }
                    }
                }

                if (X) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oAC = (ja*JC + jc)*(MA*MC);

                            RTac(Fac + oAC, &RMV);
                        }
                    }

                    for (int ja=0; ja<JA; ++ja) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oAD = (ja*JD + jd)*(MA*MD);

                            RTad(Fad + oAD, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oBC = (jb*JC + jc)*(MB*MC);

                            RTbc(Fbc + oBC, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            RTbd(Fbd + oBD, &RMV);
                        }
                    }
                }

            }

            // update fock matrix
            for (UI32 k=0; k<DPC; ++k) {
                if (TheBatch.TileList64[i].use[k] == 0) continue;
                UI32 ap12 = TheBatch.TileList64[i].ap12[k];
                UI32 ap34 = TheBatch.TileList64[i].ap34[k];

                UI32 ta = TheBatch.SP12[ap12].ata;
                UI32 tb = TheBatch.SP12[ap12].atb;
                UI32 tc = TheBatch.SP34[ap34].ata;
                UI32 td = TheBatch.SP34[ap34].atb;

                //coulomb
                if (J) {
                    double * dFab = Jcont(ta,tb) + offAB;
                    double * dFcd = Jcont(tc,td) + offCD;
                    for (int mm=0; mm<sAB; ++mm) dFab[mm] += Fab[mm](k);
                    for (int mm=0; mm<sCD; ++mm) dFcd[mm] += Fcd[mm](k);
                }

                //exchange
                if (X) {
                    double * dFac = Xcont(ta,tc) + offAC;
                    double * dFad = Xcont(ta,td) + offAD;
                    double * dFbc = Xcont(tb,tc) + offBC;
                    double * dFbd = Xcont(tb,td) + offBD;
                    for (int mm=0; mm<sAC; ++mm) dFac[mm] -= Fac[mm](k);
                    for (int mm=0; mm<sAD; ++mm) dFad[mm] -= Fad[mm](k);
                    for (int mm=0; mm<sBC; ++mm) dFbc[mm] -= Fbc[mm](k);
                    for (int mm=0; mm<sBD; ++mm) dFbd[mm] -= Fbd[mm](k);
                }

            }

        }
    }

    {
        JXD32 Jabcd, Jcdab, Xacbd, Xadbc, Xbdac, Xbcad;
        JXD32 Jabab, Xaabb, Xabba, Xbbaa;

        JXrot32 Rab,  Rac,  Rad,  Rbc,  Rbd,  Rcd;
        JXrot32 RTab, RTac, RTad, RTbc, RTbd, RTcd;

        // select functions
        {
            {
                Jabab = Echidna::dmd.J32abab[la][lb];
                Xaabb = Echidna::dmd.X32aabb[la][lb];
                Xabba = Echidna::dmd.X32abba[la][lb];
                Xbbaa = Echidna::dmd.X32bbaa[la][lb];
            }

            {
                Jabcd = Echidna::dmd.J32abcd[la][lb][lc][ld];
                Jcdab = Echidna::dmd.J32cdab[la][lb][lc][ld];
                Xacbd = Echidna::dmd.X32acbd[la][lb][lc][ld];
                Xadbc = Echidna::dmd.X32adbc[la][lb][lc][ld];
                Xbdac = Echidna::dmd.X32bdac[la][lb][lc][ld];
                Xbcad = Echidna::dmd.X32bcad[la][lb][lc][ld];
            }

            Rab = Echidna::dmd.R32[la][lb];
            Rac = Echidna::dmd.R32[la][lc];
            Rad = Echidna::dmd.R32[la][ld];
            Rbc = Echidna::dmd.R32[lb][lc];
            Rbd = Echidna::dmd.R32[lb][ld];
            Rcd = Echidna::dmd.R32[lc][ld];

            RTab = Echidna::dmd.RT32[la][lb];
            RTac = Echidna::dmd.RT32[la][lc];
            RTad = Echidna::dmd.RT32[la][ld];
            RTbc = Echidna::dmd.RT32[lb][lc];
            RTbd = Echidna::dmd.RT32[lb][ld];
            RTcd = Echidna::dmd.RT32[lc][ld];
        }


        cacheline32        * t8  = b4block.T8s;
        RotationMatrices32 * rm8 = b4block.RM8s;

        cacheline32 *Dab, *Dac, *Dad, *Dbc, *Dbd, *Dcd;
        cacheline32 *Fab, *Fac, *Fad, *Fbc, *Fbd, *Fcd;

        Dab = (cacheline32*)b4block.Dab;
        Dcd = (cacheline32*)b4block.Dcd;
        Dac = (cacheline32*)b4block.Dac;
        Dad = (cacheline32*)b4block.Dad;
        Dbc = (cacheline32*)b4block.Dbc;
        Dbd = (cacheline32*)b4block.Dbd;

        Fab = (cacheline32*)b4block.Fab;
        Fcd = (cacheline32*)b4block.Fcd;
        Fac = (cacheline32*)b4block.Fac;
        Fad = (cacheline32*)b4block.Fad;
        Fbc = (cacheline32*)b4block.Fbc;
        Fbd = (cacheline32*)b4block.Fbd;


        for (UI32 i=0; i<TheBatch.Ntiles32; ++i) {

            RotationMatrices32 & RMV = rm8[i];

            //compute offsets, digest
            UI32 offAB = TheBatch.TileList32[i].offAB;
            UI32 offCD = TheBatch.TileList32[i].offCD;

            UI32 offAC = TheBatch.TileList32[i].offAC;
            UI32 offAD = TheBatch.TileList32[i].offAD;
            UI32 offBC = TheBatch.TileList32[i].offBC;
            UI32 offBD = TheBatch.TileList32[i].offBD;

            if (J) {
                memset(Fab, 0, sAB*sizeof(cacheline32));
                memset(Fcd, 0, sCD*sizeof(cacheline32));
            }

            if (X) {
                memset(Fac, 0, sAC*sizeof(cacheline32));
                memset(Fad, 0, sAD*sizeof(cacheline32));
                memset(Fbc, 0, sBC*sizeof(cacheline32));
                memset(Fbd, 0, sBD*sizeof(cacheline32));
            }

            // set density matrix blocks
            for (UI32 k=0; k<FPC; ++k) {
                if (TheBatch.TileList32[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k];

                UI32 ta = TheBatch.SP12[ap12].ata;
                UI32 tb = TheBatch.SP12[ap12].atb;
                UI32 tc = TheBatch.SP34[ap34].ata;
                UI32 td = TheBatch.SP34[ap34].atb;

                if (J) {
                    const double * dDab = Ds(ta,tb) + offAB;
                    const double * dDcd = Ds(tc,td) + offCD;
                    for (int mm=0; mm<sAB; ++mm) Dab[mm](k) = dDab[mm];
                    for (int mm=0; mm<sCD; ++mm) Dcd[mm](k) = dDcd[mm];
                }

                if (X) {
                    const double * dDac = Ds(ta,tc) + offAC;
                    const double * dDad = Ds(ta,td) + offAD;
                    const double * dDbc = Ds(tb,tc) + offBC;
                    const double * dDbd = Ds(tb,td) + offBD;
                    for (int mm=0; mm<sAC; ++mm) Dac[mm](k) = dDac[mm];
                    for (int mm=0; mm<sAD; ++mm) Dad[mm](k) = dDad[mm];
                    for (int mm=0; mm<sBC; ++mm) Dbc[mm](k) = dDbc[mm];
                    for (int mm=0; mm<sBD; ++mm) Dbd[mm](k) = dDbd[mm];
                }

            }

            // rotate density matrices
            {
                if (J) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jb=0; jb<JB; ++jb) {
                            UI32 oAB = (ja*JB + jb)*(MA*MB);

                            Rab(Dab + oAB, &RMV);
                        }
                    }

                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            Rcd(Dcd + oCD, &RMV);
                        }
                    }
                }

                if (X) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oAC = (ja*JC + jc)*(MA*MC);

                            Rac(Dac + oAC, &RMV);
                        }
                    }

                    for (int ja=0; ja<JA; ++ja) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oAD = (ja*JD + jd)*(MA*MD);

                            Rad(Dad + oAD, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oBC = (jb*JC + jc)*(MB*MC);

                            Rbc(Dbc + oBC, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            Rbd(Dbd + oBD, &RMV);
                        }
                    }
                }

            }

            // digestion (contraction)
            // ***********************

            for (int ja=0; ja<JA; ++ja) {
                for (int jb=0; jb<JB; ++jb) {
                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {

                            UI32 oAB = (ja*JB + jb)*(MA*MB);
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            UI32 oAC = (ja*JC + jc)*(MA*MC);
                            UI32 oAD = (ja*JD + jd)*(MA*MD);
                            UI32 oBC = (jb*JC + jc)*(MB*MC);
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            if (ja*JB+jb > jc*JD+jd) {

                                if (J) {
                                    Jabcd(Fab + oAB, t8, Dcd + oCD);
                                    Jcdab(Fcd + oCD, t8, Dab + oAB);
                                }
                                if (X) {
                                    Xacbd(Fac + oAC, t8, Dbd + oBD);
                                    Xadbc(Fad + oAD, t8, Dbc + oBC);
                                    Xbdac(Fbd + oBD, t8, Dac + oAC);
                                    Xbcad(Fbc + oBC, t8, Dad + oAD);
                                }
                            }

                            else if (ja==jc && jb==jd) {

                                if (J) {
                                    Jabab(Fab + oAB, t8, Dcd + oCD);
                                }

                                if (X) {
                                    Xaabb(Fac + oAC, t8, Dbd + oBD);
                                    Xbbaa(Fbd + oBD, t8, Dac + oAC);

                                    Xabba(Fad + oAD, t8, Dad + oAD);
                                }
                            }

                            t8 += TheBatch.wsize;
                        }
                    }
                }
            }

            // rotate fock matrices
            {
                if (J) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jb=0; jb<JB; ++jb) {
                            UI32 oAB = (ja*JB + jb)*(MA*MB);

                            RTab(Fab + oAB, &RMV);
                        }
                    }

                    for (int jc=0; jc<JC; ++jc) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oCD = (jc*JD + jd)*(MC*MD);

                            RTcd(Fcd + oCD, &RMV);
                        }
                    }
                }

                if (X) {
                    for (int ja=0; ja<JA; ++ja) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oAC = (ja*JC + jc)*(MA*MC);

                            RTac(Fac + oAC, &RMV);
                        }
                    }

                    for (int ja=0; ja<JA; ++ja) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oAD = (ja*JD + jd)*(MA*MD);

                            RTad(Fad + oAD, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jc=0; jc<JC; ++jc) {
                            UI32 oBC = (jb*JC + jc)*(MB*MC);

                            RTbc(Fbc + oBC, &RMV);
                        }
                    }

                    for (int jb=0; jb<JB; ++jb) {
                        for (int jd=0; jd<JD; ++jd) {
                            UI32 oBD = (jb*JD + jd)*(MB*MD);

                            RTbd(Fbd + oBD, &RMV);
                        }
                    }
                }

            }

            // update fock matrix
            for (UI32 k=0; k<FPC; ++k) {
                if (TheBatch.TileList32[i].use[k] == 0) continue;

                UI32 ap12 = TheBatch.TileList32[i].ap12[k];
                UI32 ap34 = TheBatch.TileList32[i].ap34[k];

                UI32 ta = TheBatch.SP12[ap12].ata;
                UI32 tb = TheBatch.SP12[ap12].atb;
                UI32 tc = TheBatch.SP34[ap34].ata;
                UI32 td = TheBatch.SP34[ap34].atb;

                //coulomb
                if (J) {
                    double * dFab = Jcont(ta,tb) + offAB;
                    double * dFcd = Jcont(tc,td) + offCD;
                    for (int mm=0; mm<sAB; ++mm) dFab[mm] += Fab[mm](k);
                    for (int mm=0; mm<sCD; ++mm) dFcd[mm] += Fcd[mm](k);
                }

                //exchange
                if (X) {
                    double * dFac = Xcont(ta,tc) + offAC;
                    double * dFad = Xcont(ta,td) + offAD;
                    double * dFbc = Xcont(tb,tc) + offBC;
                    double * dFbd = Xcont(tb,td) + offBD;
                    for (int mm=0; mm<sAC; ++mm) dFac[mm] -= Fac[mm](k);
                    for (int mm=0; mm<sAD; ++mm) dFad[mm] -= Fad[mm](k);
                    for (int mm=0; mm<sBC; ++mm) dFbc[mm] -= Fbc[mm](k);
                    for (int mm=0; mm<sBD; ++mm) dFbd[mm] -= Fbd[mm](k);
                }
            }
        }
    }

    chronoC.Stop();

}




#include "../math/angular.hpp"
#include "../math/eigen.hpp"
#include "../math/tensors.hpp"
#include "../integrals/fock2e.hpp"

//find the spectral norm of the 2e-tensor block
cacheline64 BatchEvaluator::CalcCauchySchwarz(const ERIBatch & TheBatch, const ShellPair & ABs, const ERIgeometries64 & eriG8, const LibQuimera::ERITile64 & ET, UI8 lenk) {

    TheBatch.ERIalgorithm->K4(eriG8, ET, *TheBatch.ABp, *TheBatch.ABp, b4block.I8d, b4ERIs, false);

    cacheline64 * i8 = b4block.I8d;
    cacheline64 * t8 = b4block.T8d;

    const int JA = TheBatch.ABp->Ja;
    const int JB = TheBatch.ABp->Jb;

    for (int ja=0; ja<JA; ++ja) {
        for (int jb=0; jb<JB; ++jb) {
            for (int jc=0; jc<JA; ++jc) {
                for (int jd=0; jd<JB; ++jd) {
                    TheBatch.ERIalgorithm->MIRROR(i8, eriG8, t8, b4ERIs);
                    i8 += TheBatch.msize;
                    t8 += TheBatch.wsize;
                }
            }
        }
    }


    int n  = nmS[ABs.getLa()] * nmS[ABs.getLb()];
    int nn = JA * JB * n;

    cacheline64 maxv;


    //find the maximum eigenvalue
    tensor2 M;
    M.setsize(nn);

    double vv[nn];

    int nmLa = nmS[ABs.getLa()];
    int nmLb = nmS[ABs.getLb()];

    for (int t=0; t<lenk; ++t) {
        //fill the matrix

        cacheline64 * t8 = b4block.T8d;

        M.zeroize();

        for (int ja=0; ja<JA; ++ja) {
            for (int jb=0; jb<JB; ++jb) {
                for (int jc=0; jc<JA; ++jc) {
                    for (int jd=0; jd<JB; ++jd) {

                        int Jab = (ja*JB + jb) * nmLa * nmLb;
                        int Jcd = (jc*JB + jd) * nmLa * nmLb;

                        if (ja==jc && jb==jd)
                        for (int i=0; i<nmLa; ++i) {
                            for (int j=0; j<nmLb; ++j) {
                                for (int k=0; k<nmLa; ++k) {
                                    for (int l=0; l<nmLb; ++l) {

                                        int ij = i*nmLb + j;
                                        int kl = k*nmLb + l;

                                        int ijkl = (((l*nmLa+k)*nmLb+j)*nmLa+i);

                                        M(Jab+ij, Jcd+kl) = t8[ijkl](t);

                                    }
                                }
                            }
                        }

                        t8 += TheBatch.wsize;
                    }
                }
            }
        }

        DiagonalizeE(M, vv);

        double sum = 0;
        bool npd = false;
        for (int i=0; i<nn; ++i) {
            if (vv[i]<=-1e-12) npd = true;
            sum += vv[i];
        }

        if (npd) {
            Echidna::EDebugger.precision(16);
            Echidna::EDebugger << "Error! ABAB block ERI matrix is not positive-definite!!!" << std::endl;
        }
        Echidna::EDebugger << int(TheBatch.la) << " " << int(TheBatch.lb) << " " << int(TheBatch.lc) << " " << int(TheBatch.ld) << std::endl;
        Echidna::EDebugger << int(JA) << " " << int(JB) << std::endl;
        Echidna::EDebugger << t << std::endl;
        Echidna::EDebugger << " Eigenvalues: ";
        for (int i=0; i<nn; ++i) Echidna::EDebugger << " " << vv[i] << " "; Echidna::EDebugger << std::endl << std::endl;

        maxv(t) = sum;
    }

    return maxv;
}

double BatchEvaluator::CalcCauchySchwarz(const ERIBatch & TheBatch, const ShellPair & ABs, const ERIgeometry & eriG) {

    double * I0 = (double*)b4block.I8d;
    double * T0 = (double*)b4block.T8d;

    TheBatch.ERIalgorithm->K4(*TheBatch.ABp, *TheBatch.CDp, I0, b4ERIs, false);

    const int JA = TheBatch.ABp->Ja;
    const int JB = TheBatch.ABp->Jb;

    {
        double * i0 = I0;
        double * t0 = T0;

        for (int ja=0; ja<JA; ++ja) {
            for (int jb=0; jb<JB; ++jb) {
                for (int jc=0; jc<JA; ++jc) {
                    for (int jd=0; jd<JB; ++jd) {
                        TheBatch.ERIalgorithm->MIRROR(i0, t0, b4ERIs);
                        i0 += TheBatch.msize;
                        t0 += TheBatch.wsize;
                    }
                }
            }
        }
    }


    int n  = nmS[ABs.getLa()] * nmS[ABs.getLb()];
    int nn = JA * JB * n;
    int nmLa = nmS[ABs.getLa()];
    int nmLb = nmS[ABs.getLb()];

    double max;

    tensor2 M;
    M.setsize(nn);


    double vv[nn];

    double * t0 = T0;

    M.zeroize();


    for (int ja=0; ja<JA; ++ja) {
        for (int jb=0; jb<JB; ++jb) {
            for (int jc=0; jc<JA; ++jc) {
                for (int jd=0; jd<JB; ++jd) {

                    int Jab = (ja*JB + jb) * nmLa * nmLb;
                    int Jcd = (jc*JB + jd) * nmLa * nmLb;

                    if (ja==jc && jb==jd)
                    for (int i=0; i<nmLa; ++i) {
                        for (int j=0; j<nmLb; ++j) {
                            for (int k=0; k<nmLa; ++k) {
                                for (int l=0; l<nmLb; ++l) {

                                    int ij = i*nmLb + j;
                                    int kl = k*nmLb + l;

                                    int ijkl = (((l*nmLa+k)*nmLb+j)*nmLa+i);

                                    M(Jab+ij, Jcd+kl) = t0[ijkl];

                                }
                            }
                        }
                    }

                    t0 += TheBatch.wsize;
                }
            }
        }
    }

    DiagonalizeE(M, vv);

    double sum = 0;
    bool npd = false;
    for (int i=0; i<nn; ++i) {
        if (vv[i]<=-1e-12) npd = true;
        sum += vv[i];
    }

    if (npd) {
        Echidna::EDebugger.precision(16);
        Echidna::EDebugger << "Error! AAAA block ERI matrix is not positive-definite!!!" << std::endl;
    }
    Echidna::EDebugger << int(TheBatch.la) << " " << int(TheBatch.lb) << " " << int(TheBatch.lc) << " " << int(TheBatch.ld) << std::endl;
    Echidna::EDebugger << int(JA) << " " << int(JB) << std::endl;
    Echidna::EDebugger << " Eigenvalues: ";
    for (int i=0; i<nn; ++i) Echidna::EDebugger << " " << vv[i] << " "; Echidna::EDebugger << std::endl << std::endl;

    max = sum;

    return max;
}

void BatchEvaluator::ConstructW(const ERIBatch & TheBatch, const ShellPair & ABs, const ShellPair & CDs, const ERIgeometry & eriG, double * WWW, UI32 wa, bool OnlyJ) {

    double * I0 = (double*)b4block.I8d;
    double * T0 = (double*)b4block.T8d;

    TheBatch.ERIalgorithm->K4(*TheBatch.ABp, *TheBatch.CDp, I0, b4ERIs, OnlyJ);

    double * II = I0;
    double * TT = T0;

    const int JA = TheBatch.ABp->Ja;
    const int JB = TheBatch.ABp->Jb;
    const int JC = TheBatch.CDp->Ja;
    const int JD = TheBatch.CDp->Jb;


    for (int ja=0; ja<JA; ++ja) {
        for (int jb=0; jb<JB; ++jb) {
            for (int jc=0; jc<JC; ++jc) {
                for (int jd=0; jd<JD; ++jd) {
                    TheBatch.ERIalgorithm->MIRROR(II, TT, b4ERIs);

                    //cout << "II: ";for (int mm=0; mm<TheBatch.msize; ++mm) cout << " " << II[mm]; cout << endl;
                    //cout << "TT: ";for (int mm=0; mm<TheBatch.wsize; ++mm) cout << " " << TT[mm]; cout << endl;

                    II += TheBatch.msize;
                    TT += TheBatch.wsize;
                }
            }
        }
    }
    ConstructBigW(*TheBatch.ABp, *TheBatch.CDp, wa, T0, WWW);
}

