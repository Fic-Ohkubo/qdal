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


#ifndef __ERIblocks__
#define __ERIblocks__
#include <set>
#include <map>
#include <list>

#include "../defs.hpp"
#include "../basis/SPprototype.hpp"
#include "../integrals/atomprod.hpp"
#include "../integrals/fock2e.hpp"
#include "../libquimera/libquimera.hpp"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif



class ShellPair;
class AtomProdPrototypes;
class APlists;
class symtensor;
class Interactions;
class ERIgeometry;
class ERIBatch;
class BlockBuffer;
class ERIbatchBuffer;


// CONCURRENT, MULTISTEP PRODUCERS/CONSUMERS PROBLEM

template<class T> class TaskQueue {
  private:
    std::deque<T*> TaskFIFOQueue;
    #ifdef _OPENMP
        omp_lock_t Lock;
    #endif

  public:

    TaskQueue() {
        #ifdef _OPENMP
        omp_init_lock (&Lock);
        #endif
    }

    ~TaskQueue() {
        #ifdef _OPENMP
        omp_destroy_lock (&Lock);
        #endif
    }

    void push(T * rhs) {
        #ifdef _OPENMP
        omp_set_lock (&Lock);
        #endif
        {
            TaskFIFOQueue.push_back(rhs);
        }
        #ifdef _OPENMP
        omp_unset_lock (&Lock);
        #endif
    }

    void push_top(T * rhs) {
        #ifdef _OPENMP
        omp_set_lock (&Lock);
        #endif
        {
            TaskFIFOQueue.push_front(rhs);
        }
        #ifdef _OPENMP
        omp_unset_lock (&Lock);
        #endif
    }

    T * pop() {
        T * ret;
        #ifdef _OPENMP
        omp_set_lock (&Lock);
        #endif
        {
            if (TaskFIFOQueue.empty())
                ret = NULL;
            else {
                ret = TaskFIFOQueue.front();
                TaskFIFOQueue.pop_front();
            }
        }
        #ifdef _OPENMP
        omp_unset_lock (&Lock);
        #endif

        return ret;
    }

    int length() {
        return TaskFIFOQueue.size();
    }
};

//semi-useless classes
//********************

class CasePrescreener {
    public:

    static const int NLogsShort = 1024;
    double LogShort[NLogsShort];
    double maxR2;


    double LogShortCoulomb (double R2);
    void   InitLogsShort ();

    double LogsShort (double R2);
};

class SPdata {
    public:

    point c;
    UI32  n; //number in the SPlist
};

extern CasePrescreener CPEB;

//actual classes
//**************


struct ShellPairPair {
    UI32 n12;
    UI32 n34;

    ShellPairPair() {};

    ShellPairPair(UI32 a, UI32 b) {
        n12 = a;
        n34 = b;
    }
};

//to classify shellpairs according to total contraction degree
struct SPlist {
    UI32  Nsps;
    UI32  Isps;
    UI32  vList;
    UI32  rList; //not used for now

    SPlist() {
        Nsps  = Isps = 0;
        vList = 0;
        rList = 0;
    }
};

struct SPPblock {
    UI64            Nints; //LOL @ UI16 bug
    UI64            Iints;
    ShellPairPair * vList;


    SPPblock() {
        Nints = Iints = 0;
        vList = NULL;
    }

    void reset() {
        Nints = Iints = 0;
        vList = NULL;
    }

    ShellPairPair pop() {
        if (Iints>0) {
            --Iints;
            return vList[Iints];
            //vList[Iints].n12 = vList[Iints].n34 = 0; //unnecessary
        }
        else {

        }
    }

    const ShellPairPair & top() const {
        return vList[Iints-1];
    }

    void push(const ShellPairPair & rhs) {
        vList[Iints] = rhs;
        ++Iints;
    }

    int size() const {
        return Iints;
    }

    int capacity() const {
        return Nints;
    }

    void reserve() {
        ++Nints;
    }

};

struct SPJ {
    float cs;
    UI32   n;

    inline bool operator< (const SPJ & rhs) const {return (cs<rhs.cs);}
    inline bool operator==(const SPJ & rhs) const {return (cs<rhs.cs);}
    inline bool operator> (const SPJ & rhs) const {return (cs>rhs.cs);}
    inline bool operator>=(const SPJ & rhs) const {return (cs>=rhs.cs);}
    inline bool operator<=(const SPJ & rhs) const {return (cs<=rhs.cs);}
};

//pointers to functions and t_sizes uniquely dependant on geometry and shell pair 'prototypes'
class BatchInfo {
    friend class BlockBuffer;
    friend class ERIblock;
    friend class Fock2e;

  protected:
    const ShellPairPrototype * ABp;
    const ShellPairPrototype * CDp;

    bool SkipSame;
    bool SameShell;
    bool SameList;

    GEOM geometry;

    UI8 la;
    UI8 lb;
    UI8 lc;
    UI8 ld;
    UI8 Lt;
    UI8 Lmax;

    //size of the ERI batch
    UI16 J4;
    UI32 fsize;

    UI16 msize; //size of the kernel list
    UI32 msize4; //size of the kernel list (times number of ERI batches)

    UI16 wsize;
    UI32 wsize4;

    //funciones necesarias
    //********************
    const LibQuimera::Qalgorithm * ERIalgorithm;


    //window arguments
    //****************

    //atom lists
    const Array<AtomProd> * AP12list;
    const Array<AtomProd> * AP34list;

    const ShellPair * SP12;
    const ShellPair * SP34;

    UI32 tmax12;
    UI32 tmax34;


    //to resume
    //*********
    struct {
        UI32 storeAB;
        UI32 storeCD;
        UI16 storeABb;
        UI16 storeCDb;
        bool storeBlock;
        bool store3center;
    } State;


//functions to set geometry, etc.
  public:
    void Set(GEOM geom, const ShellPairPrototype & AB, const ShellPairPrototype & CD, bool samef);
    void SetLists(const ShellPair * ABs, const ShellPair * CDs, const Array<AtomProd> * ABap, const Array<AtomProd> * CDap, UI32 tmax12, UI32 tmax34, bool skipsame);
    void SetLists(const ShellPair * ABs,                         const Array<AtomProd> * ABap                              , UI32 tmax12,              bool skipsame);


    bool operator<(const BatchInfo & rhs) const {
        UI8 gl, gr;
        if      (geometry==ABCD) gl = 1;
        else if (geometry==AACD) gl = 2;
        else if (geometry==AACC) gl = 3;
        else if (geometry==ABAB) gl = 4;
        else if (geometry==AAAA) gl = 5;
        else  gl = 0;

        if      (rhs.geometry==ABCD) gr = 1;
        else if (rhs.geometry==AACD) gr = 2;
        else if (rhs.geometry==AACC) gr = 3;
        else if (rhs.geometry==ABAB) gr = 4;
        else if (rhs.geometry==AAAA) gr = 5;
        else  gr = 0;


        if (gl!=gr) return (gl<gr);

        if (ABp->l1!=rhs.ABp->l1) return (ABp->l1<rhs.ABp->l1);
        if (ABp->l2!=rhs.ABp->l2) return (ABp->l2<rhs.ABp->l2);
        if (CDp->l1!=rhs.CDp->l1) return (CDp->l1<rhs.CDp->l1);
        if (CDp->l2!=rhs.CDp->l2) return (CDp->l2<rhs.CDp->l2);

        if (ABp->ka!=rhs.ABp->ka) return (ABp->ka<rhs.ABp->ka);
        if (ABp->kb!=rhs.ABp->kb) return (ABp->kb<rhs.ABp->kb);
        if (CDp->ka!=rhs.CDp->ka) return (CDp->ka<rhs.CDp->ka);
        if (CDp->kb!=rhs.CDp->kb) return (CDp->kb<rhs.CDp->kb);

        return (ABp<CDp);
    }

};



//blox will be dispatched to any available hardware
class ERIblock:public BatchInfo {
    friend class ERIBatch;
    friend class Fock2e;

    //memory buffers and pointers for ERI evaluation
  private:
    ERIbatchBuffer * EBB;

    size_t ElementSize;
    ShellPairPair * SPPlist;
    double       * R2;


    //atom blocks
    UI32 min12;
    UI32 max12;
    UI32 * min34l;
    UI32 * max34l;

    //the transposed ones
    UI32 min34;
    UI32 max34;
    UI32 * min12l;
    UI32 * max12l;

    r1tensor <SPlist>   SortedABs;
    r1tensor <SPlist>   SortedCDs;
    r2tensor <SPPblock> SortedERIs;


    //number of GDO buckets for each elemental basis function pair
    UI16 K2max12;
    UI16 K2max34;

    SPJ * ABlogs;
    SPJ * CDlogs;

    static const double maxX = 1e100;


    template <class T> struct Xrow {
        std::map <UI32,T> r;
        //int n;
        float minf;

        void add (UI32 atb, const T & v) {
            r[atb] = v;
        }

        int size() const {
            return r.size();
        }

        T & operator[] (UI32 atb) {
            return r[atb]; //r[atb];
        }

        T operator[] (UI32 atb) const {
            if (r.count(atb)) return r.at(atb);
            else              return maxX;
        }

        T GetMin() {
            minf = maxX;

            std::map <UI32,float>::const_iterator it;

            for (it=r.begin(); it!=r.end(); ++it) minf = std::min (minf, it->second);

            return minf;
        }
    };

    template <class T> struct Xsparse {
        Xrow<T> * v;
        int n;
        float minf;

        Xsparse() {
            n = 0;
            v = NULL;
        }

        void set(int nn) {
            n = nn;
            v = new Xrow<T>[n];
        }

        void clean() {
            delete[] v;
            v = NULL;
            n = 0;
        }

        void add(int a, int b, float vv) {
            v[a][b] = vv;
        }

        ~Xsparse() {
            delete[] v;
        }

        const Xrow<T> & operator[] (int i) const {
            return v[i];
        }

        Xrow<T> & operator[] (int i) {
            return v[i];
        }


        T operator() (int i, int j) const {
            if (v[i].r.count(j)==0) return maxX;
            return v[i].r.at(j);
        }

        T & operator() (int i, int j) {
            return v[i].r[j];
        }

        float GetMin() {
            minf = v[0].GetMin();

            for (int i=1; i<n; ++i) minf = std::min (minf, v[i].GetMin());

            return minf;
        }
    };

    struct {
        int natoms;

        //make a list of ata, atb, atc, atd
        std::set<UI32> AtomsA, AtomsB, AtomsC, AtomsD;

        Xsparse<UI32>    AtomsAB, AtomsCD;
        Xsparse<float>   AtomsAC, AtomsBD, AtomsAD, AtomsBC;
    } ONXtensors;


    void InitONX(const r2tensor< r2tensor<float> > & logFs, const r2tensor< r2tensor<float> > & logDs);
    void CleanONX();


    //FUNCTIONS FOR GENERATING LISTS, PRESCREENING AND PACKING ERIs
    //*************************************************************

    //mode: 0 -> regular  1 -> only count 3-center ABAD integrals   2 -> discard 3-center ABAD integrals
    template<bool G, bool J, bool X, bool R, int mode> UI64 PrescreenedList (float logD, float logS, const r2tensor< r2tensor<float> > & logAD, const r2tensor< r2tensor<float> > & logDs, UI64 NINTS);
    template<bool G, bool J, bool X, bool R>           UI64 PrescreenedListS(float logD, float logS, const r2tensor< r2tensor<float> > & logAD, const r2tensor< r2tensor<float> > & logDs, UI64 NINTS);


    void GenerateCasePrescreenedList(float logD, float logS, const r2tensor< r2tensor<float> > & logDs);

    template<bool G, int mode> UI64 ONXPrescreenedList (float logD, float logS,  const r2tensor< r2tensor<float> > & logFs, const r2tensor< r2tensor<float> > & logDs, UI64 NINTS);
    template<bool G, int mode> UI64 JPrescreenedList   (float logD, float logS,  const r2tensor< r2tensor<float> > & logFs, const r2tensor< r2tensor<float> > & logDs, UI64 NINTS);

    void SortPrescreenedLists();
    void CaseSortPrescreenedLists();
    void ClearBlocks();
    void InitLists();
    void InitListsJ(const r2tensor< r2tensor<float> > & logDs);
    template <int LEN> void InitBlocks();
    template <int LEN> void FuseIncompleteTiles();


    //functions for prescreening
    template <int LEN> UI64 GetNTiles() const;
    UI64 MakeAllTiles(LibQuimera::ERITile64 * pTile);
    UI64 MakeAllTiles(LibQuimera::ERITile32 * pTile);


    void SetInitialList       (float lmo);
    void SetInitialList_ABCD  (float lmo);
    void SetInitialList_ABCDs (float lmo);
    void SetInitialList_ABAB  (float lmo);
    void SetInitialList_AACD  (float lmo);
    void SetInitialList_AACC  (float lmo);
    void SetInitialList_AACCs (float lmo);

    bool Check3Center();

    //add several blocks to a list
    // J: add only Coulomb, X: only exchange, mode: 1: ABAB integrals, 2: ABCD integrals with potential center degeneracy, 0: all other kinds
    template<bool J, bool X, int mode> bool MakeBlocks      (float logD, float logS, const r2tensor< r2tensor<float> > & logADs, const r2tensor< r2tensor<float> > & logDs, TaskQueue<ERIBatch> & TaskList);

    bool include(bool J, bool X, bool R, bool useJD, bool useJS, bool useXD, bool useXS);

//functions and variables to be accessed publicly
  public:

    //chronometers
    //************
    UI64 nJs;
    UI64 nXs;

    UI64 mJs;
    UI64 mXs;

    double tJs;
    double tXs;

    double pJs;
    double pXs;


    ERIblock();
    ~ERIblock();

    void Init(size_t max_at_int, ERIbatchBuffer * ERIBB);              //very first initialization
    void Reset();
    const ERIblock & operator=(const BatchInfo & rhs); //sets the information for the next shellpair pair batch of integrals


    //FUNCTIONS FOR BLOCK EVALUATION OF ERIs
    //**************************************

    bool MakeCoulomb  (float logD, float logS,  const r2tensor< r2tensor<float> > & logADs, const r2tensor< r2tensor<float> > & logDs, TaskQueue<ERIBatch> & TaskList);
    bool MakeExchange (float logD, float logS,  const r2tensor< r2tensor<float> > & logADs, const r2tensor< r2tensor<float> > & logDs, TaskQueue<ERIBatch> & TaskList);
    bool MakeAll      (float logD, float logS,  const r2tensor< r2tensor<float> > & logADs, const r2tensor< r2tensor<float> > & logDs, TaskQueue<ERIBatch> & TaskList);




    //OTHER FUNCTION EVALUATIONS
    //**************************
    cacheline64  CalcCauchySchwarz(const ShellPair & ABs, const ERIgeometries64 & eriG8, const LibQuimera::ERITile64 & ET, UI8 lenk);
    float        CalcCauchySchwarz(const ShellPair & ABs, const ERIgeometry & eriG);

    void      ConstructW       (           const ShellPair & ABs, const ShellPair & CDs, const ERIgeometry & eriG, double * WWW, UI32 wa, bool UseCase = false);
    void      ConstructWS      (double w, const ShellPair & ABs, const ShellPair & CDs, const ERIgeometry & eriG, double * WWW, UI32 wa);
    void      ConstructWL      (double w, const ShellPair & ABs, const ShellPair & CDs, const ERIgeometry & eriG, double * WWW, UI32 wa);
} __attribute__((aligned(64))); //don't wverlap in cache!



#endif
