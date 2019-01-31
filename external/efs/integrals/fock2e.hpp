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



#ifndef __Fock2eS__
#define __Fock2eS__
#include <set>
#include <list>
#include <queue>

#include "../defs.hpp"
#include "../integrals/atomprod.hpp"
#include "../integrals/ERIblocks.hpp"
#include "../math/newsparse.hpp"
#include "ERIblocks.hpp"
#include "ERIbatch.hpp"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

const bool PrintStatistics = true;

void ConstructBigW(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, int wa, const double * w, double * W);

class ShellPair;
class AtomProdPrototypes;
class APlists;
class symtensor;
class Interactions;
class ERIgeometry;
class BatchInfo;
class ERIblock;
class BatchEvaluator;
class ERIbatchBuffer;



//classes for managing concurrent algorithms

class FockJobQueue {
    friend class Fock2e;

  private:
    bool Initialized;
    std::set<BatchInfo> JobQueue;

  public:

    FockJobQueue() {
        Initialized = false;
    }

    void push(BatchInfo & Info) {
        JobQueue.insert(Info);
    }

    void Done() {
        Initialized = true;
    }

    void Clear() {
        if (Initialized) {
            Initialized = false;
            JobQueue.clear();
        }
    }

};

template<class T> class TaskQueue;

template<class T> class ConstTaskQueue;


//this class handles the 3electron contribution to the Fock matrix
//****************************************************************
class Fock2e {
    friend class ERIblock;

    static size_t MEMPERBATCH;

  private:


    UI64 max_at_int;  //maximum number of atom-atom interactions to be considered
    UI16 OGMT;        //number of threads

    //pointers to required data
    const APlists            * APL;

    Sparse D2sparse;
    Sparse AD2sparse;

    //input and output
    r2tensor< r2tensor<float> > logDs;
    r2tensor< r2tensor<float> > logAD;


    double XS;
    float  logmover;            //-log of the prescreen threshold
    bool   GeneralContraction;
    bool   firstcycle;          //some lazy initializers required in the first call to the method


    //multithraded classes for handling each of the steps of the generation of the Fock matrix
    FockJobQueue TwoElectronJobs;
    ERIbatchBuffer * EbatchBuffer;

    ERIblock       * Blocks;
    BatchEvaluator * BatchEval;

    void InitERIroutines();
    void InitJobQueue(); //just once per molecule should suffice

    void CalcDensityPrescreeners();

    void CalcSparsicity(Sparse & SS);

    //void ExecuteJobs          (float logD, float logS);
    //void ExecuteJobsAlterning (float logD, float logS);
    template <bool J, bool X>  void ExecuteJobs          (float logD, float logS);


    void AddCoulomb(Sparse & F, int at, int wa, const double * const T, const Sparse & D);
    void AddXchange(Sparse & F, int at, int wa, const double * const T, const Sparse & D);

    void Trim(Sparse & F);

  public:

    Fock2e();

    void Init (UI32 nInteractingAtoms, const APlists & pAPL, const sparsetensorpattern & SparseTensorPattern); //iniciador
    void SetXS(double Xscaling);
    void CalcCauchySwarz();

    // on entrance, F contains a copy of the total density matrix, which along with the density matrix increment is used for efficient
    // prescreening of interactions;
    // on output, F contains the 2-electron Fock matrix update
    // AD can be either the total density matrix or just the density matrix update of the cycle
    void FockUpdate(tensor2 & F, const tensor2 & AD, double Dthresh, double Fthresh, bool SplitJX = false);

    ~Fock2e();
};


#endif
