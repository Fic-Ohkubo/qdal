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


#include "../2eints/quimera.hpp"

#include "../basis/SPprototype.hpp"
#include "../math/angular.hpp"
#include "../low/cache64.hpp"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif


using namespace std;
using namespace LibQuimera;


K4benchmark::K4benchmark() {
    ncallsK4 = ncallsMIRROR = 0;
    deltaK4  = deltaMIRROR  = 0;
    tflopsK4 = tflopsMIRROR = 0;

    aK4 = aKa = aKb = aKc = aKd = 0;
}

void K4benchmark::AddCall(const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp)  {

    const PrimitiveSet & PSab = ABp.Psets[ET.nKab];
    const PrimitiveSet & PScd = CDp.Psets[ET.nKcd];

    int K1 = PScd.nKb;
    int K2 = ET.nKcd;

    int K3 = PSab.nKb;
    int K4 = ET.nKab;

    K3*=K2;
    K4*=K2;

    int Ka = PSab.nKa[0];
    int Kb = PSab.nKb;
    int Kc = PScd.nKa[0];
    int Kd = PScd.nKb;

    const int JA = ABp.Ja;
    const int JB = ABp.Jb;
    const int JC = CDp.Ja;
    const int JD = CDp.Jb;

    UI32 NFLOPS = (5 + nK4J1*JA) * K4 + (2 + nK3J1*JA + nK3J2*JA*JB) * K3 + (6 + nK2J2*JA*JB + nK2J3*JA*JB*JC) * K2 + (3 + nK1J3*JA*JB*JC + nK1J4*JA*JB*JC*JD) * K1 + (nK0J4*JA*JB*JC*JD + 17);
    NFLOPS  *= 8;


    #pragma omp atomic
    tflopsK4 += NFLOPS;

    #pragma omp atomic
    aK4     += K4;
    #pragma omp atomic
    aKa     += Ka;
    #pragma omp atomic
    aKb     += Kb;
    #pragma omp atomic
    aKc     += Kc;
    #pragma omp atomic
    aKd     += Kd;


    #pragma omp atomic
    ++ncallsK4;
}

double K4benchmark::Statistics() const {

    if (ncallsK4==0 && ncallsMIRROR==0) return 0;

    LibQuimera::Quimera::QBenchmarker << "ERI batch: ";
    if      (type.geometry==ABCD) LibQuimera::Quimera::QBenchmarker << "ABCD  ";
    else if (type.geometry==AACD) LibQuimera::Quimera::QBenchmarker << "AACD  ";
    else if (type.geometry==AACC) LibQuimera::Quimera::QBenchmarker << "AACC  ";

    LibQuimera::Quimera::QBenchmarker << int(type.La) << " " << int(type.Lb) << " " << int(type.Lc) << " " << int(type.Ld) << " :" << endl;

    double Tt=0;

    {
        LibQuimera::Quimera::QBenchmarker << " K4 contraction code: ";
        if (type.isCDR) LibQuimera::Quimera::QBenchmarker << "CDR  ";
        else            LibQuimera::Quimera::QBenchmarker << "NO CDR  ";

        //if (mC8!=NULL) LibQuimera::Quimera::QBenchmarker << " static" << endl;
        //else           LibQuimera::Quimera::QBenchmarker << " interpreted" << endl;
        LibQuimera::Quimera::QBenchmarker << endl;
        LibQuimera::Quimera::QBenchmarker << "  Calls:       " << ncallsK4 << endl;

        double At = deltaK4;
        double MFlops = double(tflopsK4)/(At*1000000);
        double time   = (1000000*At)/double(DOUBLES_PER_CACHE_LINE*ncallsK4);

        double aFLOPs = double(tflopsK4)/double(DOUBLES_PER_CACHE_LINE*ncallsK4);

        LibQuimera::Quimera::QBenchmarker << "   Total time:     " << At << " sec" << endl;
        LibQuimera::Quimera::QBenchmarker << "   Total FLOPS:    " << tflopsK4 << endl;
        LibQuimera::Quimera::QBenchmarker << "   Performance:    " << MFlops << " MFflops" << endl;
        LibQuimera::Quimera::QBenchmarker << "   FLOP breakdown: " << nK4J0 << " K4J0 + " << nK4J1 << " K4J1 + " << nK3J1 << " K3J1 + " << nK3J2 << " K3J2 + " << nK2J2 << " K2J2 + " << nK2J3 << " K2J3 + " << nK1J3 << " K1J3 + " << nK1J4 << " K1J4 + " << nK0J4 << " K0J4 FLOPS" << endl;

        LibQuimera::Quimera::QBenchmarker << "   Average FLOPs:  " << aFLOPs << " FLOPs" << endl;

        LibQuimera::Quimera::QBenchmarker << "   Average K4:     " << double(aK4)/double(ncallsK4) << endl;
        LibQuimera::Quimera::QBenchmarker << "   Average Ks:     " << double(aKa)/double(ncallsK4) << ", " << double(aKb)/double(ncallsK4) << ", " << double(aKc)/double(ncallsK4) << ", " << double(aKd)/double(ncallsK4) << endl;
        LibQuimera::Quimera::QBenchmarker << "   Average time:   " << time   << " usec" << endl;

        Tt += At;
    }


    {
        LibQuimera::Quimera::QBenchmarker << " MIRROR transformation code: ";
        //if (mC8!=NULL) LibQuimera::Quimera::QBenchmarker << " static" << endl;
        //else
        LibQuimera::Quimera::QBenchmarker << " interpreted" << endl;

        LibQuimera::Quimera::QBenchmarker << "  Calls:       " << ncallsMIRROR << endl;

        double At = deltaMIRROR;
        UI64 tflopsMIRROR = (DOUBLES_PER_CACHE_LINE*ncallsMIRROR)*nMIRROR;
        double MFlops = double(tflopsMIRROR)/(At*1000000);
        double time   = (1000000*At)/double(DOUBLES_PER_CACHE_LINE*ncallsMIRROR);

        //LibQuimera::Quimera::QBenchmarker << " Kernel transformations:" << endl;
        LibQuimera::Quimera::QBenchmarker << "   Total time:     " << At << " sec" << endl;
        LibQuimera::Quimera::QBenchmarker << "   Total FLOPs:    " << tflopsMIRROR << endl;
        LibQuimera::Quimera::QBenchmarker << "   Performance:    " << MFlops << " MFflops" << endl;
        LibQuimera::Quimera::QBenchmarker << "   FLOP count:     " << nMIRROR << " FLOPs" << endl;
        LibQuimera::Quimera::QBenchmarker << "   Average time:   " << time   << " usec" << endl;

        Tt += At;
    }

    LibQuimera::Quimera::QBenchmarker << endl;

    return Tt;
}

void K4benchmark::Reset() {
    ncallsK4 = ncallsMIRROR = 0;
    deltaK4  = deltaMIRROR  = 0;
    tflopsK4 = tflopsMIRROR = 0;

    aK4 = aKa = aKb = aKc = aKd = 0;
}


size_t p_ERIbuffer::SetBuffers(cacheline * mem, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp,  const p_Qalgorithm & algo) {

    const int JA = ABp.Ja;
    const int JB = ABp.Jb;
    const int JC = CDp.Ja;
    const int JD = CDp.Jb;

    const int KA = ABp.Ka;
    const int KB = ABp.Kb;
    const int KC = CDp.Ka;
    const int KD = CDp.Kb;


    int maxK4 = KA * KB * KC * KD;

    int mF0  = 2 * maxK4;
    int mF0e = algo.memF0e;
    int mF0f = algo.memF0f;

    int mJ1e = algo.memK4J1e * JA;
    int mJ1f = algo.memK4J1f * JA;

    int mJ2e = algo.memK3J2e * JA*JB;
    int mJ2f = algo.memK3J2f * JA*JB;

    int mJ3e = algo.memK2J3e * JA*JB*JC;
    int mJ3f = algo.memK2J3f * JA*JB*JC;

    int mJ4e = algo.memK1J4e * JA*JB*JC*JD;
    int mJ4f = algo.memK1J4f * JA*JB*JC*JD;



    size_t offset = 0;

    //for cacheline64

    //bufferK4      = mem + offset; offset += algo.K4Mem;

    KAB      = mem + offset; offset += KA * KB;
    KCD      = mem + offset; offset += KC * KD;



    F0    = mem + offset; offset += mF0;

    F0e   = mem + offset; offset += mF0e;
    F0f   = mem + offset; offset += mF0f;

    K4J1e = mem + offset; offset += mJ1e;
    K4J1f = mem + offset; offset += mJ1f;

    K3J2e = mem + offset; offset += mJ2e;
    K3J2f = mem + offset; offset += mJ2f;

    K2J3e = mem + offset; offset += mJ3e;
    K2J3f = mem + offset; offset += mJ3f;

    K1J4e = mem + offset; offset += mJ4e;
    K1J4f = mem + offset; offset += mJ4f;

    bufferMIRROR  = mem + offset; offset += algo.MaxMem; //to store intermediate values in transformations

    // for one-center integrals

    //bufferk4 = (double*)bufferK4;
    f0    = (double*) F0;
    f0e   = (double*) F0e;
    f0f   = (double*) F0f;
    k4j1e = (double*) K4J1e;
    k4j1f = (double*) K4J1f;
    k3j2e = (double*) K3J2e;
    k3j2f = (double*) K3J2f;
    k2j3e = (double*) K2J3e;
    k2j3f = (double*) K2J3f;
    k1j4e = (double*) K1J4e;
    k1j4f = (double*) K1J4f;

    buffer   = (double*)bufferMIRROR;


    return (offset * sizeof(cacheline));
}


//generic ERI algorithm
p_Qalgorithm::p_Qalgorithm() {
    IC = NULL;
    SC = NULL;
    //bench = NULL;
    IsDynamic = false;
}

p_Qalgorithm::~p_Qalgorithm() {
    delete IC;
    delete SC;
    //delete bench;
}


UI32 p_Qalgorithm::GetNKernels() const {
    return nKernels;
}

size_t p_Qalgorithm::MemSize(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp) const {


    const int JA = ABp.Ja;
    const int JB = ABp.Jb;
    const int JC = CDp.Ja;
    const int JD = CDp.Jb;

    const int KA = ABp.Ka;
    const int KB = ABp.Kb;
    const int KC = CDp.Ka;
    const int KD = CDp.Kb;

    /*
    int mF0 = algo.memF0;
    int mJ1 = algo.memJ1*JA;
    int mJ2 = algo.memJ2*JA*JB;
    int mJ3 = algo.memJ3*JA*JB*JC;
    int mJ4 = algo.memJ4*JA*JB*JC*JD;
    */


    int mF0  = memF0;
    int mF0e = memF0e;
    int mF0f = memF0f;

    int mJ1e = memK4J1e * JA;
    int mJ1f = memK4J1f * JA;

    int mJ2e = memK3J2e * JA*JB;
    int mJ2f = memK3J2f * JA*JB;

    int mJ3e = memK2J3e * JA*JB*JC;
    int mJ3f = memK2J3f * JA*JB*JC;

    int mJ4e = memK1J4e * JA*JB*JC*JD;
    int mJ4f = memK1J4f * JA*JB*JC*JD;


    int maxK4 = KA * KB * KC * KD;


    size_t offset = 0;

    //for cacheline64
    offset += MaxMem; //to store intermediate values in transformations
    //offset += K4Mem;

    offset += KA * KB;
    offset += KC * KD;



    offset += mF0*maxK4;

    offset += mF0e;
    offset += mF0f;

    offset += mJ1e;
    offset += mJ1f;

    offset += mJ2e;
    offset += mJ2f;

    offset += mJ3e;
    offset += mJ3f;

    offset += mJ4e;
    offset += mJ4f;


    return (offset * sizeof(cacheline));
}



p_Qalgorithm::p_Qalgorithm(K4routine & routine) {
    IC = NULL;
    SC = &routine;
    IsDynamic = false;
}

void p_Qalgorithm::SetStatic(ERItype   & type, K4routine & routine) {
    if (IC   !=NULL) delete IC;
    //if (bench!=NULL) delete bench;

    IC = NULL;
    SC = &routine;
    //bench = NULL;

    IsDynamic = false;

    eritype = type;
}

void p_Qalgorithm::SetIC     (ERItype   & type) {
    if (IC   !=NULL) delete IC;
    //if (bench!=NULL) delete bench;

    IC = NULL;
    SC = NULL;
    //bench = NULL;

    IsDynamic = true;

    eritype = type;
}


void p_Qalgorithm::Initialize() {

    //initialized already

    //if (IsDynamic)
    {
        //complete initialization
        IC = new ERIroutine;

        IC->Set(eritype.La, eritype.Lb, eritype.Lc, eritype.Ld, eritype.geometry, eritype.isCDR);
        IC->Init();


        //memory, etc.
        memF0    = IC->memF0;
        memF0e   = IC->memF0e;
        memF0f   = IC->memF0f;
        memK4J1e = IC->memK4J1e;
        memK4J1f = IC->memK4J1f;
        memK3J2e = IC->memK3J2e;
        memK3J2f = IC->memK3J2f;
        memK2J3e = IC->memK2J3e;
        memK2J3f = IC->memK2J3f;
        memK1J4e = IC->memK1J4e;
        memK1J4f = IC->memK1J4f;

        K4Mem    = IC->K4Mem;
        MaxMem   = IC->MaxMem;
        nKernels = IC->nKernels;

        //performance
        /*
        bench->nK4J0   = IC->nK4J0;
        bench->nK4J1   = IC->nK4J1;
        bench->nK3J1   = IC->nK3J1;
        bench->nK3J2   = IC->nK3J2;
        bench->nK2J2   = IC->nK2J2;
        bench->nK2J3   = IC->nK2J3;
        bench->nK1J3   = IC->nK1J3;
        bench->nK1J4   = IC->nK1J4;
        bench->nK0J4   = IC->nK0J4;
        bench->nMIRROR = IC->NFLOPS;
        */
    }

}

void p_Qalgorithm::Clear() {
    delete IC;
    delete SC;

    IC    = NULL;
    SC    = NULL;
}




void p_Qalgorithm::K4(const ERIgeometries32 & vars16, const ERITile32 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline32 * uv_m_st16, p_ERIbuffer & buffer, bool OnlyJ) const {

    IC->CalcGammas     (vars16, ET, ABp, CDp, buffer, OnlyJ);
    IC->ContractCDR_GC (vars16, ET, ABp, CDp, uv_m_st16, buffer);

    //if (IsDynamic) IC->ContractCDR_GC(vars8,ET,ABp,CDp,uv_m_st8,buffer);
    //else          (*SC)(buffer,vars8, ET, ABp, CDp, uv_m_st8);
}

void p_Qalgorithm::K4(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, double * __restrict__ uv_m_st, p_ERIbuffer & buffer, bool OnlyJ) const {
    IC->ContractCDR_GC(ABp, CDp, uv_m_st, buffer, OnlyJ);
}



void p_Qalgorithm::MIRRORbench(const cacheline64 * __restrict__  mem, const ERIgeometries64 & vars8, cacheline64 * __restrict__ ERI8, p_ERIbuffer & buffer) const {

    Chronometer chronoT;
    chronoT.Start();


    if      (eritype.geometry==ABCD) IC->TransformABCD(mem, vars8, ERI8, buffer);
    else if (eritype.geometry==AACD) IC->TransformAACD(mem, vars8, ERI8, buffer);
    else if (eritype.geometry==AACC) IC->TransformAACC(mem, vars8, ERI8, buffer);
    else if (eritype.geometry==ABAB) IC->TransformABCD(mem, vars8, ERI8, buffer);
    else                             IC->TransformABCD(mem, vars8, ERI8, buffer);


    chronoT.Stop();

    double At = chronoT.GetTotalTime();

    K4benchmark * bench = p_Q.Benchmarks[eritype];

    #pragma omp atomic
    ++bench->ncallsMIRROR;

    #pragma omp atomic
    bench->deltaMIRROR += At;
}


void p_Qalgorithm::MIRROR(const cacheline64 * __restrict__  mem, const ERIgeometries64 & vars8, cacheline64 * __restrict__ ERI8, p_ERIbuffer & buffer) const {
    if      (eritype.geometry==ABCD) IC->TransformABCD(mem, vars8, ERI8, buffer);
    else if (eritype.geometry==AACD) IC->TransformAACD(mem, vars8, ERI8, buffer);
    else if (eritype.geometry==AACC) IC->TransformAACC(mem, vars8, ERI8, buffer);
    else if (eritype.geometry==ABAB) IC->TransformABAB(mem, vars8, ERI8, buffer);
    else                             IC->TransformABCD(mem, vars8, ERI8, buffer);
}

void p_Qalgorithm::MIRROR(const cacheline32 * __restrict__  mem, const ERIgeometries32 & vars16, cacheline32 * __restrict__ ERI16, p_ERIbuffer & buffer) const {
    if      (eritype.geometry==ABCD) IC->TransformABCD(mem, vars16, ERI16, buffer);
    else if (eritype.geometry==AACD) IC->TransformAACD(mem, vars16, ERI16, buffer);
    else if (eritype.geometry==AACC) IC->TransformAACC(mem, vars16, ERI16, buffer);
    else if (eritype.geometry==ABAB) IC->TransformABAB(mem, vars16, ERI16, buffer);
    else                             IC->TransformABCD(mem, vars16, ERI16, buffer);
}


void p_Qalgorithm::MIRROR(const double * __restrict__ uv_m_st, double * W, p_ERIbuffer & buffer) const {
    IC->TransformAAAA(uv_m_st, W, buffer);
}




p_Quimera p_Q;

p_Quimera::p_Quimera() {

    //link all static code
    LinkStatic(true);
}

p_Quimera::~p_Quimera() {

    map<ERItype, Qalgorithm*>::iterator it;

    for (it = Solvers.begin(); it!=Solvers.end(); ++it) {
        //remember to destroy the pointers
        delete it->second;
        it->second = NULL;
    }

    map<ERItype, K4benchmark*>::iterator itb;

    for (itb = Benchmarks.begin(); itb!=Benchmarks.end(); ++itb) {
        //remember to destroy the pointers
        delete itb->second;
        itb->second = NULL;
    }
}

void p_Quimera::SetStatic(GEOM geometry, UI8 la, UI8 lb, UI8 lc, UI8 ld,    K4routine routine, bool overwrite=false) {

    ERItype type;
    type(geometry, la, lb, lc, ld, UseGC, false, UseCDR);

    //check whether an algorithm already exists
    if (!overwrite && Solvers.count(type) == 1) return;

    //clear if already existing
    if (Solvers.count(type) == 1) Solvers[type]->pQalgorithm->Clear();

    //set the algorithm as a static link to the routine
    Solvers[type]->pQalgorithm->SetStatic(type, routine);
}

void p_Quimera::LinkStatic(bool overwrite = false) {


    //some static libraries
    //ERItype GC_CDR;
    //GC_CDR(ABCD, 0,0,0,0, false, false, true);
    //SetStatic ( GC_CDR , CDRGC_ABCD_SSSS, overwrite);

    /*
    SetStatic ( ABCD, 0,0,0,0, false, false, true , K4GC_ABCD_SSSS, overwrite);
    SetStatic ( ABCD, 1,0,0,0, false, false, true , K4GC_ABCD_PSSS, overwrite);
    SetStatic ( ABCD, 1,0,1,0, false, false, true , K4GC_ABCD_PSPS, overwrite);
    SetStatic ( ABCD, 1,1,0,0, false, false, true , K4GC_ABCD_PPSS, overwrite);
    SetStatic ( ABCD, 1,1,1,0, false, false, true , K4GC_ABCD_PPPS, overwrite);
    SetStatic ( ABCD, 1,1,1,1, false, false, true , K4GC_ABCD_PPPP, overwrite);
    */

    /*
    SetStatic ( ABCD, 2,0,0,0, false, false, true , K4GC_ABCD_DSSS, overwrite);
    SetStatic ( ABCD, 2,0,1,0, false, false, true , K4GC_ABCD_DSPS, overwrite);
    SetStatic ( ABCD, 2,0,1,1, false, false, true , K4GC_ABCD_DSPP, overwrite);
    SetStatic ( ABCD, 2,0,2,0, false, false, true , K4GC_ABCD_DSDS, overwrite);
    SetStatic ( ABCD, 2,1,0,0, false, false, true , K4GC_ABCD_DPSS, overwrite);
    SetStatic ( ABCD, 2,1,1,0, false, false, true , K4GC_ABCD_DPPS, overwrite);
    SetStatic ( ABCD, 2,1,1,1, false, false, true , K4GC_ABCD_DPPP, overwrite);
    SetStatic ( ABCD, 2,1,2,0, false, false, true , K4GC_ABCD_DPDS, overwrite);
    SetStatic ( ABCD, 2,1,2,1, false, false, true , K4GC_ABCD_DPDP, overwrite);
    SetStatic ( ABCD, 2,2,0,0, false, false, true , K4GC_ABCD_DDSS, overwrite);
    SetStatic ( ABCD, 2,2,1,0, false, false, true , K4GC_ABCD_DDPS, overwrite);
    SetStatic ( ABCD, 2,2,1,1, false, false, true , K4GC_ABCD_DDPP, overwrite);
    SetStatic ( ABCD, 2,2,2,0, false, false, true , K4GC_ABCD_DDDS, overwrite);
    SetStatic ( ABCD, 2,2,2,1, false, false, true , K4GC_ABCD_DDDP, overwrite);
    SetStatic ( ABCD, 2,2,2,2, false, false, true , K4GC_ABCD_DDDD, overwrite);
    */
}


void p_Quimera::ListNeeded  (GEOM geometry, UI8 la, UI8 lb, UI8 lc, UI8 ld) {
    ERItype type;
    type(geometry, la, lb, lc, ld, UseGC, false, UseCDR);

    SetLater.push(type);
}

void p_Quimera::GenerateNeeded(bool overwrite) {


    //fill the map with the needed keys for later initialization
    while (SetLater.size() != 0) {
        ERItype type;
        type = SetLater.top();
        SetLater.pop();

        //marks the algorithm for later initialization
        if (overwrite || Solvers.count(type)==0) {
            Solvers[type] = new Qalgorithm; //->pQalgorithm->SetIC(type);
            Solvers[type]->pQalgorithm->SetIC(type);
        }
    }




    Quimera::QMessenger << "Initializing interpreted routines";
    Quimera::QMessenger.Push(); {
        map<ERItype, Qalgorithm*>::iterator it = Solvers.begin();
        int n = Solvers.size();

        #pragma omp parallel for schedule(dynamic)
        for (UI32 i=0; i<n; ++i) {

            Qalgorithm * algo;

            #pragma omp critical
            {
                algo = (it->second);
                ++it;
            }

            algo->pQalgorithm->Initialize();

            //QMessenger.Out("*");
        }
    } Quimera::QMessenger.Pop(true);



    std::map <ERItype, Qalgorithm*>::iterator it;

    for (it=Solvers.begin(); it!=Solvers.end(); ++it) {
        ERItype type = it->first;

        Benchmarks[type] = new K4benchmark;
        Benchmarks[type]->type = type;

        //performance
        Benchmarks[type]->nK4J0   = it->second->pQalgorithm->IC->nK4J0;
        Benchmarks[type]->nK4J1   = it->second->pQalgorithm->IC->nK4J1;
        Benchmarks[type]->nK3J1   = it->second->pQalgorithm->IC->nK3J1;
        Benchmarks[type]->nK3J2   = it->second->pQalgorithm->IC->nK3J2;
        Benchmarks[type]->nK2J2   = it->second->pQalgorithm->IC->nK2J2;
        Benchmarks[type]->nK2J3   = it->second->pQalgorithm->IC->nK2J3;
        Benchmarks[type]->nK1J3   = it->second->pQalgorithm->IC->nK1J3;
        Benchmarks[type]->nK1J4   = it->second->pQalgorithm->IC->nK1J4;
        Benchmarks[type]->nK0J4   = it->second->pQalgorithm->IC->nK0J4;
        Benchmarks[type]->nMIRROR = it->second->pQalgorithm->IC->NFLOPS;
    }

}

void p_Quimera::Statistics() const {

    LibQuimera::Quimera::QBenchmarker << "ERI statistics:" << endl;

    map<ERItype, K4benchmark*>::const_iterator it;

    double AT = 0;

    for (it = Benchmarks.begin(); it!=Benchmarks.end(); ++it)
        AT += it->second->Statistics();

    LibQuimera::Quimera::QBenchmarker << endl;
    LibQuimera::Quimera::QBenchmarker << "Total time spent in contractions + tranformations: " << AT << endl;
}

const Qalgorithm * p_Quimera::SelectAlgorithm (GEOM geometry, UI8 la, UI8 lb, UI8 lc, UI8 ld) {
    ERItype type;
    type(geometry, la, lb, lc, ld, UseGC, false, UseCDR);

    //if there is no algorithm in the pool
    if (Solvers.count(type) == 0) {
        //lazy initialization
        {
            Solvers[type] = new Qalgorithm;
            Solvers[type]->pQalgorithm->SetIC(type);
            Solvers[type]->pQalgorithm->Initialize();
        }

        {
            Benchmarks[type] = new K4benchmark;
            Benchmarks[type]->type = type;

            //performance
            Benchmarks[type]->nK4J0   = Solvers[type]->pQalgorithm->IC->nK4J0;
            Benchmarks[type]->nK4J1   = Solvers[type]->pQalgorithm->IC->nK4J1;
            Benchmarks[type]->nK3J1   = Solvers[type]->pQalgorithm->IC->nK3J1;
            Benchmarks[type]->nK3J2   = Solvers[type]->pQalgorithm->IC->nK3J2;
            Benchmarks[type]->nK2J2   = Solvers[type]->pQalgorithm->IC->nK2J2;
            Benchmarks[type]->nK2J3   = Solvers[type]->pQalgorithm->IC->nK2J3;
            Benchmarks[type]->nK1J3   = Solvers[type]->pQalgorithm->IC->nK1J3;
            Benchmarks[type]->nK1J4   = Solvers[type]->pQalgorithm->IC->nK1J4;
            Benchmarks[type]->nK0J4   = Solvers[type]->pQalgorithm->IC->nK0J4;
            Benchmarks[type]->nMIRROR = Solvers[type]->pQalgorithm->IC->NFLOPS;
        }
    }

    return (Solvers[type]);
}






