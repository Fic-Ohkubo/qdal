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


#include <cstdlib>
#include <fstream>
#include <cmath>
#include <malloc.h>

#include "../defs.hpp"

#include "../low/chrono.hpp"
#include "../low/cache64.hpp"
#include "../math/angular.hpp"
#include "../math/gamma.hpp"

#include "../integrals/fock2e.hpp"
#include "../integrals/rotations.hpp"
#include "../integrals/atomprod.hpp"
#include "../integrals/ERIblocks.hpp"
#include "../integrals/newcontractions.hpp"

#include "../basis/shellpair.hpp"
#include "../basis/SPprototype.hpp"

#include "../libquimera/libquimera.hpp"
#include "../libechidna/libechidna.hpp"

using namespace std;
using namespace LibIGamma;


#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

#define __KC__
#define __TR__
#define __JXC__

using namespace LibAngular;
using namespace LibQuimera;

const UI16 DPC = DOUBLES_PER_CACHE_LINE;
const UI16 FPC = FLOATS_PER_CACHE_LINE;

const UI64 MAX_INTS_PER_BLOCK  = FPC*MAX_TILES_PER_BLOCK;
const UI64 EXTRA_INTS_PER_BLOCK = maxK2*maxK2*FPC;         //this is an absolute upper bound; in practice only a few are needed

static const int MINR2 = 1e-20;

//find number of interactions in O(ln(N))
static UI32 maxSP(float lmo, const ShellPair * CDs, UI32 max) {
    UI32 min = 0;

    while (max!=min) {
        UI32 med = (min+max)/2;

        if (CDs[med].logCS > lmo) max = med;
        else                      min = med+1;
    }

    return max;
}


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

//sets the general properties of the ERI block from the geometry and shell pair prototypes
void BatchInfo::Set(GEOM geom, const ShellPairPrototype & AB, const ShellPairPrototype & CD, bool samef) {
    //set prototypes and such
    ABp = &AB;
    CDp = &CD;

    la = ABp->l1;
    lb = ABp->l2;
    lc = CDp->l1;
    ld = CDp->l2;

    Lt = la + lb + lc + ld;
    UI8 maxab = la>lb?la:lb;
    UI8 maxcd = lc>ld?lc:ld;
    Lmax = maxab>maxcd?maxab:maxcd;

    SameShell = samef;


    //set geometry
    geometry = geom;

    if      (geometry==ABCD) ERIalgorithm = Q.SelectAlgorithm(ABCD,la,lb,lc,ld); // || geometry==ABAB
    else if (geometry==AACD) ERIalgorithm = Q.SelectAlgorithm(AACD,la,lb,lc,ld);
    else if (geometry==AACC) ERIalgorithm = Q.SelectAlgorithm(AACC,la,lb,lc,ld);
    else if (geometry==ABAB) ERIalgorithm = Q.SelectAlgorithm(ABAB,la,lb,lc,ld);
    else if (geometry==AAAA) ERIalgorithm = Q.SelectAlgorithm(AAAA,la,lb,lc,ld);


    //set values required for computing memory requirements
    J4     = ABp->Ja * ABp->Jb * CDp->Ja * CDp->Jb;
    Lt     = ABp->l1 + ABp->l2 + CDp->l1 + CDp->l2;
    fsize  = (ABp->Ka * ABp->Kb *CDp->Ka * CDp->Kb)*(Lt+1);
    msize  = ERIalgorithm->GetNKernels();
    wsize  = nmS[la]*nmS[lb]*nmS[lc]*nmS[ld];
    msize4 = J4*msize;
    wsize4 = J4*wsize;

    //SetJXC();

    State.storeABb  = 0;
    State.storeCDb  = 0;
    State.storeAB   = 0;
    State.storeCD   = 0;
    State.storeBlock   = false;
    State.store3center = false;
}

void BatchInfo::SetLists(const ShellPair * ABs, const ShellPair * CDs, const Array<AtomProd> * ABap, const Array<AtomProd> * CDap, UI32 max12, UI32 max34, bool skipsame) {
    AP12list = ABap;
    AP34list = CDap;
    SP12 = ABs;
    SP34 = CDs;
    SkipSame = skipsame;
    SameList = false;

    tmax12 = max12;
    tmax34 = max34;
}

void BatchInfo::SetLists(const ShellPair * ABs, const Array<AtomProd> * ABap, UI32 max12, bool skipsame) {
    AP12list = ABap;
    AP34list = ABap;
    SP12 = ABs;
    SP34 = ABs;
    SkipSame = skipsame;
    SameList = true;

    tmax12 = max12;
    tmax34 = max12;
}



ERIblock::ERIblock() {
    min12l = NULL;
    max12l = NULL;
    min34l = NULL;
    max34l = NULL;

    SPPlist = NULL;
    R2 = NULL;

    ABlogs = NULL;
    CDlogs = NULL;

    nJs = 0;
    nXs = 0;
}

ERIblock::~ERIblock() {
    delete[] min12l;
    delete[] max12l;
    delete[] min34l;
    delete[] max34l;

    delete[] SPPlist;
    delete[] R2;

    delete[] ABlogs;
    delete[] CDlogs;
}

void ERIblock::Init(size_t max_at_int, ERIbatchBuffer * ERIBB) {

    SPPlist = new ShellPairPair[MAX_TILES_PER_BLOCK*FPC*FPC]; // if every integral is put in a different box, DPC * NINT space is needed

    R2      = new double[MAX_TILES_PER_BLOCK*FPC]; // maximum possible number of SPP in the same block

    min34l  = new UI32[max_at_int];
    max34l  = new UI32[max_at_int];

    min12l  = new UI32[max_at_int];
    max12l  = new UI32[max_at_int];

    ABlogs  = new SPJ[max_at_int];
    CDlogs  = new SPJ[max_at_int];

    EBB = ERIBB;
}

void ERIblock::Reset() {
    nJs = 0;
    nXs = 0;
    mJs = 0;
    mXs = 0;
    tJs = 0;
    tXs = 0;
    pJs = 0;
    pXs = 0;
}

//default copy constructor plus initialization of some stuff
const ERIblock & ERIblock::operator=(const BatchInfo & rhs) {


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

    State = rhs.State;

    return *this;
}



//finds minimum distance between the segments AB and CD
double R2distance(const AtomProd & ABp, const AtomProd & CDp) {
    const point & A = ABp.A;
    const point & B = ABp.B;
    const point & C = CDp.A;
    const point & D = CDp.B;

    vector3 AB = B-A;
    vector3 AC = C-A;
    vector3 CD = D-C;

    double ab2  = AB*AB;
    double abcd = AB*CD;
    double cd2  = CD*CD;

    double acab = AC*AB;
    double accd = AC*CD;

    double dis = ab2*cd2 - abcd*abcd;

    double v = cd2*acab + abcd*accd;
    double w = ab2*accd + abcd*acab;

    point ABmin, CDmin;

    if      (v<=0)   ABmin = A;
    else if (v>=dis) ABmin = B;
    else             ABmin = A + (v/dis) * AB;

    if      (w<=0)   CDmin = C;
    else if (w>=dis) CDmin = D;
    else             CDmin = C + (w/dis) * CD;

    vector3 vmin = CDmin - ABmin;
    double R2 = vmin*vmin;
    if (R2<MINR2) R2=0;
    return (R2);
}


CasePrescreener CPEB;

double CasePrescreener::LogShortCoulomb (double R2) {

    double maxR2 = (vg_max-vg_step)/case_w2; //absolute maximum (any more gives an infinite logarithm)

    double logSC;

    double ** gamma = ::IncompleteGamma.gamma_table;

    //close enough to interact
    if (R2<maxR2) {
        double w2R2 = case_w2 * R2;

        double p = ivg_step*(w2R2-vg_min);
        int pos = int(p+0.5);
        double x0 = vg_min + vg_step*double(pos);
        double Ax1 = x0-w2R2;
        double Ax2 = 0.5 * Ax1*Ax1;
        double Ax3 = 0.33333333333333333333 * Ax1*Ax2;

        double G = (gamma[pos][1] + gamma[pos][2] * Ax1 + gamma[pos][3] * Ax2 + gamma[pos][4] * Ax3);

        G *= 2*sqrt(case_w2/PI);

        double S = 1/sqrt(R2) - G;

        logSC = -log(S);
    }
    else
        logSC = 100;

    return logSC;
}

void   CasePrescreener::InitLogsShort () {
     //absolute maximum (any more gives an infinite logarithm)

    double AR2 = maxR2 / double(NLogsShort);

    for (int i=0; i<NLogsShort; ++i) {
        double R2 = AR2 * NLogsShort;

        LogShort[i] = LogShortCoulomb(R2);
    }
}

double CasePrescreener::LogsShort (double R2) {
    if (R2>=maxR2) return 100;

    int p = (R2/maxR2)*NLogsShort;

    return LogShort[p]; // fast and cheap
}

//improve coulomb CASE prescreening
//use alternative CASE prescreening


void ERIblock::GenerateCasePrescreenedList(float logD, float logS,  const r2tensor< r2tensor<float> > & logDs) {
/*
    double iwmin = 1./case_w2 + 1./ABp->kmin + 1./CDp->kmin;
    double Dmax2 = ShortCoulomb.CalcR2(1./iwmin, case_err);
    double Dmax = 2* sqrt(Dmax2);

    PWBoxing TheBoxingScheme;

    //compute position
    PW * SPD = new PW [tmax12];

    {
        double k1    = ABp->ka[0];
        double k2    = ABp->kb[0];
        double k12   = k1+k2;
        double ik12  = 1/k12;

        for (UI32 ab=0; ab<tmax12; ++ab) {
            const AtomProd & AP12 = (*AP12list)[ab];

            SPD[ab].P.x = (k1*AP12.A.x + k2*AP12.B.x)*ik12;
            SPD[ab].P.y = (k1*AP12.A.y + k2*AP12.B.y)*ik12;
            SPD[ab].P.z = (k1*AP12.A.z + k2*AP12.B.z)*ik12;

            SPD[ab].w.n[0] = ab;
        }
    }

    //register everything
    TheBoxingScheme.InitFrom(SPD, tmax12, Dmax);

    delete[] SPD;



    //now compute interactions
    UI64 Nprescreened = 0;

    {
        double k3    = CDp->ka[0];
        double k4    = CDp->kb[0];
        double k34   = k3+k4;
        double ik34  = 1/k34;

        const PWBox * IntBoxes[8];

        for (UI32 cd=0; cd<tmax34; ++cd) {
            const AtomProd & AP34 = (*AP34list)[cd];

            point Q;

            Q.x = (k3*AP34.A.x + k4*AP34.B.x)*ik34;
            Q.y = (k3*AP34.A.y + k4*AP34.B.y)*ik34;
            Q.z = (k3*AP34.A.z + k4*AP34.B.z)*ik34;

            //retrieve the (maximally) 8 interacting boxes
            TheBoxingScheme.RetrieveBoxes(Q, IntBoxes);

            for (int nb=0; nb<8; ++nb) {
                if (IntBoxes[nb] == NULL) continue;

                //loop over the box' elements
                list<PW>::const_iterator it;

                for (it=IntBoxes[nb]->elements.begin(); it!=IntBoxes[nb]->elements.end(); ++it) {
                    PW P = *it;

                    //check whether it's necessary to compute
                    int ab = P.w.n[0];

                    if (cd>=max34l[ab])     continue;
                    if (SkipSame && ab==cd) continue;



                    UI32 ata = SP12[ab].ata;
                    UI32 atb = SP12[ab].atb;
                    UI32 atc = SP34[cd].ata;
                    UI32 atd = SP34[cd].atb;

                    UI16  fa = ABp->nb1;
                    UI16  fb = ABp->nb2;
                    UI16  fc = CDp->nb1;
                    UI16  fd = CDp->nb2;

                    //simple check (could be done better taking distance into account)
                    double wabcd = SP12[ab].logCS + SP34[cd].logCS;

                    double logJ;
                    logJ  = logDs[ata][atb][fa][fb] + logDs[atc][atd][fc][fd] + wabcd;

                    if (logJ<lmo)  {
                        ++Nprescreened;

                        UI16 i = K2max12-SP12[ab].nK2; //>0
                        UI16 j = K2max34-SP34[cd].nK2; //>0

                        ++SortedERIs(i,j).Nints;
                    }

                }
            }
        }

        InitBlocks();

        for (UI32 cd=0; cd<tmax34; ++cd) {
            const AtomProd & AP34 = (*AP34list)[cd];

            point Q;

            Q.x = (k3*AP34.A.x + k4*AP34.B.x)*ik34;
            Q.y = (k3*AP34.A.y + k4*AP34.B.y)*ik34;
            Q.z = (k3*AP34.A.z + k4*AP34.B.z)*ik34;

            //retrieve the (maximally) 8 interacting boxes
            TheBoxingScheme.RetrieveBoxes(Q, IntBoxes);

            for (int nb=0; nb<8; ++nb) {
                if (IntBoxes[nb] == NULL) continue;

                //loop over the box' elements
                list<PW>::const_iterator it;

                for (it=IntBoxes[nb]->elements.begin(); it!=IntBoxes[nb]->elements.end(); ++it) {
                    PW P = *it;

                    //check whether it's necessary to compute
                    int ab = P.w.n[0];

                    if (cd>=max34l[ab])     continue;
                    if (SkipSame && ab==cd) continue;



                    UI32 ata = SP12[ab].ata;
                    UI32 atb = SP12[ab].atb;
                    UI32 atc = SP34[cd].ata;
                    UI32 atd = SP34[cd].atb;

                    UI16  fa = ABp->nb1;
                    UI16  fb = ABp->nb2;
                    UI16  fc = CDp->nb1;
                    UI16  fd = CDp->nb2;

                    //simple check (could be done better taking distance into account)
                    double wabcd = SP12[ab].logCS + SP34[cd].logCS;

                    double logJ;
                    logJ  = logDs[ata][atb][fa][fb] + logDs[atc][atd][fc][fd] + wabcd;

                    if ( logJ<lmo )  {
                        UI16 i = K2max12-SP12[ab].nK2; //>0
                        UI16 j = K2max34-SP34[cd].nK2; //>0

                        UI64 t = SortedERIs(i,j).Iints;

                        SortedERIs(i,j).vList[t].n12 = ab;
                        SortedERIs(i,j).vList[t].n34 = cd;

                        ++SortedERIs(i,j).Iints;
                    }

                }
            }
        }
    }
*/
}

inline bool ERIblock::include(bool J, bool X, bool R, bool useJD, bool useJS, bool useXD, bool useXS) {
    bool useD = (J && useJD || X && useXD);
    bool useS = (J && useJS || X && useXS);


    return ( !R && useD ) || (R && useS && !useD); // || (R && (J && useJS || X && useXS) && !(J && useJD || X && useXD) ) );
    //return !R;
}


// G: generate integral list
// J: include Coulomb
// X: include Exchange
// R: single/double precission
// mode: 1:     only count ABAD integrals in the ABCD block
//       2:     skip ABAD-type integrals in ABCD block
//       other: count all integrals

template<bool G, bool J, bool X, bool R>           UI64 ERIblock::PrescreenedListS(float logD, float logS,  const r2tensor< r2tensor<float> > & logFs, const r2tensor< r2tensor<float> > & logDs, UI64 NINTS) {

    UI64 LPL = 0;

    for (UI32 ab=State.storeAB; ab<tmax12; ++ab) {
        if (!ThisNode.CheckMod(ab)) continue;

        UI32 ata = SP12[ab].ata;
        UI32 atb = SP12[ab].atb;
        UI32 atc = SP34[ab].ata;
        UI32 atd = SP34[ab].atb;

        UI16  fa = ABp->nb1;
        UI16  fb = ABp->nb2;
        UI16  fc = CDp->nb1;
        UI16  fd = CDp->nb2;


        float wabcd = SP12[ab].logCS + SP34[ab].logCS;

        float logJa;
        float logJb;

        float logX1a;
        float logX2a;
        float logX1b;
        float logX2b;

        bool useJD, useJS, useXD, useXS;
        useJD = useJS = useXD = useXS = false;

        if (J) {
            logJa  = logFs[ata][atb][fa][fb] + logDs[atc][atd][fc][fd] + wabcd;
            logJb  = logDs[ata][atb][fa][fb] + logFs[atc][atd][fc][fd] + wabcd;

            useJD = (logJa<logD || logJb<logD);
            useJS = (logJa<logS || logJb<logS);
        }
        if (X) {
            logX1a = logFs[ata][atc][fa][fc] + logDs[atb][atd][fb][fd] + wabcd;
            logX2a = logFs[ata][atd][fa][fd] + logDs[atb][atc][fb][fc] + wabcd;
            logX1b = logDs[ata][atc][fa][fc] + logFs[atb][atd][fb][fd] + wabcd;
            logX2b = logDs[ata][atd][fa][fd] + logFs[atb][atc][fb][fc] + wabcd;

            useXD = (logX1a<logD || logX2a<logD || logX1b<logD || logX2b<logD);
            useXS = (logX1a<logS || logX2a<logS || logX1b<logS || logX2b<logS);
        }


        //generate the list
        if (G) {
            if (include(J,X,R,useJD,useJS,useXD,useXS)) {
                UI16 i = K2max12-SP12[ab].nK2; //>0
                UI16 j = K2max34-SP34[ab].nK2; //>0

                ShellPairPair SPP(ab,ab);
                SortedERIs(i,j).push(SPP);

                ++LPL;

                if (LPL==NINTS) {
                    State.storeAB = ab+1;
                    return LPL; //interrupt the thread
                }
            }
        }
        //only count integrals
        else {
            if (include(J,X,R,useJD,useJS,useXD,useXS)) {
                ++LPL;

                UI16 i = K2max12-SP12[ab].nK2; //>0
                UI16 j = K2max34-SP34[ab].nK2; //>0

                SortedERIs(i,j).reserve();

                if (LPL==NINTS) {
                    return LPL;
                }
            }

        }


    }

    if(G) State.storeAB = tmax12;

    return LPL;
}

template<bool G, bool J, bool X, bool R, int mode> UI64 ERIblock::PrescreenedList(float logD, float logS,  const r2tensor< r2tensor<float> > & logFs, const r2tensor< r2tensor<float> > & logDs, UI64 NINTS) {

    UI64 LPL = 0;


    UI16 minnab = State.storeABb;
    UI16 maxnab = K2max12;

    for (UI16 nab=minnab; nab<maxnab; ++nab) {
        UI16 minncd = (nab==State.storeABb)?State.storeCDb:0;
        UI16 maxncd = K2max34;

        for (UI16 ncd=minncd; ncd<maxncd; ++ncd) {

            UI32 Aab = SortedABs[nab].vList;
            UI32 Acd = SortedCDs[ncd].vList;

            UI32 minab = (nab==State.storeABb && ncd==State.storeCDb)?State.storeAB:0;
            UI32 maxab = SortedABs[nab].Nsps;

            UI16 i = nab; //K2max12-SP12[ab].nK2; //>0
            UI16 j = ncd; //K2max34-SP34[cd].nK2; //>0

            //skips in none of the shellpair pairs survives CS prescreening
            //if (Acd >= max34l[Aab+maxab-1]) continue;

            for (UI32 abp=minab; abp<maxab; ++abp) {
                UI32 ab = Aab + abp;

                if (!ThisNode.CheckMod(ab)) continue;

                UI32 mincd = (nab==State.storeABb && ncd==State.storeCDb && abp==State.storeAB)?State.storeCD:0;
                UI32 maxcd = SortedCDs[ncd].Nsps;

                if (Acd >= max34l[ab]) continue;


                for (UI32 cdp=mincd; cdp<maxcd; ++cdp) {
                    UI32 cd = Acd + cdp;

                    if (cd >= max34l[ab]) continue;
                    if (SkipSame && ab==cd) continue;

                    UI32 ata = SP12[ab].ata;
                    UI32 atb = SP12[ab].atb;
                    UI32 atc = SP34[cd].ata;
                    UI32 atd = SP34[cd].atb;

                    UI16  fa = ABp->nb1;
                    UI16  fb = ABp->nb2;
                    UI16  fc = CDp->nb1;
                    UI16  fd = CDp->nb2;

                    if (ata==atc || ata==atd || atb==atc || atb==atd) {
                        // skip it if computing strictly 4-center integrals
                        if (mode==2) continue;
                    }
                    else {
                        // skip it if computing strictly 3-center integrals
                        if (mode==1) continue;
                    }


                    float wabcd = SP12[ab].logCS + SP34[cd].logCS;


                    float logJa;
                    float logJb;

                    float logX1a;
                    float logX2a;
                    float logX1b;
                    float logX2b;

                    bool useJD, useJS, useXD, useXS;
                    useJD = useJS = useXD = useXS = false;

                    if (J) {
                        logJa  = logFs[ata][atb][fa][fb] + logDs[atc][atd][fc][fd] + wabcd;
                        logJb  = logDs[ata][atb][fa][fb] + logFs[atc][atd][fc][fd] + wabcd;

                        useJD = (logJa<logD || logJb<logD);
                        useJS = (logJa<logS || logJb<logS);
                    }
                    if (X) {
                        logX1a = logFs[ata][atc][fa][fc] + logDs[atb][atd][fb][fd] + wabcd;
                        logX2a = logFs[ata][atd][fa][fd] + logDs[atb][atc][fb][fc] + wabcd;
                        logX1b = logDs[ata][atc][fa][fc] + logFs[atb][atd][fb][fd] + wabcd;
                        logX2b = logDs[ata][atd][fa][fd] + logFs[atb][atc][fb][fc] + wabcd;

                        useXD = (logX1a<logD || logX2a<logD || logX1b<logD || logX2b<logD);
                        useXS = (logX1a<logS || logX2a<logS || logX1b<logS || logX2b<logS);
                    }

                    //if generate
                    if (G) {
                        if (include(J,X,R,useJD,useJS,useXD,useXS)) {
                            ShellPairPair SPP(ab,cd);
                            SortedERIs(i,j).push(SPP);


                            ++LPL;

                            if (LPL==NINTS) {
                                State.storeABb = nab;
                                State.storeCDb = ncd;
                                State.storeAB = abp; //save the state for the next iteration
                                State.storeCD = cdp+1;
                                return LPL;
                            }
                        }
                    }
                    // else only count
                    else {
                        if (include(J,X,R,useJD,useJS,useXD,useXS)) {
                            ++LPL;
                            SortedERIs(i,j).reserve();

                            if (LPL==NINTS) {
                                return LPL;
                            }
                        }
                    }
                }
            }
        }
    }

    if (G) {
        State.storeABb = K2max12;
        State.storeCDb = K2max34;
        State.storeAB = tmax12;
        State.storeCD = tmax34;
    }

    return LPL;
}



void ERIblock::InitONX(const r2tensor< r2tensor<float> > & logFs, const r2tensor< r2tensor<float> > & logDs) {

    int natoms = ONXtensors.natoms = logDs.N();

    //make a list of ata, atb, atc, atd
    for (UI32 ab=0; ab<tmax12; ++ab) {
       if (!ThisNode.CheckMod(ab)) continue;

        UI32 ata = SP12[ab].ata;
        UI32 atb = SP12[ab].atb;

        ONXtensors.AtomsA.insert(ata);
        ONXtensors.AtomsB.insert(atb);
    }

    for (UI32 cd=0; cd<tmax34; ++cd) {
        UI32 atc = SP34[cd].ata;
        UI32 atd = SP34[cd].atb;

        ONXtensors.AtomsC.insert(atc);
        ONXtensors.AtomsD.insert(atd);
    }


    //AtomsAB.set(natoms);
    ONXtensors.AtomsCD.set(natoms);
    ONXtensors.AtomsAC.set(natoms);
    ONXtensors.AtomsAD.set(natoms);
    ONXtensors.AtomsBC.set(natoms);
    ONXtensors.AtomsBD.set(natoms);

    set<UI32>::iterator it1, it2;

    UI16  fa = ABp->nb1;
    UI16  fb = ABp->nb2;
    UI16  fc = CDp->nb1;
    UI16  fd = CDp->nb2;

    float maxx = log(1e99)/2;

/*
    for (UI32 ab=0; ab<tmax12; ++ab) {
        UI32 ata = SP12[ab].ata;
        UI32 atb = SP12[ab].atb;

        AtomsAB(ata,atb) = ab; //SP12[ab].logCS;
    }
*/
    for (UI32 cd=0; cd<tmax34; ++cd) {
        UI32 atc = SP34[cd].ata;
        UI32 atd = SP34[cd].atb;

        ONXtensors.AtomsCD(atc,atd) = cd; // SP34[cd].logCS;
    }

    for (it1=ONXtensors.AtomsA.begin(); it1!=ONXtensors.AtomsA.end(); ++it1) {
        for (it2=ONXtensors.AtomsC.begin(); it2!=ONXtensors.AtomsC.end(); ++it2) {
            UI32 ata = *it1;
            UI32 atc = *it2;
            if (logDs[ata][atc][fa][fc] < maxx)
                ONXtensors.AtomsAC(ata,atc) = logDs[ata][atc][fa][fc];
        }
    }

    for (it1=ONXtensors.AtomsB.begin(); it1!=ONXtensors.AtomsB.end(); ++it1) {
        for (it2=ONXtensors.AtomsD.begin(); it2!=ONXtensors.AtomsD.end(); ++it2) {
            UI32 atb = *it1;
            UI32 atd = *it2;
            if (logDs[atb][atd][fb][fd] < maxx)
                ONXtensors.AtomsBD(atb,atd) = logDs[atb][atd][fb][fd];
        }
    }

    for (it1=ONXtensors.AtomsA.begin(); it1!=ONXtensors.AtomsA.end(); ++it1) {
        for (it2=ONXtensors.AtomsD.begin(); it2!=ONXtensors.AtomsD.end(); ++it2) {
            UI32 ata = *it1;
            UI32 atd = *it2;
            if (logDs[ata][atd][fa][fd] < maxx)
                ONXtensors.AtomsAD(ata,atd) = logDs[ata][atd][fa][fd];
        }
    }

    for (it1=ONXtensors.AtomsB.begin(); it1!=ONXtensors.AtomsB.end(); ++it1) {
        for (it2=ONXtensors.AtomsC.begin(); it2!=ONXtensors.AtomsC.end(); ++it2) {
            UI32 atb = *it1;
            UI32 atc = *it2;
            if (logDs[atb][atc][fb][fc] < maxx)
                ONXtensors.AtomsBC(atb,atc) = logDs[atb][atc][fb][fc];
        }
    }


    float minBC =  ONXtensors.AtomsBC.GetMin();
    float minAD =  ONXtensors.AtomsAD.GetMin();
    float minAC =  ONXtensors.AtomsAC.GetMin();
    float minBD =  ONXtensors.AtomsBD.GetMin();

    //cout << "1" << endl; cout.flush();
}

void ERIblock::CleanONX() {
    ONXtensors.AtomsA.clear();
    ONXtensors.AtomsB.clear();
    ONXtensors.AtomsC.clear();
    ONXtensors.AtomsD.clear();
    ONXtensors.AtomsAB.clean();
    ONXtensors.AtomsCD.clean();
    ONXtensors.AtomsAC.clean();
    ONXtensors.AtomsBD.clean();
    ONXtensors.AtomsAD.clean();
    ONXtensors.AtomsBC.clean();
}

template<bool G, int mode> UI64 ERIblock::ONXPrescreenedList(float logD, float logS,  const r2tensor< r2tensor<float> > & logFs, const r2tensor< r2tensor<float> > & logDs, UI64 NINTS) {

    //cout << "2" << endl; cout.flush();

    map<UI32,UI32>::iterator itab, itcd;
    map<UI32,float>::iterator ita, itb, itc, itd;


    float minAB =  SP12[0].logCS;
    float minCD =  SP34[0].logCS;

    UI64 LPL = 0;


    for (UI32 ab=State.storeAB; ab<tmax12; ++ab) {
        if (!ThisNode.CheckMod(ab)) continue;

        UI32 ata = SP12[ab].ata;
        UI32 atb = SP12[ab].atb;

        float lAB = SP12[ab].logCS;

        float minBC = ONXtensors.AtomsBC[atb].minf;
        float minAD = ONXtensors.AtomsAD[ata].minf;

        float minAC = ONXtensors.AtomsAC[ata].minf;
        float minBD = ONXtensors.AtomsBD[atb].minf;

        if (lAB+minBC+minCD+minAD>=logD) break;

        //eliminate duplicates in SortedERIs and keeps it sorted
        set<UI32> cds;

        //make a list

        for (itc=ONXtensors.AtomsBC[atb].r.begin(); itc!=ONXtensors.AtomsBC[atb].r.end(); ++itc) {
            UI32 atc = itc->first;

            if (mode==2) if (ata==atc || atb==atc) continue; // skip it if computing strictly 4-center integrals

            float lBC = itc->second;

            if (lAB+lBC+minCD+minAD>=logD) continue;

            for (itcd=ONXtensors.AtomsCD[atc].r.begin(); itcd!=ONXtensors.AtomsCD[atc].r.end(); ++itcd) {
                UI32 atd = itcd->first;

                if (mode==2) if   (ata==atd || atb==atd) continue; // skip it if computing strictly 4-center integrals
                if (mode==1) if (!(ata==atc || ata==atd || atb==atc || atb==atd)) continue; // skip it if computing strictly 3-center integrals

                UI32 cd = itcd->second;

                if (!ThisNode.CheckMod(ab+cd)) continue;
                if (cd >= max34l[ab]) continue;
                if (SkipSame && ab==cd) continue;

                float lCD = SP34[cd].logCS;
                float lAD = ONXtensors.AtomsAD[ata][atd];

                if (lAB+lBC+lCD+lAD < logD) {
                    cds.insert(cd);
                }
            }
        }


        for (itc=ONXtensors.AtomsAC[ata].r.begin(); itc!=ONXtensors.AtomsAC[ata].r.end(); ++itc) {
            UI32 atc = itc->first;

            if (mode==2) if (ata==atc || atb==atc) continue; // skip it if computing strictly 4-center integrals

            float lAC = itc->second;

            if (lAB+lAC+minCD+minBD>=logD) continue;

            for (itcd=ONXtensors.AtomsCD[atc].r.begin(); itcd!=ONXtensors.AtomsCD[atc].r.end(); ++itcd) {
                UI32 atd = itcd->first;

                if (mode==2) if   (ata==atd || atb==atd) continue; // skip it if computing strictly 4-center integrals
                if (mode==1) if (!(ata==atc || ata==atd || atb==atc || atb==atd)) continue; // skip it if computing strictly 3-center integrals

                UI32 cd = itcd->second;

                if (!ThisNode.CheckMod(ab+cd)) continue;
                if (cd >= max34l[ab]) continue;
                if (SkipSame && ab==cd) continue;

                float lCD = SP34[cd].logCS;
                float lBD = ONXtensors.AtomsBD[atb][atd];

                if (lAB+lAC+lCD+lBD < logD) {
                    cds.insert(cd);
                }
            }
        }


        set<UI32>::iterator it = cds.begin();

        UI32 sCD = (ab==State.storeAB)?State.storeCD:0;

        // skip the elements that have been computed
        for (int ii=0; ii<sCD; ++ii) ++it;


        for (; it!=cds.end(); ++it, ++sCD) {

            UI32 cd = *it;

            UI16 i = K2max12 - SP12[ab].nK2;
            UI16 j = K2max34 - SP34[cd].nK2;


            //if generate
            if (G) {
                //if (include(J,X,R,useJD,useJS,useXD,useXS))
                {
                    ShellPairPair SPP(ab,cd);
                    SortedERIs(i,j).push(SPP);
                    ++LPL;

                    if (LPL==NINTS) {
                        State.storeAB = ab;
                        State.storeCD = sCD +1;
                        return LPL;
                    }
                }
            }

            // else only count
            else {
                //if (include(J,X,R,useJD,useJS,useXD,useXS))
                {
                    ++LPL;
                    SortedERIs(i,j).reserve();

                    if (LPL==NINTS) {
                        return LPL;
                    }
                }
            }
        }
    }

    if (G) {
        State.storeABb = K2max12;
        State.storeCDb = K2max34;
        State.storeAB  = tmax12;
        State.storeCD  = tmax34;
    }

    return LPL;
}

template<bool G, int mode> UI64 ERIblock::JPrescreenedList(float logD, float logS,  const r2tensor< r2tensor<float> > & logFs, const r2tensor< r2tensor<float> > & logDs, UI64 NINTS) {

    UI64 LPL = 0;


    UI16 minnab = State.storeABb;
    UI16 maxnab = K2max12;

    for (UI16 nab=minnab; nab<maxnab; ++nab) {
        UI16 minncd = (nab==State.storeABb)?State.storeCDb:0;
        UI16 maxncd = K2max34;

        for (UI16 ncd=minncd; ncd<maxncd; ++ncd) {

            UI32 Aab = SortedABs[nab].vList;
            UI32 Acd = SortedCDs[ncd].vList;

            UI32 minab = (nab==State.storeABb && ncd==State.storeCDb)?State.storeAB:0;
            UI32 maxab = SortedABs[nab].Nsps;

            UI16 i = nab; //K2max12-SP12[ab].nK2; //>0
            UI16 j = ncd; //K2max34-SP34[cd].nK2; //>0

            //skips in none of the shellpair pairs survives CS prescreening
            //if (Acd >= max34l[Aab+maxab-1]) continue;


            for (UI32 abp=minab; abp<maxab; ++abp) {
                UI32 ab = Aab + abp;

                if (!ThisNode.CheckMod(ab)) continue;

                UI32 mincd = (nab==State.storeABb && ncd==State.storeCDb && abp==State.storeAB)?State.storeCD:0;
                UI32 maxcd = SortedCDs[ncd].Nsps;

                if (Acd >= max34l[ab]) continue;

                UI32 ata = SP12[ab].ata;
                UI32 atb = SP12[ab].atb;
                UI16  fa = ABp->nb1;
                UI16  fb = ABp->nb2;

                float logAB = SP12[ab].logCS + logDs[ata][atb][fa][fb];


                for (UI32 cdp=mincd; cdp<maxcd; ++cdp) {
                    UI32 cd = Acd + cdp;

                    if (cd >= max34l[ab]) continue;
                    if (SkipSame && ab==cd) continue;

                    UI32 atc = SP34[cd].ata;
                    UI32 atd = SP34[cd].atb;

                    if (ata==atc || ata==atd || atb==atc || atb==atd) {
                        // skip it if computing strictly 4-center integrals
                        if (mode==2) continue;
                    }
                    else {
                        // skip it if computing strictly 3-center integrals
                        if (mode==1) continue;
                    }

                    UI16  fc = CDp->nb1;
                    UI16  fd = CDp->nb2;

                    float logJ = logAB + logDs[atc][atd][fc][fd] +  SP34[cd].logCS;

                    //if generate
                    if (G) {
                        if (logJ<logD) {
                            ShellPairPair SPP(ab,cd);
                            SortedERIs(i,j).push(SPP);


                            ++LPL;

                            if (LPL==NINTS) {
                                State.storeABb = nab;
                                State.storeCDb = ncd;
                                State.storeAB = abp; //save the state for the next iteration
                                State.storeCD = cdp+1;
                                return LPL;
                            }
                        }
                    }
                    // else only count
                    else {
                        if (logJ<logD) {
                            ++LPL;
                            SortedERIs(i,j).reserve();

                            if (LPL==NINTS) {
                                return LPL;
                            }
                        }
                    }
                }
            }
        }
    }

    if (G) {
        State.storeABb = K2max12;
        State.storeCDb = K2max34;
        State.storeAB = tmax12;
        State.storeCD = tmax34;
    }

    return LPL;
}


//this is just a quicksort implementation
template <class T, class K> void FlashSort(T * keys, K * elements, int N) {

    //for short lists default to selection sort
    if (N<32) {
        for (int i=0; i<N; ++i)  {
            T best = keys[i];
            int   k   = i;

            for (int j=i+1; j<N; ++j)
                if (keys[j]<best) {best = keys[j]; k = j;}

            swap(keys[i]    , keys[k]);
            swap(elements[i], elements[k]);
        }
        return;
    }

    T & k1 = keys[0];
    T & k2 = keys[N/2];
    T & k3 = keys[N-1];
    T kpiv;

    if      (k1<=k2 && k2<=k3) kpiv = k2;
    else if (k3<=k2 && k2<=k1) kpiv = k2;
    else if (k1<=k3 && k3<=k2) kpiv = k3;
    else if (k2<=k3 && k3<=k1) kpiv = k3;
    else                       kpiv = k1;
    //else if (k2<k1 && k1<k3) kpiv = k1;
    //else if (k3<k1 && k1<k2) kpiv = k1;


    //count number of elements lower than the pivot
    int nl=0;
    int ng=0;
    for (int i=0; i<N; ++i) {
        if      (keys[i] < kpiv) ++nl;
        else if (keys[i] > kpiv) ++ng;
    }

    int i,j,k;


    i=0;
    j=nl;

    while (1) {
        //find the first useful place in the first half of the list
        while (i<nl && keys[i]<kpiv) ++i;
        if (i>=nl) break;
        while (j<N  && keys[j]>=kpiv) ++j;
        if (j>=N) break;

        swap(keys[i]    , keys[j]);
        swap(elements[i], elements[j]);
    }

    j = nl;
    k = N-ng;

    while (1) {
        //find the first useful place in the first half of the list
        while (j<N-ng && keys[j]==kpiv) ++j;
        if (j>=N-ng) break;
        while (k<N  && keys[k]>kpiv) ++k;
        if (k>=N) break;

        swap(keys[j]    , keys[k]);
        swap(elements[j], elements[k]);
    }


    //FlashSort both sublists (skipping the pivot)
    FlashSort(keys         ,  elements         ,  nl);

    FlashSort(keys + (N-ng),  elements + (N-ng),  ng);
}

template <class T> void FlashSort(T * keys, int N) {

    //for short lists default to selection sort
    if (N<32) {
        for (int i=0; i<N; ++i)  {
            T best = keys[i];
            int   k   = i;

            for (int j=i+1; j<N; ++j)
                if (keys[j]<best) {best = keys[j]; k = j;}

            swap(keys[i]    , keys[k]);
        }
        return;
    }

    T & k1 = keys[0];
    T & k2 = keys[N/2];
    T & k3 = keys[N-1];
    T kpiv;

    if      (k1<=k2 && k2<=k3) kpiv = k2;
    else if (k3<=k2 && k2<=k1) kpiv = k2;
    else if (k1<=k3 && k3<=k2) kpiv = k3;
    else if (k2<=k3 && k3<=k1) kpiv = k3;
    else                       kpiv = k1;
    //else if (k2<k1 && k1<k3) kpiv = k1;
    //else if (k3<k1 && k1<k2) kpiv = k1;


    //count number of elements lower than the pivot
    int nl=0;
    int ng=0;
    for (int i=0; i<N; ++i) {
        if      (keys[i] < kpiv) ++nl;
        else if (keys[i] > kpiv) ++ng;
    }

    int i,j,k;


    i=0;
    j=nl;

    while (1) {
        //find the first useful place in the first half of the list
        while (i<nl && keys[i]<kpiv) ++i;
        if (i>=nl) break;
        while (j<N  && keys[j]>=kpiv) ++j;
        if (j>=N) break;

        swap(keys[i]    , keys[j]);
    }

    j = nl;
    k = N-ng;

    while (1) {
        //find the first useful place in the first half of the list
        while (j<N-ng && keys[j]==kpiv) ++j;
        if (j>=N-ng) break;
        while (k<N  && keys[k]>kpiv) ++k;
        if (k>=N) break;

        swap(keys[j]    , keys[k]);
    }


    //FlashSort both sublists (skipping the pivot)
    FlashSort(keys         ,  nl);
    FlashSort(keys + (N-ng),  ng);
}


void ERIblock::SortPrescreenedLists() {

    // sort the shell quadruplets according to distance
    // ************************************************

    for (UI32 i=0; i<K2max12; ++i) {
        for (UI32 j=0; j<K2max34; ++j) {

            int N = SortedERIs(i,j).size();

            if (N<2) continue;

            ShellPairPair * SPPblock = SortedERIs(i,j).vList;


            //compute the (squared) distances
            for (int n=0; n<N; ++n) {
                UI32 ab = SPPblock[n].n12;
                UI32 cd = SPPblock[n].n34;

                const AtomProd & AP12 = (*AP12list)[ ab ];
                const AtomProd & AP34 = (*AP34list)[ cd ];

                R2[n] = R2distance(AP12, AP34);
            }

            FlashSort(R2, SPPblock, N);
        }
    }

}

void ERIblock::CaseSortPrescreenedLists() {
    // sort the shell quadruplets according to distance; discard elements far apart
    // ****************************************************************************

    const int NVIG = 1024;
    const double vg_min  = 0;
    const double vg_max  = 32;
    const double vg_step  = (vg_max-vg_min)/double(NVIG);


    //correct, but extremely conservative
    double maxR2 = (vg_max-vg_step)*(1./ABp->kmin + 1./CDp->kmin + 1./case_w2);

    UI64 N0 = 0;
    UI64 NF = 0;


    for (UI32 i=0; i<K2max12; ++i) {
        for (UI32 j=0; j<K2max34; ++j) {

            UI64 N = SortedERIs(i,j).Nints;

            N0 += N;
            //if (N<2) continue;

            ShellPairPair * SPPblock = SortedERIs(i,j).vList;

            double * R2 = new double[N];

            //compute the (squared) distances
            for (int n=0; n<N; ++n) {
                UI32 ab = SPPblock[n].n12;
                UI32 cd = SPPblock[n].n34;

                const AtomProd & AP12 = (*AP12list)[ ab ];
                const AtomProd & AP34 = (*AP34list)[ cd ];

                R2[n] = R2distance(AP12, AP34);
            }

            FlashSort(R2, SPPblock, N);

            //find the last allowed element
            int last=0;
            for (; last<N; ++last) if (R2[last]>maxR2) break;

            //set the list shorter
            SortedERIs(i,j).Nints = last;
            NF += last;


            delete[] R2;
        }
    }

    //Benchmarker << "Initial CASE Coulomb list contained " << N0 << " elements" << endl;
    //Benchmarker << " prescreened list contains "  << NF << " elements" << endl;
    //Benchmarker << " cutoff distance for CASE with present pair of prototypes is " << sqrt(maxR2)*BOHR << " Angstroms" << endl << endl;
}

template <int LEN> void ERIblock::InitBlocks() {

    //set list pointers
    UI64 Icount = 0;

    for (UI32 i=0; i<K2max12; ++i) {
        for (UI32 j=0; j<K2max34; ++j) {
            SortedERIs(i,j).vList = SPPlist + LEN*Icount;
            int nblocks = (SortedERIs(i,j).Nints + (LEN-1)) / LEN;
            SortedERIs(i,j).Nints = nblocks*LEN;
            Icount += nblocks;
        }
    }
}

void ERIblock::InitLists() {

    for (UI32 ab=0; ab<tmax12; ++ab) {
        UI16 i = K2max12-SP12[ab].nK2; //>0
        ++SortedABs[i].Nsps;
    }

    for (UI32 cd=0; cd<tmax34; ++cd) {
        UI16 j = K2max34-SP34[cd].nK2; //>0
        ++SortedCDs[j].Nsps;
    }

    UI32 Icount12 = 0;
    UI32 Icount34 = 0;

    for (UI16 p=0; p<K2max12; ++p) {
        SortedABs[p].vList = Icount12;
        Icount12 += SortedABs[p].Nsps;
    }

    for (UI16 p=0; p<K2max34; ++p) {
        SortedCDs[p].vList = Icount34;
        Icount34 += SortedCDs[p].Nsps;
    }
}


void ERIblock::InitListsJ(const r2tensor< r2tensor<float> > & logDs) {

    UI16  fa = ABp->nb1;
    UI16  fb = ABp->nb2;

    for (UI32 ab=0; ab<tmax12; ++ab) {
        UI32 ata = SP12[ab].ata;
        UI32 atb = SP12[ab].atb;
        ABlogs[ab].cs = SP12[ab].logCS + logDs[ata][atb][fa][fb];
        ABlogs[ab].n  = ab;
    }

    UI16  fc = CDp->nb1;
    UI16  fd = CDp->nb2;

    for (UI32 cd=0; cd<tmax34; ++cd) {
        UI32 atc = SP34[cd].ata;
        UI32 atd = SP34[cd].atb;
        CDlogs[cd].cs = SP34[cd].logCS + logDs[atc][atd][fc][fd];
        CDlogs[cd].n  = cd;
    }

    //sort both lists
    FlashSort(ABlogs, tmax12);
    FlashSort(CDlogs, tmax34);

    for (UI32 ab=0; ab<tmax12; ++ab) {
        UI16 i = K2max12-SP12[ab].nK2; //>0
        ++SortedABs[i].Nsps;
    }

    for (UI32 cd=0; cd<tmax34; ++cd) {
        UI16 j = K2max34-SP34[cd].nK2; //>0
        ++SortedCDs[j].Nsps;
    }

    UI32 Icount12 = 0;
    UI32 Icount34 = 0;

    // these will now be offset, but no shellpair will demote

    for (UI16 p=0; p<K2max12; ++p) {
        SortedABs[p].vList = Icount12;
        Icount12 += SortedABs[p].Nsps;
    }

    for (UI16 p=0; p<K2max34; ++p) {
        SortedCDs[p].vList = Icount34;
        Icount34 += SortedCDs[p].Nsps;
    }
}


//efficient integral vector packing algorithm
template <int LEN> void ERIblock::FuseIncompleteTiles() {
    //fuse incomplete tiles
    //*********************

    for (int e=1; e<LEN; e*=2) {

        //a better algorithm should eliminate with the closest element in the table class, not just the next in the row/column
        //also: just check active rows and columns

        //eliminate by rows
        for (int p=0; p<SortedERIs.N(); ++p) {
            int qq;
            bool first=true;

            for (int q=0; q<SortedERIs.M(); ++q) {

                if (SortedERIs(p,q).size()&e) {
                    if (first) {
                        if (SortedERIs(p,q).size() + e <= SortedERIs(p,q).capacity()) {
                            qq = q;
                            first = false;
                        }
                    }
                    else {
                        //upgrade 'e' elements to a higher class and update
                        for (int ee=0; ee<e; ++ee) {
                            ShellPairPair SP = SortedERIs(p,q).pop();
                            SortedERIs(p,qq).push(SP);
                        }
                        first = true;
                    }
                }
            }
        }

        //eliminate by columns
        for (int q=0; q<SortedERIs.M(); ++q) {
            int pp;
            bool first=true;

            for (int p=0; p<SortedERIs.N(); ++p) {
                if (SortedERIs(p,q).size()&e) {
                    if (first) {
                        //adding 'e' integrals to this class should not overflow it
                        if (SortedERIs(p,q).size() + e <= SortedERIs(p,q).capacity()) {
                            pp = p;
                            first = false;
                        }
                    }
                    else {
                        //upgrade 'e' elements to a higher class and update
                        for (int ee=0; ee<e; ++ee) {
                            ShellPairPair SP = SortedERIs(p,q).pop();
                            SortedERIs(pp,q).push(SP);
                        }
                        first = true;
                    }
                }
            }
        }



        //eliminate permitted diagonals (without generating new tiles)
        UI32 LastInRow[maxK2];
        UI32 LastAccepting[maxK2];

        for (int p=0; p<SortedERIs.N(); ++p) {
            int last = -1;
            for (int q=0; q<SortedERIs.M(); ++q)
                if (SortedERIs(p,q).size()&e) last = q;

            LastInRow[p] = last;

            if ( last>-1 && SortedERIs(p,last).size()+e<=SortedERIs(p,last).capacity() )
                LastAccepting[p] = last;
            else
                LastAccepting[p] = -1;
        }

        for (int p=SortedERIs.N()-1; p>=0; --p) {
            if (LastInRow[p]==-1) continue;

            int pp;
            for (pp=p-1; pp>=0; --pp) {
                if (LastAccepting[pp]>-1 && LastAccepting[pp]<LastInRow[p]) break;
            }

            if (pp==-1) continue;

            //merge
            int q  = LastInRow[p];
            int qq = LastAccepting[pp];

            for (int ee=0; ee<e; ++ee) {
                ShellPairPair SP = SortedERIs(p,q).pop();
                SortedERIs(pp,qq).push(SP);
            }

            LastInRow[p] = LastInRow[pp] = LastAccepting[p] = LastAccepting[pp] = -1;
        }

    }


    // more sophisticated alternatives are possible, but they come at a cost
    // simply move all orphan integrals together to the bucket which minimally
    // satisfies the required conditions

    int cup  = 0;
    int minp = SortedERIs.N();
    int minq = SortedERIs.M();

    for (int p=0; p<SortedERIs.N(); ++p) {
        for (int q=0; q<SortedERIs.M(); ++q) {

            if (SortedERIs(p,q).size()%LEN != 0) {
                cup += SortedERIs(p,q).size()%LEN;
                minp = min(minp, p);
                minq = min(minq, q);
            }
        }
    }

    if (cup==0) return; // no leftovers; return with glory

    cup -= SortedERIs(minp,minq).size()%LEN; //don't count it twice!!

    //merge the list and the leftovers in a unique array
    int totmin = cup + SortedERIs(minp,minq).size();

    ShellPairPair leftovers[totmin];
    int ncup  = 0;

    while (SortedERIs(minp,minq).size() > 0) {
        leftovers[ncup] = SortedERIs(minp,minq).pop();
        ++ncup;
    }

    for (int p=0; p<SortedERIs.N(); ++p) {
        for (int q=0; q<SortedERIs.M(); ++q) {
            //if (p==minp && q==minq) continue; //doesn't matter

            while (SortedERIs(p,q).size()%LEN != 0) {
                leftovers[ncup] =  SortedERIs(p,q).pop();
                ++ncup;
            }
        }
    }



    // compact everything else!!!
    ShellPairPair * sppl = SPPlist;

    for (int p=0; p<SortedERIs.N(); ++p) {
        for (int q=0; q<SortedERIs.M(); ++q) {
            if (p==minp && q==minq) continue; //leave for last

            for (int n=0; n<SortedERIs(p,q).size(); ++n)
                sppl[n] = SortedERIs(p,q).vList[n];

            SortedERIs(p,q).vList = sppl;

            sppl += SortedERIs(p,q).size();
        }
    }

    // append the merged list at the end of the array
    for (int n=0; n<totmin; ++n) sppl[n] = leftovers[n];

    SortedERIs(minp,minq).Iints = totmin;
    SortedERIs(minp,minq).Nints = totmin;
    SortedERIs(minp,minq).vList = sppl;

}

// make the tiles
// **************

template <int LEN> UI64 ERIblock::GetNTiles() const {

    UI64 NNtiles = 0;

    for (UI16 i=0; i<SortedERIs.N(); ++i) {
        for (UI16 j=0; j<SortedERIs.M(); ++j) {

            UI64 Nints = SortedERIs(i,j).size();

            //set same contraction header for all tiles
            NNtiles += Nints/LEN;

            //fill remnant
            if (Nints%LEN) ++NNtiles;
        }
    }

    return NNtiles;
}

UI64 ERIblock::MakeAllTiles(LibQuimera::ERITile64 * pTile) {

    //compute the offsets for later digestion step
    UI32 ta = SP12->ata;
    UI32 tb = SP12->atb;
    UI32 tc = SP34->ata;
    UI32 td = SP34->atb;


    UI32 offAB = sparsetensorpattern::GetOffset(ta, tb, ABp->nb1, ABp->nb2);
    UI32 offCD = sparsetensorpattern::GetOffset(tc, td, CDp->nb1, CDp->nb2);

    UI32 offAC = sparsetensorpattern::GetOffset(ta, tc, ABp->nb1, CDp->nb1);
    UI32 offAD = sparsetensorpattern::GetOffset(ta, td, ABp->nb1, CDp->nb2);
    UI32 offBC = sparsetensorpattern::GetOffset(tb, tc, ABp->nb2, CDp->nb1);
    UI32 offBD = sparsetensorpattern::GetOffset(tb, td, ABp->nb2, CDp->nb2);


    UI64 Ntiles = 0;

    for (UI16 i=0; i<SortedERIs.N(); ++i) {
        for (UI16 j=0; j<SortedERIs.M(); ++j) {

            UI64 Nints = SortedERIs(i,j).size();

            //set same contraction header for all tiles
            UI64 t;
            for (t=0; t<Nints/DPC; ++t) {

                pTile[Ntiles].nKab = SortedERIs.N()-i-1;
                pTile[Ntiles].nKcd = SortedERIs.M()-j-1;
                pTile[Ntiles].used = DPC;

                for (int k=0; k<DPC; ++k) {
                    pTile[Ntiles].ap12[k] = SortedERIs(i,j).vList[DPC*t+k].n12;
                    pTile[Ntiles].ap34[k] = SortedERIs(i,j).vList[DPC*t+k].n34;
                    pTile[Ntiles].use[k] = 1;
                }

                ++Ntiles;
            }

            //fill remnant
            if (Nints%DPC) {

                pTile[Ntiles].nKab = SortedERIs.N()-i-1;
                pTile[Ntiles].nKcd = SortedERIs.M()-j-1;
                pTile[Ntiles].used = Nints%DPC;

                for (int k=0; k<Nints%DPC; ++k) {
                    pTile[Ntiles].ap12[k] = SortedERIs(i,j).vList[DPC*t+k].n12;
                    pTile[Ntiles].ap34[k] = SortedERIs(i,j).vList[DPC*t+k].n34;

                    pTile[Ntiles].use[k] = 1;
                }

                for (int k=Nints%DPC; k<DPC; ++k) {
                    pTile[Ntiles].ap12[k] = SortedERIs(i,j).vList[DPC*t].n12;
                    pTile[Ntiles].ap34[k] = SortedERIs(i,j).vList[DPC*t].n34;

                    pTile[Ntiles].use[k] = 0;
                }

                ++Ntiles;
            }
        }
    }

    for (UI64 n=0; n<Ntiles; ++n) {
        pTile[n].wSP12 = SP12->WW;
        pTile[n].wSP34 = SP34->WW;
        pTile[n].nK2ab = ABp->Ka * ABp->Kb;
        pTile[n].nK2cd = CDp->Ka * CDp->Kb;

        pTile[n].offAB = offAB;
        pTile[n].offCD = offCD;
        pTile[n].offAC = offAC;
        pTile[n].offAD = offAD;
        pTile[n].offBC = offBC;
        pTile[n].offBD = offBD;

        pTile[n].lowgeom = geometry;
    }

    return Ntiles;
}

UI64 ERIblock::MakeAllTiles(LibQuimera::ERITile32 * pTile) {

    //compute the offsets for later digestion step
    UI32 ta = SP12->ata;
    UI32 tb = SP12->atb;
    UI32 tc = SP34->ata;
    UI32 td = SP34->atb;


    UI32 offAB = sparsetensorpattern::GetOffset(ta, tb, ABp->nb1, ABp->nb2);
    UI32 offCD = sparsetensorpattern::GetOffset(tc, td, CDp->nb1, CDp->nb2);

    UI32 offAC = sparsetensorpattern::GetOffset(ta, tc, ABp->nb1, CDp->nb1);
    UI32 offAD = sparsetensorpattern::GetOffset(ta, td, ABp->nb1, CDp->nb2);
    UI32 offBC = sparsetensorpattern::GetOffset(tb, tc, ABp->nb2, CDp->nb1);
    UI32 offBD = sparsetensorpattern::GetOffset(tb, td, ABp->nb2, CDp->nb2);


    UI64 Ntiles = 0;

    for (UI16 i=0; i<SortedERIs.N(); ++i) {
        for (UI16 j=0; j<SortedERIs.M(); ++j) {

            UI64 Nints = SortedERIs(i,j).size();

            //set same contraction header for all tiles
            UI64 t;
            for (t=0; t<Nints/FPC; ++t) {

                pTile[Ntiles].nKab = SortedERIs.N()-i-1;
                pTile[Ntiles].nKcd = SortedERIs.M()-j-1;
                pTile[Ntiles].used = FPC;

                for (int k=0; k<FPC; ++k) {
                    pTile[Ntiles].ap12[k] = SortedERIs(i,j).vList[FPC*t+k].n12;
                    pTile[Ntiles].ap34[k] = SortedERIs(i,j).vList[FPC*t+k].n34;
                    pTile[Ntiles].use[k] = 1;
                }

                ++Ntiles;
            }

            //fill remnant
            if (Nints%FPC) {

                pTile[Ntiles].nKab = SortedERIs.N()-i-1;
                pTile[Ntiles].nKcd = SortedERIs.M()-j-1;
                pTile[Ntiles].used = Nints%FPC;

                for (int k=0; k<Nints%FPC; ++k) {
                    pTile[Ntiles].ap12[k] = SortedERIs(i,j).vList[FPC*t+k].n12;
                    pTile[Ntiles].ap34[k] = SortedERIs(i,j).vList[FPC*t+k].n34;

                    pTile[Ntiles].use[k] = 1;
                }

                for (int k=Nints%FPC; k<FPC; ++k) {
                    pTile[Ntiles].ap12[k] = SortedERIs(i,j).vList[FPC*t].n12;
                    pTile[Ntiles].ap34[k] = SortedERIs(i,j).vList[FPC*t].n34;

                    pTile[Ntiles].use[k] = 0;
                }

                ++Ntiles;
            }
        }
    }

    for (UI64 n=0; n<Ntiles; ++n) {
        pTile[n].wSP12 = SP12->WW;
        pTile[n].wSP34 = SP34->WW;
        pTile[n].nK2ab = ABp->Ka * ABp->Kb;
        pTile[n].nK2cd = CDp->Ka * CDp->Kb;
        pTile[n].offAB = offAB;
        pTile[n].offCD = offCD;
        pTile[n].offAC = offAC;
        pTile[n].offAD = offAD;
        pTile[n].offBC = offBC;
        pTile[n].offBD = offBD;
    }

    return Ntiles;
}




//evaluate blocks of integrals
//****************************

void ERIblock::SetInitialList_ABCD (float lmo) {

    UI64 NElements = 0;

    for (max12=0; max12<tmax12; ++max12) {
        float maxCS = lmo - SP12[max12].logCS;
        UI32 max34 = maxSP( maxCS, SP34, tmax34);
        max34l[max12] = max34;

        NElements += max34;
    }

}

void ERIblock::SetInitialList_ABCDs(float lmo) {
    UI64 NElements = 0;

    for (max12=0; max12<tmax12; ++max12) {
        float maxCS = lmo - SP12[max12].logCS;
        UI32 max34 = maxSP( maxCS, SP34, tmax34);
        max34 = min(max12, max34);
        max34l[max12] = max34;

        NElements += max34;
    }
}

void ERIblock::SetInitialList_ABAB (float lmo) {
    //new method
    UI64 NElements = 0;

    for (UI32 ap12=0; ap12<tmax12; ++ap12) {
        const ShellPair & AB = SP12[ap12];
        const ShellPair & CD = SP34[ap12];
        if ((AB.logCS + CD.logCS > lmo)) {tmax12 = ap12; break;}
    }
    NElements = tmax12;

}

void ERIblock::SetInitialList_AACD (float lmo) {

    UI64 NElements = 0;

    for (max12=0; max12<tmax12; ++max12) {
        float maxCS = lmo - SP12[max12].logCS;
        UI32 max34 = maxSP( maxCS, SP34, tmax34);
        max34l[max12] = max34;

        NElements += max34;
    }

}

void ERIblock::SetInitialList_AACC (float lmo) {
    UI64 NElements = 0;

    for (max12=0; max12<tmax12; ++max12) {
        min34l[max12] = 0;
        max34l[max12] = tmax34;

        NElements += tmax34;
    }

}

void ERIblock::SetInitialList_AACCs(float lmo) {
    UI64 NElements = 0;

    for (max12=0; max12<tmax12; ++max12) {
        min34l[max12] = 0;
        max34l[max12] = max12;

        NElements += max12;
    }

}


void ERIblock::SetInitialList      (float lmo) {

    if       (geometry==ABCD) {
        if (!SameList) SetInitialList_ABCD  (lmo);
        else           SetInitialList_ABCDs (lmo);
    }
    else if (geometry==AACD) {
                        SetInitialList_AACD  (lmo);
    }
    else if (geometry==AACC) {
        if (!SameList) SetInitialList_AACC  (lmo);
        else           SetInitialList_AACCs (lmo);
    }
    else if (geometry==ABAB) {
                       SetInitialList_ABAB  (lmo);
    }

    K2max12 = 0;
    K2max34 = 0;

    for (UI32 ab=0; ab<tmax12; ++ab) K2max12 = max (K2max12, SP12[ab].nK2);
    for (UI32 cd=0; cd<tmax34; ++cd) K2max34 = max (K2max34, SP34[cd].nK2);

    SortedABs.set(K2max12);
    SortedCDs.set(K2max34);

    SortedERIs.set(K2max12, K2max34);
}

//empty everything
void ERIblock::ClearBlocks() {

    for (UI16 i=0; i<K2max12; ++i) {
        for (UI16 j=0; j<K2max34; ++j) {
            SortedERIs(i,j).reset();
        }
    }
}

#include "ERIbatch.hpp"

enum accuracy {LOW, MEDIUM, HIGH};

// S = 3 ABAB two-center integrals
// S = 2 ABCD integrals with some repeated element
// S = 1 ABCD integrals with some repeated element which become 3-center ABAD integrals
// S = 0 rest of integrals
template<bool J, bool X, int S> bool ERIblock::MakeBlocks      (float logD, float logS, const r2tensor< r2tensor<float> > & logAD,  const r2tensor< r2tensor<float> > & logDs, TaskQueue<ERIBatch> & TaskList) {

    // switch between 3-center and 4-center algorithm for CDR
    // necessary to use CDR on the regular ABCD integrals
    if (S==1) ERIalgorithm = Q.SelectAlgorithm(ABAD,la,lb,lc,ld);
    if (S==2) ERIalgorithm = Q.SelectAlgorithm(ABCD,la,lb,lc,ld);


    accuracy acc = HIGH;

    if      (acc==HIGH) {
        //always use double precission
        logD = logS;
    }
    else if (acc==MEDIUM) {
        //always use double precision for 1,2 and 3 center integrals
        if (geometry!=ABCD || (S==1) ) logD = logS;
    }
    else if (acc==LOW) {
        //always use double precision for 1,2 and 3 center integrals
        if (geometry!=ABCD || (S==1) ) logD = logS;
        else {
            logD = -100; //never use double for ABCD integrals
        }
    }



    size_t mem   = ERIalgorithm->MemSize(*ABp, *CDp);
    size_t block = Fock2e::MEMPERBATCH;

    int fa = ABp->Ja*(2*la+1);
    int fb = ABp->Jb*(2*lb+1);
    int fc = CDp->Ja*(2*lc+1);
    int fd = CDp->Jb*(2*ld+1);

    //memory for contractions (constant)
    mem += 2 * (fa*fb + fa*fc + fa*fd + fc*fd + fb*fc + fb*fd) * sizeof(cacheline64);

    if (mem>block) {
        Echidna::EMessenger << "Error trying to allocate more memory for ERI buffers than available; increase program memory or reduce umber of threads" << endl;
        throw (1523);
    }

    //do double packed
    if (!State.storeBlock) {

        SetInitialList(logD);
        InitLists();

        size_t size8  = (wsize4 + msize4) * sizeof(cacheline64) + DPC * sizeof(ERIgeometry) + sizeof(ERIgeometries64) + sizeof(RotationMatrices64); //maximum size (assuming a float tile)

        UI64 max_tiles = min ( UI64( (block-mem) / size8), MAX_TILES_PER_BLOCK);

        // if only exchange is being computed, initialize tensors for ONX
        if (!J && X) InitONX(logAD, logDs);

        //generate batches of ERIs until the buffer mnemory is exhausted or the whole class is finished
        while (1) {

            UI16 sABb = State.storeABb;
            UI16 sCDb = State.storeCDb;
            UI32 sAB  = State.storeAB;
            UI32 sCD  = State.storeCD;

            UI64 Nprescreened;

            //counts the number of shell quartets in each GDO prescreening bucket
            ClearBlocks();

            if (!J && X) {
                if (S==3) Nprescreened = PrescreenedListS    <false, J,X, false  > (logD, logS, logAD, logDs, max_tiles*DPC);
                else      Nprescreened = ONXPrescreenedList  <false, S>            (logD, logS, logAD, logDs, max_tiles*DPC);
            }
            else {
                if (S==3) Nprescreened = PrescreenedListS <false, J,X, false  > (logD, logS, logAD, logDs, max_tiles*DPC);
                else      Nprescreened = PrescreenedList  <false, J,X, false,S> (logD, logS, logAD, logDs, max_tiles*DPC);
            }

            if (Nprescreened==0) break; //no quadruplets left

            UI64 nTiles2Get = (Nprescreened + (DPC-1))/DPC;

            UI64 size = nTiles2Get * sizeof(ERITile64);
            ERITile64 * tList = (ERITile64*)EBB->Allocate(size);

            // no memory at all!!!
            if (tList==NULL) return false;

            UI64 nTilesWeGot = size/sizeof(ERITile64);


            ERIBatch * newbatch  = new ERIBatch;
            *newbatch                 = *this;
            newbatch->TileList64      = tList;
            newbatch->memAllocated    = size;
            newbatch->J               = J;
            newbatch->X               = X;

            // must compute again, with fewer tiles
            if (nTiles2Get!=nTilesWeGot) {
                ClearBlocks();

                if (!J && X) {
                    if (S==3)                    Nprescreened = PrescreenedListS <false, J,X, false  > (logD, logS, logAD, logDs, nTilesWeGot*DPC);
                    else                         Nprescreened = ONXPrescreenedList  <false, S> (logD, logS, logAD, logDs, nTilesWeGot*DPC);
                }
                else {
                    if (S==3)                    Nprescreened = PrescreenedListS <false, J,X, false  > (logD, logS, logAD, logDs, nTilesWeGot*DPC);
                    else                         Nprescreened = PrescreenedList  <false, J,X, false,S> (logD, logS, logAD, logDs, nTilesWeGot*DPC);
                }
            }

            //fills the GDO pescreening buckets with integrals
            InitBlocks<DPC>();

            if (!J && X) {
                if (S==3)                    PrescreenedListS   <true, J,X, false  > (logD, logS, logAD, logDs, Nprescreened);
                else                         ONXPrescreenedList  <true, S> (logD, logS, logAD, logDs, Nprescreened);
            }
            else {
                if (S==3)                    PrescreenedListS <true, J,X, false  > (logD, logS, logAD, logDs, Nprescreened);
                else                         PrescreenedList  <true, J,X, false,S> (logD, logS, logAD, logDs, Nprescreened);
            }


            //merges all incomplete tiles to maximize vectorization and sorts each GDO bucket quadruplet list by increasing distance (to later use long range approximations)
            FuseIncompleteTiles<DPC>();
            SortPrescreenedLists();

            newbatch->Ntiles64        = MakeAllTiles(newbatch->TileList64);

            if (newbatch->Ntiles64 != nTilesWeGot) {
                cout << "problem!!!    " << newbatch->Ntiles64 << " " << nTiles2Get << " " << nTilesWeGot << " " << GetNTiles<DPC>() << "   " << Nprescreened << endl;
            }

            TaskList.push(newbatch);
        }

        State.storeABb = 0;
        State.storeCDb = 0;
        State.storeAB  = 0;
        State.storeCD  = 0;

        //free the memory
        if (!J && X) CleanONX();

        if (S==1) {
            State.store3center = true;
            return false;
        }
        else if (geometry==ABCD) {
            State.storeBlock   = true;
            return false;
        }
        else return true;
    }

    //do float packed
    if (State.storeBlock) {

        SetInitialList(logS);
        InitLists();

        size_t size8 = (wsize4 + msize4) * sizeof(cacheline32) + FPC * sizeof(ERIgeometry) + sizeof(ERIgeometries32) + sizeof(RotationMatrices32); //maximum size (assuming a float tile)
        UI64 max_tiles = min ( UI64( (block-mem) / size8), MAX_TILES_PER_BLOCK);

        while (1) {

            UI16 sABb = State.storeABb;
            UI16 sCDb = State.storeCDb;
            UI32 sAB  = State.storeAB;
            UI32 sCD  = State.storeCD;

            UI64 Nprescreened;

            //counts the number of shell quartets in each GDO prescreening bucket
            ClearBlocks();

            if (S==3)                    Nprescreened = PrescreenedListS <false, J,X, true  > (logD, logS, logAD, logDs, max_tiles*FPC);
            else                         Nprescreened = PrescreenedList  <false, J,X, true,S> (logD, logS, logAD, logDs, max_tiles*FPC);


            if (Nprescreened==0) break;

            UI64 nTiles2Get = (Nprescreened + (FPC-1))/FPC; //at least; let's think about it l8r

            UI64 size = nTiles2Get * sizeof(ERITile32);
            ERITile32 * tList = (ERITile32*)EBB->Allocate(size);

            // no memory at all!!!
            if (tList==NULL) return false;

            UI64 nTilesWeGot = size/sizeof(ERITile32);

            ERIBatch * newbatch  = new ERIBatch;
            *newbatch                 = *this;
            newbatch->TileList32      = tList;
            newbatch->memAllocated    = size;
            newbatch->J               = J;
            newbatch->X               = X;



            // must compute again, with fewer tiles
            if (nTiles2Get!=nTilesWeGot) {
                ClearBlocks();

                if (S==3)                    Nprescreened = PrescreenedListS <false, J,X, true  > (logD, logS, logAD, logDs, nTilesWeGot*FPC);
                else                         Nprescreened = PrescreenedList  <false, J,X, true,S> (logD, logS, logAD, logDs, nTilesWeGot*FPC);
            }

            //fills the GDO pescreening buckets with integrals
            InitBlocks<FPC>();

            if (S==3)                    PrescreenedListS <true, J,X, true  > (logD, logS, logAD, logDs, Nprescreened);
            else                         PrescreenedList  <true, J,X, true,S> (logD, logS, logAD, logDs, Nprescreened);


            //merges all incomplete tiles to maximize vectorization and sorts each GDO bucket quadruplet list by increasing distance (to later use long range approximations)
            FuseIncompleteTiles<FPC>();
            SortPrescreenedLists();

            if (GetNTiles<FPC>() != nTilesWeGot) {
                cout << "problem!!!    " << nTiles2Get << " " << nTilesWeGot << " " << GetNTiles<FPC>() << "   " << Nprescreened << endl;;
            }

            newbatch->Ntiles32        = MakeAllTiles(newbatch->TileList32);

            TaskList.push(newbatch);
        }

        return true;
    }

    return true;
}



void PrintEvaluationInfo(GEOM geometry, const ShellPairPrototype * ABp, const ShellPairPrototype * CDp) {
    Echidna::EBenchmarker << "Thread number " << omp_get_thread_num() << " beginning ERI block evaluation" << endl;
    Echidna::EBenchmarker << "   Geometry ";

    if      (geometry==ABCD) Echidna::EBenchmarker << "ABCD";
    else if (geometry==ABAB) Echidna::EBenchmarker << "ABAB";
    else if (geometry==AACD) Echidna::EBenchmarker << "AACD";
    else if (geometry==AACC) Echidna::EBenchmarker << "AACC";

    Echidna::EBenchmarker << "   Ls: " << int(ABp->l1) << " " << int(ABp->l2) << " " << int(CDp->l1) << " " << int(CDp->l2) << endl;
    Echidna::EBenchmarker << "   Elements: " << ABp->atype1 << " " << ABp->atype2 << " " << CDp->atype1 << " "  << CDp->atype2 << endl;
    Echidna::EBenchmarker << "   BF: " << int(ABp->nb1) << " " << int(ABp->nb2) << " " << int(CDp->nb1) << " " << int(CDp->nb2) << endl;

}

bool ERIblock::Check3Center() {
    bool check3C = false;
    if (geometry==ABCD) {
        //check if it is possible to have overlapping centers in ABCD integrals
        if (ABp->atype1==CDp->atype1 || ABp->atype1==CDp->atype2 || ABp->atype2==CDp->atype1 || ABp->atype2==CDp->atype2)
            check3C = true;
    }

    return check3C;
}

bool ERIblock::MakeCoulomb (float logD, float logS,  const r2tensor< r2tensor<float> > & logAD, const r2tensor< r2tensor<float> > & logDs, TaskQueue<ERIBatch> & TaskList) {

    //PrintEvaluationInfo(geometry, ABp, CDp);
    //Echidna::EBenchmarker << " Coulomb only" << endl;

    if      (geometry==ABAB)     return MakeBlocks <true, false, 3> (logD, logS, logAD, logDs, TaskList);
    else if (Check3Center()) {
        if (!State.store3center) return MakeBlocks <true, false, 1> (logD, logS, logAD, logDs, TaskList);
        else                     return MakeBlocks <true, false, 2> (logD, logS, logAD, logDs, TaskList);
    }
    else                         return MakeBlocks <true, false, 0> (logD, logS, logAD, logDs, TaskList);
}

bool ERIblock::MakeExchange(float logD, float logS,  const r2tensor< r2tensor<float> > & logAD, const r2tensor< r2tensor<float> > & logDs, TaskQueue<ERIBatch> & TaskList) {

    //PrintEvaluationInfo(geometry, ABp, CDp);
    //Echidna::EBenchmarker << " Exchange only" << endl;

    if (geometry==ABAB)      return MakeBlocks <false, true, 3> (logD, logS, logAD, logDs, TaskList);
    else if (Check3Center()) {
        if (!State.store3center) return MakeBlocks <false, true, 1> (logD, logS, logAD, logDs, TaskList);
        else                     return MakeBlocks <false, true, 2> (logD, logS, logAD, logDs, TaskList);
    }
    else                     return MakeBlocks <false, true, 0> (logD, logS, logAD, logDs, TaskList);
}

bool ERIblock::MakeAll     (float logD, float logS,  const r2tensor< r2tensor<float> > & logAD, const r2tensor< r2tensor<float> > & logDs, TaskQueue<ERIBatch> & TaskList) {

    //PrintEvaluationInfo(geometry, ABp, CDp);
    //Echidna::EBenchmarker << " Coulomb and Exchange simultaneously" << endl;

    if (geometry==ABAB)      return MakeBlocks <true, true, 3> (logD, logS, logAD, logDs, TaskList);
    else if (Check3Center()) {
        if (!State.store3center) return MakeBlocks <true, true, 1> (logD, logS, logAD, logDs, TaskList);
        else                     return MakeBlocks <true, true, 2> (logD, logS, logAD, logDs, TaskList);
    }
    else                     return MakeBlocks <true, true, 0> (logD, logS, logAD, logDs, TaskList);
}


