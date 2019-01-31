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
#include "../math/angular.hpp"
#include "../math/tensors.hpp"
#include "../math/eigen.hpp"

#include "../integrals/fock2e.hpp"
#include "../integrals/rotations.hpp"
#include "../integrals/atomprod.hpp"

#include "../basis/shellpair.hpp"
#include "../basis/SPprototype.hpp"

#include "../libechidna/libechidna.hpp"


#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
    #define omp_get_num_threads() 1
#endif

using namespace std;
using namespace LibAngular;
using namespace LibQuimera;


size_t Fock2e::MEMPERBATCH;


static inline UI8 lmax(UI8 l1, UI8 l2) {
    l1 = (l1==LSP)?1:l1;
    l2 = (l2==LSP)?1:l2;
    return max(l1,l2);
}

static inline int min(int a, int b) {
    return (a<b)?a:b;
}

static inline int Tpos (int M1, int M2, int M3, int M4, int v1, int v2, int v3, int v4) {
    return v4 * M3*M2*M1 + v3*M2*M1 + v2*M1 + v1;
}


//very inefficient
//the point is that the matrix is constructed only once for each type of atom
void ConstructBigW(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, int wa, const double * w, double * W) {

    const int maxLa = nmS[ABp.l1];
    const int maxLb = nmS[ABp.l2];
    const int maxLc = nmS[CDp.l1];
    const int maxLd = nmS[CDp.l2];

    int wsize = maxLa * maxLb * maxLc * maxLd;

    const int Ja = ABp.Ja;
    const int Jb = ABp.Jb;
    const int Jc = CDp.Ja;
    const int Jd = CDp.Jb;

    const double * TT = w;

    int fa = ABp.f1;
    for (int ja=0; ja<Ja; ++ja) {
        int fb = ABp.f2;
        for (int jb=0; jb<Jb; ++jb) {
            int fc = CDp.f1;
            for (int jc=0; jc<Jc; ++jc) {
                int fd = CDp.f2;
                for (int jd=0; jd<Jd; ++jd) {

                    for (int ma=0; ma<maxLa; ++ma) {
                        for (int mb=0; mb<maxLb; ++mb) {
                            for (int mc=0; mc<maxLc; ++mc) {
                                for (int md=0; md<maxLd; ++md) {
                                    int posw = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);

                                    int faa, fbb, fcc, fdd;

                                    if (ABp.inverted) {faa = fa+ma; fbb = fb+mb;}
                                    else              {faa = fb+mb; fbb = fa+ma;}

                                    if (CDp.inverted) {fcc = fd+md; fdd = fc+mc;}
                                    else              {fcc = fc+mc; fdd = fd+md;}

                                    int posW1 = Tpos (wa, wa, wa, wa, faa, fbb, fcc, fdd);
                                    int posW2 = Tpos (wa, wa, wa, wa, fcc, fdd, faa, fbb);

                                    int posW3 = Tpos (wa, wa, wa, wa, fbb, faa, fcc, fdd);
                                    int posW4 = Tpos (wa, wa, wa, wa, fcc, fdd, fbb, faa);

                                    int posW5 = Tpos (wa, wa, wa, wa, faa, fbb, fdd, fcc);
                                    int posW6 = Tpos (wa, wa, wa, wa, fdd, fcc, faa, fbb);

                                    int posW7 = Tpos (wa, wa, wa, wa, fbb, faa, fdd, fcc);
                                    int posW8 = Tpos (wa, wa, wa, wa, fdd, fcc, fbb, faa);


                                    W[posW1] = TT[posw];
                                    W[posW2] = TT[posw];
                                    W[posW3] = TT[posw];
                                    W[posW4] = TT[posw];
                                    W[posW5] = TT[posw];
                                    W[posW6] = TT[posw];
                                    W[posW7] = TT[posw];
                                    W[posW8] = TT[posw];

                                }
                            }
                        }
                    }

                    TT += wsize;

                    fd += maxLd;
                }
                fc += maxLc;
            }
            fb += maxLb;
        }
        fa += maxLa;
    }
}

void Fock2e::AddCoulomb(Sparse & F, int at, int wa, const double * const T, const Sparse & D) {

    UI16 id = D.ids[at];

    //loop over shells
    for (int b1=0, ma=0; b1<D.nfs[id]; ++b1) {
        UI32 mj1 = D.js[id][b1];
        UI32 mm1 = D.ms[id][b1];

        for (UI8 j1=0; j1<mj1; ++j1) {
            for (UI8 m1=0; m1<mm1; ++m1, ++ma) {


                for (int b2=0, mb=0; b2<D.nfs[id]; ++b2) {
                    UI32 mj2 = D.js[id][b2];
                    UI32 mm2 = D.ms[id][b2];

                    for (UI8 j2=0; j2<mj2; ++j2) {
                        for (UI8 m2=0; m2<mm2; ++m2, ++mb) {
                            if (mb>ma) continue;

                            double sum = 0;
                            int pos = (ma*wa+mb)*wa*wa;


                            for (int b3=0, mc=0; b3<D.nfs[id]; ++b3) {
                                UI32 mj3 = D.js[id][b3];
                                UI32 mm3 = D.ms[id][b3];

                                for (UI8 j3=0; j3<mj3; ++j3) {
                                    for (UI8 m3=0; m3<mm3; ++m3, ++mc) {


                                        for (int b4=0, md=0; b4<D.nfs[id]; ++b4) {
                                            UI32 mj4 = D.js[id][b4];
                                            UI32 mm4 = D.ms[id][b4];

                                            for (UI8 j4=0; j4<mj4; ++j4) {
                                                for (UI8 m4=0; m4<mm4; ++m4, ++md) {

                                                    sum += T[pos] * D(at,at, b3,b4, j3,j4, m3,m4);

                                                    ++pos;
                                                }
                                            }
                                        }

                                    }
                                }
                            }

                            if (ma==mb) sum*= 0.5;

                            F(at,at, b1,b2, j1,j2, m1,m2) += 2*sum;
                        }
                    }
                }

            }
        }
    }

}

void Fock2e::AddXchange(Sparse & F, int at, int wa, const double * const T, const Sparse & D) {


    UI16 id = D.ids[at];

    //loop over shells
    int ma = 0;
    for (int b1=0; b1<D.nfs[id]; ++b1) {
        UI32 mj1 = D.js[id][b1];
        UI32 mm1 = D.ms[id][b1];

        for (UI8 j1=0; j1<mj1; ++j1) {
            for (UI8 m1=0; m1<mm1; ++m1, ++ma) {

                int mc = 0;
                for (int b3=0; b3<D.nfs[id]; ++b3) {
                    UI32 mj3 = D.js[id][b3];
                    UI32 mm3 = D.ms[id][b3];

                    for (UI8 j3=0; j3<mj3; ++j3) {
                        for (UI8 m3=0; m3<mm3; ++m3, ++mc) {

                            if (mc>ma) continue;


                            double sum = 0;
                            int pos = (ma*wa*wa + mc)*wa;

                            int mb = 0;

                            for (int b2=0; b2<D.nfs[id]; ++b2) {
                                UI32 mj2 = D.js[id][b2];
                                UI32 mm2 = D.ms[id][b2];

                                for (UI8 j2=0; j2<mj2; ++j2) {
                                    for (UI8 m2=0; m2<mm2; ++m2, ++mb) {

                                        int pos2 = pos + mb*wa*wa;

                                        int md = 0;

                                        for (int b4=0; b4<D.nfs[id]; ++b4) {
                                            UI32 mj4 = D.js[id][b4];
                                            UI32 mm4 = D.ms[id][b4];

                                            for (UI8 j4=0; j4<mj4; ++j4) {
                                                for (UI8 m4=0; m4<mm4; ++m4, ++md) {

                                                    sum += T[pos2] * D(at,at, b2,b4, j2,j4, m2,m4);

                                                    ++pos2;
                                                }
                                            }
                                        }

                                    }
                                }
                            }

                            if (ma==mc) sum*= 0.5;

                            F(at,at, b1,b3, j1,j3, m1,m3) -= sum;
                        }
                    }
                }

            }
        }
    }

}

void Fock2e::Trim(Sparse & F) {

    for (int at=0; at<F.natoms; ++at) {
        UI16 id = F.ids[at];

        //loop over shells
        int ma = 0;
        for (int b1=0; b1<F.nfs[id]; ++b1) {
            UI32 mj1 = F.js[id][b1];
            UI32 mm1 = F.ms[id][b1];

            for (UI8 j1=0; j1<mj1; ++j1) {
                for (UI8 m1=0; m1<mm1; ++m1, ++ma) {

                    int mc = 0;
                    for (int b3=0; b3<F.nfs[id]; ++b3) {
                        UI32 mj3 = F.js[id][b3];
                        UI32 mm3 = F.ms[id][b3];

                        for (UI8 j3=0; j3<mj3; ++j3) {
                            for (UI8 m3=0; m3<mm3; ++m3, ++mc) {

                                if (mc>ma) F(at,at, b1,b3, j1,j3, m1,m3) = 0;

                                //F(at,at, b1,b3, j1,j3, m1,m3) -= sum;
                            }
                        }
                    }

                }
            }
        }

    }

}

static inline bool InvertPrototypePair(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp) {
    UI8 ma = ABp.l1==LSP?4:2*ABp.l1+1;
    UI8 mb = ABp.l2==LSP?4:2*ABp.l2+1;
    UI8 mc = CDp.l1==LSP?4:2*CDp.l1+1;
    UI8 md = CDp.l2==LSP?4:2*CDp.l2+1;

    return (ma>mc) || (ma==mc && mb>=md);
}



Fock2e::Fock2e() {
    Blocks    = NULL;
    BatchEval = NULL;

    //APP = NULL;
    APL = NULL;
}

Fock2e::~Fock2e() {
    delete[] Blocks;
    delete[] BatchEval;
}


void Fock2e::Init(UI32 nInteractingAtoms, const APlists & pAPL, const sparsetensorpattern & SparseTensorPattern) {

    XS = 1.;

    //these pointers hold the information about the shellpairs, etc.
    APL = &pAPL;

    //all the following sizes for buffers and arrays should be calculated properly
    OGMT = omp_get_max_threads();

    //MAXMEM       = 64  * Gword; //1*Gword;
    //MEMBUFFER    = 256 * Mword; //8*MAX_TILES_PER_BLOCK*sizeof(ERITile64) ; //2*Gword;
    MEMPERBATCH  = Echidna::MAXMEM / OGMT;


    //maximum number of atom interactions
    max_at_int = nInteractingAtoms; //AtomInteractions.GetTotN();

    EbatchBuffer = new ERIbatchBuffer;
    EbatchBuffer->Init(Echidna::MEMBUFFER);

    Blocks = new ERIblock[OGMT];
    for (int p=0; p<OGMT; ++p)
        Blocks[p].Init(max_at_int, EbatchBuffer);

    BatchEval = new BatchEvaluator[OGMT];
    for (int p=0; p<OGMT; ++p)
        BatchEval[p].Init(MEMPERBATCH);

    for (int p=0; p<OGMT; ++p) {
        //BatchEval[p].JX.Set();
        BatchEval[p].Jcont.Set();
        BatchEval[p].Xcont.Set();
    }

    //initialize density prescreen matrix
    UI32 natoms = SparseTensorPattern.natoms;

    r1tensor<int> DPlengths(natoms);

    for (UI32 at1=0; at1<natoms; ++at1) {
        UI16 id1 = SparseTensorPattern.ids[at1];
        UI16 nb1 = SparseTensorPattern.nfs[id1];

        DPlengths[at1] = nb1;
    }

    logDs.set(DPlengths, DPlengths);
    logAD.set(DPlengths, DPlengths);


    InitERIroutines(); //initializes the ERI evaluation routines needed for the system
    InitJobQueue();    //initialized the job queue
}


void Fock2e::SetXS(double Xscaling) {
    XS = Xscaling;
}

void Fock2e::InitERIroutines() {

    //sets the geometry and angular momenta only for the required routines
    //********************************************************************

    //all four atoms different (4-center ERIs) (k² N²)
    for (int n12=0; n12<APL->AtomPairBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomPairBatch[n12];
        if (APB12.APlist.n==0) continue;

        for (int n34=0; n34<=n12; ++n34) {
            const APbatch & APB34 = APL->AtomPairBatch[n34];
            if (APB34.APlist.n==0) continue;

            //loop over shell pair prototypes
            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];
                for (int cd=0; cd<APB34.APP->nGPt; ++cd) {
                    const ShellPairPrototype & CDp = APB34.APP->Gs[cd];
                    if (APB12.nAPS[ab]==0) continue;
                    if (APB34.nAPS[cd]==0) continue;
                    bool rev = InvertPrototypePair(ABp, CDp);
                    if (rev) Q.ListNeeded(ABCD, ABp.l1, ABp.l2, CDp.l1, CDp.l2);
                    else     Q.ListNeeded(ABCD, CDp.l1, CDp.l2, ABp.l1, ABp.l2);
                }
            }

            // for integrals with at least one atom center shared

            // needed because CDR becomes numerically unstable when trying to apply
            // the general, 4-center integral algorithm to degenerate 3-center and 2-center
            // cases in conjunction with basis sets containing very diffuse primitives
            if ( (APB12.APP->Atype1 == APB34.APP->Atype1) ||
                 (APB12.APP->Atype1 == APB34.APP->Atype2) ||
                 (APB12.APP->Atype2 == APB34.APP->Atype1) ||
                 (APB12.APP->Atype2 == APB34.APP->Atype2) ) {

                //loop over shell pair prototypes
                for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                    const ShellPairPrototype & ABp = APB12.APP->Gs[ab];
                    for (int cd=0; cd<APB34.APP->nGPt; ++cd) {
                        const ShellPairPrototype & CDp = APB34.APP->Gs[cd];
                        if (APB12.nAPS[ab]==0) continue;
                        if (APB34.nAPS[cd]==0) continue;
                        bool rev = InvertPrototypePair(ABp, CDp);

                        if (rev) Q.ListNeeded(ABAD, ABp.l1, ABp.l2, CDp.l1, CDp.l2);
                        else     Q.ListNeeded(ABAD, CDp.l1, CDp.l2, ABp.l1, ABp.l2);
                    }
                }

            }


        }
    }


    //same atom pair (2-center ABAB-type ERIs) (k N)
    for (int n12=0; n12<APL->AtomPairBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomPairBatch[n12];

        if (APB12.APlist.n==0) continue;

        for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
            const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

            for (int cd=0; cd<=ab; ++cd) {
                const ShellPairPrototype & CDp = APB12.APP->Gs[cd];

                if (APB12.nAPS[ab]==0) continue;
                if (APB12.nAPS[cd]==0) continue;

                bool rev = InvertPrototypePair(ABp, CDp);

                if (rev) Q.ListNeeded(ABAB, ABp.l1, ABp.l2, CDp.l1, CDp.l2);
                else     Q.ListNeeded(ABAB, CDp.l1, CDp.l2, ABp.l1, ABp.l2);
            }
        }
    }


    //two centers + same atom (3-center ERIs) (k N²)
    for (int n12=0; n12<APL->AtomSameBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomSameBatch[n12];
        if (APB12.APlist.n==0) continue;

        for (int n34=0; n34<APL->AtomPairBatch.n; ++n34) {
            const APbatch & APB34 = APL->AtomPairBatch[n34];
            if (APB34.APlist.n==0) continue;

            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];
                for (int cd=0; cd<APB34.APP->nGPt; ++cd) {
                    const ShellPairPrototype & CDp = APB34.APP->Gs[cd];
                    if (APB12.nAPS[ab]==0) continue;
                    if (APB34.nAPS[cd]==0) continue;
                    Q.ListNeeded(AACD, ABp.l1, ABp.l2, CDp.l1, CDp.l2);
                }
            }
        }
    }

    //two different atoms (2-center ERIs) (N²)
    for (int n12=0; n12<APL->AtomSameBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomSameBatch[n12];
        if (APB12.APlist.n==0) continue;

        for (int n34=0; n34<=n12; ++n34) {
            const APbatch & APB34 = APL->AtomSameBatch[n34];
            if (APB34.APlist.n==0) continue;
            if (n12==n34 && APB34.APlist.n==1) continue; //only one atom!

            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                for (int cd=0; cd<APB34.APP->nGPt; ++cd) {
                    const ShellPairPrototype & CDp = APB34.APP->Gs[cd];
                    bool rev = InvertPrototypePair(ABp, CDp);
                    if (rev) Q.ListNeeded(AACC, ABp.l1, ABp.l2, CDp.l1, CDp.l2);
                    else     Q.ListNeeded(AACC, CDp.l1, CDp.l2, ABp.l1, ABp.l2);
                }
            }
        }
    }

    //one atom only (1-center)
    for (int n12=0; n12<APL->AtomSameBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomSameBatch[n12];
        const AtomProdPrototype * APP = APB12.APP;

        //calcula las integrales para SOLO un atomo de cada tipo
        for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
            const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

            for (int cd=0; cd<=ab; ++cd) {
                const ShellPairPrototype & CDp = APB12.APP->Gs[cd];

                if ((ABp.l1+ABp.l2+CDp.l1+CDp.l2)%2) continue;
                if ((ABp.l2==0) && (CDp.l2==0) && (ABp.l1!=CDp.l1)) continue;

                bool rev = InvertPrototypePair(ABp,CDp);

                if (rev) Q.ListNeeded(AAAA, ABp.l1, ABp.l2, CDp.l1, CDp.l2);
                else     Q.ListNeeded(AAAA, CDp.l1, CDp.l2, ABp.l1, ABp.l2);
            }
        }
    }

    //initializes (computes) the required routines
    Q.GenerateNeeded(false);
}

//put in the queue all the calls to different ERI routines, except for the 1-nuclear
void Fock2e::InitJobQueue() {

    TwoElectronJobs.Clear();

    // ERI calc
    // ********

    //4-center integrals
    // ****************

    //all four atoms different (4-center ERIs) (k² N²)
    for (int n12=0; n12<APL->AtomPairBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomPairBatch[n12];

        if (APB12.APlist.n==0) continue;

        for (int n34=0; n34<n12; ++n34) {
            const APbatch & APB34 = APL->AtomPairBatch[n34];

            if (APB34.APlist.n==0) continue;

            //loop over shell pair prototypes
            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                for (int cd=0; cd<APB34.APP->nGPt; ++cd) {
                    const ShellPairPrototype & CDp = APB34.APP->Gs[cd];

                    if (APB12.nAPS[ab]==0) continue;
                    if (APB34.nAPS[cd]==0) continue;

                    bool rev = InvertPrototypePair(ABp, CDp);

                    BatchInfo Info;

                    //sets the general properties of the ERI block
                    if (rev) {
                        Info.Set(ABCD, ABp, CDp, false);
                        Info.SetLists(APB12.SP[ab],APB34.SP[cd],&APB12.APlist,&APB34.APlist, APB12.nAPS[ab], APB34.nAPS[cd],false);
                    }
                    else {
                        Info.Set(ABCD, CDp, ABp, false);
                        Info.SetLists(APB34.SP[cd],APB12.SP[ab],&APB34.APlist,&APB12.APlist, APB34.nAPS[cd], APB12.nAPS[ab],false);
                    }

                    TwoElectronJobs.push(Info);
                }
            }
        }
    }

    //n12==n34
    for (int n12=0; n12<APL->AtomPairBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomPairBatch[n12];

        if (APB12.APlist.n==0) continue;

        {
            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                for (int cd=0; cd<ab; ++cd) {
                    const ShellPairPrototype & CDp = APB12.APP->Gs[cd];

                    if (APB12.nAPS[ab]==0) continue;
                    if (APB12.nAPS[cd]==0) continue;

                    bool rev = InvertPrototypePair(ABp, CDp);

                    BatchInfo Info;

                    if (rev) {
                        Info.Set(ABCD, ABp, CDp, false);
                        Info.SetLists(APB12.SP[ab],APB12.SP[cd],&APB12.APlist,&APB12.APlist, APB12.nAPS[ab], APB12.nAPS[cd],true);
                    }
                    else {
                        Info.Set(ABCD, CDp, ABp, false);
                        Info.SetLists(APB12.SP[cd],APB12.SP[ab],&APB12.APlist,&APB12.APlist, APB12.nAPS[cd], APB12.nAPS[ab],true);
                    }

                    TwoElectronJobs.push(Info);
                }

                //ab==cd
                {
                    BatchInfo Info;

                    Info.Set(ABCD, ABp, ABp, false);
                    Info.SetLists(APB12.SP[ab],&APB12.APlist, APB12.nAPS[ab],false);

                    TwoElectronJobs.push(Info);
                }
            }
        }
    }

    //same atom pair (2-center ERIs, spread integrals) (k N)
    for (int n12=0; n12<APL->AtomPairBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomPairBatch[n12];

        if (APB12.APlist.n==0) continue;


        for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
            const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

            for (int cd=0; cd<ab; ++cd) {
                const ShellPairPrototype & CDp = APB12.APP->Gs[cd];

                if (APB12.nAPS[ab]==0) continue;
                if (APB12.nAPS[cd]==0) continue;

                bool rev = InvertPrototypePair(ABp, CDp);
                UI32 tmax = min(APB12.nAPS[ab], APB12.nAPS[cd]);

                BatchInfo Info;

                if (rev) {
                    Info.Set(ABAB, ABp, CDp, false);
                    Info.SetLists(APB12.SP[ab],APB12.SP[cd],&APB12.APlist,&APB12.APlist,tmax,tmax,false);
                }
                else {
                    Info.Set(ABAB, CDp, ABp, false);
                    Info.SetLists(APB12.SP[cd],APB12.SP[ab],&APB12.APlist,&APB12.APlist,tmax,tmax,false);
                }

                TwoElectronJobs.push(Info);
            }

            //ab==cd
            {
                UI32 tmax = APB12.nAPS[ab];

                BatchInfo Info;

                Info.Set(ABAB, ABp, ABp, true);
                Info.SetLists(APB12.SP[ab],&APB12.APlist,tmax,false);

                TwoElectronJobs.push(Info);
            }
        }
    }


    //3-center integrals
    // ****************

    //two centers + same atom (3-center ERIs) (k N²)
    for (int n12=0; n12<APL->AtomSameBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomSameBatch[n12];

        if (APB12.APlist.n==0) continue;

        for (int n34=0; n34<APL->AtomPairBatch.n; ++n34) {
            const APbatch & APB34 = APL->AtomPairBatch[n34];

            if (APB34.APlist.n==0) continue;

            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                for (int cd=0; cd<APB34.APP->nGPt; ++cd) {
                    const ShellPairPrototype & CDp = APB34.APP->Gs[cd];

                    if (APB34.nAPS[cd]==0) continue;

                    BatchInfo Info;

                    Info.Set(AACD, ABp, CDp, false);
                    Info.SetLists(APB12.SP[ab],APB34.SP[cd],&APB12.APlist,&APB34.APlist, APB12.nAPS[ab], APB34.nAPS[cd],false);

                    TwoElectronJobs.push(Info);
                }
            }
        }
    }

    //2-center integrals
    // ****************

    //two different atoms (2-center ERIs) (N²)
    for (int n12=0; n12<APL->AtomSameBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomSameBatch[n12];

        if (APB12.APlist.n==0) continue;

        for (int n34=0; n34<n12; ++n34) {
            const APbatch & APB34 = APL->AtomSameBatch[n34];

            if (APB34.APlist.n==0) continue;

            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                for (int cd=0; cd<APB34.APP->nGPt; ++cd) {
                    const ShellPairPrototype & CDp = APB34.APP->Gs[cd];

                    bool rev = InvertPrototypePair(ABp, CDp);

                    BatchInfo Info;

                    if (rev) {
                        Info.Set(AACC, ABp, CDp, false);
                        Info.SetLists(APB12.SP[ab],APB34.SP[cd],&APB12.APlist,&APB34.APlist, APB12.APlist.n, APB34.APlist.n,false);
                    }
                    else {
                        Info.Set(AACC, CDp, ABp, false);
                        Info.SetLists(APB34.SP[cd],APB12.SP[ab],&APB34.APlist,&APB12.APlist, APB34.APlist.n, APB12.APlist.n,false);
                    }

                    TwoElectronJobs.push(Info);
                }
            }
        }
    }

    //n12==n34
    for (int n12=0; n12<APL->AtomSameBatch.n; ++n12) {
        const APbatch & APB12 = APL->AtomSameBatch[n12];

        if (APB12.APlist.n<2) continue;

        {
            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                for (int cd=0; cd<ab; ++cd) {
                    const ShellPairPrototype & CDp = APB12.APP->Gs[cd];

                    bool rev = InvertPrototypePair(ABp, CDp);

                    BatchInfo Info;

                    if (rev) {
                        Info.Set(AACC, ABp, CDp, false);
                        Info.SetLists(APB12.SP[ab],APB12.SP[cd],&APB12.APlist,&APB12.APlist, APB12.APlist.n, APB12.APlist.n,true);
                    }
                    else {
                        Info.Set(AACC, CDp, ABp, false);
                        Info.SetLists(APB12.SP[cd],APB12.SP[ab],&APB12.APlist,&APB12.APlist, APB12.APlist.n, APB12.APlist.n,true);
                    }

                    TwoElectronJobs.push(Info);
                }

                {
                    BatchInfo Info;

                    Info.Set(AACC, ABp,ABp, false);
                    Info.SetLists(APB12.SP[ab],&APB12.APlist, APB12.APlist.n,false);

                    TwoElectronJobs.push(Info);
                }
            }
        }
    }

    TwoElectronJobs.Done();
}

//compute the frobenius norm of every block and extract the logarithm
void Fock2e::CalcDensityPrescreeners() {

    UI32 natoms = D2sparse.natoms;

    for (UI32 at1=0; at1<natoms; ++at1) {
        for (UI32 at2=0; at2<natoms; ++at2) {
            UI16 id1 = D2sparse.ids[at1];
            UI16 id2 = D2sparse.ids[at2];

            for (UI16 bf1=0; bf1<D2sparse.nfs[id1]; ++bf1) {
                for (UI16 bf2=0; bf2<D2sparse.nfs[id2]; ++bf2) {
                    UI32 len = D2sparse.flen(id1,bf1) * D2sparse.flen(id2,bf2);

                    const double * vals = D2sparse(at1, at2, bf1, bf2);

                    double sum2 = 0;
                    for (UI32 it=0; it<len; ++it) {
                        sum2 += vals[it]*vals[it];
                    }

                    sum2 = max(sum2, 1.e-100);
                    sum2 = -log(sum2)/2;

                    logDs[at1][at2][bf1][bf2] = sum2;
                }
            }
        }
    }

    for (UI32 at1=0; at1<natoms; ++at1) {
        for (UI32 at2=0; at2<natoms; ++at2) {
            UI16 id1 = AD2sparse.ids[at1];
            UI16 id2 = AD2sparse.ids[at2];

            for (UI16 bf1=0; bf1<AD2sparse.nfs[id1]; ++bf1) {
                for (UI16 bf2=0; bf2<AD2sparse.nfs[id2]; ++bf2) {
                    UI32 len = AD2sparse.flen(id1,bf1) * AD2sparse.flen(id2,bf2);

                    const double * vals = AD2sparse(at1, at2, bf1, bf2);

                    double sum2 = 0;
                    for (UI32 it=0; it<len; ++it) {
                        sum2 += vals[it]*vals[it];
                    }

                    sum2 = max(sum2, 1.e-100);
                    sum2 = -log(sum2)/2;

                    logAD[at1][at2][bf1][bf2] = sum2;
                }
            }
        }
    }

    /*
    for (UI32 at1=0; at1<natoms; ++at1) {
        UI16 id1 = D2sparse.ids[at1];

        for (UI16 bf1=0; bf1<D2sparse.nfs[id1]; ++bf1) {

            for (UI32 at2=0; at2<natoms; ++at2) {
                UI16 id2 = D2sparse.ids[at2];

                for (UI16 bf2=0; bf2<D2sparse.nfs[id2]; ++bf2) {
                    cout << logDs[at1][at2][bf1][bf2] << " ";
                }
                cout << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    cout << endl;
    cout << endl;

    for (UI32 at1=0; at1<natoms; ++at1) {
        UI16 id1 = D2sparse.ids[at1];

        for (UI16 bf1=0; bf1<D2sparse.nfs[id1]; ++bf1) {

            for (UI32 at2=0; at2<natoms; ++at2) {
                UI16 id2 = D2sparse.ids[at2];

                for (UI16 bf2=0; bf2<D2sparse.nfs[id2]; ++bf2) {
                    cout << logAD[at1][at2][bf1][bf2] << " ";
                }
                cout << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    char a;
    cin >> a;
    */
}

void Fock2e::CalcSparsicity(Sparse & SS) {

    UI32 natoms = SS.natoms;

    for (UI32 at1=0; at1<natoms; ++at1) {
        for (UI32 at2=0; at2<natoms; ++at2) {
            UI16 id1 = SS.ids[at1];
            UI16 id2 = SS.ids[at2];

            for (UI16 bf1=0; bf1<SS.nfs[id1]; ++bf1) {
                for (UI16 bf2=0; bf2<SS.nfs[id2]; ++bf2) {
                    UI32 len = SS.flen(id1,bf1) * SS.flen(id2,bf2);

                    const double * vals = SS(at1, at2, bf1, bf2);

                    bool filled = false;

                    for (UI32 it=0; it<len; ++it)
                        filled |= (vals[it]!=0);

                    logDs[at1][at2][bf1][bf2] = (filled?1.:0.);
                }
            }
        }
    }

}

//evaluates Cauchy-Swarz bounds for ERI
void Fock2e::CalcCauchySwarz() {

    Echidna::EMessenger << "Calculating Cauchy-Swarz bounds for integral screening";
    Echidna::EMessenger.Push(); {

        //two different atoms
        for (int n12=0; n12<APL->AtomPairBatch.n; ++n12) {
            const APbatch & APB12 = APL->AtomPairBatch[n12];

            if (APB12.APlist.n==0) continue; //no overlaps between these two specific elements

            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                if (APB12.nAPS[ab]==0) continue; //no overlaps between these two shell pair prototypes

                int OGMT = omp_get_max_threads();

                ERIBatch CSbatch;
                CSbatch.Set(ABAB, ABp, ABp, true);
                CSbatch.Ntiles64 = 1;
                for (int i=0; i<OGMT; ++i) BatchEval[i].SetBuffers(CSbatch);

                //buffers, sizes and other stuff needed for the calculation

                //ERI
                const AtomProd & AP12 = APB12.APlist[0];
                ShellPair * ABs = AP12.pSPs[ab];

                #pragma omp parallel for schedule(dynamic)
                for (int ap12=0; ap12<APB12.nAPS[ab]; ap12+=DOUBLES_PER_CACHE_LINE) {

                    ERIgeometry   eriG[DOUBLES_PER_CACHE_LINE];
                    ERIgeometries64  eriG8;

                    ERITile64 ET;

                    ET.nKab = ABp.Ka*ABp.Kb -1;
                    ET.nKcd = ABp.Ka*ABp.Kb -1;

                    int lenk = min (DOUBLES_PER_CACHE_LINE, APB12.nAPS[ab]-ap12);

                    for (int k=0; k<lenk; ++k) {
                        ET.ap12[k] = ap12+k;
                        ET.ap34[k] = ap12+k;
                    }
                    for (int k=lenk; k<DOUBLES_PER_CACHE_LINE; ++k) {
                        ET.ap12[k] = ap12; //+lenk-1;
                        ET.ap34[k] = ap12; //+lenk-1;
                    }


                    ET.wSP12 = ABs->WW;
                    ET.wSP34 = ABs->WW;
                    ET.nK2ab = ABp.Ka * ABp.Kb;
                    ET.nK2cd = ABp.Ka * ABp.Kb;

                    ET.lowgeom = ABAB;

                    for (int k=0; k<DOUBLES_PER_CACHE_LINE; ++k) {
                        const AtomProd & AP12 = APB12.APlist[ ET.ap12[k] ];
                        eriG[k].MakeRotation2S (AP12);
                        eriG[k].Adjust         (ABp.inverted, ABp.inverted);
                    }

                    PackGeometries           (eriG, eriG8);

                    cacheline64 maxv;
                    maxv = BatchEval[omp_get_thread_num()].CalcCauchySchwarz(CSbatch, *ABs, eriG8, ET, lenk);

                    //set the Cauchy-Schwarz parameter
                    for (int k=0; k<lenk; ++k) {
                        double logCS = -log(maxv(k))/2;
                        ABs[ap12+k].logCS = logCS;
                        //cout << logCS << " ";
                    }
                }

                //cout << endl;
            }
        }

        //same atom
        for (int n12=0; n12<APL->AtomSameBatch.n; ++n12) {
            const APbatch & APB12 = APL->AtomSameBatch[n12];

            ERIgeometry eriG;
            point A(0,0,0);
            eriG.MakeRotation1(A);

            //#pragma omp parallel for schedule(dynamic)
            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];
                const ShellPair & ABs = APB12.SP[ab][0];

                ERIBatch CSbatch;
                CSbatch.Ntiles64 = 1;
                CSbatch.Set(AAAA, ABp, ABp, true);
                BatchEval[0].SetBuffers(CSbatch);

                double max = BatchEval[0].CalcCauchySchwarz(CSbatch, ABs,eriG);
                double logCS = -log(max)/2;

                //set the cauchy-schwarz parameter to ALL atom pairs
                //**************************************************

                for (int ap12=0; ap12<APB12.APlist.n; ++ap12)  {
                    const AtomProd & AP12 = APB12.APlist[ap12];
                    ShellPair & ABs = *AP12.pSPs[ab];
                    ABs.logCS = logCS;
                }
            }
        }

    } Echidna::EMessenger.Pop();
}


void Fock2e::FockUpdate(tensor2 & F, const tensor2 & D2, double Dthresh, double Sthresh, bool SplitJX) {

    //make a copy of the D2 tensor to a format which exploits data locality
    D2sparse  = D2;
    AD2sparse = F;

    CalcDensityPrescreeners();
    float logD = -log(Dthresh);
    float logS = -log(Sthresh);

    //set the partial JX matrices to 0
    for (int p=0; p<OGMT; ++p) BatchEval[p].Jcont.zeroize();
    if (XS>0.)
        for (int p=0; p<OGMT; ++p) BatchEval[p].Xcont.zeroize();

    for (int p=0; p<OGMT; ++p) Blocks[p].Reset();


    Echidna::EMessenger << "Computing two-electron contributions to the new Fock matrix";
    Echidna::EMessenger.Push();

    Chronometer ERIchrono;
    ERIchrono.Start();


    if (ThisNode.IsMaster()) {

        //1-center integrals
        // ****************
        //this solution is a bit unfortunate, but the proper solution involves
        //coding the many permutational symmetries and special cases

        for (int n12=0; n12<APL->AtomSameBatch.n; ++n12) {
            const APbatch & APB12 = APL->AtomSameBatch[n12];
            const AtomProdPrototype * APP = APB12.APP;

            UI64 wa = APB12.APP->wa1;
            UI64 wa4;
            wa4 = wa * wa * wa * wa;

            double * WWW = new double[wa4];

            ERIgeometry eriG;
            point A(0,0,0);
            eriG.MakeRotation1(A);


            //COULOMB
            //*******

            for (int i=0; i<wa4; ++i) WWW[i]=0;

            //calcula las integrales para SOLO un atomo de cada tipo
            #pragma omp parallel for schedule(dynamic)
            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPair & ABs = APB12.SP[ab][0];
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                for (int cd=0; cd<=ab; ++cd) {
                    const ShellPair & CDs = APB12.SP[cd][0];
                    const ShellPairPrototype & CDp = APB12.APP->Gs[cd];

                    if ((ABp.l1+ABp.l2+CDp.l1+CDp.l2)%2) continue;
                    if ((ABp.l2==0) && (CDp.l2==0) && (ABp.l1!=CDp.l1)) continue;

                    bool rev = InvertPrototypePair(ABp,CDp);

                    int ogtn = omp_get_thread_num();

                    if (rev) {
                        ERIBatch C1batch;
                        C1batch.Set(AAAA, ABp, CDp, (ab==cd));
                        C1batch.Ntiles64 = 1;

                        BatchEval[ogtn].SetBuffers(C1batch);
                        BatchEval[ogtn].ConstructW(C1batch, ABs,CDs,eriG,WWW,wa, true); //OnlyJ
                    }
                    else {
                        ERIBatch C1batch;
                        C1batch.Set(AAAA, CDp, ABp, (ab==cd));
                        C1batch.Ntiles64 = 1;

                        BatchEval[ogtn].SetBuffers(C1batch);
                        BatchEval[ogtn].ConstructW(C1batch, CDs,ABs,eriG,WWW,wa, true); //OnlyJ
                    }

                }
            }

            #pragma omp parallel for schedule(dynamic)
            for (int ap12=0; ap12<APB12.APlist.n; ++ap12) {
                const AtomProd & AP12 = APB12.APlist[ap12];
                int at1 = AP12.at1;
                Sparse & Jc = BatchEval[omp_get_thread_num()].Jcont;

                AddCoulomb(Jc, at1,  wa, WWW, D2sparse);
            }


            if (XS == 0.) {delete[] WWW; continue;} //skip exchange if not needed

            //EXCHANGE
            //********

            for (int i=0; i<wa4; ++i) WWW[i]=0;

            //calcula las integrales para SOLO un atomo de cada tipo
            #pragma omp parallel for schedule(dynamic)
            for (int ab=0; ab<APB12.APP->nGPt; ++ab) {
                const ShellPair & ABs = APB12.SP[ab][0];
                const ShellPairPrototype & ABp = APB12.APP->Gs[ab];

                for (int cd=0; cd<=ab; ++cd) {
                    const ShellPair & CDs = APB12.SP[cd][0];
                    const ShellPairPrototype & CDp = APB12.APP->Gs[cd];

                    if ((ABp.l1+ABp.l2+CDp.l1+CDp.l2)%2) continue;
                    if ((ABp.l2==0) && (CDp.l2==0) && (ABp.l1!=CDp.l1)) continue;

                    bool rev = InvertPrototypePair(ABp,CDp);

                    int ogtn = omp_get_thread_num();

                    if (rev) {
                        ERIBatch C1batch;
                        C1batch.Set(AAAA, ABp, CDp, (ab==cd));
                        C1batch.Ntiles64 = 1;

                        BatchEval[ogtn].SetBuffers(C1batch);
                        BatchEval[ogtn].ConstructW(C1batch, ABs,CDs,eriG,WWW,wa, false);
                    }
                    else {
                        ERIBatch C1batch;
                        C1batch.Set(AAAA, CDp, ABp, (ab==cd));
                        C1batch.Ntiles64 = 1;

                        BatchEval[ogtn].SetBuffers(C1batch);
                        BatchEval[ogtn].ConstructW(C1batch, CDs,ABs,eriG,WWW,wa, false);
                    }

                }
            }


            #pragma omp parallel for schedule(dynamic)
            for (int ap12=0; ap12<APB12.APlist.n; ++ap12) {
                const AtomProd & AP12 = APB12.APlist[ap12];
                int at1 = AP12.at1;
                Sparse & Xc = BatchEval[omp_get_thread_num()].Xcont;

                AddXchange(Xc, at1,  wa, WWW, D2sparse);
            }

            delete[] WWW;
        }
    }


    //evaluate all n-center integrals, n>1
    //************************************
    if    (XS>0.) ExecuteJobs<true, true > (logD, logS); // for DFT with no exact exchange contribution
    else          ExecuteJobs<true, false> (logD, logS); // for DFT with exact exchange contribution or HF

    ERIchrono.Stop();

    if (PrintStatistics) {

        Q.Statistics();

        double Gsec = 0;
        double Fsec = 0;
        double Ksec = 0;
        double Tsec = 0;
        double Csec = 0;

        UI64 NJs = 0;
        UI64 NXs = 0;

        for (int i=0; i<omp_get_max_threads(); ++i) {
            Gsec += BatchEval[i].chronoG.GetTotalTime();
            Fsec += BatchEval[i].chronoF.GetTotalTime();
            Ksec += BatchEval[i].chronoK.GetTotalTime();
            Tsec += BatchEval[i].chronoT.GetTotalTime();
            Csec += BatchEval[i].chronoC.GetTotalTime();

            NJs  += Blocks[i].nJs;
            NXs  += Blocks[i].nXs;
        }

        Echidna::EBenchmarker << "          Cumulative time in geometry     : " << Gsec << endl;
        Echidna::EBenchmarker << "          Cumulative time in gamma        : " << Fsec << endl;
        Echidna::EBenchmarker << "          Cumulative time in kernels      : " << Ksec << endl;
        Echidna::EBenchmarker << "          Cumulative time in transforms   : " << Tsec << endl;
        Echidna::EBenchmarker << "          Cumulative time in contractions : " << Csec << endl;

        Echidna::EBenchmarker << "          Number of Coulomb block integrals  : " << NJs << endl;
        Echidna::EBenchmarker << "          Number of Exchange block integrals : " << NXs << endl;
    }

    //print more statistics
    if (0) {
        double tJs, pJs, tXs, pXs;
        tJs = pJs = tXs = pXs = 0;

        UI64 nJs, mJs, nXs, mXs;
        nJs = mJs = nXs = mXs = 0;

        for (int i=0; i<omp_get_max_threads(); ++i) {
            tJs  += Blocks[i].tJs;
            pJs  += Blocks[i].pJs;
            tXs  += Blocks[i].tXs;
            pXs  += Blocks[i].pXs;

            nJs  += Blocks[i].nJs;
            mJs  += Blocks[i].mJs;
            nXs  += Blocks[i].nXs;
            mXs  += Blocks[i].mXs;
        }

        Echidna::EMessenger << "Total       J " << tJs << endl;
        Echidna::EMessenger << "Partial     J " << pJs << endl;
        Echidna::EMessenger << "Fraction of J " << pJs/tJs << endl;

        Echidna::EMessenger << "Total       X " << tXs << endl;
        Echidna::EMessenger << "Partial     X " << pXs << endl;
        Echidna::EMessenger << "Fraction of X " << pXs/tXs << endl;

        Echidna::EMessenger << "Total blocks   J " << mJs << endl;
        Echidna::EMessenger << "Partial blocks J " << nJs << endl;

        Echidna::EMessenger << "Total blocks   X " << mXs << endl;
        Echidna::EMessenger << "Partial blocks X " << nXs << endl;
    }

    Echidna::EMessenger.Pop();

    //add each thread's contributions to the final fock matrix
    {

        //reduce the coulomb matrix
        {
            //potential race condition avoided
            for (int e=1; e<OGMT; ++e) BatchEval[0].Jcont += BatchEval[e].Jcont;

            //reduce in a tree-like fashion
            /*
            #pragma omp parallel
            {
                int id  = omp_get_thread_num();

                for (int e=1; e<OGMT; e*=2) {
                    if (id%(2*e) == 0 && id+e<OGMT) {
                        BatchEval[id].Jcont += BatchEval[id+e].Jcont;
                    }

                    #pragma omp barrier
                }
            }
            */
        }



        //add exchange contributions (if needed)
        if (XS>0.) {
            
            for (int e=1; e<OGMT; ++e) BatchEval[0].Xcont += BatchEval[e].Xcont;

            /*
            #pragma omp parallel
            {
                //exchange matrix
                int id  = omp_get_thread_num();

                for (int e=1; e<OGMT; e*=2) {
                    if (id%(2*e) == 0 && id+e<OGMT) {
                        BatchEval[id].Xcont += BatchEval[id+e].Xcont;
                    }

                    #pragma omp barrier
                }
            }
            */
        }

        //add exchange contributions (if needed)
        if (XS>0.) {
            //scale exchange
            if (XS<1.) BatchEval[0].Xcont *= XS;
            //add coulomb and exchange
            BatchEval[0].Jcont += BatchEval[0].Xcont;
        }


        //sums JX, reflecting the upper right half over the diagonal and adding it
        //it also doubles the values at the diagonal
        {
            F = BatchEval[0].Jcont;
            F.Symmetrize();
        }


        //reduce - broadcast at the MPI level
        Tensor2Reduce(F);
    }

}


template <bool J, bool X> void Fock2e::ExecuteJobs         (float logD, float logS) {

    enum JobDescription{JOB_NONE, JOB_MAKE_BATCH, JOB_COMPUTE_BATCH};

    if (UseCASE) Echidna::EDebugger << "Warning: using CASE operator without splitting of Coulomb and Exchange contributions!" << endl;


    if (J && X)
        Echidna::EMessenger << "Computing 2 electron Fock contributions";
    else if (J && !X)
        Echidna::EMessenger << "Computing Coulomb interaction";
    else if (!J && X)
        Echidna::EMessenger << "Computing exchange interaction";

    Echidna::EMessenger.Push();


    //copy the JobQueue (which is a set) to a FIFO queue
    TaskQueue<BatchInfo> BlockList;
    {
        set<BatchInfo>::const_iterator it;

        for (it=TwoElectronJobs.JobQueue.begin(); it!=TwoElectronJobs.JobQueue.end(); ++it) {
            BatchInfo * BI =  new BatchInfo;
            *BI = *it;
            //const BatchInfo * BI = &(*it);
            BlockList.push(BI);
        }
    }

    TaskQueue<ERIBatch> OCLDList;

    //begin
    volatile int workingthreads = omp_get_num_threads();

    #pragma omp parallel shared(workingthreads)
    {
        UI32 tn = omp_get_thread_num();
        bool lastfailed = false;

        ERIBatch  * OCLD;
        BatchInfo * Block;

        while (workingthreads>0 || !lastfailed ) {

            JobDescription JobType = JOB_NONE;

            //first try to retrieve a batch
            OCLD = OCLDList.pop();

            if (OCLD!=NULL) JobType = JOB_COMPUTE_BATCH;
            //if it doesn't succeed, try to retrieve a block from the queue
            else {
                Block = BlockList.pop();
                if (Block!=NULL) JobType = JOB_MAKE_BATCH;
            }

            //never actually executed
            if (lastfailed && JobType!=JOB_NONE) {
                lastfailed = false;
                #pragma omp atomic
                ++workingthreads;
            }

            //OCLD PATHJOB_MAKE_BATCH
            if      (JobType==JOB_COMPUTE_BATCH) {
                //cout << "Compute!" << endl;
                BatchEval[tn].Evaluate (*OCLD, D2sparse);
                OCLD->clear(EbatchBuffer);
                delete OCLD; //recycle ?
            }
            //
            else if (JobType==JOB_MAKE_BATCH) {
                //cout << "Generate!" << endl;
                bool finish;
                Blocks[tn] = *Block;
                if (J &&  X) finish = Blocks[tn].MakeAll      (logD, logS, logAD, logDs, OCLDList);
                if (J && !X) finish = Blocks[tn].MakeCoulomb  (logD, logS, logAD, logDs, OCLDList);
                if (!J && X) finish = Blocks[tn].MakeExchange (logD, logS, logAD, logDs, OCLDList);
                if (finish) delete Block;
                else        {
                    Block->State    = Blocks[tn].State;
                    BlockList.push_top(Block);
                }
            }

            else {
                if (!lastfailed) {
                    lastfailed = true;
                    #pragma omp atomic
                    --workingthreads;
                }
            }
        }
    }

    if (BlockList.length() != 0) Echidna::EDebugger << "Warning! not all ERI blocks were successfully evaluated!" << endl;
    if (OCLDList.length()  != 0) Echidna::EDebugger << "Warning! not all ERI blocks were successfully evaluated!" << endl;

    EbatchBuffer->Check();

    Echidna::EMessenger.Pop();
}
