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


#include <iostream>
#include <set>

#include "../defs.hpp"

#include "IICgen.hpp"
#include "SKS.hpp"

#include "../2eints/IICinit.hpp"
#include "../math/angular.hpp"


using namespace std;
using namespace LibAngular;


static const int PF = 8; //prefetch diostabnce


void ERIroutine::Set(UI8 La, UI8 Lb, UI8 Lc, UI8 Ld, GEOM geom, bool cdr=false) {
    geometry = geom;
    la = La;
    lb = Lb;
    lc = Lc;
    ld = Ld;
    IsSet = true;
    useCDR = cdr;

    for (int i=0; i<32; ++i) ninstrK4[i] = 0;
    for (int i=0; i<32; ++i) ninstr[i] = 0;

    InnerContractionRoutine = NULL;
}

void ERIroutine::Write() {

    if (IsInitialized) return;
    IsInitialized = true;

    Lt = la+lb+lc+ld;
    Am = la+lb+lc+ld+1;


    bool ABx, ABy, ABz;
    bool CDx, CDy, CDz;
    bool ACx, ACy, ACz;

    bool ABt, CDt, ACt;

    //initialize

    switch(geometry) {
        case ABCD:
            ABx = ABy = CDx = false;
            ABz = CDy = CDz = true;
            ACx = ACy = ACz = true;
        break;

        case ABAD:
            ABx = ABy = CDx = false;
            ABz = CDy = CDz = true;
            ACx = ACy = ACz = true;
        /*
            ABx = ABy = CDx = false;
            ABz = CDy = CDz = true;
            ACx = ACy = ACz = false;
            */
        break;

        case AACD:
            ABx = ABy = ABz = CDx = CDy = false;
            CDz = true;
            ACx = false;
            ACy = ACz = true;
        break;

        case AACC:
            ABx = ABy = ABz = CDx = CDy = CDz = false;
            ACx = ACy = false;
            ACz = true;
        break;

        case ABAB:
            ABx = ABy = CDx = CDy = false;
            ABz = CDz = true;
            ACx = ACy = false;
            ACz = true;
        break;

        case AAAD:
            ABx = ABy = ABz = CDx = CDy = false;
            CDz = true;
            ACx = ACy = ACz = false;
        break;

        case AAAA:
            ABx = ABy = ABz = CDx = CDy = CDz = false;
            ACx = ACy = ACz = false;
        break;


        default:
            ABx = ABy = ABz = CDx = CDy = CDz = false;
            ACx = ACy = ACz = false;
        break;
    }

    ABt = ABx or ABy or ABz;
    CDt = CDx or CDy or CDz;
    ACt = ACx or ACy or ACz;

    //don't use CDR for one-center integrals (will crash) or for center-degenerate ABCD integrals
    //(produces numerical inaccuracy in some cases)
    useCDR =  false; // (geometry==ABCD); // (geometry==AACD) || (geometry==AACC) || (geometry==AACD) || (geometry==ABCD); //(geometry!=AAAA) && (geometry!=ABAB)  && (geometry!=ABAD) && (geometry!=ABCD);
    useGC = true;


    EvaluationScheme Scheme;

    {
        //populates the list of needed integrals
        {
            set<integralCart> setERI;
            set<ikernel> UVmST;

            for (UI8 a=0; a<nmS[la]; ++a) {
                for (UI8 b=0; b<nmS[lb]; ++b) {
                    for (UI8 c=0; c<nmS[lc]; ++c) {
                        for (UI8 d=0; d<nmS[ld]; ++d) {
                            integralCart ERI(a,b,c,d);
                            ERI.m = -1; //mark it
                            setERI.insert(ERI);
                        }
                    }
                }
            }

            Scheme.Set(la,lb,lc,ld,geometry,useGC);

            //MIRROR
            {
                set<integralCart> setSSSS;
                set<integralCart> setSSSC;
                set<integralCart> setSSCC;
                set<integralCart> setSCCC;
                set<integralCart> setCCCC;

                set<integralCart> setHRRkx;
                set<integralCart> setHRRbx;
                set<integralCart> setHRRky;
                set<integralCart> setHRRby;
                set<integralCart> setHRRkz;
                set<integralCart> setHRRbz;

                set<integralCart> setCTEkx;
                set<integralCart> setCTEbx;
                set<integralCart> setCTEky;
                set<integralCart> setCTEby;
                set<integralCart> setCTEkz;
                set<integralCart> setCTEbz;

                set<integralCart> setRx;
                set<integralCart> setRy;
                set<integralCart> setRz;


                Scheme.addInitial(setSSSS);

                Scheme.addSphD (setSSSS , setSSSC);
                Scheme.addSphC (setSSSC , setSSCC);

                Scheme.addHRRkx(setSSCC , setHRRkx, CDx);
                Scheme.addHRRky(setHRRkx, setHRRky, CDy);
                Scheme.addHRRkz(setHRRky, setHRRkz, CDz);

                Scheme.addSphB (setHRRkz, setSCCC);
                Scheme.addSphA (setSCCC , setCCCC);

                Scheme.addHRRbx(setCCCC , setHRRbx, ABx);
                Scheme.addHRRby(setHRRbx, setHRRby, ABy);
                Scheme.addHRRbz(setHRRby, setHRRbz, ABz);

                Scheme.addCTEbx(setHRRbz, setCTEbx, ABx);
                Scheme.addCTEby(setCTEbx, setCTEby, ABy);
                Scheme.addCTEkx(setCTEby, setCTEkx, CDx);
                Scheme.addCTEky(setCTEkx, setCTEky, CDy);

                Scheme.addRxcc (setCTEky, setRx   , ABx, CDx, ACx);
                Scheme.addRycc (setRx   , setRy   , ABy, CDy, ACy);
                Scheme.addCTEkz(setRy   , setCTEkz, CDz);
                Scheme.addCTEbz(setCTEkz, setCTEbz, ABz);
                Scheme.addRzcc (setCTEbz, setRz   , ABz,  CDz, ACz);

                Scheme.MakeKernels(setRz, UVmST);
            }


            //generate the set of spherical kernels
            set<ikernel> sphKs;
            TestSKS(sphKs, geometry, la, lb, lc, ld);


            //K4
            if (!useCDR) {
                set<ikernel> UVmT;
                set<ikernel> UVm;
                set<ikernel> Vm;
                set<ikernel> m;
                set<ikernel> e;

                Scheme.addS(sphKs, UVmT);
                Scheme.addT(UVmT, UVm);
                Scheme.addU(UVm, Vm);
                Scheme.addV(Vm, m);

                if (geometry!=AAAA) {
                    Scheme.UseCDR4 (m, e, Lt);
                    Scheme.UseAERR4(e, Lt);
                }
                else
                    Scheme.add00(Lt);
            }
            else {
                set<ikernel> K1fe;
                set<ikernel> K2fe;


                set<ikernel> K0c, K0b, K0f, K0e;
                set<ikernel> K1c, K1b, K1f, K1e;
                set<ikernel> K2c, K2b, K2f, K2e;
                set<ikernel> K3c, K3b, K3f, K3e;
                set<ikernel> K4f, K4e;

                Scheme.UseCDR0  (sphKs, K0e, K0c, ABt, CDt, ACt);
                Scheme.UseAERR0 (K0e, K0b);

                Scheme.addS(K0c,  K1f);
                Scheme.addS(K0b,  K1e);

                Scheme.UseCDR1  (K1f, K1e, K1c, ABt, CDt, ACt);
                Scheme.UseAERR1 (K1e, K1b);

                Scheme.addT(K1c,  K2f);
                Scheme.addT(K1b,  K2e);

                Scheme.UseCDR2  (K2f, K2e, K2c, ABt, CDt, ACt);
                Scheme.UseAERR2 (K2e, K2b);

                Scheme.addU (K2c, K3f);
                Scheme.addU (K2b, K3e);

                Scheme.UseCDR3  (K3f, K3e, K3c, ABt, CDt, ACt);
                Scheme.UseAERR3 (K3e, K3b);

                Scheme.addV (K3c, K4f);
                Scheme.addV (K3b, K4e);

                if (geometry!=AAAA) {
                    Scheme.UseCDR4 (K4f, K4e, Lt);
                    Scheme.UseAERR4(K4e, Lt);
                }
                else
                    Scheme.add00(Lt);
            }


            Scheme.LinkVarsK4(UVmST);
            Scheme.LinkVars(setERI);

            Scheme.LinkBackK4();
            Scheme.LinkBack();

            Scheme.Simplify();
            Scheme.AssignMemory();

            Scheme.GenerateCode();
        }

        //K4
        {
            UI32 totalK4 = 0;
            for (UI8 i=BOYS; i<NADA; ++i) ninstrK4[i] = Scheme.ninstrK4[i];
            ninstrK4[NADA] = 0;

            for (UI8 iblock=BOYS; iblock<NADA; ++iblock) totalK4 += ninstrK4[iblock];

            eseqK4 = new Op_K4[totalK4 + PF]; //add the prefetching space for the last

            int p=0;


            for (UI8 iblock=BOYS; iblock<=NADA; ++iblock) {
                nseqK4[iblock] = eseqK4 + p;

                for (int i=0; i<ninstrK4[iblock]; ++i) {
                    nseqK4[iblock][i].dest  = Scheme.seqsK4[iblock][i].dest;
                    nseqK4[iblock][i].op1   = Scheme.seqsK4[iblock][i].op1;
                    nseqK4[iblock][i].op2   = Scheme.seqsK4[iblock][i].op2;
                    nseqK4[iblock][i].op3   = Scheme.seqsK4[iblock][i].op3;
                    nseqK4[iblock][i].op4   = Scheme.seqsK4[iblock][i].op4;
                    nseqK4[iblock][i].op5   = Scheme.seqsK4[iblock][i].op5;
                    nseqK4[iblock][i].op6   = Scheme.seqsK4[iblock][i].op6;
                    nseqK4[iblock][i].ope   = Scheme.seqsK4[iblock][i].ope;
                    nseqK4[iblock][i].aux   = Scheme.seqsK4[iblock][i].aux;
                }
                p += ninstrK4[iblock];
            }
        }

        //MIRROR
        {
            for (UI8 i=KERNELS; i<=REORDER; ++i) ninstr[i] = Scheme.ninstr[i];

            UI32 total = 0;
            for (UI8 iblock=KERNELS; iblock<=REORDER; ++iblock) total += ninstr[iblock];

            eseq = new Op_MIRROR[total + PF]; //add the prefetching space for the last

            int p=0;

            for (UI8 iblock=KERNELS; iblock<=REORDER; ++iblock) {
                for (int i=0; i<ninstr[iblock]; ++i) {
                    eseq[p].dest  = Scheme.seqs[iblock][i].dest;
                    eseq[p].op1   = Scheme.seqs[iblock][i].op1;
                    eseq[p].op2   = Scheme.seqs[iblock][i].op2;
                    eseq[p].op3   = Scheme.seqs[iblock][i].op3;
                    eseq[p].op4   = Scheme.seqs[iblock][i].op4;
                    eseq[p].op5   = Scheme.seqs[iblock][i].op5;
                    eseq[p].op6   = Scheme.seqs[iblock][i].op6;
                    eseq[p].aux   = Scheme.seqs[iblock][i].aux;
                    ++p;
                }
            }

        }

        MaxMem  = Scheme.MaxMem;
        NFLOPS  = Scheme.NFLOPS;

        WriteIC2File();
    }

    //COPY EVERYTHING
    //===============

    {
        nKernels =  ninstr[KERNELS];

        // set pointers; observe that K4 and MIRROR behave differently regarding the instruction pointer
        int p=0, q=0;

        for (UI8 iblock=BOYS; iblock<=NADA; ++iblock) {
            nseqK4[iblock] = eseqK4 + p;
            p += ninstrK4[iblock];
        }



        for (UI8 iblock=KERNELS; iblock<=REORDER; ++iblock) {
            q += ninstr[iblock];
            nseq[iblock] = eseq + q;
        }


        // find the largest powers of the gaussian exponents needed to compute the kernels

        maxV = 0;
        for (int i=0; i<ninstrK4[K4E]; ++i) maxV = max(maxV, nseqK4[K4E][i].aux);
        for (int i=0; i<ninstrK4[K4F]; ++i) maxV = max(maxV, nseqK4[K4F][i].aux);
        maxU = 0;
        for (int i=0; i<ninstrK4[K3E]; ++i) maxU = max(maxU, nseqK4[K3E][i].aux);
        for (int i=0; i<ninstrK4[K3F]; ++i) maxU = max(maxU, nseqK4[K3F][i].aux);
        maxT = 0;
        for (int i=0; i<ninstrK4[K2E]; ++i) maxT = max(maxT, nseqK4[K2E][i].aux);
        for (int i=0; i<ninstrK4[K2F]; ++i) maxT = max(maxT, nseqK4[K2F][i].aux);
        maxS = 0;
        for (int i=0; i<ninstrK4[K1E]; ++i) maxS = max(maxS, nseqK4[K1E][i].aux);
        for (int i=0; i<ninstrK4[K1F]; ++i) maxS = max(maxS, nseqK4[K1F][i].aux);


/*

        // make a local copy of the Spherical Harmonic coefficients

        if (la>1) {
            for (UI8 m=0; m<2*la+1; ++m) {
                for (int n=0; n<6; ++n) Ca[m][n] = 0;
                for (int n=0; n<SHList[la][m].nps; ++n) Ca[m][n] = SHList[la][m].T[n].cN;
            }

            for (UI8 m=0; m<2*la+1; ++m) {
                for (int n=0; n<6; ++n) fCa[m][n] = Ca[m][n];
            }
        }

        if (lb>1) {
            for (UI8 m=0; m<2*lb+1; ++m) {
                for (UI8 n=0; n<6; ++n) Cb[m][n] = 0;
                for (UI8 n=0; n<SHList[lb][m].nps; ++n) Cb[m][n] = SHList[lb][m].T[n].cN;
            }

            for (UI8 m=0; m<2*lb+1; ++m) {
                for (int n=0; n<6; ++n) fCb[m][n] = Cb[m][n];
            }
        }

        if (lc>1) {
            for (UI8 m=0; m<2*lc+1; ++m) {
                for (UI8 n=0; n<6; ++n) Cc[m][n] = 0;
                for (UI8 n=0; n<SHList[lc][m].nps; ++n) Cc[m][n] = SHList[lc][m].T[n].cN;
            }

            for (UI8 m=0; m<2*lc+1; ++m) {
                for (int n=0; n<6; ++n) fCc[m][n] = Cc[m][n];
            }
        }

        if (ld>1) {
            for (UI8 m=0; m<2*ld+1; ++m) {
                for (UI8 n=0; n<6; ++n) Cd[m][n] = 0;
                for (UI8 n=0; n<SHList[ld][m].nps; ++n) Cd[m][n] = SHList[ld][m].T[n].cN;
            }

            for (UI8 m=0; m<2*ld+1; ++m) {
                for (int n=0; n<6; ++n) fCd[m][n] = Cd[m][n];
            }
        }
*/

        // pad the IC array so that the last data prefecthes won't attempt to load garbage into the cache

        K4Mem   = 0;
        for (UI8 i=BOYS; i<NADA; ++i) K4Mem += ninstrK4[i];

        for (int i=0; i<PF; ++ i) {
            UI32 p = K4Mem + i;
            eseqK4[p].dest = 0;
            eseqK4[p].op1  = 0;
            eseqK4[p].op2  = 0;
            eseqK4[p].op3  = 0;
            eseqK4[p].op4  = 0;
            eseqK4[p].op5  = 0;
            eseqK4[p].op6  = 0;
            eseqK4[p].ope  = 0;
            eseqK4[p].aux  = 0;
        }


        UI32 total = 0;
        for (UI8 iblock=KERNELS; iblock<=REORDER; ++iblock)
            total += ninstr[iblock];

        for (int i=0; i<PF; ++ i) {
            UI32 p = total + i;
            eseq[p].dest = 0;
            eseq[p].op1  = 0;
            eseq[p].op2  = 0;
            eseq[p].op3  = 0;
            eseq[p].op4  = 0;
            eseq[p].op5  = 0;
            eseq[p].op6  = 0;
            eseq[p].aux  = 0;
        }


        // number of FLOPs for each block

        nK4J0 =   ninstrK4[AERR4] +  3*ninstrK4[CDR4];
        nK4J1 = 2*ninstrK4[K4E]   +  2*ninstrK4[K4F];
        nK3J1 = 2*ninstrK4[AERR3] +  7*ninstrK4[CDR3] - ninstrK4[K4E] - ninstrK4[K4F];
        nK3J2 = 2*ninstrK4[K3E]   +  2*ninstrK4[K3F];
        nK2J2 = 2*ninstrK4[AERR2] +  7*ninstrK4[CDR2] - ninstrK4[K3E] - ninstrK4[K3F];
        nK2J3 = 2*ninstrK4[K2E]   +  2*ninstrK4[K2F];
        nK1J3 =   ninstrK4[AERR1] + 13*ninstrK4[CDR1] - ninstrK4[K2E] - ninstrK4[K2F];
        nK1J4 = 2*ninstrK4[K1E]   +  2*ninstrK4[K1F];
        nK0J4 =   ninstrK4[AERR0] + 13*ninstrK4[CDR0] - ninstrK4[K1E] - ninstrK4[K1F];


        nK4J0e =   ninstrK4[AERR4]-1 + ninstrK4[CDR4]-1;
        nK4J1e = 2*ninstrK4[K4E];
        nK3J1e = 2*ninstrK4[AERR3] + ninstrK4[CDR3] - ninstrK4[K4E];
        nK3J2e = 2*ninstrK4[K3E];
        nK2J2e = 2*ninstrK4[AERR2] + ninstrK4[CDR2] - ninstrK4[K3E];
        nK2J3e = 2*ninstrK4[K2E];
        nK1J3e =   ninstrK4[AERR1] + ninstrK4[CDR1] - ninstrK4[K2E];
        nK1J4e = 2*ninstrK4[K1E];
        nK0J4e =   ninstrK4[AERR0] + ninstrK4[CDR0] - ninstrK4[K1E];

        nK4J0f =  2*ninstrK4[CDR4]-2;
        nK4J1f =  2*ninstrK4[K4F];
        nK3J1f =  6*ninstrK4[CDR3] - ninstrK4[K4F];
        nK3J2f =  2*ninstrK4[K3F];
        nK2J2f =  6*ninstrK4[CDR2] - ninstrK4[K3F];
        nK2J3f =  2*ninstrK4[K2F];
        nK1J3f = 12*ninstrK4[CDR1] - ninstrK4[K2F];
        nK1J4f =  2*ninstrK4[K1F];
        nK0J4f = 12*ninstrK4[CDR0] - ninstrK4[K1F];


        // memory required for the K4 contraction
        memF0  =     2; //ninstrK4[BOYS];
        memF0e =  Lt+1; //ninstrK4[AERR4];
        memF0f =  Lt+1; //ninstrK4[CDR4];

        memK4J1e =  ninstrK4[K4E]  + ninstrK4[AERR3];
        memK4J1f =  ninstrK4[K4F]  + ninstrK4[CDR3];
        memK3J2e =  ninstrK4[K3E]  + ninstrK4[AERR2];
        memK3J2f =  ninstrK4[K3F]  + ninstrK4[CDR2];
        memK2J3e =  ninstrK4[K2E]  + ninstrK4[AERR1];
        memK2J3f =  ninstrK4[K2F]  + ninstrK4[CDR1];
        memK1J4e =  ninstrK4[K1E]  + ninstrK4[AERR0];
        memK1J4f =  ninstrK4[K1F]  + ninstrK4[CDR0];
    }

    //K4VILcode.Init(ninstrK4, nseqK4);

    // do this only for 4-center integrals
    if (geometry==ABCD) WriteInnerContraction();
}



string ERIroutine::IdString(bool fuse) const {

    string fname;

    string geo;

    switch(geometry) {
        case ABCD:
            geo = "ABCD";
        break;

        case ABAD:
            geo = "ABAD";
        break;

        case AACD:
            geo = "AACD";
        break;

        case AACC:
            geo = "AACC";
        break;

        case ABAB:
            geo = "ABAB";
        break;

        case AAAA:
            geo = "AAAA";
        break;

        case NONE:
            geo = "NONE";
        break;

        default:
        break;
    }

    fname = geo + "_";

    fname += L2S(la);
    fname += L2S(lb);

    if (fuse)
        fname += L2S(lc+ld);
    else {
        fname += L2S(lc);
        fname += L2S(ld);
    }

    return fname;
}

#include "../ICgen/QICgen.hpp"

void ERIroutine::WriteIC2File() const {

    string fname = "IC_" + IdString() + ".qic";
    string fdir = QIC_DIR;

    if (fdir=="") return; //don't write
    fname = fdir + "/" + fname;
    //write to file
    ofstream file;
    file.open(fname.c_str(), ios::out | ios::binary);

    // write info about the routine

    file.write ((char*)(&la), sizeof (UI8));
    file.write ((char*)(&lb), sizeof (UI8));
    file.write ((char*)(&lc), sizeof (UI8));
    file.write ((char*)(&ld), sizeof (UI8));
    file.write ((char*)&geometry, sizeof (GEOM));
    file.write ((char*)&useCDR,   sizeof (bool));
    file.write ((char*)&useGC,    sizeof (bool));


    //K4
    {
        //number of instructions in each step
        for (UI8 i=0; i<32; ++i)     file.write ((char*)&ninstrK4[i], sizeof (UI32));

        int totalK4 = 0;
        for (UI8 iblock=BOYS; iblock<NADA; ++iblock) totalK4 += ninstrK4[iblock];

        //actual instructions
        for (UI32 i=0; i<totalK4; ++i) file.write ((char*)&eseqK4[i], sizeof (Op_K4));
    }

    //MIRROR
    {
        //number of instructions in each step
        for (UI8 i=0; i<32; ++i)     file.write ((char*)&ninstr[i], sizeof (UI32));

        int total = 0;
        for (UI8 iblock=KERNELS; iblock<=REORDER; ++iblock)  total += ninstr[iblock];

        for (UI32 i=0; i<total; ++i) file.write ((char*)&eseq[i], sizeof (Op_MIRROR));
    }

    file.write ((char*)&MaxMem, sizeof (UI32));
    file.write ((char*)&NFLOPS, sizeof (UI32));

    file.close();
}


void ERIroutine::WriteLoopsContraction2(ofstream & file, const string & D, const string & S, const Op_K4 * IICbegin, const Op_K4 * IICend, bool firstloop) {

    for (const Op_K4 * s2=IICbegin; s2<IICend; ++s2) {
        //const cacheline64 * op1 = S + s2->op1;
        //cacheline64       * dest = D + s2->dest;
        //UI16 aux                 =     s2->aux;

        //__builtin_prefetch(D + s2[PF].dest, 1, 1);
        if (firstloop)
            file << "        store("<<D<<"+"<<s2->dest<<", load("<<S<<"+"<<s2->op1<<") * Ws["<<s2->aux<<"]);" << endl; // " << s2->K << endl;
        else
            file << "        store("<<D<<"+"<<s2->dest<<", load("<<D<<"+"<<s2->dest<<") + load("<<S<<"+"<<s2->op1<<") * Ws["<<s2->aux<<"]);" << endl; // " << s2->K << endl;
    }

}


void ERIroutine::WriteLoopsContraction(ofstream & file, const string & D, const string & S,
                                int memOffset,
                                const string & W, int maxW,
                                const string & NN,
                                const string & JJ,
                                const Op_K4 * IICbegin,
                                const Op_K4 * IICend,
                                bool firstloop
                                ) {

    //for (int j=0; j<min(J, k+1); ++j)

    file << "      for (int j=0; j<"<<JJ<<"; ++j) {" << endl;
    file << "        double Ws["<<maxW+1<<"] __attribute__((aligned("<<CACHE_LINE_SIZE<<")));" << endl;

//        for (int w=0; w<=maxW; ++w)
//    file << "        Ws["<<w<<"] = "<<W<<"["<<w<<"] * "<<NN<<"[j];" << endl;

    file << "        Ws[0] = "<<NN<<"[j];" << endl;

        for (int w=1; w<=maxW; ++w)
    file << "        Ws["<<w<<"] = Ws["<<w-1<<"] * "<<W<<";" << endl;

        WriteLoopsContraction2(file, D, S, IICbegin, IICend, firstloop);

    file << "        "<<D<<" += "<<memOffset<<";" << endl;
    file << "      }" << endl;
}


void ERIroutine::WriteInnerContraction() {

    string ICRname = "K2C_" + IdString(true);
    string EFSdir = K2C_DIR;
    if (EFSdir=="") return; //skip

    string K2file = EFSdir + "/" + ICRname + ".cpp";

    ofstream file;
    file.open(K2file.c_str());
    file.precision(16);

    file << "#include <string.h>" << endl;
    file << "#include \"low/cache64.hpp\"" << endl;
    file << "#include \"integrals/rotations.hpp\"" << endl;
    file << "#include \"basis/SPprototype.hpp\"" << endl;
    file << "#include \"2eints/quimera.hpp\"" << endl;

    file << endl;

    file << "static const double im2[32] = {";
    for (int m=0; m<31; ++m)
        file << "1./" << 2*m+1 << "., ";
    file << "1./63.};" << endl;

    file << endl;



    file << "void " << ICRname << "(cacheline64 * (&F0), const ERIgeometries64 & vars8, const PrimitiveSet & PSab, const ShellPairPrototype & ABp, double ikcd, double rcd, p_ERIbuffer & buffer, cacheline64 & AB2, cacheline64 & X2) { " << endl;

    file << "  const cacheline64 & ABz = vars8.ABz;" << endl;
    file << "  const cacheline64 & CDy = vars8.CDy;" << endl;
    file << "  const cacheline64 & CDz = vars8.CDz;" << endl;

    file << "  const cacheline64 & ACx = vars8.ACx;" << endl;
    file << "  const cacheline64 & ACy = vars8.ACy;" << endl;
    file << "  const cacheline64 & ACz = vars8.ACz;" << endl;

    file << "  const int JA = ABp.Ja;" << endl;
    file << "  const int JB = ABp.Jb;" << endl;

    file << "  const int nJ1 = JA;" << endl;
    file << "  const int nJ2 = JA*JB;" << endl;

    file << "  const int mK3f = nJ2*"<<memK3J2f<<";" << endl;

    file << "  const int mK4f = nJ1*"<<memK4J1f<<";" << endl;

    file << "  const int mK3e = nJ2*"<<memK3J2e<<";" << endl;
    file << "  const int mK4e = nJ1*"<<memK4J1e<<";" << endl;

        if (geometry==ABCD || geometry==ABAD) {
    file << "  cacheline64 AQy;" << endl;
    file << "  cacheline64 AQz;" << endl;

    file << "  cacheline64 X2Y2;" << endl;
    file << "  cacheline64 AQ2;" << endl;

    file << "  cacheline64 AQAB;" << endl;

    file << "  AQy = ACy + CDy * rcd;" << endl;
    file << "  AQz = ACz + CDz * rcd;" << endl;

    file << "  X2Y2 = X2   + AQy*AQy;" << endl;
    file << "  AQ2  = X2Y2 + AQz*AQz;" << endl;

    file << "  AQAB = (AQz*ABz); AQAB += AQAB;" << endl;
        }
        else if (geometry==AACD) {
    file << "  cacheline64 AQz;" << endl;
    file << "  cacheline64 AQ2;" << endl;

    file << "  AQz = ACz + CDz * rcd;" << endl;
    file << "  AQ2  = X2 + AQz*AQz;" << endl;
        }



    file << "  memset(buffer.K3J2e, 0, mK3e*sizeof(cacheline64));" << endl;
    file << "  memset(buffer.K3J2f, 0, mK3f*sizeof(cacheline64));" << endl;


    file << "  for (int b=0; b<PSab.nKb; ++b) {" << endl; {
    file << "    {" << endl; {

        if (geometry==ABCD || geometry==ABAD) {
    file << "      double rab  = ABp.BCP.r[b][0];" << endl;
    file << "      cacheline64 PQz;  PQz = AQz - ABz*rab;" << endl;
    file << "      cacheline64 R2;   R2 = X2Y2 + PQz*PQz;" << endl;
                }

    file << "      double ikab = ABp.BCP.k[b][0];" << endl;
    file << "      double iKpq  = (ikab + ikcd);" << endl;

    file << "      cacheline64 * F0e = (cacheline64*)buffer.F0e;" << endl;
    file << "      cacheline64 * F0f = (cacheline64*)buffer.F0f;" << endl;

                //copy the gamma of higher L and the contracted exponential
    file << "      F0e[0] = F0[0];" << endl;
    file << "      F0f[0] = F0[1];" << endl;

                //AERR K4 (uncontracted)
                {
                    Op_K4 * s2=nseqK4[AERR4]; ++s2; //skip the first exponential
                    for (; s2<nseqK4[AERR4+1]; ++s2)
    file << "      F0e["<<s2->dest<<"] = F0e["<<s2->op1<<"] * iKpq;" << endl;// " << s2->K << endl;
                }

                //CDR K4 (uncontracted downward recursion)
                {
                    Op_K4 * s2=nseqK4[CDR4]; ++s2; //skip the first element, which is the gamma function
                    for (; s2<nseqK4[CDR4+1]; ++s2) {
                        if (geometry==ABCD || geometry==ABAD)
    file << "      F0f["<<s2->dest<<"] = (R2 *  F0f["<<s2->op1<<"] + F0e["<<s2->ope<<"]) * im2["<<s2->aux<<"];" << endl; // " << s2->K << endl;
                        else if (geometry==AACD)
    file << "      F0f["<<s2->dest<<"] = (AQ2 * F0f["<<s2->op1<<"] + F0e["<<s2->ope<<"]) * im2["<<s2->aux<<"];" << endl; // " << s2->K << endl;
                        else if (geometry==AACC)
    file << "      F0f["<<s2->dest<<"] = (X2  * F0f["<<s2->op1<<"] + F0e["<<s2->ope<<"]) * im2["<<s2->aux<<"];" << endl; // " << s2->K << endl;
                    }
                }


                //K4J1
                {
    file << "      cacheline64 * Je = (cacheline64*)buffer.K4J1e;" << endl;
    file << "      cacheline64 * Jf = (cacheline64*)buffer.K4J1f;" << endl;

    //for (int j=0; j<min(J, k+1); ++j)
    file << "      for (int j=0; j<JA; ++j) {" << endl;
    file << "        double Ws["<<maxV+1<<"] __attribute__((aligned("<<CACHE_LINE_SIZE<<")));" << endl;
    file << "        Ws[0] = ABp.BCP.Na[0][j];" << endl;
        for (int w=1; w<=maxV; ++w)
    file << "        Ws["<<w<<"] = Ws["<<w-1<<"] * ABp.BCP.k[b][0];" << endl;
        WriteLoopsContraction2(file, "Je", "F0e", nseqK4[K4E], nseqK4[K4E+1], true);
        WriteLoopsContraction2(file, "Jf", "F0f", nseqK4[K4F], nseqK4[K4F+1], true);
    file << "        Je += "<<memK4J1e<<";" << endl;
    file << "        Jf += "<<memK4J1f<<";" << endl;
    file << "      }" << endl;

                }

    file << "      F0 += 2;" << endl;
    file << "    }" << endl; }


    file << "    for (int a=1; a<PSab.nKa[b]; ++a) {" << endl; {

        if (geometry==ABCD || geometry==ABAD) {
    file << "      double rab  = ABp.BCP.r[b][a];" << endl;
    file << "      cacheline64 PQz;  PQz = AQz - ABz*rab;" << endl;
    file << "      cacheline64 R2;   R2 = X2Y2 + PQz*PQz;" << endl;
                }

    file << "      double ikab = ABp.BCP.k[b][a];" << endl;
    file << "      double iKpq  = (ikab + ikcd);" << endl;

    file << "      cacheline64 * F0e = (cacheline64*)buffer.F0e;" << endl;
    file << "      cacheline64 * F0f = (cacheline64*)buffer.F0f;" << endl;

                //copy the gamma of higher L and the contracted exponential
    file << "      F0e[0] = F0[0];" << endl;
    file << "      F0f[0] = F0[1];" << endl;


                //AERR K4 (uncontracted)
                {
                    Op_K4 * s2=nseqK4[AERR4]; ++s2; //skip the first exponential
                    for (; s2<nseqK4[AERR4+1]; ++s2)
    file << "      F0e["<<s2->dest<<"] = F0e["<<s2->op1<<"] * iKpq;" << endl; // " << s2->K << endl;
                }

                //CDR K4 (uncontracted downward recursion)
                {
                    Op_K4 * s2=nseqK4[CDR4]; ++s2; //skip the first element, which is the gamma function
                    for (; s2<nseqK4[CDR4+1]; ++s2) {
                        if (geometry==ABCD || geometry==ABAD)
    file << "      F0f["<<s2->dest<<"] = (R2 * F0f["<<s2->op1<<"] + F0e["<<s2->ope<<"]) * im2["<<s2->aux<<"];" << endl; // " << s2->K << endl;
                        else if (geometry==AACD)
    file << "      F0f["<<s2->dest<<"] = (AQ2 * F0f["<<s2->op1<<"] + F0e["<<s2->ope<<"]) * im2["<<s2->aux<<"];" << endl; // " << s2->K << endl;
                        else if (geometry==AACC)
    file << "      F0f["<<s2->dest<<"] = (X2  * F0f["<<s2->op1<<"] + F0e["<<s2->ope<<"]) * im2["<<s2->aux<<"];" << endl; // " << s2->K << endl;
                    }
                }




                {
    file << "      cacheline64 * Je = (cacheline64*)buffer.K4J1e;" << endl;
    file << "      cacheline64 * Jf = (cacheline64*)buffer.K4J1f;" << endl;

    file << "      for (int j=0; j<JA; ++j) {" << endl;
    file << "        double Ws["<<maxV+1<<"] __attribute__((aligned("<<CACHE_LINE_SIZE<<")));" << endl;
    file << "        Ws[0] = ABp.BCP.Na[a][j];" << endl;
        for (int w=1; w<=maxV; ++w)
    file << "        Ws["<<w<<"] = Ws["<<w-1<<"] * ABp.BCP.k[b][a];" << endl;
        WriteLoopsContraction2(file, "Je", "F0e", nseqK4[K4E], nseqK4[K4E+1], false);
        WriteLoopsContraction2(file, "Jf", "F0f", nseqK4[K4F], nseqK4[K4F+1], false);
    file << "        Je += "<<memK4J1e<<";" << endl;
    file << "        Jf += "<<memK4J1f<<";" << endl;
    file << "      }" << endl;

                }

    file << "      F0 += 2;" << endl;
    file << "    }" << endl; }



            //K3J1 //AERR K3
    file << "    {" << endl; {
    file << "      cacheline64 * Jf = (cacheline64*)buffer.K4J1f;" << endl;
    file << "      cacheline64 * Je = (cacheline64*)buffer.K4J1e;" << endl;

    file << "      __m128d ikcd_v = _mm_load1_pd(&ikcd);" << endl;
                //u_atom ikcd_v(ikcd);

    file << "      for (int nj1=0; nj1<nJ1; ++nj1) {" << endl;
                   for (Op_K4 * s2=nseqK4[AERR3]; s2<nseqK4[AERR3+1]; ++s2) {
    file << "        store(Je+"<<s2->dest<<", load(Je+"<<s2->op1<<") + load(Je+"<<s2->op2<<") * ikcd_v);" << endl; // " << s2->K << endl;
                    }

    file << "        Je += "<<memK4J1e<<";" << endl;
    file << "      }" << endl;
    file << "    }" << endl; }

            //CDR3
    file << "    {" << endl; {

                if (geometry==ABCD || geometry==ABAD) {
    file << "      cacheline64 AQABb; AQABb = AQAB * ABp.BCP.b1[b];" << endl;
    file << "      cacheline64 AB2b2; AB2b2 = AB2  * ABp.BCP.b2[b];" << endl;
                }

    file << "      cacheline64 * Jf = (cacheline64*)buffer.K4J1f;" << endl;
    file << "      cacheline64 * Je = (cacheline64*)buffer.K4J1e;" << endl;

    file << "      for (int nj1=0; nj1<nJ1; ++nj1) {" << endl;
                    for (Op_K4 * s2=nseqK4[CDR3]; s2<nseqK4[CDR3+1]; ++s2) {
                      if (geometry==ABCD || geometry==ABAD) {
    file << "        store(Jf+"<<s2->dest<<", ("
                      <<    " AQ2*load(Jf+"<<s2->op1<<")"
                     << " - AQABb*load(Jf+"<<s2->op2<<")"
                     << " + AB2b2*load(Jf+"<<s2->op3<<")"
                           << " + load(Je+"<<s2->ope<<")"
                           << ") * im2["<<s2->aux<<"] );" << endl; // " << s2->K << endl;
                      }
                      else if (geometry==AACD) {
    file << "        store(Jf+"<<s2->dest<<", ("
                      <<    " AQ2*load(Jf+"<<s2->op1<<")"
                           << " + load(Je+"<<s2->ope<<")"
                           << ") * im2["<<s2->aux<<"] );" << endl; // " << s2->K << endl;
                      }
                      else if (geometry==AACC) {
    file << "        store(Jf+"<<s2->dest<<", ("
                      <<    "  X2*load(Jf+"<<s2->op1<<")"
                           << " + load(Je+"<<s2->ope<<")"
                           << ") * im2["<<s2->aux<<"] );" << endl; // " << s2->K << endl;
                      }
                    }

    file << "        Jf += "<<memK4J1f<<";" << endl;
    file << "        Je += "<<memK4J1e<<";" << endl;
    file << "      }" << endl;
    file << "    }" << endl; }

            //K3E
    file << "    {" << endl; {
    file << "      const cacheline64 * Je = (cacheline64*)buffer.K4J1e;" << endl;
    file << "      const cacheline64 * Jf = (cacheline64*)buffer.K4J1f;" << endl;
    file << "      cacheline64       * J2e = (cacheline64*)buffer.K3J2e;" << endl;
    file << "      cacheline64       * J2f = (cacheline64*)buffer.K3J2f;" << endl;

    file << "      double Uu["<<int(maxJ)<<"]["<<maxU+1<<"] __attribute__((aligned("<<CACHE_LINE_SIZE<<")));" << endl;

    file << "      for (int j=0; j<JB; ++j) {" << endl;
    file << "        Uu[j][0] = ABp.BCP.Nb[b][j];" << endl;

    file << "        for (int u=1; u<="<<maxU<<"; ++u)" << endl;
    file << "            Uu[j][u] = Uu[j][u-1] * ABp.BCP.b1[b];" << endl;
    file << "      }" << endl;


    file << "      for (int nj1=0; nj1<nJ1; ++nj1) {" << endl;

    file << "        for (int j=0; j<JB; ++j) {" << endl;

    for (const Op_K4 * s2=nseqK4[K3E]; s2<nseqK4[K3E+1]; ++s2)
            file << "          store(J2e + "<<s2->dest<<"+j*"<<memK3J2e<<", load(J2e + "<<s2->dest<<"+j*"<<memK3J2e<<") + load(Je + "<<s2->op1<<") * Uu[j]["<<s2->aux<<"]);" << endl;

    for (const Op_K4 * s2=nseqK4[K3F]; s2<nseqK4[K3F+1]; ++s2)
            file << "          store(J2f + "<<s2->dest<<"+j*"<<memK3J2f<<", load(J2f + "<<s2->dest<<"+j*"<<memK3J2f<<") + load(Jf + "<<s2->op1<<") * Uu[j]["<<s2->aux<<"]);" << endl;

    file << "        }" << endl;

    file << "        J2e += JB*"<<memK3J2e<<";" << endl;
    file << "        J2f += JB*"<<memK3J2f<<";" << endl;
    file << "        Je += "<<memK4J1e<<";" << endl;
    file << "        Jf += "<<memK4J1f<<";" << endl;

    file << "      }" << endl;
    file << "    }" << endl; }


    file << "  }" << endl; } //loop over b


    //AERR K2
    file << "  {" << endl; {
    file << "    cacheline64 * J2 = (cacheline64*)buffer.K3J2e;" << endl;

    file << "    __m128d ikcd_v = _mm_load1_pd(&ikcd);" << endl;
            //u_atom ikcd_v(ikcd);

            //K2J2
    file << "    for (int nj2=0; nj2<nJ2; ++nj2) {" << endl;
                for (Op_K4 * s2=nseqK4[AERR2]; s2<nseqK4[AERR2+1]; ++s2) {
    file << "        store(J2+"<<s2->dest<<", load(J2+"<<s2->op1<<") + load(J2+"<<s2->op2<<") * ikcd_v);" << endl; // " << s2->K << endl;

                }
    file << "      J2 += "<<memK3J2e<<";" << endl;
    file << "    }" << endl;
    file << "  }" << endl; }

    //CDR K2
    file << "  {" << endl; {
    file << "    cacheline64 * J2 = (cacheline64*)buffer.K3J2f;" << endl;
    file << "    cacheline64 * Je = (cacheline64*)buffer.K3J2e;" << endl;

    file << "    for (int nj2=0; nj2<nJ2; ++nj2) {" << endl;
                for (Op_K4*s2 = nseqK4[CDR2]; s2<nseqK4[CDR2+1]; ++s2) {


                      if (geometry==ABCD || geometry==ABAD) {
    file << "        store(J2+"<<s2->dest<<", ("
                      <<    " AQ2*load(J2+"<<s2->op1<<")"
                      << " - AQAB*load(J2+"<<s2->op2<<")"
                       << " + AB2*load(J2+"<<s2->op3<<")"
                           << " + load(Je+"<<s2->ope<<")"
                           << ") * im2["<<s2->aux<<"] );" << endl; // " << s2->K << endl;
                      }
                      else if (geometry==AACD) {
    file << "        store(J2+"<<s2->dest<<", ("
                      <<    " AQ2*load(J2+"<<s2->op1<<")"
                           << " + load(Je+"<<s2->ope<<")"
                           << ") * im2["<<s2->aux<<"] );" << endl; // " << s2->K << endl;
                      }
                      else if (geometry==AACC) {
    file << "        store(J2+"<<s2->dest<<", ("
                      <<    " X2*load(J2+"<<s2->op1<<")"
                           << " + load(Je+"<<s2->ope<<")"
                           << ") * im2["<<s2->aux<<"] );" << endl; // " << s2->K << endl;
                      }
                }
    file << "      J2 += "<<memK3J2f<<";" << endl;
    file << "      Je += "<<memK3J2e<<";" << endl;
    file << "    }" << endl;
    file << "  }" << endl; }



    file << "}" << endl;


    file.close();
}




