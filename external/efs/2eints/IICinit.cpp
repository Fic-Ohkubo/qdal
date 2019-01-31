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
#include "../libquimera/libquimera.hpp"

#include "../math/angular.hpp"
#include "../2eints/IICinit.hpp"


#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

using namespace std;
using namespace LibQuimera;
using namespace LibAngular;

static const UI8 PF = 8; //prefetch distance


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

void ERIroutine::Init() {

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

    //don't use CDR for now
    useCDR =  false;
    useGC = true;


    {
        ReadICfromFile();
    }

    //COPY EVERYTHING
    //===============

    //grand finale
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

    //generate compact, variable-length instructions (improves cache performance)
    K4VILcode.Init(ninstrK4, nseqK4);

    //load the inner contraction specialized routines
    LoadK2C();
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

void ERIroutine::ReadICfromFile() {

    string fname, fdir;
    fname = "IC_" + IdString() + ".qic";
    fdir = IC_DIR;

    if (fdir=="") {
        cout << "Error: environment variable QIC_DIR not set" << endl;
        throw 1442; // random error code
    }
    fname = fdir + "/" + fname;
    //cout << fname << endl;


    //write to file
    ifstream file;
    file.open(fname.c_str(), ios::in | ios::binary);


    // read info about the routine
    UI8 LA, LB, LC, LD;
    GEOM geometry2;
    bool useCDRr, useGCr;

    file.read ((char*)(&LA), sizeof (UI8));
    file.read ((char*)(&LB), sizeof (UI8));
    file.read ((char*)(&LC), sizeof (UI8));
    file.read ((char*)(&LD), sizeof (UI8));
    file.read ((char*)&geometry2, sizeof (GEOM));
    file.read ((char*)&useCDRr,   sizeof (bool));
    file.read ((char*)&useGCr,    sizeof (bool));

    //check everything OK
    if (la!=LA || lb!=LB || lc!=LC || ld!=LD || geometry!=geometry2) {
        //print something
        //Quimera::QMessenger << "Error: IC file's angular momenta/geometry does not coincide" << endl;
        //Quimera::QMessenger << int(la) << " " << int(lb) << " " << int(lc) << " " << int(ld) << "  " <<  geometry << endl;
        //Quimera::QMessenger << int(LA) << " " << int(LB) << " " << int(LC) << " " << int(LD) << "  " <<  geometry2 << endl;

        cout << "Error: IC file's angular momenta/geometry does not coincide" << endl;
        cout << int(la) << " " << int(lb) << " " << int(lc) << " " << int(ld) << "  " <<  geometry << endl;
        cout << int(LA) << " " << int(LB) << " " << int(LC) << " " << int(LD) << "  " <<  geometry2 << endl;


        //exit gracefully?
        file.close();
        throw 1234;
    }

    //K4
    {
        //number of instructions in each step
        for (UI8 i=0; i<32; ++i)     file.read ((char*)&ninstrK4[i], sizeof (UI32));

        int totalK4 = 0;
        for (UI8 iblock=BOYS; iblock<NADA; ++iblock) totalK4 += ninstrK4[iblock];

        eseqK4 = new Op_K4[totalK4 + PF];

        //actual instructions
        for (UI32 i=0; i<totalK4; ++i) file.read ((char*)&eseqK4[i], sizeof (Op_K4));
    }

    //MIRROR
    {
        //number of instructions in each step
        for (UI8 i=0; i<32; ++i)     file.read ((char*)&ninstr[i], sizeof (UI32));

        int total = 0;
        for (UI8 iblock=KERNELS; iblock<=REORDER; ++iblock)  total += ninstr[iblock];

        eseq   = new Op_MIRROR[total + PF];

        for (UI32 i=0; i<total; ++i) file.read ((char*)&eseq[i], sizeof (Op_MIRROR));
    }

    file.read ((char*)&MaxMem, sizeof (UI32));
    file.read ((char*)&NFLOPS, sizeof (UI32));

    file.close();
}


