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

#include <algorithm>
#include "../basis/SPprototype.hpp"
#include "../math/angular.hpp"
#include "../basis/GTO.hpp"
using namespace LibAngular;


class primitivepair {
  public:
    double ik;
    UI8 ka;
    UI8 kb;
};

ShellPairPrototype::ShellPairPrototype() {
    Kt = 0;
    Ka = Kb = 0;
}

void ShellPairPrototype::BasisProd (const GTO & gA, const GTO & gB) {
    inverted = ((gA.l==LSP)?1:gA.l) < ((gB.l==LSP)?1:gB.l);
    //inverted = false;

    ll1 = ((gA.l==LSP)?1:gA.l) + ((gB.l==LSP)?1:gB.l);

    const GTO & A = inverted?gB:gA;
    const GTO & B = inverted?gA:gB;

    {
        l1 = A.l;
        l2 = B.l;

        SP1 = A.SP;
        SP2 = B.SP;

        Ka = A.K;
        Kb = B.K;
        Ja = A.J;
        Jb = B.J;

        for (int k=0; k<Ka; ++k) ka[k] = A.k[k];
        for (int j=0; j<Ja; ++j) for (int k=0; k<Ka; ++k) Na[j][k]  = A.N[j][k];
        if (SP1)
        for (int j=0; j<Ja; ++j) for (int k=0; k<Ka; ++k) Nap[j][k] = A.Np[j][k];

        for (int k=0; k<Kb; ++k) kb[k] = B.k[k];
        for (int j=0; j<Jb; ++j) for (int k=0; k<Kb; ++k) Nb[j][k]  = B.N[j][k];
        if (SP2)
        for (int j=0; j<Jb; ++j) for (int k=0; k<Kb; ++k) Nbp[j][k] = B.Np[j][k];
    }

    // *************************************************************************** //


    for (int p=0; p<Ka; ++p) {
        for (int q=0; q<Kb; ++q) {
            kr[p][q] = (ka[p]*kb[q])/(ka[p]+kb[q]);
        }
    }

    kmin = ka[0] + kb[0];

    wb1 = nmS[l1];
    wb2 = nmS[l2];

    Kt  = Ka * Kb;

    samef = false;


    // *************************************************************************** //

    //generate primitive sets
    {
        //generate list of primitives
        primitivepair pps[maxK2];

        int Kab=0;

        for (int p=0; p<Ka; ++p) {
            for (int q=0; q<Kb; ++q) {
                pps[Kab].ik = 1/kr[p][q];
                pps[Kab].ka  = p;
                pps[Kab].kb  = q;
                ++Kab;
            }
        }

        //sort them
        for (int i=0; i<Kab; ++i) {
            double kmax = pps[i].ik;
            int    imax = i;

            for (int j=i+1; j<Kab; ++j) {
                if (pps[j].ik > kmax) {
                    imax = j;
                    kmax = pps[j].ik;
                }
            }

            std::swap(pps[i], pps[imax]);
        }

        //fill primitive set
        UI8 nKa[maxK];
        UI8 nKb = 0;
        for (int p=0; p<maxK; ++p) nKa[p] = 0;

        //loop over all primitive sets
        for (int i=0; i<Kab; ++i) {
            //add a primitive
            UI8 ka = pps[i].ka;
            UI8 kb = pps[i].kb;

            if (nKa[kb] == 0) ++nKb;
            ++nKa[kb];

            Psets[i].nK2 = i+1;
            Psets[i].maxD2 = pps[i].ik;
            for (int p=0; p<maxK; ++p) Psets[i].nKa[p]   = nKa[p];
            Psets[i].nKb   = nKb;
        }

    }


    BCP.From(*this);
}

//same shell, same atom
void ShellPairPrototype::BasisProd (const GTO & A) {
    inverted = false;

    l1 = l2 = A.l;

    ll1 = 2*((A.l==LSP)?1:A.l);

    SP1 = SP2 = A.SP;

    Ka = Kb = A.K;
    Ja = Jb = A.J;

    // *************************************************************************** //

    for (int k=0; k<Ka; ++k) ka[k] = A.k[k];
    for (int j=0; j<Ja; ++j) for (int k=0; k<Ka; ++k) Na[j][k]  = A.N[j][k];
    if (SP1)
    for (int j=0; j<Ja; ++j) for (int k=0; k<Ka; ++k) Nap[j][k] = A.Np[j][k];

    for (int k=0; k<Kb; ++k) kb[k] = A.k[k];
    for (int j=0; j<Jb; ++j) for (int k=0; k<Kb; ++k) Nb[j][k]  = A.N[j][k];
    if (SP2)
    for (int j=0; j<Jb; ++j) for (int k=0; k<Kb; ++k) Nbp[j][k] = A.Np[j][k];


    for (int p=0; p<Ka; ++p) {
        for (int q=0; q<Kb; ++q) {
            kr[p][q] = (ka[p]*kb[q])/(ka[p]+kb[q]);
        }
    }

    kmin = ka[0] + kb[0];

    wb1 = nmS[l1];
    wb2 = nmS[l2];

    Kt = (Ka*Ka + Ka) / 2;

    samef = true;


    //generate primitive sets
    {
        //generate list of primitives
        primitivepair pps[maxK2];

        int Kab=0;

        for (int p=0; p<Ka; ++p) {
            for (int q=0; q<Kb; ++q) {
                pps[Kab].ik = 1/kr[p][q];
                pps[Kab].ka  = p;
                pps[Kab].kb  = q;
                ++Kab;
            }
        }

        //sort them
        for (int i=0; i<Kab; ++i) {
            double kmax = pps[i].ik;
            int    imax = i;

            for (int j=i+1; j<Kab; ++j) {
                if (pps[j].ik > kmax) {
                    imax = j;
                    kmax = pps[j].ik;
                }
            }

           std::swap(pps[i], pps[imax]);
        }

        //fill primitive set
        UI8 nKa[maxK];
        UI8 nKb = 0;
        for (int p=0; p<maxK; ++p) nKa[p] = 0;

        //loop over all primitive sets
        for (int i=0; i<Kab; ++i) {
            //add a primitive
            UI8 ka = pps[i].ka;
            UI8 kb = pps[i].kb;

            if (nKa[kb] == 0) ++nKb;
            ++nKa[kb];

            Psets[i].nK2 = i+1;
            Psets[i].maxD2 = pps[i].ik;
            for (int p=0; p<maxK; ++p) Psets[i].nKa[p]   = nKa[p];
            Psets[i].nKb   = nKb;
        }


        /*
        for (int i=0; i<Kab; ++i) {
            cout << Psets[i].nK2 << "  " << Psets[i].maxD2 << endl;
            cout << "  " << int(Psets[i].nKb) << "    ";
            for (int p=0; p<maxK; ++p) cout  << " " << int(Psets[i].nKa[p]);
            cout << endl;
        }
        char a;
        cin >> a;
        */


    }




    BCP.From(*this);
}

//this holds the bool UseGC for now
#include "../2eints/quimera.hpp"


BatchConstantsPair::BatchConstantsPair() {
}

BatchConstantsPair::~BatchConstantsPair() {
}

//data for applying CTEs during ERI evaluation phase
//for use later on
void BatchConstantsPair::From(const ShellPairPrototype & AB) {
    const bool shorten = true;
    const double NormS = SHList[0][0].T[0].cN;
    const double NormP = SHList[1][0].T[0].cN;

    UI8 La = AB.l1;
    UI8 Lb = AB.l2;

    UI8 Ka = AB.Ka;
    UI8 Kb = AB.Kb;

    UI8 Ja = AB.Ja;
    UI8 Jb = AB.Jb;

    for (UI8 ja=0; ja<Ja; ++ja) {
        for (UI8 ka=0; ka<Ka; ++ka) {
            Na [ka][ja] = AB.Na [ja][ka];

            //absorb angular normalization for L=s,l,p into coefficients
            if (shorten) {
                if (La==0 || La==LSP) Na [ka][ja] *= NormS;
                if (La==1)            Na [ka][ja] *= NormP;
                //if (La==LSP)          Nap[ka][ja] *= NormP;
            }
        }
    }

    for (UI8 jb=0; jb<Jb; ++jb) {
        for (UI8 kb=0; kb<Kb; ++kb) {
            Nb [kb][jb] = AB.Nb [jb][kb];

            //absorb angular normalization for L=s,l,p into coefficients
            if (shorten) {
                if (Lb==0 || Lb==LSP) Nb [kb][jb] *= NormS;
                if (Lb==1)            Nb [kb][jb] *= NormP;
                //if (Lb==LSP)          Nbp[kb][jb] *= NormP;
            }
        }
    }

    for (UI8 b=0; b<AB.Kb; ++b) {
        double kb = 2 * AB.kb[b];

        //double Kbn[LM8]; Kbn[0] = 1; for (int u=1; u<LM4; ++u) Kbn[u] = kb * Kbn[u-1];


        double Ns = 1;
        //double Np = 1;


        if (!UseGC) {
            Ns *= Nb [b][0];
            //Np *= AB.Nbp[0][b];
        }

         //exact length of the array is (lab + lcd)
        //for (int u=0; u<LM4; ++u) U[b][u] = Ns * Kbn[u];
    }


    for (UI8 b=0; b<AB.Kb; ++b) {
        for (UI8 a=0; a<AB.Ka; ++a) {
            double ke = 1 / (2 * (AB.ka[a] + AB.kb[b]));

            double kke = 1 / (sqrt(AB.ka[a] + AB.kb[b]) * (AB.ka[a] + AB.kb[b]));

            //double Ken[LM6]; Ken[0] = 1; for (int v=1; v<LM6; ++v) Ken[v] = ke * Ken[v-1];

            double Ns = 1;
            //double Np = 1;

            if (!UseGC) {
                Ns *= Na [a][0];
                //Np *= AB.Nap[0][a];
            }

            //exact length of the array is (2*lab + lcd)
            //for (int v=0; v<LM6; ++v) V[b][a][v] = Ns * Ken[v];


            K[b][a] = kke;
        }
    }

    for (UI8 b=0; b<AB.Kb; ++b) ikb[b] = 1/AB.kb[b];
    for (UI8 a=0; a<AB.Ka; ++a) ika[a] = 1/AB.ka[a];

    //internuclear distance fractions
    //*******************************
    for (UI8 a=0; a<AB.Ka; ++a) {
        for (UI8 b=0; b<AB.Kb; ++b) {
            r[b][a] = AB.kb[b] / (AB.ka[a] + AB.kb[b]);
        }
    }

    //in order to calculate the gaussian exponent and AERR
    //****************************************************
    for (UI8 a=0; a<AB.Ka; ++a) {
        for (UI8 b=0; b<AB.Kb; ++b) {
            k[b][a] = 0.5/(AB.ka[a] + AB.kb[b]);
        }
    }


    //for center-degenerate, long range integrals
    //*******************************************

    /*
    for (int i=0; i<LM6; ++i) UV[i] = 0;

    for (int b=0; b<AB.Kb; ++b) {
        double VV[LM6]; for (int i=0; i<LM6; ++i) VV[i] = 0;
        for (int a=0; a<AB.Ka; ++a) {
            for (int i=0; i<LM6; ++i) VV[i] += V[b][a][i];
        }
        for (int i=0; i<LM6; ++i) UV[i] += U[b][0] * VV[i];
    }
    */

    //for CDR in K1 and K3 loops
    //**************************

    for (UI8 b=0; b<AB.Kb; ++b) {
        double kb = 2 * AB.kb[b];

        b1[b] = kb;
        b2[b] = kb*kb;
    }

}

