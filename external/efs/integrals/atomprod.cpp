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
#include "../basis/GTO.hpp"
#include "../basis/shellpair.hpp"
#include "../basis/SPprototype.hpp"
#include "../math/angular.hpp"
#include "../integrals/atomprod.hpp"
#include "../libechidna/libechidna.hpp"

using namespace std;
using namespace LibAngular;

AtomProdPrototypes::AtomProdPrototypes() {
    AP = NULL;
    APsame = NULL;
    APpool = NULL;
}

void AtomProdPrototypes::GenerateFrom (const r1tensor<rAtomPrototype*> & AtomTypes) {

     //inicia el array de atomos
    {
        const int n = AtomTypes.N();
        //copia la lista de tipos atomicos
        //ATypes = AtomTypes;

        nA = n;
        nAP = ((n+1)*n)/2;

        //inicia tensores
        APpool = new AtomProdPrototype[nAP+n];
        AP = new AtomProdPrototype*[n];

        int it = 0;
        for (int i=0; i<n; ++i) {
            it += i;
            AP[i] = APpool + it;
        }

        it += n;
        APsame = APpool + it;
    }

    // initializes all shell pair prototypes
    {
        int t=0;

        // atom pairs
        for (int i=0; i<AtomTypes.N(); ++i) {
            for (int j=0; j<=i; ++j) {
                t += AtomTypes[i]->basis.N() * AtomTypes[j]->basis.N();
            }
        }

        // same atom
        for (int i=0; i<AtomTypes.N(); ++i) {
            t += (AtomTypes[i]->basis.N()+1) * AtomTypes[i]->basis.N() / 2;
        }

        GPall = new ShellPairPrototype[t];
    }


    //calcula los productos para productos de atomos distintos
    int n=0;
    short int num=0;

    for (int i=0; i<AtomTypes.N(); ++i) {
        int id1 = AtomTypes[i]->DId; //AtomTypes[i];
        int nb1 = AtomTypes[i]->basis.N(); //BasisSet[id1].n;

        for (int j=0; j<=i; ++j) {
            int id2 = AtomTypes[j]->DId; //AtomTypes[j];
            int nb2 = AtomTypes[j]->basis.N(); //BasisSet[id2].n;

            AtomProdPrototype & APP = AP[i][j];
            APP.Atype1 = id1;
            APP.Atype2 = id2;
            APP.nGP1   = nb1;
            APP.nGP2   = nb2;
            APP.nGPt   = nb1*nb2;
            APP.num    = num;
            APP.same   = false;
            APP.Gs = GPall + n;
            APP.Rs = new float[APP.nGPt];

            int p=0;

            int offset1, offset2;
            offset1 = 0;

            for (int b1=0; b1<nb1; ++b1) {
                const GTO & g1 = AtomTypes[i]->basis[b1]; //BasisSet[id1][b1];

                offset2 = 0;

                for (int b2=0; b2<nb2; ++b2) {
                    const GTO & g2 = AtomTypes[j]->basis[b2]; //BasisSet[id2][b2];

                    APP.Gs[p].BasisProd(g1, g2);

                    if (!APP.Gs[p].inverted) {
                        APP.Gs[p].atype1 = id1; //AtomTypes[i];
                        APP.Gs[p].atype2 = id2; //AtomTypes[j];
                        APP.Gs[p].nb1    = b1;
                        APP.Gs[p].nb2    = b2;
                        APP.Gs[p].f1 = offset1;
                        APP.Gs[p].f2 = offset2;
                    }
                    else {
                        APP.Gs[p].atype1 = id2; //AtomTypes[j];
                        APP.Gs[p].atype2 = id1; //AtomTypes[i];
                        APP.Gs[p].nb1    = b2;
                        APP.Gs[p].nb2    = b1;
                        APP.Gs[p].f1 = offset2;
                        APP.Gs[p].f2 = offset1;
                    }

                    APP.Gs[p].sameatom = false;
                    APP.Rs[p] = 1./g1.k[0] + 1./g2.k[0];

                    ++p;

                    offset2 += g2.J * nmS[g2.l];
                }
                offset1 += g1.J * nmS[g1.l];
            }

            //int offset2 = 0;
            //for (int b2=0; b2<nb2; ++b2) offset2 += nmS[BasisSet[id2][b2].l];
            APP.wa1 = offset1;
            APP.wa2 = offset2;

            //genera y ordena las Rs de mayor a menor
            for (int p=0; p<APP.nGPt; ++p) {
                int    pmax = p;
                float  rmax = APP.Rs[p];

                for (int q=p+1; q<APP.nGPt; ++q) {
                    if (APP.Rs[q]>rmax) {
                        pmax = q;
                        rmax = APP.Rs[q];
                    }
                }

                //for some reason valgrind keeps bitchin about these ones
                if (p!=pmax) {
                    swap(APP.Rs[p], APP.Rs[pmax]);
                    swap(APP.Gs[p], APP.Gs[pmax]);
                }
            }

            //assigns the corresponding numbers in the tensor
            for (int p=0; p<APP.nGPt; ++p)
                APP.Gs[p].num = n + p;

            n += APP.nGPt;

            //sets highest l
            APP.maxL = 0;

            for (int b1=0; b1<nb1; ++b1) {
                const GTO & g1 = AtomTypes[i]->basis[b1]; //BasisSet[id1][b1];
                int l = (g1.l==LSP)?1:g1.l;
                if (l>APP.maxL) APP.maxL = l;
            }

            for (int b2=0; b2<nb2; ++b2) {
                const GTO & g2 = AtomTypes[j]->basis[b2]; //BasisSet[id2][b2];
                int l = (g2.l==LSP)?1:g2.l;
                if (l>APP.maxL) APP.maxL = l;
            }

            ++num;
        }
    }

    nGP = n;

    //calcula los productos para el mismo atomo
    for (int i=0; i<AtomTypes.N(); ++i) {
        int id1 = AtomTypes[i]->DId; //AtomTypes[i];
        int nb1 = AtomTypes[i]->basis.N(); //BasisSet[id1].n;

        //i==j
        {
            AtomProdPrototype & APP = APsame[i];

            APP.Atype1 = id1;
            APP.Atype2 = id1;
            APP.nGP1   = nb1;
            APP.nGP2   = nb1;
            APP.nGPt   = (nb1+1)*nb1/2;
            APP.num    = num;
            APP.same   = true;
            APP.Rs = new float[APP.nGPt];
            APP.Gs = GPall + n;
            int p=0;

            int offset1 = 0;
            for (int b1=0; b1<nb1; ++b1) {
                const GTO & g1 = AtomTypes[i]->basis[b1]; //BasisSet[id1][b1];

                int offset2 = 0;
                for (int b2=0; b2<=b1; ++b2) {
                    const GTO & g2 = AtomTypes[i]->basis[b2]; //BasisSet[id1][b2];

                    if (b1>b2) APP.Gs[p].BasisProd(g1, g2);
                    else       APP.Gs[p].BasisProd(g1);

                    if (!APP.Gs[p].inverted) {
                        APP.Gs[p].nb1    = b1;
                        APP.Gs[p].nb2    = b2;
                        APP.Gs[p].f1 = offset1;
                        APP.Gs[p].f2 = offset2;

                    }
                    else {
                        APP.Gs[p].nb1    = b2;
                        APP.Gs[p].nb2    = b1;
                        APP.Gs[p].f1 = offset2;
                        APP.Gs[p].f2 = offset1;
                    }

                    APP.Gs[p].atype1 = id1; //AtomTypes[i];
                    APP.Gs[p].atype2 = id1; //AtomTypes[i];
                    APP.Gs[p].sameatom = true;
                    APP.Rs[p] = 1./g1.k[0] + 1./g2.k[0];

                    ++p;

                    offset2 += g2.J * nmS[g2.l];
                }
                offset1 += g1.J * nmS[g1.l];
            }

            APP.wa1 = offset1;
            APP.wa2 = offset1;


            //genera y ordena las Rs de mayor a menor
            for (int p=0; p<APP.nGPt; ++p) {
                int    pmax = p;
                float  rmax = APP.Rs[p];

                for (int q=p+1; q<APP.nGPt; ++q) {
                    if (APP.Rs[q]>rmax) {
                        pmax = q;
                        rmax = APP.Rs[q];
                    }
                }
                if (p!=pmax) {
                    swap(APP.Rs[p], APP.Rs[pmax]);
                    swap(APP.Gs[p], APP.Gs[pmax]);
                }
            }


            //assigns the correspondent numbers in the tensor
            for (int p=0; p<APP.nGPt; ++p)
                APP.Gs[p].num = n + p;

            n += APP.nGPt;

            //sets highest l
            APP.maxL = 0;

            for (int b1=0; b1<nb1; ++b1) {
                const GTO & g1 = AtomTypes[i]->basis[b1]; //BasisSet[id1][b1];
                int l = (g1.l==LSP)?1:g1.l;
                if (l>APP.maxL) APP.maxL = l;
            }

            ++num;
        }
    }


    nGPsame = n - nGP;
}

AtomProdPrototypes::~AtomProdPrototypes() {
    delete[] APpool;
    delete[] AP;
}


APlists::APlists() {
    SPpool = NULL;
    WWpool = NULL;

    nGTOps = 0;
}

APlists::~APlists() {
    delete[] SPpool;
    delete[] WWpool;
}

static inline int Id12(int id1, int id2) {
    if (id1>=id2) return id1*(id1+1)/2 + id2;
    else          return id2*(id2+1)/2 + id1;
}

static void FlashSort(AtomProd * prods, int N) {

    //for short lists default to selection sort
    if (N<16) {
        for (int i=0; i<N; ++i)  {
            double best = prods[i].r2;
            int   k   = i;

            for (int j=i+1; j<N; ++j)
                if (prods[j].r2<best) {best = prods[j].r2; k = j;}

            swap(prods[i]    , prods[k]);
        }
        return;
    }

    double kmin, kmax;
    kmin = kmax = prods[0].r2;

    for (int i=0; i<N; ++i) {
        kmin = min(kmin, prods[i].r2);
        kmax = max(kmax, prods[i].r2);
    }

    double kmed = (kmin+kmax)*0.5;

    if ( kmed==kmin || kmax==kmed ) return;

    //count number of elements lower than the pivot
    int nl=0;
    for (int i=0; i<N; ++i)
        if (prods[i].r2 < kmed) ++nl;

    //count number of elements higher than the pivot in the lower part of the list
    int sw=0;
    for (int i=0; i<nl; ++i)
        if (prods[i].r2 >= kmed) ++sw;

    {
        int l=0;
        int h=nl;

        for (int nns=0; nns<sw; ++nns) {
            //find one value larger than the pivot
            while (prods[l].r2 < kmed) ++l;

            //find low in the higher part
            while (prods[h].r2 >= kmed) ++h;

            //swap
            swap(prods[l]    , prods[h]);
        }
    }


    //FlashSort both lists
    FlashSort(prods    ,   nl);
    FlashSort(prods+nl , N-nl);
}



//symtensor & Sp, symtensor & Tp, symtensor & Xp, symtensor & Yp, symtensor & Zp) {
void APlists::From(const r1tensor<rAtomPrototype*> & AtomTypes, const AtomProdPrototypes & APprototypes, const r2tensor<int> & IL, const rAtom * Atoms, int nAtoms, double logGDO) {

    UI16 nA = AtomTypes.N();
    UI32 nAP = nA*(nA+1)/2;

    //generates sorted lists of interacting atom pairs
    //************************************************
    {
        AtomSameBatch.set(nA);
        AtomPairBatch.set(nAP);

        UI32 * na  = new UI32[nA+nAP];
        UI32 * nap = na + nA;

        for (int a=0; a<nA;  ++a) na[a]  = 0;
        for (int a=0; a<nAP; ++a) nap[a] = 0;

        for (int at1=0; at1<IL.N(); ++at1) {
            int id = Atoms[at1].rAP->RId;
            ++na[id];
        }

        for (int at1=0; at1<IL.N(); ++at1) {
            int id1 = Atoms[at1].rAP->RId;

            //at1 > at2
            for (int j=0; j<IL.M(at1); ++j) {
                int at2 = IL[at1][j];
                if (at2>=at1) break;

                int id2  = Atoms[at2].rAP->RId;

                int n12 = Id12(id1, id2);

                ++nap[n12];
            }
        }

        //sets the prototypes for each list
        for (int id1=0; id1<nA; ++id1) {
            AtomProdPrototype * APP = &APprototypes.APsame[id1];
            AtomSameBatch[id1].APP = APP;
        }

        for (int id1=0; id1<nA; ++id1) {
            for (int id2=0; id2<=id1; ++id2) {
                int n12 = Id12(id1, id2);
                AtomProdPrototype * APP = &(APprototypes.AP[id1][id2]);
                AtomPairBatch[n12].APP = APP;
            }
        }

        for (int a=0; a<nA;  ++a) AtomSameBatch[a].APlist.set(na[a]);
        for (int a=0; a<nAP; ++a) AtomPairBatch[a].APlist.set(nap[a]);

        //populates the lists
        for (int a=0; a<nA;  ++a) na[a]  = 0;
        for (int a=0; a<nAP; ++a) nap[a] = 0;

        //same atom
        for (int at1=0; at1<nAtoms; ++at1) {
            int id = Atoms[at1].rAP->RId;

            AtomProd & AP = AtomSameBatch[id].APlist[na[id]];
            AtomProdPrototype * APP = &APprototypes.APsame[id];
            AP.APprototype       = APP;
            AP.reverse_prototype = false;
            AP.nAP               = APP->num;
            AP.r2                = 0;
            AP.norm              = 0;
            AP.at1               = at1;
            AP.at2               = at1;
            AP.A                 = Atoms[at1].c;
            AP.B                 = Atoms[at1].c;

            ++na[id];
        }

        //different atoms
        for (int at1=0; at1<IL.N(); ++at1) {
            int id1 = Atoms[at1].rAP->RId;

            //at1 > at2
            for (int j=0; j<IL.M(at1); ++j) {
                int at2 = IL[at1][j];
                if (at2>=at1) break;

                int id2  = Atoms[at2].rAP->RId;

                //distance, squared
                double r2; {
                    vector3 c = Atoms[at1].c - Atoms[at2].c;
                    r2 = c*c;
                }

                int n12 = Id12(id1, id2);
                AtomProdPrototype * APP;

                if (id1>=id2) {
                    APP = &(APprototypes.AP[id1][id2]);

                    AtomProd & AP = AtomPairBatch[n12].APlist[nap[n12]];
                    AP.APprototype       = APP;
                    AP.nAP               = APP->num;
                    AP.r2                = r2;
                    AP.norm              = sqrt(r2);
                    AP.at1               = at1;
                    AP.at2               = at2;
                    AP.A                 = Atoms[at1].c;
                    AP.B                 = Atoms[at2].c;
                    AP.reverse_prototype = false;
                    AP.MakeRotation();
                }
                else {
                    APP = &(APprototypes.AP[id2][id1]);

                    AtomProd & AP = AtomPairBatch[n12].APlist[nap[n12]];
                    AP.APprototype       = APP;
                    AP.nAP               = APP->num;
                    AP.r2                = r2;
                    AP.norm              = sqrt(r2);
                    AP.at1               = at2;
                    AP.at2               = at1;
                    AP.A                 = Atoms[at2].c;
                    AP.B                 = Atoms[at1].c;
                    AP.reverse_prototype = true;
                    AP.MakeRotation();
                }

                ++nap[n12];
            }
        }

        //sorts the atom pairs by increasing distance (could be done in O(N))
        // flashsort this!!!!!!
        for (int n12=0; n12<AtomPairBatch.n; ++n12) {
            Array<AtomProd> & APlist = AtomPairBatch[n12].APlist;

            FlashSort(APlist.keys, APlist.n);

            /*
            //quick and dirt, O(N^2) sorting algorithm; just for testing the idea
            for (int i=0; i<APlist.n; ++i) {
                AtomProd & APi = APlist[i];
                double rmin = APi.r2;
                int    jmin= i;

                for (int j=i+1; j<APlist.n; ++j) {
                    AtomProd & APj = APlist[j];
                    if (APj.r2<rmin) {rmin = APj.r2; jmin = j;}
                }
                if (i!=jmin) swap(APlist[i], APlist[jmin]);
            }
            */
        }

        delete[] na;
    }


    //calculates the memory required for the shell pairs and the expansions
    //*********************************************************************
    {
        nGTOps = 0;

        for (int id=0; id<nA; ++id) {
            AtomProdPrototype & APP = APprototypes.APsame[id];

            int maxGP = APP.nGPt;
            AtomSameBatch[id].nAPS.set(maxGP);
            AtomSameBatch[id].SP.set(maxGP);

            nGTOps += AtomSameBatch[id].APlist.n * APprototypes.APsame[id].nGPt;

            for (int i=0; i<AtomSameBatch[id].APlist.n; ++i)
                AtomSameBatch[id].APlist[i].nGTOps = APprototypes.APsame[id].nGPt;

            for (int i=0; i<maxGP; ++i) AtomSameBatch[id].nAPS[i] = AtomSameBatch[id].APlist.n;
        }



        for (int id=0; id<nAP; ++id) {
            AtomProdPrototype & APP = APprototypes.APpool[id];

            AtomPairBatch[id].nAPS.set(APP.nGPt);
            AtomPairBatch[id].SP.set(APP.nGPt);

            int maxGP = APP.nGPt;
            // for each atom pair, discard all shell pairs which dont survive GDO screening
            for (int i=0; i<AtomPairBatch[id].APlist.n; ++i) {

                while(maxGP>0 && AtomPairBatch[id].APlist[i].r2 > logGDO*APP.Rs[maxGP-1]) {
                    AtomPairBatch[id].nAPS[maxGP-1] = i;
                    --maxGP;
                }

                AtomPairBatch[id].APlist[i].nGTOps = maxGP;
                nGTOps += maxGP;
            }

            // for a given elemental shell pair, discard all atom pairs with large enough distance
            for (int i=0; i<APP.nGPt; ++i) {

                int jmax;
                for (jmax=0; jmax<AtomPairBatch[id].APlist.n; ++jmax) {
                    if (AtomPairBatch[id].APlist[jmax].r2 * APP.Gs[i].kr[0][0] > logGDO)
                        break;
                }

                AtomPairBatch[id].nAPS[i] = jmax;
            }
        }
    }


    //initializes and calculates the array of GTO products
    //****************************************************
    {
        ////////////////////////

        for (int id=0; id<nA; ++id) {
            AtomProdPrototype & APP = APprototypes.APsame[id];
            for (int i=0; i<AtomSameBatch[id].APlist.n; ++i) {
                AtomSameBatch[id].APlist[i].pSPs = new ShellPair*[AtomSameBatch[id].APlist[i].nGTOps];
            }
        }

        for (int id=0; id<nAP; ++id) {
            AtomProdPrototype & APP = APprototypes.APpool[id];
            for (int i=0; i<AtomPairBatch[id].APlist.n; ++i) {
                AtomPairBatch[id].APlist[i].pSPs = new ShellPair*[AtomPairBatch[id].APlist[i].nGTOps];
            }
        }



        UI64 Wmem = 0;

        for (int id=0; id<nA; ++id) {
            AtomProdPrototype & APP = APprototypes.APsame[id];
            for (int b12=0; b12<APP.nGPt; ++b12) {
                int s2 = APP.Gs[b12].Ka * APP.Gs[b12].Kb;
                Wmem += AtomSameBatch[id].APlist.n * s2;
            }
        }

        for (int id=0; id<nAP; ++id) {
            AtomProdPrototype & APP = APprototypes.APpool[id];
            for (int b12=0; b12<APP.nGPt; ++b12) {
                int s2 = APP.Gs[b12].Ka * APP.Gs[b12].Kb;
                Wmem += AtomPairBatch[id].nAPS[b12] * s2;
            }
        }


        //
        int nSP = 0;
        SPpool = new ShellPair[nGTOps];
        UI64 wSP = 0;
        WWpool = new double[Wmem];


        for (int id=0; id<nA; ++id) {
            AtomProdPrototype & APP = APprototypes.APsame[id];

            for (int b12=0; b12<APP.nGPt; ++b12) {
                int s2 = APP.Gs[b12].Ka * APP.Gs[b12].Kb;

                AtomSameBatch[id].SP[b12] = SPpool + nSP;

                for (int i=0; i<AtomSameBatch[id].APlist.n; ++i) {
                    AtomSameBatch[id].APlist[i].pSPs[b12] = &(SPpool[nSP]);
                    int at1 = AtomSameBatch[id].APlist[i].at1;

                    SPpool[nSP].ata = at1;
                    SPpool[nSP].atb = at1;
                    SPpool[nSP].Form (APP.Gs[b12], Atoms[at1].c, WWpool + wSP);

                    wSP += s2;
                    ++nSP;
                }
            }
        }

        for (int id=0; id<nAP; ++id) {
            AtomProdPrototype & APP = APprototypes.APpool[id];

            for (int b12=0; b12<APP.nGPt; ++b12) {
                int s2 = APP.Gs[b12].Ka * APP.Gs[b12].Kb;

                AtomPairBatch[id].SP[b12] = SPpool + nSP;

                for (int i=0; i<AtomPairBatch[id].nAPS[b12]; ++i) {
                    AtomPairBatch[id].APlist[i].pSPs[b12] = &(SPpool[nSP]);
                    int at1 = AtomPairBatch[id].APlist[i].at1;
                    int at2 = AtomPairBatch[id].APlist[i].at2;

                    if (!APP.Gs[b12].inverted) {
                        SPpool[nSP].ata = at1;
                        SPpool[nSP].atb = at2;
                    }
                    else {
                        SPpool[nSP].ata = at2;
                        SPpool[nSP].atb = at1;
                    }

                    SPpool[nSP].Form (APP.Gs[b12], Atoms[at1].c, Atoms[at2].c, WWpool + wSP, logGDO);

                    wSP += s2;
                    ++nSP;
                }
            }
        }

    }

}



