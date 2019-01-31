#include "../math/newsparse.hpp"
#include "../math/tensors.hpp"
#include "../math/angular.hpp"
#include "../libechidna/libechidna.hpp"
using namespace LibAngular;

//THIS CLASS IS UBER CONFUSING FOR THE USER
//TAKE THE TIME TO ENCAPSULATE IT PROPERLY AND FIX THE CODE
//ALSO: ADD SPARSICITY
//*********************************************************

r2tensor<UI32> sparsetensorpattern::a2pos; //posicion del bloque para cada atomo

UI16 * sparsetensorpattern::alen; //length of the atom block
UI16 * sparsetensorpattern::nfs;  //number of functions per atom type

//second level, GC function blocks
r2tensor< r2tensor<UI32> > sparsetensorpattern::f2pos;
r2tensor<UI16> sparsetensorpattern::flen; //length of the function block
r2tensor<UI8> sparsetensorpattern::js;
r2tensor<UI8> sparsetensorpattern::ms;

//third level, individual function blocks
r2tensor< r2tensor<r2tensor<UI32> > > sparsetensorpattern::j2pos;


//one has to store which 'element' every atom is so to avoid wasting too much space
UI16 sparsetensorpattern::nids;
UI16 * sparsetensorpattern::ids;

UI32 sparsetensorpattern::natoms;
UI32 sparsetensorpattern::nbasis;
UI32 sparsetensorpattern::ntotal;
UI32 sparsetensorpattern::ntotal2;


//first level, atom blocks
UI32 * sparsetensorpattern::a1pos; //posiciones en el array lineal para cada atomo

//second level, GC function blocks
r2tensor<UI32> sparsetensorpattern::f1pos; //posiciones en el array lineal para cada atomo

//third level, individual function blocks
r2tensor<r2tensor<UI32> > sparsetensorpattern::j1pos;


void sparsetensorpattern::SetBlocks(const rAtom * Atoms, int nAtoms, const r1tensor<rAtomPrototype*> & AtomTypeList) {
    natoms = nAtoms;


    //count how many different 'elements' really are; set a correspondence between elements and a local identifier
    ids = new UI16[natoms];

    //map<UI16, UI16> mapIds;
    //UI16 nid = 0;
    for (UI32 at=0; at<natoms; ++at) {
        UI16 id = Atoms[at].rAP->DId;
        ids[at] = id;
        /*
        if (mapIds.count(id)==0) {
            mapIds[id] = nid;
            nid++;
        }
        */
    }
    nids = AtomTypeList.N(); //mapIds.size();

    //hace el mapa inverso

    /*
    UI16 * types = new UI16[nids];
    for (map<UI16, UI16>::iterator it=mapIds.begin(); it!=mapIds.end(); ++it) {
        types[it->second] = it->first;
    }
    */


    //ATOM BLOCK LEVEL
    //****************

    //compute the width/length of the element block
    alen  = new UI16 [nids];

    for (UI16 id=0; id<nids; ++id) {
        UI16 w = 0;
        //UI16 type = types[id];

        //sum all basis functions
        for (UI8 b=0; b<AtomTypeList[id]->basis.N(); ++b) {
            UI8 l = AtomTypeList[id]->basis[b].l;
            w += (2*l+1) * AtomTypeList[id]->basis[b].J;
        }
        /*
        for (UI8 b=0; b<BasisSet[type].n; ++b) {
            UI8 l = BasisSet[type][b].l;
            w += (2*l+1) * BasisSet[type][b].J;
        }
        */
        alen[id] = w;
    }


    //offsets from the pointer to the linear array
    a2pos.set(natoms, natoms);
    a1pos = new UI32[natoms];

    UI32 Apos2 = 0; //total accumulated offset
    UI32 Apos1 = 0; //total accumulated offset

    for (UI32 ati=0; ati<natoms; ++ati) {
        a1pos[ati] = Apos1;
        UI16 id1   = Atoms[ati].rAP->RId;  //mapIds[ids[ati]];
        Apos1 += alen[id1];

        for (UI32 atj=0; atj<natoms; ++atj) {
            a2pos[ati][atj] = Apos2;
            UI16 id2   = Atoms[atj].rAP->RId; //mapIds[ids[atj]];
            Apos2 += alen[id1]*alen[id2];
        }
    }


    //GC FUNCTION BLOCK LEVEL
    //***********************

    //MemAllocator palloc;

    //first initialize the function block lengths' array
    nfs  = new UI16 [nids];
    for (UI32 id=0; id<nids; ++id) nfs[id] = AtomTypeList[id]->basis.N(); //BasisSet[types[id]].n;

    flen.set(nfs, nids);
    js.set(nfs, nids);
    ms.set(nfs, nids);

    for (UI16 id=0; id<nids; ++id) {

        //sum all basis functions
        for (UI8 b=0; b<AtomTypeList[id]->basis.N(); ++b) {
            UI8 l =AtomTypeList[id]->basis[b].l;
            flen[id][b] = nmS[l] * AtomTypeList[id]->basis[b].J;
            js[id][b] = AtomTypeList[id]->basis[b].J;
            ms[id][b] = nmS[l];
        }

        /*
        UI16 type = types[id];

        //sum all basis functions
        for (UI8 b=0; b<BasisSet[type].n; ++b) {
            UI8 l = BasisSet[type][b].l;
            flen[id][b] = nmS[l] * BasisSet[type][b].J;
            js[id][b] = BasisSet[type][b].J;
            ms[id][b] = nmS[l];
        }
        */
    }

    f2pos.set(nids, nids);

    //for each element pair combination
    for (UI16 id1=0; id1<nids; ++id1) {
        for (UI16 id2=0; id2<nids; ++id2) {
            f2pos[id1][id2].set(nfs[id1], nfs[id2]);

            UI32 ww = 0;

            //loop over both atom's function blocks
            for (int b1=0; b1<nfs[id1]; ++b1) {
                for (int b2=0; b2<nfs[id2]; ++b2) {
                    f2pos[id1][id2][b1][b2] = ww;
                    ww += flen[id1][b1]*flen[id2][b2];
                }
            }
        }
    }


    f1pos.set(nfs, nids);
    for (UI16 id1=0; id1<nids; ++id1) {
        UI32 ww = 0;
        for (int b1=0; b1<nfs[id1]; ++b1) {
            f1pos[id1][b1] = ww;
            ww += flen[id1][b1];
        }
    }



    //change the identifiers  for the inner definition
    for (UI32 at=0; at<natoms; ++at)
        ids[at] = Atoms[at].rAP->RId; //mapIds[ids[at]];


    //calcula el numero final de filas/columnas y valores que se deben almacenar
    ntotal = 0;
    for (UI32 at=0; at<natoms; ++at)
        ntotal += alen[ids[at]];

    ntotal2 = ntotal*ntotal;

    //delete[] types;
}

UI32 sparsetensorpattern::GetOffset(UI32 at1, UI32 at2) {
    return a2pos(at1, at2);
}

UI32 sparsetensorpattern::GetOffset(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2) {
    UI16 id1 = ids[at1];
    UI16 id2 = ids[at2];

    return f2pos(id1,id2)(nf1,nf2);
}


UI32 sparsetensorpattern::GetOffset(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2, UI16 nj1, UI16 nj2) {
    UI16 id1 = ids[at1];
    UI16 id2 = ids[at2];

    return f2pos(id1,id2)(nf1,nf2);
}





static inline UI32 Tpos (UI32 M1, UI32 M2, UI32 v1, UI32 v2) {
    return v2*M1 + v1;
}

Sparse & Sparse::operator=(const tensor2 & rhs) {
    if (values==NULL) Set();


    //foreach atom
    for (UI32 ati=0; ati<natoms; ++ati) {
        for (UI32 atj=0; atj<natoms; ++atj) {

            UI16 id1   = ids[ati];
            UI16 id2   = ids[atj];

            double * offset1 = values + a2pos[ati][atj];

            UI32 p1 = a1pos[ati];
            UI32 q1 = a1pos[atj];

            for (int b1=0; b1<nfs[id1]; ++b1) {
                for (int b2=0; b2<nfs[id2]; ++b2) {
                    double * offset2 = offset1 + f2pos[id1][id2][b1][b2];

                    UI32 p2 = p1 + f1pos[id1][b1];
                    UI32 q2 = q1 + f1pos[id2][b2];

                    UI32 mj1 = js[id1][b1];
                    UI32 mj2 = js[id2][b2];
                    UI32 mm1 = ms[id1][b1];
                    UI32 mm2 = ms[id2][b2];

                    for (UI8 j1=0; j1<mj1; ++j1) {
                        for (UI8 j2=0; j2<mj2; ++j2) {
                            double * offset3 = offset2  +  (j1*mj2+j2)*(mm1*mm2);

                            UI32 p3 = p2 + j1*mm1;
                            UI32 q3 = q2 + j2*mm2;


                            for (UI8 m1=0; m1<mm1; ++m1) {
                                for (UI8 m2=0; m2<mm2; ++m2) {
                                    UI32 Ap = Tpos(mm1,mm2, m1,m2);

                                    UI32 p4 = p3 + m1;
                                    UI32 q4 = q3 + m2;

                                    (*this)(ati,atj, b1,b2, j1,j2)[Ap] = rhs(p4, q4);
                                }
                            }

                        }
                    }

                }
            }


        }
    }

    return *this;
}

Sparse & Sparse::operator=(const Sparse & rhs) {
    if (values==NULL) Set();
    for (int i=0; i<ntotal2; ++i) values[i] = rhs.values[i];
    return *this;
}

Sparse & Sparse::operator+=(const Sparse & rhs) {
    if (values==NULL) {Set(); zeroize();}
    for (int i=0; i<ntotal2; ++i) values[i] += rhs.values[i];
    return *this;
}

Sparse & Sparse::operator-=(const Sparse & rhs) {
    if (values==NULL) {Set(); zeroize();}
    for (int i=0; i<ntotal2; ++i) values[i] -= rhs.values[i];
    return *this;
}

Sparse & Sparse::operator*=(double rhs) {
    if (values==NULL) {Set(); zeroize();}
    for (int i=0; i<ntotal2; ++i) values[i] *= rhs;
    return *this;
}

Sparse & Sparse::operator/=(double rhs) {
    if (values==NULL) {Set(); zeroize();}
    double irhs = 1./rhs;
    for (int i=0; i<ntotal2; ++i) values[i] *= irhs;
    return *this;
}



