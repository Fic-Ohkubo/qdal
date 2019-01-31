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


#include <fstream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <cstdlib>
#include <string>
#include <queue>

#include "IICgen.hpp"
#include "../math/angular.hpp"
#include "../low/chrono.hpp"
#include "../low/cache.hpp"

using namespace std;
using namespace LibAngular;


/* *********************************** */

EvaluationScheme::EvaluationScheme() {
    useGC  = false;

    NK4 = 0;
    sequenceK4 = NULL;
    nblocksK4 = 0;

    for (int i=0; i<32; ++i) ninstrK4[i] = 0;
    for (int i=0; i<32; ++i) seqsK4[i]   = NULL;


    RRSK4[nblocksK4] = BOYS;  ++nblocksK4;

    RRSK4[nblocksK4] = ADRR4; ++nblocksK4;
    RRSK4[nblocksK4] = AERR4; ++nblocksK4;
    RRSK4[nblocksK4] = CDR4 ; ++nblocksK4;
    RRSK4[nblocksK4] = K4D  ; ++nblocksK4;
    RRSK4[nblocksK4] = K4E  ; ++nblocksK4;
    RRSK4[nblocksK4] = K4F  ; ++nblocksK4;

    RRSK4[nblocksK4] = ADRR3; ++nblocksK4;
    RRSK4[nblocksK4] = AERR3; ++nblocksK4;
    RRSK4[nblocksK4] = CDR3 ; ++nblocksK4;
    RRSK4[nblocksK4] = K3D  ; ++nblocksK4;
    RRSK4[nblocksK4] = K3E  ; ++nblocksK4;
    RRSK4[nblocksK4] = K3F  ; ++nblocksK4;

    RRSK4[nblocksK4] = ADRR2; ++nblocksK4;
    RRSK4[nblocksK4] = AERR2; ++nblocksK4;
    RRSK4[nblocksK4] = CDR2 ; ++nblocksK4;
    RRSK4[nblocksK4] = K2D  ; ++nblocksK4;
    RRSK4[nblocksK4] = K2E  ; ++nblocksK4;
    RRSK4[nblocksK4] = K2F  ; ++nblocksK4;

    RRSK4[nblocksK4] = ADRR1; ++nblocksK4;
    RRSK4[nblocksK4] = AERR1; ++nblocksK4;
    RRSK4[nblocksK4] = CDR1 ; ++nblocksK4;
    RRSK4[nblocksK4] = K1D  ; ++nblocksK4;
    RRSK4[nblocksK4] = K1E  ; ++nblocksK4;
    RRSK4[nblocksK4] = K1F  ; ++nblocksK4;

    RRSK4[nblocksK4] = ADRR0; ++nblocksK4;
    RRSK4[nblocksK4] = AERR0; ++nblocksK4;
    RRSK4[nblocksK4] = CDR0 ; ++nblocksK4;


    RDep = NULL;
    N   = 0;

    sequence = NULL;
    nblocks = 0;

    for (int i=0; i<32; ++i) ninstr[i] = 0;
    for (int i=0; i<32; ++i) seqs[i]   = NULL;


    RRS[nblocks] = KERNELS; ++nblocks;
    RRS[nblocks] = MMDZ;    ++nblocks;
    RRS[nblocks] = CTEBZ;   ++nblocks;
    RRS[nblocks] = CTEKZ;   ++nblocks;
    RRS[nblocks] = MMDY;    ++nblocks;
    RRS[nblocks] = MMDX;    ++nblocks;
    RRS[nblocks] = CTEKY;   ++nblocks;
    RRS[nblocks] = CTEKX;   ++nblocks;
    RRS[nblocks] = CTEBY;   ++nblocks;
    RRS[nblocks] = CTEBX;   ++nblocks;
    RRS[nblocks] = HRRBZ;   ++nblocks;
    RRS[nblocks] = HRRBY;   ++nblocks;
    RRS[nblocks] = HRRBX;   ++nblocks;
    RRS[nblocks] = SPHA;    ++nblocks;
    RRS[nblocks] = SPHB;    ++nblocks;
    RRS[nblocks] = HRRKZ;   ++nblocks;
    RRS[nblocks] = HRRKY;   ++nblocks;
    RRS[nblocks] = HRRKX;   ++nblocks;
    RRS[nblocks] = SPHC;    ++nblocks;
    RRS[nblocks] = SPHD;    ++nblocks;
    RRS[nblocks] = REORDER; ++nblocks;
}

EvaluationScheme::~EvaluationScheme() {
    LastDep = NULL;
    delete[] RDep;
    delete[] sequence;
    EvaluationOrder.clear();
    Info.clear();
    IntIdMap.clear();
    IdIntMap.clear();

    K4EvaluationOrder.clear();
    K4Info.clear();
    K4IdMap.clear();
    IdK4Map.clear();
}

UI32 EvaluationScheme::GetId(const integralCart & i) const {
    std::map <integralCart, UI32>::const_iterator it;
    it = IntIdMap.find(i);
    if (it==IntIdMap.end()) return NOTIN;
    else return it->second;
}

void EvaluationScheme::Append(const IntVar & IV) {
    EVList.push_front(IV);
}

void EvaluationScheme::AppendK4(const K4Var & KV) {
    EVListK4.push_front(KV);
}

int EvaluationScheme::CalcLife(list<IntVar> * glist, set<integralCart> & setf) {
    map<integralCart, UI32> Life;

    list<IntVar>::iterator it;
    for (it=glist->begin(); it!=glist->end(); ++it) {
        Life[it->I] = 0;
    }
    integralCart null;
    Life[null] = 0;

    //set lifelihood of variables
    UI32 pos = 0;
    for (it=glist->begin(); it!=glist->end(); ++it) {
        ++pos;
        Life[it->V.VD1] = max(Life[it->V.VD1], pos);
        Life[it->V.VD2] = max(Life[it->V.VD2], pos);
        Life[it->V.VD3] = max(Life[it->V.VD3], pos);
        Life[it->V.VD4] = max(Life[it->V.VD4], pos);
        Life[it->V.VD5] = max(Life[it->V.VD5], pos);
        Life[it->V.VD6] = max(Life[it->V.VD6], pos);
    }

    //the final set, of course is alive after the evaluation
    ++pos;
    set<integralCart>::iterator sit;
    for(sit=setf.begin(); sit!=setf.end(); ++sit) {
        Life[*sit] = pos;
    }

    //compute minimum graph diameter
    int nn = glist->size();
    char * An = new char[nn+2];
    for (int t=1; t<=nn+1; ++t) An[t] = 0;

    map<integralCart, UI32>::iterator iit;
    for(iit=Life.begin(); iit!=Life.end(); ++iit) {
        ++An[iit->second];
    }

    int maxalive = 0;
    int totvars = 1;

    for (int t=2; t<=nn; ++t) {
        totvars -= An[t-1]; //variable muertas tras la ecuacion anterior
        ++totvars; //nueva variable declarada
        maxalive = max(maxalive, totvars);
    }
    delete[] An;
    Life.clear();

    return maxalive;
}

void EvaluationScheme::AddGroup(map<integralCart,VarDep> (&imap)[1+4*LMAX+4], set<integralCart> & setf) {
    //find AN (difference of variables) and AV (maximum number of intermediate variables)

    //compute life
    list<IntVar> * glist = NULL;
    int maxlife1 = 1024*1024*1024;

    //expensive and almost useless
    if (0) {
        glist = new list<IntVar>;

        map<integralCart,VarDep> gmap;
        map<integralCart, UI32> nrefs;
        set<integralCart> unassigned;
        set<integralCart> assigned;


        for (int v=4*LMAX+4; v>=0; --v) {
            map<integralCart,VarDep>::iterator it4;
            for (it4=imap[v].begin(); it4!=imap[v].end(); ++it4) {
                gmap[it4->first] = it4->second;
            }
        }

        //count how many times some variable is referenced
        const integralCart null;
        assigned.insert(null); //insert a null integral, so it will appear as if null reference point  top an existing, assigned variable
        nrefs[null] = 2*1024*1024*1024; //enough not to be the deciding option


        for (int v=4*LMAX+4; v>=0; --v)  {
            map<integralCart,VarDep>::iterator it4;
            for (it4=imap[v].begin(); it4!=imap[v].end(); ++it4) {
                nrefs[it4->first] = 0; //initialize number of references
                if (v>0) unassigned.insert(it4->first);
                else     assigned.insert(it4->first);
            }
        }

        //increment number of dependencies
        for (int v=4*LMAX+4; v>0; --v)  {
            map<integralCart,VarDep>::iterator it4;
            for (it4=imap[v].begin(); it4!=imap[v].end(); ++it4) {
                ++nrefs[it4->second.VD1];
                ++nrefs[it4->second.VD2];
                ++nrefs[it4->second.VD3];
                ++nrefs[it4->second.VD4];
                ++nrefs[it4->second.VD5];
                ++nrefs[it4->second.VD6];
            }
        }

        //add the assigned variables to the list
        set<integralCart>::iterator sit;

        for (sit = assigned.begin(); sit != assigned.end(); ++sit) {
            if (*sit != null) {
                IntVar IV;
                IV.I = *sit;
                IV.V = gmap[IV.I];
                glist->push_back(IV);
            }
        }

        //assign variables and add them to the list
        while (unassigned.size() > 0) {
            set<integralCart>::iterator it;

            integralCart next = null;
            unsigned short int  ncancel = 0;
            unsigned short int  nforw   = 1024*1024*1024;
            unsigned short int  nbforw  = 1024*1024*1024;

            for (it = unassigned.begin(); it!=unassigned.end(); ++it) {
                //check dependencies are already available
                integralCart IC = *it;
                integralCart V1 = gmap[IC].VD1;
                integralCart V2 = gmap[IC].VD2;
                integralCart V3 = gmap[IC].VD3;
                integralCart V4 = gmap[IC].VD4;
                integralCart V5 = gmap[IC].VD5;
                integralCart V6 = gmap[IC].VD6;

                //if is valid
                if ((assigned.count(V1) == 1) && (assigned.count(V2) == 1) && (assigned.count(V3) == 1) && (assigned.count(V4) == 1)  && (assigned.count(V5) == 1)  && (assigned.count(V6) == 1)) {
                    //order of evaluation:
                    // FIRST:   cancel as many variables as possible
                    // SECOND:  make it so the variable will be cancelled soon (count self forward dependencies)
                    // THIRD:   leave as little forward dependencies as possible

                    unsigned short int c=0;
                    if (nrefs[V1]==1) ++c;
                    if (nrefs[V2]==1) ++c;
                    if (nrefs[V3]==1) ++c;
                    if (nrefs[V4]==1) ++c;
                    if (nrefs[V5]==1) ++c;
                    if (nrefs[V6]==1) ++c;

                    unsigned short int  f = min(min(min(nrefs[V1],nrefs[V2]), min(nrefs[V3],nrefs[V4])), min(nrefs[V5],nrefs[V6]));
                    unsigned short int  b = nrefs[IC];

                    //compare
                    bool isbetter = false;
                    if (c>ncancel) isbetter = true;
                    else if (c == ncancel) {
                        if (f<nforw) isbetter = true;
                        else if (f==nforw) {
                            if (b<nbforw) isbetter = true;
                        }
                    }

                    if (isbetter || next==null) {
                        next = IC;
                        ncancel = c;
                        nforw   = f;
                        nbforw  = b;
                    }
                }
            }


            //decrement the number of references
            integralCart V1 = gmap[next].VD1;
            integralCart V2 = gmap[next].VD2;
            integralCart V3 = gmap[next].VD3;
            integralCart V4 = gmap[next].VD4;
            integralCart V5 = gmap[next].VD5;
            integralCart V6 = gmap[next].VD6;
            --nrefs[V1];
            --nrefs[V2];
            --nrefs[V3];
            --nrefs[V4];
            --nrefs[V5];
            --nrefs[V6];

            //move from unassigned to assigned
            unassigned.erase(next);
            assigned.insert(next);

            //push the integral in teh evaluation list
            IntVar IV;
            IV.I = next;
            IV.V = gmap[IV.I];
            glist->push_back(IV);
        }

        nrefs.clear();
        gmap.clear();
        unassigned.clear();
        assigned.clear();

        int maxlife1  = CalcLife(glist, setf);
    }


    list<IntVar> * nlist = NULL;
    int maxlife2 = 1024*1024*1024;

    if (1) {
        nlist = new list<IntVar>;

        for (int v=4*LMAX+4; v>0; --v) {
            map<integralCart,VarDep>::iterator it4;
            for (it4=imap[v].begin(); it4!=imap[v].end(); ++it4) {
                IntVar IV;
                IV.I = it4->first;
                IV.V = it4->second;
                nlist->push_front(IV);
            }
        }

        /*
        {
            map<integralCart,VarDep>::iterator it4;
            for (it4=imap[0].begin(); it4!=imap[0].end(); ++it4) {
                IntVar IV;
                IV.I = it4->first;
                IV.V = it4->second;
                IV.V.ETYPE = OTHER;
                nlist->push_front(IV);
            }
        }
        */
        int maxlife2  = CalcLife(nlist, setf);
    }

    int increment =  int(setf.size()) - int(imap[0].size());


    list<IntVar> * flist;

    if (maxlife1<maxlife2) {
        flist = glist;
        delete nlist;
        cout << "Variable number increment:  " << increment << endl;
        cout << "Alive vars 1:               " << maxlife1 << endl;
        cout << "Alive vars 2:               " << maxlife2 << endl;
    }
    else {
        flist = nlist;
        delete glist;
    }

    Var2 NV;
    NV.AN = increment;
    NV.AV = maxlife1;

    groups.insert ( pair<Var2, list<IntVar>* > (NV, flist) );
}

void EvaluationScheme::AddGroups() {
    //sort groups and append to the list
    multimap<Var2, list<IntVar>* >::reverse_iterator rit;

    for (rit = groups.rbegin(); rit != groups.rend(); ++rit) {
        list<IntVar> * group = rit->second;
        list<IntVar>::reverse_iterator rit2;

        for (rit2 = group->rbegin(); rit2 != group->rend(); ++rit2) {
            Append(*rit2);
        }
        delete rit->second;
    }
    groups.clear();
}


void EvaluationScheme::LinkVars(const set<integralCart> & set0) {
    std::list<IntVar>::iterator itiv;
    UI32 Id;


    //generate list of RR blocks before any var is eliminated
    //*******************************************************
    Id = 0;
    for (itiv = EVList.begin(); itiv!=EVList.end();) {
        //delete duplicated integrals from list &
        if (IntIdMap.count(itiv->I)==1) {
            itiv = EVList.erase(itiv);
        }
        //assign a unique identifier to every integral
        else {
            ++Id;
            IdIntMap[Id] = itiv->I;
            IntIdMap[itiv->I] = Id;
            ++itiv;
        }
    }

    int pos = 0;

    //assign variables for each integral
    for (itiv = EVList.begin(); itiv!=EVList.end();) {
        ++pos;

        UI32 Id = IntIdMap[itiv->I];

        VarInfo & info = Info[Id];

        info.position = pos;
        info.final    = false;

        if (itiv->V.ETYPE==KERNELS) {
            ikernel K;
            K = itiv->V.VD1;
            info.aux = (K4IdMap.count(K) == 1)?K4IdMap[K]:NOTIN;
        }
        else {
            info.VD1 = (IntIdMap.count(itiv->V.VD1) == 1)?IntIdMap[itiv->V.VD1]:NOTIN;
            info.VD2 = (IntIdMap.count(itiv->V.VD2) == 1)?IntIdMap[itiv->V.VD2]:NOTIN;
            info.VD3 = (IntIdMap.count(itiv->V.VD3) == 1)?IntIdMap[itiv->V.VD3]:NOTIN;
            info.VD4 = (IntIdMap.count(itiv->V.VD4) == 1)?IntIdMap[itiv->V.VD4]:NOTIN;
            info.VD5 = (IntIdMap.count(itiv->V.VD5) == 1)?IntIdMap[itiv->V.VD5]:NOTIN;
            info.VD6 = (IntIdMap.count(itiv->V.VD6) == 1)?IntIdMap[itiv->V.VD6]:NOTIN;
            info.aux = itiv->V.aux;
        }


        if (info.VD1 != NOTIN) ++info.nBDep;
        if (info.VD2 != NOTIN) ++info.nBDep;
        if (info.VD3 != NOTIN) ++info.nBDep;
        if (info.VD4 != NOTIN) ++info.nBDep;
        if (info.VD5 != NOTIN) ++info.nBDep;
        if (info.VD6 != NOTIN) ++info.nBDep;

        if (info.nBDep==0) info.zero = true;


        info.ETYPE = itiv->V.ETYPE;

        EvaluationOrder.push_back(Id);

        itiv = EVList.erase(itiv);
    }

    N = IntIdMap.size();


    //marks output integrals
    std::set <integralCart>::const_iterator it3;
    for (it3=set0.begin(); it3!=set0.end(); ++it3) {
        UI32 Id = (IntIdMap.count(*it3) == 1)?IntIdMap[*it3]:NOTIN;

        if (Id!=-1) {
            Info[Id].final = true;
            Info[Id].lifelihood = N+1; //must be alive at the end
        }
    }
}

void EvaluationScheme::LinkVarsK4(const set<ikernel> & setk) {
    std::list<K4Var>::iterator itiv;
    UI32 Id = 0;

    //generate list of RR blocks before any var is eliminated
    //*******************************************************


    for (itiv = EVListK4.begin(); itiv!=EVListK4.end();) {

        //delete duplicated integrals from list &
        if (K4IdMap.count(itiv->K)==1) {
            itiv = EVListK4.erase(itiv);
        }
        //assign a unique identifier to every integral
        else {
            ++Id;
            IdK4Map[Id] = itiv->K;
            K4IdMap[itiv->K] = Id;
            ++itiv;
        }
    }

    UI32 pos = 0;

    //assign variables for each integral
    for (itiv = EVListK4.begin(); itiv!=EVListK4.end();) {
        ++pos;

        UI32 Id = K4IdMap[itiv->K];

        VarInfoK4 & info = K4Info[Id];

        info.position = pos;
        info.final    = false;

        info.VD1 = (K4IdMap.count(itiv->V.VD1) == 1)?K4IdMap[itiv->V.VD1]:NOTIN;
        info.VD2 = (K4IdMap.count(itiv->V.VD2) == 1)?K4IdMap[itiv->V.VD2]:NOTIN;
        info.VD3 = (K4IdMap.count(itiv->V.VD3) == 1)?K4IdMap[itiv->V.VD3]:NOTIN;
        info.VD4 = (K4IdMap.count(itiv->V.VD4) == 1)?K4IdMap[itiv->V.VD4]:NOTIN;
        info.VD5 = (K4IdMap.count(itiv->V.VD5) == 1)?K4IdMap[itiv->V.VD5]:NOTIN;
        info.VD6 = (K4IdMap.count(itiv->V.VD6) == 1)?K4IdMap[itiv->V.VD6]:NOTIN;
        info.VDE = (K4IdMap.count(itiv->V.VDE) == 1)?K4IdMap[itiv->V.VDE]:NOTIN;

        if (info.VD1 != NOTIN) ++info.nBDep;
        if (info.VD2 != NOTIN) ++info.nBDep;
        if (info.VD3 != NOTIN) ++info.nBDep;
        if (info.VD4 != NOTIN) ++info.nBDep;
        if (info.VD5 != NOTIN) ++info.nBDep;
        if (info.VD6 != NOTIN) ++info.nBDep;
        if (info.VDE != NOTIN) ++info.nBDep;

        if (info.nBDep==0) info.zero = true;

        info.aux   = itiv->V.aux;
        info.ETYPE = itiv->V.ETYPE;

        K4EvaluationOrder.push_back(Id);

        itiv = EVListK4.erase(itiv);
    }

    NK4 = K4IdMap.size();

    //marks output integrals
    std::set <ikernel>::const_iterator it3;
    for (it3=setk.begin(); it3!=setk.end(); ++it3) {
        UI32 Id = (K4IdMap.count(*it3) == 1)?K4IdMap[*it3]:NOTIN;

        if (Id!=-1) {
            K4Info[Id].final = true;
            K4Info[Id].lifelihood = NK4+1; //must be alive at the end
        }
    }
}


void EvaluationScheme::LinkBack() {
    UI32 ndeps = 0;

    map <UI32, VarInfo>::iterator it;

    //assign variables for each integral
    for (it = Info.begin(); it != Info.end(); ++it) {
        it->second.nDep = 0;
    }

    for (it = Info.begin(); it != Info.end(); ++it) {
        //if (it->second.ETYPE==KERNELS) continue;
        UI32 R1 = it->second.VD1;
        UI32 R2 = it->second.VD2;
        UI32 R3 = it->second.VD3;
        UI32 R4 = it->second.VD4;
        UI32 R5 = it->second.VD5;
        UI32 R6 = it->second.VD6;
        if (R1!=NOTIN) {++Info[R1].nDep; ++ndeps;}
        if (R2!=NOTIN) {++Info[R2].nDep; ++ndeps;}
        if (R3!=NOTIN) {++Info[R3].nDep; ++ndeps;}
        if (R4!=NOTIN) {++Info[R4].nDep; ++ndeps;}
        if (R5!=NOTIN) {++Info[R5].nDep; ++ndeps;}
        if (R6!=NOTIN) {++Info[R6].nDep; ++ndeps;}
    }

    NEVS = ndeps; //use extra variables (for later simplifying)

    //create list of integrals in order to link back
    RDep = new UI32[ndeps+NEVS];

    //assign pointers
    UI32 * piC = RDep;

    for (it = Info.begin(); it != Info.end(); ++it) {
        //if (it->second.ETYPE==KERNELS) continue;
        it->second.iDep = piC;
        piC += it->second.nDep;
        it->second.nDep = 0;
    }

    LastDep = piC; //pointer to first free var


    //copies integral to where it belongs in the list and increments counter
    for (it = Info.begin(); it != Info.end(); ++it) {
        //if (it->second.ETYPE==KERNELS) continue;

        UI32 R1 = it->second.VD1;
        UI32 R2 = it->second.VD2;
        UI32 R3 = it->second.VD3;
        UI32 R4 = it->second.VD4;
        UI32 R5 = it->second.VD5;
        UI32 R6 = it->second.VD6;

        if (R1!=NOTIN) {
            int n = Info[R1].nDep;
            Info[R1].iDep[n] = it->first;
            ++Info[R1].nDep;
        }

        if (R2!=NOTIN) {
            int n = Info[R2].nDep;
            Info[R2].iDep[n] = it->first;
            ++Info[R2].nDep;
        }

        if (R3!=NOTIN) {
            int n = Info[R3].nDep;
            Info[R3].iDep[n] = it->first;
            ++Info[R3].nDep;
        }

        if (R4!=NOTIN) {
            int n = Info[R4].nDep;
            Info[R4].iDep[n] = it->first;
            ++Info[R4].nDep;
        }

        if (R5!=NOTIN) {
            int n = Info[R5].nDep;
            Info[R5].iDep[n] = it->first;
            ++Info[R5].nDep;
        }

        if (R6!=NOTIN) {
            int n = Info[R6].nDep;
            Info[R6].iDep[n] = it->first;
            ++Info[R6].nDep;
        }
    }
}

void EvaluationScheme::LinkBackK4() {
    UI32 ndeps = 0;

    map <UI32, VarInfoK4>::iterator it;

    //assign variables for each integral
    for (it = K4Info.begin(); it != K4Info.end(); ++it) {
        it->second.nDep = 0;
    }

    for (it = K4Info.begin(); it != K4Info.end(); ++it) {
        UI32 R1 = it->second.VD1;
        UI32 R2 = it->second.VD2;
        UI32 R3 = it->second.VD3;
        UI32 R4 = it->second.VD4;
        UI32 R5 = it->second.VD5;
        UI32 R6 = it->second.VD6;
        UI32 RE = it->second.VDE;
        if (R1!=NOTIN) {++K4Info[R1].nDep; ++ndeps;}
        if (R2!=NOTIN) {++K4Info[R2].nDep; ++ndeps;}
        if (R3!=NOTIN) {++K4Info[R3].nDep; ++ndeps;}
        if (R4!=NOTIN) {++K4Info[R4].nDep; ++ndeps;}
        if (R5!=NOTIN) {++K4Info[R5].nDep; ++ndeps;}
        if (R6!=NOTIN) {++K4Info[R6].nDep; ++ndeps;}
        if (RE!=NOTIN) {++K4Info[RE].nDep; ++ndeps;}
    }

    NEVS = ndeps; //use extra variables (for later simplifying)

    //create list of integrals in order to link back
    RDep = new UI32[ndeps+NEVS];

    //assign pointers
    UI32 * piC = RDep;

    for (it = K4Info.begin(); it != K4Info.end(); ++it) {
        it->second.iDep = piC;
        piC += it->second.nDep;
        it->second.nDep = 0;
    }

    LastDep = piC; //pointer to first free var


    //copies integral to where it belongs in the list and increments counter
    for (it = K4Info.begin(); it != K4Info.end(); ++it) {
        UI32 R1 = it->second.VD1;
        UI32 R2 = it->second.VD2;
        UI32 R3 = it->second.VD3;
        UI32 R4 = it->second.VD4;
        UI32 R5 = it->second.VD5;
        UI32 R6 = it->second.VD6;
        UI32 RE = it->second.VDE;

        if (R1!=NOTIN) {
            int n = K4Info[R1].nDep;
            K4Info[R1].iDep[n] = it->first;
            ++K4Info[R1].nDep;
        }

        if (R2!=NOTIN) {
            int n = K4Info[R2].nDep;
            K4Info[R2].iDep[n] = it->first;
            ++K4Info[R2].nDep;
        }

        if (R3!=NOTIN) {
            int n = K4Info[R3].nDep;
            K4Info[R3].iDep[n] = it->first;
            ++K4Info[R3].nDep;
        }

        if (R4!=NOTIN) {
            int n = K4Info[R4].nDep;
            K4Info[R4].iDep[n] = it->first;
            ++K4Info[R4].nDep;
        }

        if (R5!=NOTIN) {
            int n = K4Info[R5].nDep;
            K4Info[R5].iDep[n] = it->first;
            ++K4Info[R5].nDep;
        }

        if (R6!=NOTIN) {
            int n = K4Info[R6].nDep;
            K4Info[R6].iDep[n] = it->first;
            ++K4Info[R6].nDep;
        }

        if (RE!=NOTIN) {
            int n = K4Info[RE].nDep;
            K4Info[RE].iDep[n] = it->first;
            ++K4Info[RE].nDep;
        }
    }
}


void EvaluationScheme::ComputeLife() {
    list<UI32>::iterator it;
    //compute lifelihood of all variables

    UI32 * Life = new UI32[N+2];

    for (int i=1; i<=N+1; ++i) Life[i] = 0;

    for (it=EvaluationOrder.begin(); it!=EvaluationOrder.end(); ++it) {
        UI32 Id = *it;

        UI32 Pos = Info[Id].position;

        UI32 Id1 = Info[Id].VD1;
        UI32 Id2 = Info[Id].VD2;
        UI32 Id3 = Info[Id].VD3;
        UI32 Id4 = Info[Id].VD4;
        UI32 Id5 = Info[Id].VD5;
        UI32 Id6 = Info[Id].VD6;

        if (Id1!=NOTIN) {
            UI32 Pos1 = Info[Id1].position;
            Life[Pos1] = max(Life[Pos1], Pos);
        }
        if (Id2!=NOTIN) {
            UI32 Pos2 = Info[Id2].position;
            Life[Pos2] = max(Life[Pos2], Pos);
        }
        if (Id3!=NOTIN) {
            UI32 Pos3 = Info[Id3].position;
            Life[Pos3] = max(Life[Pos3], Pos);
        }
        if (Id4!=NOTIN) {
            UI32 Pos4 = Info[Id4].position;
            Life[Pos4] = max(Life[Pos4], Pos);
        }
        if (Id5!=NOTIN) {
            UI32 Pos5 = Info[Id5].position;
            Life[Pos5] = max(Life[Pos5], Pos);
        }
        if (Id6!=NOTIN) {
            UI32 Pos6 = Info[Id6].position;
            Life[Pos6] = max(Life[Pos6], Pos);
        }
    }

    for (it=EvaluationOrder.begin(); it!=EvaluationOrder.end(); ++it) {
        UI32 Id = *it;
        UI32 Pos = Info[Id].position;
        //if (Info[Id].final) Info[Id].lifelihood = N+1;
        //else                Info[Id].lifelihood = Life[Pos];
        Info[Id].lifelihood = Life[Pos];
    }

    delete[] Life;
}

void EvaluationScheme::Plot() {
    list<UI32>::const_iterator it;

    cout << int(La) << " " << int(Lb) << " " << int(Lc) << " " << int(Ld) << endl;

    for (it=K4EvaluationOrder.begin(); it!=K4EvaluationOrder.end(); ++it) {
        UI32 Id = *it;
        cout << K4Info[Id].position;
        if (K4Info[Id].final) cout << " * ";
        else                cout << "   ";

        cout << IdK4Map[Id] << "  ";

        UI32 Id1 = K4Info[Id].VD1;
        UI32 Id2 = K4Info[Id].VD2;
        UI32 Id3 = K4Info[Id].VD3;
        UI32 Id4 = K4Info[Id].VD4;
        UI32 Id5 = K4Info[Id].VD5;
        UI32 Id6 = K4Info[Id].VD6;
        UI32 IdE = K4Info[Id].VDE;
        if (Id1!=NOTIN) cout << K4Info[Id1].position << " ";
        else         cout << "* ";
        if (Id2!=NOTIN) cout << K4Info[Id2].position << " ";
        else         cout << "* ";
        if (Id3!=NOTIN) cout << K4Info[Id3].position << " ";
        else         cout << "* ";
        if (Id4!=NOTIN) cout << K4Info[Id4].position << " ";
        else         cout << "* ";
        if (Id5!=NOTIN) cout << K4Info[Id5].position << " ";
        else         cout << "* ";
        if (Id6!=NOTIN) cout << K4Info[Id6].position << " ";
        else         cout << "* ";
        if (IdE!=NOTIN) cout << K4Info[IdE].position << " ";
        else         cout << "* ";

        cout << "  " << K4Info[Id].lifelihood << "  " << K4Info[Id].memory << "  " << K4Info[Id].ETYPE << endl;


        UI32 * R = K4Info[Id].iDep;
        int n = K4Info[Id].nDep;

        for (int i=0; i<n; ++i) cout  << "   " << K4Info[R[i]].position;
        cout << endl;
    }

    //char a; cin >> a;

    return;

    //just print the K4 contraction loops; the MIRROR routines are already working



    for (it=EvaluationOrder.begin(); it!=EvaluationOrder.end(); ++it) {
        UI32 Id = *it;
        cout << Info[Id].position;
        if (Info[Id].final) cout << " * ";
        else                cout << "   ";

        cout << IdIntMap[Id] << "  ";

        UI32 Id1 = Info[Id].VD1;
        UI32 Id2 = Info[Id].VD2;
        UI32 Id3 = Info[Id].VD3;
        UI32 Id4 = Info[Id].VD4;
        UI32 Id5 = Info[Id].VD5;
        UI32 Id6 = Info[Id].VD6;
        if (Id1!=NOTIN) cout << Info[Id1].position << " ";
        else         cout << "* ";
        if (Id2!=NOTIN) cout << Info[Id2].position << " ";
        else         cout << "* ";
        if (Id3!=NOTIN) cout << Info[Id3].position << " ";
        else         cout << "* ";
        if (Id4!=NOTIN) cout << Info[Id4].position << " ";
        else         cout << "* ";
        if (Id5!=NOTIN) cout << Info[Id5].position << " ";
        else         cout << "* ";
        if (Id6!=NOTIN) cout << Info[Id6].position << " ";
        else         cout << "* ";

        cout << "  " << Info[Id].lifelihood << "  " << Info[Id].memory << "  " << Info[Id].ETYPE << endl;


        UI32 * R = Info[Id].iDep;
        int n = Info[Id].nDep;

        for (int i=0; i<n; ++i) cout  << "   " << Info[R[i]].position;
        cout << endl;
    }

}

void EvaluationScheme::Simplify() {
    list<UI32>::iterator it;

    //find if expression is simplifiable
    //**********************************
    bool mod = false;

    bool compacted = true;


    while (compacted) {
        compacted = false;

        //first pass: intermediate variables
        for (it=EvaluationOrder.begin(); it!=EvaluationOrder.end();) {
            UI32 Id = *it;

            //DON'T attempt to simplify the kernels !!!!
            if (Info[Id].ETYPE==KERNELS) {
                ++it;
                continue;
            }

            //right now, it won't work with final variables
            if (Info[Id].final) {
                ++it;
                continue;
            }


            UI32 Idb = 0;

            //NULL HRR
            if ( (Info[Id].VD1!=NOTIN && Info[Id].VD2==NOTIN) &&
                 (Info[Id].ETYPE == HRRBX || Info[Id].ETYPE == HRRBY || Info[Id].ETYPE == HRRBZ || Info[Id].ETYPE == HRRKX || Info[Id].ETYPE == HRRKY || Info[Id].ETYPE == HRRKZ) ) {

                Idb = Info[Id].VD1;
            }

            //NULL CTE
            if ( (Info[Id].VD1==NOTIN && Info[Id].VD2==NOTIN && Info[Id].VD3!=NOTIN) &&
                 (Info[Id].ETYPE == CTEBX || Info[Id].ETYPE == CTEBY  || Info[Id].ETYPE == CTEBZ) ) {
                     // || Info[Id].ETYPE == CTEKX || Info[Id].ETYPE == CTEKY || Info[Id].ETYPE == CTEKZ

                Idb = Info[Id].VD3;
            }



            //SPH A COINCIDES WITH CARTESIAN
            if ( (Info[Id].VD1!=NOTIN && Info[Id].VD2==NOTIN) && (La<2) && (Info[Id].ETYPE == SPHA) ) {
                Idb = Info[Id].VD1;
            }

            //SPH B COINCIDES WITH CARTESIAN
            if ( (Info[Id].VD1!=NOTIN && Info[Id].VD2==NOTIN) && (Lb<2) && (Info[Id].ETYPE == SPHB) ) {
                Idb = Info[Id].VD1;
            }

            //SPH C COINCIDES WITH CARTESIAN
            if ( (Info[Id].VD1!=NOTIN && Info[Id].VD2==NOTIN) && (Lc<2) && (Info[Id].ETYPE == SPHC) ) {
                Idb = Info[Id].VD1;
            }

            //SPH D COINCIDES WITH CARTESIAN
            if ( (Info[Id].VD1!=NOTIN && Info[Id].VD2==NOTIN) && (Ld<2) && (Info[Id].ETYPE == SPHD) ) {
                Idb = Info[Id].VD1;
            }



            //remove Id from evaluation
            if (Idb != 0 && !Info[Idb].zero) {

                //check if dependency list is large enough

                int ndep = Info[Idb].nDep + Info[Id].nDep -1;

                if (ndep>=NEVS) {
                    //compact the list
                    compacted = true;
                    break;
                }

                //substitute all forward dependencies
                for (int n=0; n<Info[Id].nDep; ++n) {
                    UI32 Idf = Info[Id].iDep[n];

                    if (Info[Idf].VD1 == Id) Info[Idf].VD1 = Idb;
                    if (Info[Idf].VD2 == Id) Info[Idf].VD2 = Idb;
                    if (Info[Idf].VD3 == Id) Info[Idf].VD3 = Idb;
                    if (Info[Idf].VD4 == Id) Info[Idf].VD4 = Idb;
                    if (Info[Idf].VD5 == Id) Info[Idf].VD5 = Idb;
                    if (Info[Idf].VD6 == Id) Info[Idf].VD6 = Idb;
                }

                //substitute all backward dependencies
                //&& merge forward dependency lists

                //if the number of forward dependencies of any of the two variables is one,
                //use that place in the total list
                if (Info[Idb].nDep == 1) {
                    //NULL the variable
                    *(Info[Idb].iDep) = 0;

                    //copy
                    Info[Idb].nDep = Info[Id].nDep;
                    Info[Idb].iDep = Info[Id].iDep;
                }
                else if (Info[Id].nDep == 1) {
                    //NULL the variable
                    UI32 IdDep = *(Info[Id].iDep);
                    *(Info[Id].iDep) = 0;

                    //find the place for the variable ans substitute
                    for (int n=0; n<Info[Idb].nDep; ++n) {
                        if (Info[Idb].iDep[n] == Id)
                            Info[Idb].iDep[n] = IdDep;
                    }
                }
                //null the old dependency lists and append a the merged list to the end of the total list
                else {
                    UI32 * merge;

                    if (ndep<NEVS) {
                        merge = LastDep;
                        LastDep += ndep;
                    }
                    else {
                        merge = new UI32[ndep];
                        cout << "beware the pointer!!!" << endl; //small memory leak
                    }
                    NEVS -= ndep;

                    int m=0;
                    for (int n=0; n<Info[Idb].nDep; ++n) {
                        //don't copy the original one
                        if (Info[Idb].iDep[n] != Id) {
                            merge[m] = Info[Idb].iDep[n];
                            ++m;
                        }
                        Info[Idb].iDep[n] = 0; //NULL the variable
                    }
                    for (int n=0; n<Info[Id].nDep; ++n) {
                        merge[m] = Info[Id].iDep[n];
                        ++m;
                        Info[Id].iDep[n] = 0; //NULL the variable
                    }

                    Info[Idb].nDep = ndep;
                    Info[Idb].iDep = merge;
                }


                //erase this variable
                integralCart Int = IdIntMap[Id];
                Info[Id].iDep = NULL;
                IdIntMap.erase(Id);
                IntIdMap.erase(Int);
                Info.erase(Id);

                //don't increment iterator!!!!
                it = EvaluationOrder.erase(it);
                --N;

                mod = true;
            }
            else
                ++it;

        }

        //recompute all positions
        if (mod == true) {
            if (NEVS<0) cout << "needed: " << (-NEVS) << endl;

            int pos = 0;
            for (it=EvaluationOrder.begin(); it!=EvaluationOrder.end(); ++it) {
                ++pos;
                UI32 Id = *it;
                Info[Id].position = pos;
            }
        }

        if (compacted == true) {
            delete[] RDep;
            LinkBack();
        }
    }

    //(re)compute lifeness of vars
    ComputeLife();

    //count number of instructions for each block
    UI32 iblock = 0;

    for (it=K4EvaluationOrder.begin(); it!=K4EvaluationOrder.end(); ++it) {
        UI32 Id = *it;
        while (K4Info[Id].ETYPE != RRSK4[iblock]) ++iblock;
        ++ninstrK4[iblock];
    }


    iblock = 0;

    for (it=EvaluationOrder.begin(); it!=EvaluationOrder.end(); ++it) {
        UI32 Id = *it;
        while (Info[Id].ETYPE != RRS[iblock]) ++iblock;
        ++ninstr[iblock];
    }
}

//count number of instructions for each block
void EvaluationScheme::Anneal() {
    list<UI32>::iterator it;

    UI32 iblock = 0;

    for (int i=0; i<nblocks; ++i) ninstr[i] = 0;

    for (it=EvaluationOrder.begin(); it!=EvaluationOrder.end(); ++it) {
        UI32 Id = *it;
        while (Info[Id].ETYPE != RRS[iblock]) ++iblock;
        ++ninstr[iblock];
    }

    /*
    UI32 * weights = new UI32[N+1];

    list<UI32>::iterator it;

    for (it=EvaluationOrder.begin(); it!=EvaluationOrder.end(); ++it) {
        VarInfo & info = Info[*it];
        weights[info.position] = info.nBDep - info.nDep;
    }

    for (int iblock=1; iblock<nblocks; ++iblock) {

    }

    delete[] weights;
    */
}

//memory required for the MIRROR transformations
UI32 EvaluationScheme::RequiredMem() {
    char * An = new char[N+3];

    for (int i=0; i<=N+2; ++i) An[i] = 0;

    map<UI32, VarInfo>::const_iterator it;

    //could be done in parallel
    for(it=Info.begin(); it!=Info.end(); ++it) {
        UI32 l = it->second.lifelihood;

        //if (it->second.ETYPE >= KERNELS)
        //++An[l+1]; //free variables at next step
        ++An[l]; // a bit twisted; can reuse same memory slot if variable dies POTENTIAL CHAIN DEPENDENCY!!!!
    }

    //count number of variables in KERNELS
    UI32 ppos = 0;
    UI32 totvars=0;
    UI32 totvars2=0;

    totvars = 0; //ninstrK4[K4s];


    for (int iblock=KERNELS; iblock<REORDER; ++iblock) {
        maxalive[iblock] = totvars;

        totvars2 += ninstr[iblock];

        for (int ins=0; ins<ninstr[iblock]; ++ins) {
            ++ppos;
            totvars -= An[ppos]; //variable muertas en esta ecuacion
            ++totvars;           //nueva variable declarada

            if (totvars>=maxalive[iblock]) {
                maxalive[iblock] = totvars;
            }
        }
    }

    UI32 tmaxalive = 0;
    for (int iblock=KERNELS; iblock<REORDER; ++iblock) {
        if (maxalive[iblock]>tmaxalive) tmaxalive = maxalive[iblock];
    }
    delete[] An;

    return tmaxalive;
}

struct mempos {
    UI32 t;
    UI32 p;

    mempos () {}

    mempos (UI32 tt, UI32 pp) {
        t=tt;
        p=pp;
    }

    bool operator<(const mempos & rhs) const {
        if (t!=rhs.t) return (t>rhs.t);
        return (p>rhs.p);
    }
};

void EvaluationScheme::AssignMemory() {

    //assign memory offset to the vars in the K4 loops
    //================================================
    {
        //set all memory positions, kernels included
        list<UI32>::iterator itk4;
        itk4=K4EvaluationOrder.begin();

        //small stuff
        {
            UI32 mem = 0;

            for (int ins=0; ins<ninstrK4[BOYS]; ++ins) {
                UI32 Id = *itk4;
                K4Info[Id].memory = mem;
                ++mem;
                ++itk4;
            }
        }


        // using different memory pools so that there's an increased flexibility in the order of application of RRs
        UI32 memD = 0;
        UI32 memF = 0;
        UI32 memE = 0;

        //rest of RRs
        for (int iblock=BOYS+1; iblock<NADA; ++iblock) {
            UI32 * pmem;

            //in general contraction, use different buffers for each J4 step (since the size of each memory chunk cannot be predicted ono compiler time)
            if      (iblock==K4D || iblock==K3D || iblock==K2D || iblock==K1D) {memD=0; pmem = &memD;}
            else if (iblock==K4E || iblock==K3E || iblock==K2E || iblock==K1E) {memE=0; pmem = &memE;}
            else if (iblock==K4F || iblock==K3F || iblock==K2F || iblock==K1F) {memF=0; pmem = &memF;}

            else if (iblock==ADRR4 || iblock==ADRR3 || iblock==ADRR2 || iblock==ADRR1 || iblock==ADRR0) pmem = &memD;
            else if (iblock==AERR4 || iblock==AERR3 || iblock==AERR2 || iblock==AERR1 || iblock==AERR0) pmem = &memE;
            else if (iblock==CDR4  || iblock==CDR3  || iblock==CDR2  || iblock==CDR1  || iblock==CDR0 ) pmem = &memF;
            else pmem = NULL;

            for (int ins=0; ins<ninstrK4[iblock]; ++ins) {
                UI32 Id = *itk4;
                K4Info[Id].memory = *pmem;
                ++(*pmem);
                ++itk4;
            }

        }
    }

    /*
    {
        list<UI32>::iterator it=EvaluationOrder.begin();
        //assign separate memory for last block (since it's in another array)
        for (int ins=0; ins<ninstr[KERNELS]; ++ins) {
            UI32 Id = *it;
            Info[Id].VD1 = Info[Id].aux;
            ++it;
        }
    }
    */



    //now for the MIRROR transformations
    //==================================
    {
        UI32 t=0;      //time

        UI32 MemSize = RequiredMem();
        MaxMem = MemSize;

        //second simplest method possible
        priority_queue <mempos> WorkingMem;
        for (UI32 p=0; p<MemSize; ++p) {
            mempos mp(t, p);
            WorkingMem.push(mp);
        }


        list<UI32>::iterator it=EvaluationOrder.begin();

        //add rest of the blocks
        for (int iblock=KERNELS; iblock<REORDER; ++iblock) {
            for (int ins=0; ins<ninstr[iblock]; ++ins) {
                UI32 Id = *it;

                ++t; //increment time

                mempos pos; pos = WorkingMem.top(); WorkingMem.pop();

                if (pos.t>t) cout << "error" << endl;

                Info[Id].memory = pos.p;
                pos.t = Info[Id].lifelihood;

                WorkingMem.push(pos);
                ++it;
            }
        }

        //assign separate memory for last block (since it's in another array)
        for (int ins=0; ins<ninstr[REORDER]; ++ins) {
            UI32 Id = *it;
            Info[Id].memory = Info[Id].aux;
            ++it;
        }
    }


}


void EvaluationScheme::GenerateCode() {
    //count number of blocks, instructions, etc
    //don't skip kernels


    //MIRROR
    //======
    {

        list<UI32>::iterator itc;

        UI32 totalinstr = 0;
        for (int i=KERNELS; i<=REORDER; ++i) totalinstr += ninstr[i];


        //reserve memory for the instructions
        sequence = new MCinstr[totalinstr+16]; //to account for prefetching

        MCinstr * s2 = sequence;
        for (int iblock=KERNELS; iblock<=REORDER; ++iblock) {
            seqs[iblock] = s2;      //where the block starts
            s2 += ninstr[iblock];   //instructions
        }
        seqs[REORDER+1] = s2;

        //SKIP LAST STEP
        itc=EvaluationOrder.begin();

        //link kernels
        {

            for (int i=0; i<ninstr[KERNELS]; ++i) {
                UI32 Id = *itc;

                UI32 dest = Info[Id].memory;
                UI32 IDK  = Info[Id].aux;

                UI32 m1;
                UI32 aux;

                if (IDK==NOTIN) {
                    m1  = 0;
                    aux = 1; // load a 0 value
                    //cout << "Kernel not referenced!" << endl;
                }
                else {
                    m1  = K4Info[IDK].memory;
                    aux = 0;
                }

                seqs[KERNELS][i].dest = dest;
                seqs[KERNELS][i].op1 = m1;
                seqs[KERNELS][i].op2 = 0;
                seqs[KERNELS][i].op3 = 0;
                seqs[KERNELS][i].op4 = 0;
                seqs[KERNELS][i].op5 = 0;
                seqs[KERNELS][i].op6 = 0;
                seqs[KERNELS][i].aux = aux;
                seqs[KERNELS][i].ETYPE = KERNELS;
                ++itc;
            }
        }


        for (int iblock=KERNELS+1; iblock<=REORDER; ++iblock) {
            for (int i=0; i<ninstr[iblock]; ++i) {
                UI32 Id = *itc;

                UI32 dest = Info[Id].memory;

                UI32 ID1, ID2, ID3, ID4, ID5, ID6;
                ID1 = Info[Id].VD1;
                ID2 = Info[Id].VD2;
                ID3 = Info[Id].VD3;
                ID4 = Info[Id].VD4;
                ID5 = Info[Id].VD5;
                ID6 = Info[Id].VD6;

                UI32 m1 = 0;
                UI32 m2 = 0;
                UI32 m3 = 0;
                UI32 m4 = 0;
                UI32 m5 = 0;
                UI32 m6 = 0;

                if (ID1!=NOTIN) m1 = Info[ID1].memory;
                if (ID2!=NOTIN) m2 = Info[ID2].memory;
                if (ID3!=NOTIN) m3 = Info[ID3].memory;
                if (ID4!=NOTIN) m4 = Info[ID4].memory;
                if (ID5!=NOTIN) m5 = Info[ID5].memory;
                if (ID6!=NOTIN) m6 = Info[ID6].memory;

                seqs[iblock][i].dest = dest;
                seqs[iblock][i].op1 = m1;
                seqs[iblock][i].op2 = m2;
                seqs[iblock][i].op3 = m3;
                seqs[iblock][i].op4 = m4;
                seqs[iblock][i].op5 = m5;
                seqs[iblock][i].op6 = m6;
                seqs[iblock][i].aux = Info[Id].aux;
                seqs[iblock][i].ETYPE = Info[Id].ETYPE;
                ++itc;
            }
        }

    }


    //K4
    //==
    {
        list<UI32>::iterator itc;

        UI32 totalinstrK4 = 0;
        for (int i=BOYS; i<NADA; ++i) totalinstrK4 += ninstrK4[i];


        //reserve memory for the instructions
        sequenceK4 = new MCinstrK4[totalinstrK4+16]; //to account for prefetching

        MCinstrK4 * sk4 = sequenceK4;
        for (int iblock=BOYS; iblock<NADA; ++iblock) {
            seqsK4[iblock] = sk4;      //where the block starts
            sk4 += ninstrK4[iblock];   //instructions
        }
        seqsK4[NADA] = sk4;



        //fills the void in memory

        itc=K4EvaluationOrder.begin();

        for (int iblock=BOYS; iblock<NADA; ++iblock) {
            for (int i=0; i<ninstrK4[iblock]; ++i) {
                UI32 Id = *itc;

                UI32 dest = K4Info[Id].memory;

                UI32 ID1, ID2, ID3, ID4, ID5, ID6, IDE;
                ID1 = K4Info[Id].VD1;
                ID2 = K4Info[Id].VD2;
                ID3 = K4Info[Id].VD3;
                ID4 = K4Info[Id].VD4;
                ID5 = K4Info[Id].VD5;
                ID6 = K4Info[Id].VD6;
                IDE = K4Info[Id].VDE;

                UI32 m1 = 0;
                UI32 m2 = 0;
                UI32 m3 = 0;
                UI32 m4 = 0;
                UI32 m5 = 0;
                UI32 m6 = 0;
                UI32 me = 0;

                if (ID1!=NOTIN) m1 = K4Info[ID1].memory;
                if (ID2!=NOTIN) m2 = K4Info[ID2].memory;
                if (ID3!=NOTIN) m3 = K4Info[ID3].memory;
                if (ID4!=NOTIN) m4 = K4Info[ID4].memory;
                if (ID5!=NOTIN) m5 = K4Info[ID5].memory;
                if (ID6!=NOTIN) m6 = K4Info[ID6].memory;
                if (IDE!=NOTIN) me = K4Info[IDE].memory;

                seqsK4[iblock][i].dest = dest;
                seqsK4[iblock][i].op1 = m1;
                seqsK4[iblock][i].op2 = m2;
                seqsK4[iblock][i].op3 = m3;
                seqsK4[iblock][i].op4 = m4;
                seqsK4[iblock][i].op5 = m5;
                seqsK4[iblock][i].op6 = m6;
                seqsK4[iblock][i].ope = me;

                seqsK4[iblock][i].aux   = K4Info[Id].aux;
                seqsK4[iblock][i].ETYPE = K4Info[Id].ETYPE;

                ++itc;
            }
        }

    }


    //COUNT FLOPs
    //===========

    NFLOPS = 0;
    {
        for (int iblock=KERNELS+1; iblock<=REORDER; ++iblock) {
            for (int i=0; i<ninstr[iblock]; ++i) {

                if (iblock == MMDZ || iblock == MMDY || iblock == MMDX) {
                    if (seqs[iblock][i].op1!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op2!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op3!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op4!=0) NFLOPS += 2;
                    if (seqs[iblock][i].aux==1) --NFLOPS;
                    --NFLOPS;
                }

                if (iblock == CTEBX || iblock == CTEBY || iblock == CTEBZ || iblock == CTEKX || iblock == CTEKY || iblock == CTEKZ) {
                    if (seqs[iblock][i].op1!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op2!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op3!=0) NFLOPS += 1;
                    if (seqs[iblock][i].aux==1) --NFLOPS;
                    --NFLOPS;
                }

                if (iblock == HRRBX || iblock == HRRBY || iblock == HRRBZ || iblock == HRRKX || iblock == HRRKY || iblock == HRRKZ) {
                    if (seqs[iblock][i].op2!=0) NFLOPS += 2;
                }

                if (iblock == SPHA || iblock == SPHB || iblock == SPHC || iblock == SPHD) {
                    if (seqs[iblock][i].op1!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op2!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op3!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op4!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op5!=0) NFLOPS += 2;
                    if (seqs[iblock][i].op6!=0) NFLOPS += 2;
                    --NFLOPS;
                    if (iblock == SPHA && La<2) --NFLOPS;
                    if (iblock == SPHB && Lb<2) --NFLOPS;
                    if (iblock == SPHC && Lc<2) --NFLOPS;
                    if (iblock == SPHD && Ld<2) --NFLOPS;
                }
            }
        }

    }

}

