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

#ifndef __NEWRRS__
#define __NEWRRS__
#include <set>
#include <list>
#include <iostream>

#include "../defs.hpp"
#include "integrals.hpp"
#include "../2eints/IICinit.hpp"

static const UI32 NOTIN = -1;

//this is HUGE!!!
class MCinstr {
    public:
    UI32 dest;
    UI32 op1;
    UI32 op2;
    UI32 op3;
    UI32 op4;
    UI32 op5;
    UI32 op6;
    UI16 aux;
    RRTYPE ETYPE;
} __attribute__((aligned(32)));

class MCinstrK4 {
    public:
    UI32 dest;
    UI32 op1;
    UI32 op2;
    UI32 op3;
    UI32 op4;
    UI32 op5;
    UI32 op6;
    UI32 ope;
    UI16 aux;
    RRTYPE2 ETYPE;

} __attribute__((aligned(32)));



class VarDepK4 {

  public:
    ikernel VD1;
    ikernel VD2;
    ikernel VD3;
    ikernel VD4;
    ikernel VD5;
    ikernel VD6;

    ikernel VDE;

    RRTYPE2 ETYPE;
    UI16 aux;

    UI32 Life;
    UI32 nRefs;
};

class VarDep {

  public:
    integralCart VD1;
    integralCart VD2;
    integralCart VD3;
    integralCart VD4;
    integralCart VD5;
    integralCart VD6;

    RRTYPE ETYPE;
    UI16 aux;

    UI32 Life;
    UI32 nRefs;

    VarDep() : VD1(NOTIN),VD2(NOTIN),VD3(NOTIN),VD4(NOTIN),VD5(NOTIN),VD6(NOTIN) {}

    VarDep & operator=(const VarDepK4 & rhs);
};

class EvaluationScheme {

  private:


    class VarInfo {

      public:
        IDINT VD1;
        IDINT VD2;
        IDINT VD3;
        IDINT VD4;
        IDINT VD5;
        IDINT VD6;

        RRTYPE ETYPE;
        int aux;

        UI32 position;
        UI32 lifelihood;
        UI32 memory;

        UI8 nBDep;
        UI8 nDep;
        IDINT * iDep;

        bool final;
        bool zero;

        VarInfo() {
            position = lifelihood = nDep = nBDep = memory = aux = 0;
            VD1 = VD2 = VD3 = VD4 = VD5 = VD6 = -1;
            final = false;
            zero = false;
            iDep = NULL;
        }

        ~VarInfo() {
            iDep = NULL;
        }
    };

    class VarInfoK4 {

      public:
        IDINT VD1;
        IDINT VD2;
        IDINT VD3;
        IDINT VD4;
        IDINT VD5;
        IDINT VD6;
        IDINT VDE;

        RRTYPE2 ETYPE;
        int aux;

        UI32 position;
        UI32 lifelihood;
        UI32 memory;

        UI8 nBDep;
        UI8 nDep;
        IDINT * iDep;

        bool final;
        bool zero;

        VarInfoK4() {
            position = lifelihood = nDep = nBDep = memory = 0;
            VD1 = VD2 = VD3 = VD4 = VD5 = VD6 = VDE = -1;
            final = false;
            zero = false;
            iDep = NULL;
        }

        ~VarInfoK4() {
            iDep = NULL;
        }
    };

    class K4Var {
      public:
        ikernel      K;
        VarDepK4     V;
    };

    class IntVar {
      public:
        integralCart I;
        VarDep       V;
    };

    class Var2 {
      public:
        int AN;
        int AV;

        bool operator<(const Var2 & rhs) const {
            if (AN<0) {
                if (rhs.AN>=0) return true;
                return (AV<rhs.AV);
            }
            else if (AN>0) {
                if (rhs.AN<=0) return false;
                return (AV>rhs.AV);
            }
            else {
                if (rhs.AN<0) return false;
                else if (rhs.AN>0) return true;
                else return (AV<rhs.AV);
            }
        }
    };

    UI8 La;
    UI8 Lb;
    UI8 Lc;
    UI8 Ld;
    GEOM geometry;
    bool useGC;



    //Initial Data (to be descarded once the variables are linked)
    std::list<IntVar> EVList;
    std::multimap<Var2, std::list<IntVar> * > groups; //subgroups before sorting and appending to the evaluation list

    //Data to be used after linkage of variables
    std::map <IDINT, integralCart> IdIntMap; //maps each identifier to its integral
    std::map <integralCart, IDINT> IntIdMap;  //maps each integral to its Id
    std::map <IDINT, VarInfo> Info;        //Info for each variable
    std::list<IDINT> EvaluationOrder;      //order of evaluation



    std::list<K4Var> EVListK4;

    std::map <IDINT, ikernel> IdK4Map; //maps each identifier to its integral
    std::map <ikernel, IDINT> K4IdMap;  //maps each integral to its Id
    std::map <IDINT, VarInfoK4> K4Info;        //Info for each variable
    std::list<IDINT> K4EvaluationOrder;      //order of evaluation


    //Data for code generation and testing
    //************************************
    MCinstr * sequence; //total instruction sequence
    UI32 maxalive[32];
    RRTYPE RRS[32];     //type of RR for the block
    UI8 nblocks;       //number of blocks


    MCinstrK4 * sequenceK4; //total instruction sequence
    RRTYPE2 RRSK4[32];     //type of RR for the block
    UI8 nblocksK4;         //number of blocks


    IDINT * RDep;                          //array with the values of the link-back integrals
    IDINT * LastDep;                       //pointer to the last used element in RDep
    int NEVS;                             //number of extra variables

    UI32 N;      //total number of variables
    UI32 NK4;    //total number of variables

    UI32 GetId(const integralCart & i) const;

    int  CalcLife(std::list<IntVar> * glist, std::set<integralCart> & setf);
    void AddGroup(std::map<integralCart,VarDep> (&imap)[1+4*LMAX+4], std::set<integralCart> & setf);
    void AddGroup(std::map<ikernel,VarDepK4> (&imap), std::set<ikernel> & setf);
    void AddGroups();

    void Append  (const IntVar & IV);
    void AppendK4(const K4Var  & KV);

  public:

    // information to be accessed and copied by ERIroutine
    // ===================================================
    MCinstr * seqs[32]; //pointer to each block's sequence
    UI32 ninstr[32];    //number of intructions per block
    MCinstrK4 * seqsK4[32]; //pointer to each block's sequence
    UI32 ninstrK4[32];     //number of intructions per block


    UI32 MaxMem; //max memory needed for the evaluation
    UI32 NFLOPS;


    EvaluationScheme();
    ~EvaluationScheme();

    void Set(UI8 la, UI8 lb, UI8 lc, UI8 ld, GEOM geom, bool GC) {
        La = la;
        Lb = lb;
        Lc = lc;
        Ld = ld;
        geometry = geom;
        useGC    = GC;
    }



    void LinkVars(const std::set<integralCart> & set0);
    void LinkBack();
    void LinkBackK4();
    void ComputeLife();
    void Simplify();
    void Anneal();
    UI32 RequiredMem();
    void AssignMemory();

    void SphericalKernelSieve(std::set<integralCart> & setf);



    void LinkVarsK4(const std::set<ikernel> & setk);
    void Plot();
    void GenerateCode();

    void TestCode() const;
    void WriteCode(const std::string & cppnameC, const std::string & cppnameT, const std::string & fname, int la, int lb, int lc, int ld, GEOM geom);


    //Recurrence relations
    void UseCDR4 (const std::set<ikernel> & setF4, std::set<ikernel> & setE4, int totL);
    void UseAERR4(const std::set<ikernel> & setE4, int totL);

    void UseCDR3 (const std::set<ikernel> & setF3, std::set<ikernel> & setE3, std::set<ikernel> & setF3m, bool ABt, bool CDt, bool ACt);
    void UseAERR3(const std::set<ikernel> & setE3, std::set<ikernel> & setE3m);

    void UseCDR2 (const std::set<ikernel> & setF2, std::set<ikernel> & setE2, std::set<ikernel> & setF2m, bool ABt, bool CDt, bool ACt);
    void UseAERR2(const std::set<ikernel> & setE2, std::set<ikernel> & setE2m);

    void UseCDR1 (const std::set<ikernel> & setF1, std::set<ikernel> & setE1, std::set<ikernel> & setF1m, bool ABt, bool CDt, bool ACt);
    void UseAERR1(const std::set<ikernel> & setE1, std::set<ikernel> & setE1m);

    void UseCDR0 (const std::set<ikernel> & setF0, std::set<ikernel> & setE0, std::set<ikernel> & setF0m, bool ABt, bool CDt, bool ACt);
    void UseAERR0(const std::set<ikernel> & setE0, std::set<ikernel> & setE0m);



    void SieveK0K1(const std::set<ikernel> & setK0, std::set<ikernel> & setK1, bool ABt, bool CDt, bool ACt);
    void SieveK1K2(const std::set<ikernel> & setK1, std::set<ikernel> & setK2, bool ABt, bool CDt, bool ACt);
    void SieveK2K3(const std::set<ikernel> & setK2, std::set<ikernel> & setK3, bool ABt, bool CDt, bool ACt);
    void SieveK3K4(const std::set<ikernel> & setK3, std::set<ikernel> & setK4, bool ABt, bool CDt, bool ACt);

    void addS(const std::set<ikernel> & set0, std::set<ikernel> & setf);
    void addT(const std::set<ikernel> & set0, std::set<ikernel> & setf);
    void addU(const std::set<ikernel> & set0, std::set<ikernel> & setf);
    void addV(const std::set<ikernel> & set0, std::set<ikernel> & setf);

    void addBoys(int totL);

    void add0(int maxm); // regular
    void add00(int maxm); // for CDR
    void add000(int maxm); // for CDR with CASE operator

    void MakeKernels(const std::set<integralCart> & set0, std::set<ikernel> & kset);
    void AppendSet(const std::set<integralCart> & set0);

    void addRxcc(const std::set<integralCart> & set0 ,std::set<integralCart> & setf, bool ABx, bool CDx, bool ACx);
    void addRycc(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool ABy, bool CDy, bool ACy);
    void addRzcc(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool ABz, bool CDz, bool ACz);

    void addCTEbx(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool ABx);
    void addCTEby(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool ABy);
    void addCTEbz(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool ABz);

    void addCTEkx(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool CDx);
    void addCTEky(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool CDy);
    void addCTEkz(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool CDz);

    void addHRRbx(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool ABx);
    void addHRRby(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool ABy);
    void addHRRbz(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool ABz);

    void addHRRkx(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool CDx);
    void addHRRky(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool CDy);
    void addHRRkz(const std::set<integralCart> & set0, std::set<integralCart> & setf, bool CDz);

    void addSphA (const std::set<integralCart> & set0, std::set<integralCart> & fset);
    void addSphB (const std::set<integralCart> & set0, std::set<integralCart> & fset);
    void addSphC (const std::set<integralCart> & set0, std::set<integralCart> & fset);
    void addSphD (const std::set<integralCart> & set0, std::set<integralCart> & fset);

    void addInitial(std::set<integralCart> & setf);
};


struct LineNode {
    LineNode * prev;
    LineNode * next;
    UI16 key;

    LineNode() {
        prev = next = NULL;
        key = -1;
    }
};

class LinesInCache {
  public:
    UI16 MaxSize;
    UI16 nElements;

    LineNode * first;
    LineNode * last;

    std::map<UI32, LineNode*> Lines;

    //don't use it outside
    void del(UI32 Line) {
        LineNode * todel = Lines[Line];

        Lines.erase(Line);

        if (todel->prev !=NULL) todel->prev->next = todel->next;
        else                    first = first->next;
        if (todel->next !=NULL) todel->next->prev = todel->prev;
        else                    last = last->prev;

        delete todel;

        --nElements;
    }


    LinesInCache() {
        nElements = 0;
        MaxSize   = 512;
        first = NULL;
        last = NULL;
    }


    void push(UI32 line) {
        //if Line in cache: remove from list
        if (Lines.count(line)==1)
            del(line);

        //delete last element
        if (nElements==MaxSize)
            del(last->key);

        //push the new element
        LineNode * node = new LineNode;
        node->key = line;

        if (first==NULL) {
            first = node;
            last  = node;
        }
        else {
            first->prev = node;
            node->next = first;
            first = node;
        }

        Lines[line] = node;
        ++nElements;
    }

    bool incache(UI32 line) const {
        std::map<UI32, LineNode*>::const_iterator it;

        it = Lines.find(line);

        return (it!=Lines.end());
    }

    bool notincache(UI32 line) {
        std::map<UI32, LineNode*>::const_iterator it;

        it = Lines.find(line);

        if (it!=Lines.end()) return false;

        push(line);

        return true;
    }


};

#endif
