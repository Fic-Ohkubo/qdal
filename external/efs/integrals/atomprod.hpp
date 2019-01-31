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


#ifndef __ATOMPROD__
#define __ATOMPROD__

#include "../low/rtensors.hpp"
#include "../math/affine.hpp"
#include "rotations.hpp"


//forward declaration to avoid cyclic dependencies
class Atom;
class rAtom;
class AtomProdPrototype;
class AtomProdPrototypes;
class GTObasis;
class ShellPairPrototype;
class ShellPair;
class Nucleus;
class rAtomPrototype;


class AtomProdPrototype {
  public:
    float * Rs;
    ShellPairPrototype * Gs;

    int maxL;

    short int nGPt; //256*256
    short int num; //256*256

    short int wa1; //wsize of atom 1
    short int wa2; //wsize of atom 2

    char Atype1; //atom type (256 should be enough)
    char Atype2;

    char nGP1;   //number of primitives (256 should be enough)
    char nGP2;

    bool same;     //same atom
};

class AtomProdPrototypes {
  public:
    AtomProdPrototype ** AP;
    AtomProdPrototype  * APsame;
    AtomProdPrototype  * APpool;

    ShellPairPrototype * GPall;

    int nGP;
    int nGPsame;
    int nAP;
    int nA;

    AtomProdPrototypes();
    void GenerateFrom (const r1tensor<rAtomPrototype*> & AtomTypes);
    ~AtomProdPrototypes();
};


class AtomProd {
  public:
    AtomProdPrototype * APprototype; //pointer to the atom prototype

    ShellPair ** pSPs;     //pointers to the shell pairs in APlists

    point A;
    point B;
    vector3 ABv;

    vector3 vx;
    vector3 vy;
    vector3 vz;

    RotationMatrix RM;

    double r2; // distancia entre los atomos
    double norm;
    int at1;
    int at2;
    int nAP;    //position of the atom pair in the prototype list in APlists
    int nGTOps; //number of non-void GTO products/shell pairs of the atom pair
    bool reverse_prototype;

    void MakeRotation();
};

class APbatch {
  public:
    AtomProdPrototype * APP; //atom pair prototype for this batch
    Array<AtomProd> APlist;  //list of atom pairs
    Array<int> nAPS;         //number of atom pairs for which the i-th shell has not been prescreened
    Array<ShellPair *> SP;   //pointer to the beginning of the list of the i-th shell pairs
};

class APlists {
  public:
    Array< APbatch > AtomPairBatch;
    Array< APbatch > AtomSameBatch;

    ShellPair * SPpool;
    double   * WWpool;

    int nGTOps;

    APlists();
    void From (const r1tensor<rAtomPrototype*> & AtomTypes, const AtomProdPrototypes & APprototypes, const r2tensor<int> & IL, const rAtom * Atoms, int nAtoms, double lmo);

    ~APlists();
};



#endif
