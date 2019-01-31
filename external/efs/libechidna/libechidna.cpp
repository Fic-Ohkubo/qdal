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

#include <set>
#include "libechidna.hpp"
#include "../math/angular.hpp"
using namespace LibAngular;




DMD     Echidna::dmd;

MessagePlotter Echidna::EMessenger;
MessagePlotter Echidna::EBenchmarker;
MessagePlotter Echidna::EDebugger;

//default values
size_t Echidna::MAXMEM    =   4*Gword;
size_t Echidna::MEMBUFFER = 512*Mword;




Echidna::Echidna() {
    LibQuimera::Quimera::QMessenger.SetParent(&EMessenger);
    LibQuimera::Quimera::QBenchmarker.SetParent(&EBenchmarker);

    //EMessenger.SetModule   ("Echidna module");
    //EBenchmarker.SetModule ("Echidna module benchmarker");
    //EDebugger.SetModule    ("Echidna module debugger");
}

void Echidna::Init (const rAtom * AtomList, int nAtomsL) {

    EBenchmarker.SetFile("echidna.benchmark.out");
    EDebugger.SetFile   ("echidna.debug.out");

    // makes a hard copy of the atom list
    nAtoms = nAtomsL;
    Atoms  = new rAtom[nAtoms];

    for (int at1=0; at1<nAtoms; ++at1) {
        Atoms[at1] = AtomList[at1];
    }

    // generates an index of elements in the molecule
    std::map<UI16, rAtomPrototype*> ATList;
    std::map<UI16, rAtomPrototype*>::iterator it;

    for (int at1=0; at1<nAtoms; ++at1) {
        ATList[Atoms[at1].rAP->RId] =  Atoms[at1].rAP; //Atoms[at1].rAP->DId;
    }

    AtomTypeList.set(ATList.size());

    for (it=ATList.begin(); it!=ATList.end(); ++it) {
        UI16             rid = it->first;
        rAtomPrototype * pAt = it->second;

        AtomTypeList[rid] = pAt;
    }

    // initializes the sparsetensor block pattern
    SparseTensorPattern.SetBlocks(Atoms, nAtoms, AtomTypeList);


    // generates the atom pair prototypes
    Echidna::EMessenger << "Generating elemental basis pairs";

    LibQuimera::InitLibQuimera(); //needs LibAngular initialized

    Echidna::EMessenger.Push(); {
        APprototypes.GenerateFrom(AtomTypeList);
    } Echidna::EMessenger.Pop(true);
}

void Echidna::ComputeAtomPairs (const r2tensor<int> & AtomPairs, double logGDO) {

    // counts the total number of atom interactions
    TotAtomPairs = 0;

    for (int i=0; i<AtomPairs.N(); ++i) {
        TotAtomPairs += AtomPairs.M(i);
    }

    // generates the atom pair information
    AtomProductLists.From(AtomTypeList, APprototypes, AtomPairs, Atoms, nAtoms, logGDO);
}

void Echidna::InitFock2e (double Xscaling) {

    LibQuimera::Quimera::QBenchmarker.SetFile("quimera.benchmark.out");

    LibQuimera::InitLibQuimera();

    // initializes the Fock 2e solver
    WW.Init(TotAtomPairs, AtomProductLists, SparseTensorPattern);
    WW.SetXS(Xscaling);

    // computes the cauchy-schwarz parameters of the surviving ShellPairs
    WW.CalcCauchySwarz();
}


void Echidna::FockUpdate (tensor2 & F, const tensor2 & AD, double Dthresh, double Sthresh, bool SplitJX) {

    // computes the 2e part of the fock matrix
    WW.FockUpdate(F, AD, Dthresh, Sthresh, SplitJX);
}

int Echidna::nShellPairs () const {
    return AtomProductLists.nGTOps;
}

const ShellPair * Echidna::GetShellPairs() const {
    return AtomProductLists.SPpool;
}
