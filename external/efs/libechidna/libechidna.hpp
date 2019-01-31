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





    ECHIDNA library
    2-Electron Contributions through Hybrid-parallelized, Integral-Direct, Novel Algorithms


    2013.02.19, Jaime Axel Rosal Sandberg
    * initial version, refactored from previous code

*/

#ifndef __LIB_ECHIDNA__
#define __LIB_ECHIDNA__

#include <cstdlib>
#include "../integrals/fock2e.hpp"
#include "../integrals/newcontractions.hpp"
#include "../integrals/atomprod.hpp"
#include "../math/newsparse.hpp"
#include "../basis/GTO.hpp"

class rAtom;

class rAtomPrototype {
  public:
    std::string name; // name
    double mass;      // mass
    double Z;         // charge

    UI16 DId;           // Id in the original definition list
    UI16 RId;           // Id in the reduced element list
    r1tensor<GTO>  basis;

    rAtom *  Atoms;
    UI32    nAtoms;
};

class rAtomPrototypes {
  public:
    r1tensor<rAtomPrototype> prototypes;
};

class rAtom {
  public:
    point c;
    rAtomPrototype * rAP;
}; //256 bits


class Echidna {
    friend class BatchInfo;
    friend class BatchEvaluator;
    static DMD dmd;

  public:
    static MessagePlotter EMessenger;
    static MessagePlotter EBenchmarker;
    static MessagePlotter EDebugger;

    static size_t MAXMEM;
    static size_t MEMBUFFER;

  private:
    rAtom                   * Atoms;
    int                       nAtoms;
    Fock2e                    WW;
    sparsetensorpattern       SparseTensorPattern;
    AtomProdPrototypes        APprototypes;
    APlists                   AtomProductLists;
    r1tensor<rAtomPrototype*> AtomTypeList;
    UI64                      TotAtomPairs;

  public:
    Echidna();

    void Init                          (const rAtom * Atoms, int nAtoms);
    void ComputeAtomPairs              (const r2tensor<int> & AtomPairs, double logGDO);
    void InitFock2e                    (double Xscaling=1.);
    void FockUpdate                    (tensor2 & F, const tensor2 & AD, double Dthresh, double Sthresh, bool SplitJX = false);

    int nShellPairs () const;
    const ShellPair * GetShellPairs() const;

};


#endif
