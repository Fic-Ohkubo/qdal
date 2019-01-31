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


#ifndef __NEW_SPARSE__
#define __NEW_SPARSE__

#include "../defs.hpp"
#include "../basis/GTO.hpp"
#include "../low/rtensors.hpp"


class tensor2;
class tensor;
class rAtom;
class rAtomPrototype;

class sparsetensorpattern {
    friend class tensor2;
    friend class tensor;

  public:
    //first level, atom blocks
    static r2tensor<UI32> a2pos; //posiciones en el array lineal para cada atomo
    //second level, GC function blocks
    static r2tensor< r2tensor<UI32> > f2pos; //offsets for the block of GC functions of a given block respective to the beginning of the atomair block
    //third level, individual function blocks
    static r2tensor< r2tensor<r2tensor<UI32> > > j2pos; //offsets for a function pair respective to the beginning of the GC block block

    //first level, atom blocks
    static UI32 * a1pos; //posiciones en el array lineal para cada atomo
    //second level, GC function blocks
    static r2tensor<UI32> f1pos; //posiciones en el array lineal para cada atomo
    //third level, individual function blocks
    static r2tensor<r2tensor<UI32> > j1pos; //posiciones en el array lineal para cada atomo


    static UI16 nids;
    static UI16 * ids;
    static UI16 * nfs;          //number of functions for the given atom

    static UI16 * alen;         //total length of the atom block (total number of basis functions)
    static r2tensor<UI16> flen; //total length of a given function (at, nf)
    static r2tensor<UI8>  js;   //j's
    static r2tensor<UI8>  ms;   //ms

    //one has to store which 'element' every atom is so to avoid wasting so much space
    static UI32 natoms;
    static UI32 nbasis;
    static UI32 ntotal;
    static UI32 ntotal2;

  public:

    static void SetBlocks(const rAtom * Atoms, int nAtoms, const r1tensor<rAtomPrototype*> & AtomTypeList);
    static UI32 GetOffset(UI32 at1, UI32 at2);
    static UI32 GetOffset(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2);
    static UI32 GetOffset(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2, UI16 nj1, UI16 nj2);
};

class Sparse:private sparsetensorpattern {
    friend class tensor2;
    friend class tensor;
    friend class Fock2e;

  private:
    double * values;

  public:

    Sparse() {
        values   = NULL;
    }

    ~Sparse() {
        delete[] values;
    }

    void Set() {
        delete[] values;
        values = new double[ntotal2];
    }

    void zeroize() {
        for (int i=0; i<ntotal2; ++i) values[i] = 0;
    }

    inline double * operator()(UI32 at1, UI32 at2) {
        return values + a2pos(at1, at2);
    }

    inline const double * operator()(UI32 at1, UI32 at2) const {
        return values + a2pos(at1, at2);
    }

    inline double * operator()(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2) {
        UI16 id1 = ids[at1];
        UI16 id2 = ids[at2];

        return values + a2pos(at1, at2) + f2pos(id1,id2)(nf1,nf2);
    }

    inline const double * operator()(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2) const {
        UI16 id1 = ids[at1];
        UI16 id2 = ids[at2];

        return values + a2pos(at1, at2) + f2pos(id1,id2)(nf1,nf2);
    }

    inline double * operator()(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2, UI8 nj1, UI8 nj2) {
        UI16 id1 = ids[at1];
        UI16 id2 = ids[at2];
        UI32 mj1 = js[id1][nf1];
        UI32 mj2 = js[id2][nf2];
        UI32 mm1 = ms[id1][nf1];
        UI32 mm2 = ms[id2][nf2];

        return values + a2pos(at1, at2) + f2pos(id1,id2)(nf1,nf2) + (nj1*mj2 + nj2)*(mm1*mm2);
    }

    inline const double * operator()(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2, UI8 nj1, UI8 nj2) const {
        UI16 id1 = ids[at1];
        UI16 id2 = ids[at2];
        UI32 mj1 = js[id1][nf1];
        UI32 mj2 = js[id2][nf2];
        UI32 mm1 = ms[id1][nf1];
        UI32 mm2 = ms[id2][nf2];

        return values + a2pos(at1, at2) + f2pos(id1,id2)(nf1,nf2) + (nj1*mj2 + nj2)*(mm1*mm2);
    }

    inline double & operator()(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2, UI8 nj1, UI8 nj2, UI8 m1, UI8 m2) {
        UI16 id1 = ids[at1];
        UI16 id2 = ids[at2];
        UI32 mj1 = js[id1][nf1];
        UI32 mj2 = js[id2][nf2];
        UI32 mm1 = ms[id1][nf1];
        UI32 mm2 = ms[id2][nf2];

        return *(values + a2pos(at1, at2) + f2pos(id1,id2)(nf1,nf2) + (nj1*mj2 + nj2)*(mm1*mm2) + (m2*mm1 + m1) );
    }

    inline const double & operator()(UI32 at1, UI32 at2, UI16 nf1, UI16 nf2, UI8 nj1, UI8 nj2, UI8 m1, UI8 m2) const {
        UI16 id1 = ids[at1];
        UI16 id2 = ids[at2];
        UI32 mj1 = js[id1][nf1];
        UI32 mj2 = js[id2][nf2];
        UI32 mm1 = ms[id1][nf1];
        UI32 mm2 = ms[id2][nf2];

        return *(values + a2pos(at1, at2) + f2pos(id1,id2)(nf1,nf2) + (nj1*mj2 + nj2)*(mm1*mm2)  + (m2*mm1 + m1));
    }


    Sparse & operator=(const tensor2 & rhs);

    Sparse & operator=(const Sparse & rhs);

    Sparse & operator=(const tensor & rhs);

    Sparse & operator+=(const Sparse & rhs);

    Sparse & operator-=(const Sparse & rhs);

    Sparse & operator*=(double rhs);

    Sparse & operator/=(double rhs);

};

//extern sparsetensorpattern SparseTensorPattern;

#endif

