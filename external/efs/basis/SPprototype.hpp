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


#ifndef __ShellPairPrototype__
#define __ShellPairPrototype__

#include "../defs.hpp"
#include "../low/rtensors.hpp"

class GTO;
class GTObasis;
class ShellPairPrototype;

class BatchConstantsPair {
  public:

    double Na[maxK][maxJ];
    double Nb[maxK][maxJ];

    //double U[maxK][LM4];
    //double V[maxK][maxK][LM6];

    double r[maxK][maxK];
    double k[maxK][maxK];

    double K[maxK][maxK];

    double b1[maxK];
    double b2[maxK];

    double ika[maxK];
    double ikb[maxK];

    double UV[32]; //32>=6*LMAX+1

    BatchConstantsPair();
    ~BatchConstantsPair();

    void From(const ShellPairPrototype & AB);
};

class PrimitiveSet {
  public:
    UI8 nKa[maxK];
    double maxD2;
    UI16 nK2;
    UI8 nKb;

    UI8 align[21];
};
//64 bytes


class ShellPairPrototype {
  public:
    double Na[maxJ][maxK];  //pesos/coeficientes de la distribucion
    double Nb[maxJ][maxK];  //pesos/coeficientes de la distribucion

    double Nap[maxJ][maxK]; //pesos/coeficientes de la distribucion
    double Nbp[maxJ][maxK]; //pesos/coeficientes de la distribucion

    double ka[maxK];
    double kb[maxK];

    double kr[maxK][maxK];

    double kmin;

    PrimitiveSet Psets[maxK2];

    BatchConstantsPair BCP;

    UI16  atype1; //tipo de atomo
    UI16  atype2;
    UI8   nb1;    //numero de funcion
    UI8   nb2;
    UI32  num;

    UI8 l1;      //momento(s) angular(es) de a
    UI8 l2;      //momento(s) angular(es) de b

    UI8 ll1;      //
    UI8 ll2;      //

    UI16 f1;      //posicion de la funcion dentro del tensor (cartesianas)
    UI16 f2;      //posicion de la funcion dentro del tensor (cartesianas)

    UI16 fs1;     //posicion de la funcion dentro del tensor (esfericas)
    UI16 fs2;

    UI16 wb1;      //numero total de funciones 1
    UI16 wb2;      //numero total de funciones 2

    UI16 Kt;       //grado de contraccion total

    UI8 Ka;       //grado de contraccion a
    UI8 Kb;       //grado de contraccion b

    UI8 Ja;       //numero de primitivas a
    UI8 Jb;       //numero de primitivas b

    bool samef;    //producto de la misma funcion del msmo atomo
    bool sameatom; //producto de dos funciones del msmo atomo
    bool inverted; //whether or not the product is inverted

    bool SP1;
    bool SP2;


    ShellPairPrototype();
    void BasisProd (const GTO & A);
    void BasisProd (const GTO & A, const GTO & B);
};

class ShellPairs {
  public:
    ShellPairPrototype **** GP;
    ShellPairPrototype *** GP2;
    ShellPairPrototype ** GP3;
    ShellPairPrototype * GP4;

    ShellPairPrototype ** GPsame;
    ShellPairPrototype * GPsame2;

    int nGP;
    int nGPsame;

    Array<UI16> AType;

    void GenerateFrom (const Array<UI16> & AtomTypes, const GTObasis & basis);
};


#endif
