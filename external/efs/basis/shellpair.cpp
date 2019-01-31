/*
    shellpair

    definicion y calculo de productos de GTOs
*/


#include "../math/affine.hpp"
#include "../basis/SPprototype.hpp"
#include "../basis/shellpair.hpp"


using namespace std;


double ShellPair::getW(int b, int a) const {
    return WW[b * GP->Ka + a];
}


int ShellPair::getKb() const {
    return GP->Psets[nK2-1].nKb;
}

int ShellPair::getKa(int Kb) const {
    return GP->Psets[nK2-1].nKa[Kb];
}

int ShellPair::getLa() const {
    return GP->l1;
}

int ShellPair::getLb() const {
    return GP->l2;
}

int ShellPair::getFa() const {
    return GP->f1;
}

int ShellPair::getFb() const {
    return GP->f2;
}

enum SPflags {SAMEAT=1, SAMEF=2, INVERTED=4};

bool ShellPair::sameF() const {
    return (flags&SAMEF == SAMEF);
}

bool ShellPair::sameAtom() const {
    return (flags&SAMEAT == SAMEAT);
}

bool ShellPair::Inverted() const {
    return (flags&INVERTED == INVERTED);
}

//producto entre dos funciones base g1 y g2 centradas en c1 y c2. Los productos menores a minover se descartan
void ShellPair::Form(const ShellPairPrototype & prototype, const point & cc1, const point & cc2, double * ww, double logminover) {

    flags = (prototype.inverted?INVERTED:0);

    WW = ww;

    point A = cc1;
    point B = cc2;

    vector3 v = B-A;

    double r2 = v*v;
    double r  = sqrt(r2);

    v/=r;

    GP = &prototype;

    nK2 = 0;

    for (int b=0; b<prototype.Kb; ++b) {
        int a;
        for (a=0; a<prototype.Ka; ++a) {
            if (prototype.kr[a][b] * r2 > logminover) break; // should use overlap as a better estimator
        }

        nK2 += a;
    }

    //store all anyhow
    for (int b=0; b<prototype.Kb; ++b) {
        for (int a=0; a<prototype.Ka; ++a) {
            double Kab = exp(-prototype.kr[a][b] * r2);
            //W->W[b][a] =
            *ww = Kab;
            ++ww;
        }
    }

}

void ShellPair::Form(const ShellPairPrototype & prototype, const point & cc1, double * ww) {

    flags = SAMEAT + (prototype.samef?SAMEF:0) + (prototype.inverted?INVERTED:0);


    WW = ww;

    GP = &prototype;

    nK2 = prototype.Ka * prototype.Kb;


    for (int b=0; b<prototype.Kb; ++b) {
        for (int a=0; a<prototype.Ka; ++a) {
            *ww = 1;
            ++ww;
        }
    }
}



