#ifndef __GTO__
#define __GTO__

#include "../defs.hpp"

// GTO 'elemental' shell
// gaussian exponents must be sorted from smallest (more diffuse) to largest (more localized)
struct GTO {
    double k [maxK];        // gaussian exponents
    double N [maxJ][maxK];  // contraction weights
    double Np[maxJ][maxK];  // contraction weights (only SP batches; not used)

    //float ik;   // inverso de la constante mas baja (para prescreening)

    UI8 l;      // angular moment
    UI8 K;      // contraction degree
    UI8 J;      // number of general contraction functions

    bool SP;   // is an SP shell?

    GTO();
    void Normalize();
    GTO & operator=(const GTO & rhs);
	~GTO ();
};


#endif
