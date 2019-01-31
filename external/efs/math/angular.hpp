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


#ifndef __ANGULAR__
#define __ANGULAR__

#include "../defs.hpp"

namespace LibAngular {

    // terms in the expansion of a spherical harmonic
    struct SHTerm {
        double cN; //peso del armonico normalzado
        char   nc; //numero del termino correspondiente en cartesianas
        char   nx; //potencia en x
        char   ny; //potencia en y
        char   nz; //potencia en z
    };

    //definicion de un armónico esférico (real) en coordenadas cartesianas
    struct SphHar {
        int l; //momento angular
        int m; //m

        int nps;     //numero de primitivas cartesianas para construir el harmonico esferico
        SHTerm T[6]; //terminos individuales
                     //allow only up to 6 at the moment, which is enough for up to h-functions;
                     //this avoids unnecesary memory fragmentation

        //inicia los arrays de terminos
        void SetNPS(int n);
        void SetT(int t, double c, int x, int y, int z);
    };

    //this is used for metaprogrammed routines which require a compile-time static value
    template<int N> struct nmCT {
        enum {n = ((N+2)*(N+1))/2};
    };

    //last one is for SP shells
    const int nmS[LMAX+2] = {1,3,5,7,9,11,   4};
    const int nmC[LMAX+2] = {1,3,6,10,15,21, 4};

    extern int ExyzC[LMAX+2][LMAX+2][MMAX][MMAX];
    extern int ExyzT[LMAX+2][LMAX+2];

    extern SHTerm CartList[LMAX+2][MMAX];
    extern SphHar SHList  [LMAX+2][2*LMAX+1];

    extern double SHrot [LMAX+1][2*LMAX+1][2*LMAX+1];
    extern double SHrotI[LMAX+1][2*LMAX+1][2*LMAX+1];


    // inicia las tablas de esfericos harmonicos y las seminormaliza
    void InitSHList();
    void InitCartList();

    // capital letter representing the angular momentum
    char L2S (int  L);
    int S2L (char S);

    //gamma[n+1/2]
    const double hGamma[] =
    {1.77245385090551603, 0.886226925452758014, 1.32934038817913702,
    3.32335097044784255, 11.6317283965674489, 52.3427777845535202,
    287.885277815044361, 1871.25430579778835, 14034.4072934834126,
    119292.461994609007, 1.13327838894878557e6,
     1.18994230839622485e7, 1.36843365465565857e8,
     1.71054206831957322e9, 2.30923179223142384e10,
     3.34838609873556457e11, 5.18999845304012508e12,
     8.56349744751620639e13, 1.49861205331533612e15,
     2.77243229863337182e16, 5.40624298233507504e17,
     1.10827981137869038e19, 2.38280159446418433e20,
     5.36130358754441473e21, 1.25990634307293746e23,
     3.08677054052869678e24, 7.87126487834817680e25,
     2.08588519276226685e27, 5.73618428009623384e28,
     1.63481251982742664e30, 4.82269693349090860e31,
     1.47092256471472712e33};
}


#endif


