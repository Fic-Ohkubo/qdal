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


#ifndef __GAMMA__
#define __GAMMA__
#include "../defs.hpp"

namespace LibIGamma {

    const double PI3  = 31.0062766802998202;   // = PI^3
    const double PI54 = 34.986836655249725693; // = 2*PI^(5/2)

    const float fPI3  = 31.0062766802998202;   // = PI^3
    const float fPI54 = 34.986836655249725693; // = 2*PI^(5/2)

    // everything is made public since
    // the inlined and specialized algorithms perform better

    // since gamma is evaluated by downward recursion, it is in principle better to have
    // one GammaFunction for each Fm[z] and adjust the number of points in the grid, grid step
    // and critical gamma acordingly; this will be eventually implemented

    const int NEXP = 4;
    const int NVIG = 1024;
    const int NGD = 32; // >= 4*LMAX + 1 + NEXP + 2;
    const int NGD2 = 8;

    const double vg_min  = 0;
    const double vg_max  = 32;
    const double vg_step  = 32./1024.;  //(vg_max-vg_min)/double(NVIG);
    const double ivg_step = 1024./32.;    //double(NVIG)/(vg_max-vg_min);

    const float fvg_min  = 0;
    const float fvg_max  = 32;
    const float fvg_step  = 32./1024.;  //(vg_max-vg_min)/double(NVIG);
    const float fivg_step = 1024./32.;    //double(NVIG)/(vg_max-vg_min);



    class GammaFunction {
      public:

      public:
        double ** gamma_table;
        double *  gamma_table_vals;
      public:
        GammaFunction();
        ~GammaFunction();
        void InitIncompleteGammaTable();
        void InitIncompleteGammaTable(int L);
        void CalcGammas(double * F, int n, double x) const;
    };

    extern GammaFunction IncompleteGamma;
    extern GammaFunction IncompleteGammas[4*LMAX+3];

    void InitIncompleteGammas();
}


#endif
