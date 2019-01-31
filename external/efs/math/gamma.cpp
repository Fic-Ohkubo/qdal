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


#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "../math/gamma.hpp"
#include "../defs.hpp"
using namespace std;


static inline void LenzUpdate(double & A, double & Ap, const double a, const double b) {
    double tmp = A;
    A = b*A + a*Ap;
    Ap = tmp;
}

//computational chemists' "lower gamma function"
//Table[Integrate[t^(2n) exp[-x t^2],{x,0,1}],{n,0,...}]
static inline void calcF(double * F, int n, double x) {
    double expx = exp(-x);

    //calcula gamma(n + 1/2, x) / x^(n+1/2) y el resto por recurrencia descendente
    if (x < double(n) + 1.5) {
        double ac  = 1./(double(n)+0.5);
        double sum = ac;

        //loop until machine precission convergence
        for (int i=1; (sum+ac)!=sum; ++i) {
            ac *= (x)/(double(n+i)+0.5);
            sum +=ac;
        }

        F[n] = 0.5 * sum;

        //downward recursion
        for (int i=n-1; i>=0; --i)
            F[i] = (2*x*F[i+1] + 1.)/double(2*i+1);

        for (int i=n; i>=0; --i)
            F[i] *= expx;
    }

    //calcula F0[x] a partir de la fraccion continua y el resto por recurrencia ascendente
    else {
        double Ap = 1;
        double A  = 0; //b0=0

        double Bp = 0;
        double B  = 1;

        LenzUpdate(A,Ap,1.,x);
        LenzUpdate(B,Bp,1.,x);

        //OPT: optimal loop unrolling
        //OPT: function inlining
        for (int i=0; i<10; ++i) {
            LenzUpdate(A,Ap,double(i)+0.5,1.);
            LenzUpdate(B,Bp,double(i)+0.5,1.);

            LenzUpdate(A,Ap,double(i+1),x);
            LenzUpdate(B,Bp,double(i+1),x);
        }

        F[0] = sqrt(PI/(4*x)) -0.5 * expx * (A/B);

        //upward recursion
        for (int i=0; i<n; ++i)
            F[i+1] = (double(2*i+1)*F[i] - expx) / (2*x);
    }
}


namespace LibIGamma {

    GammaFunction::GammaFunction() {
        gamma_table = NULL;
        gamma_table_vals = NULL;
    }

    GammaFunction::~GammaFunction() {
        free(gamma_table);
        free(gamma_table_vals);
    }

    //nueva funcion gamma
    void GammaFunction::CalcGammas(double * F, int n, double x) const{
        //busca el valor mas cercano de F
        //cout << vg_max << " " << vg_step << endl;
        if (x>vg_max-vg_step) {
            F[0] = 0.5*sqrt(PI/x);

            double ix = 0.5/x;
            for (int i=1; i<=n; ++i) {
                F[i] = F[i-1]*ix*double(2*i-1);
            }
        }
        else {
            double p = ivg_step*(x-vg_min);
            int pos = int(p+0.5);
            double x0 = vg_min + vg_step*double(pos);
            double Ax = x0-x;

            double Axn[NEXP];
            Axn[0] = 1;
            for (int i=1; i<NEXP; ++i) {
                Axn[i] = (Ax * Axn[i-1])/double(i);
            }

            {
                double sum = 0;
                for (int j=0; j<NEXP; ++j)
                    sum += gamma_table[pos][n+j+1] * Axn[j];
                F[n] = sum;
            }

            if (n>0)
            {
                double expx = 1;

                for (int i=1; i<NEXP; ++i) {
                    expx += Axn[i];
                }

                expx *= gamma_table[pos][0];

                for (int i=n-1; i>=0; --i)
                    F[i] = (2*x*F[i+1] + expx)/double(2*i+1);
            }


        }
    }


    void GammaFunction::InitIncompleteGammaTable() {
        if (gamma_table!=NULL) return;

        //init arrays and array positions
        posix_memalign((void**)(&gamma_table)     , CACHE_LINE_SIZE, NVIG     *sizeof(double*));
        posix_memalign((void**)(&gamma_table_vals), CACHE_LINE_SIZE, NVIG*NGD2*sizeof(double));

        gamma_table      = new double*[NVIG];
        gamma_table_vals = new double[NGD*NVIG];

        for (int i=0; i<NVIG; ++i) gamma_table[i] = gamma_table_vals + i*NGD;

        double Av = (vg_max-vg_min) / double (NVIG);
        for (int i=0; i<NVIG; ++i) {
            double v = vg_min + i*Av;

            gamma_table[i][0] = exp(-v);
            calcF(gamma_table[i]+1, 4*LMAX+1+NEXP, v);
        }
    }

    void GammaFunction::InitIncompleteGammaTable(int L) {
        if (gamma_table!=NULL) return;

        //init arrays and array positions
        posix_memalign((void**)(&gamma_table)     , CACHE_LINE_SIZE, NVIG     *sizeof(double*));
        posix_memalign((void**)(&gamma_table_vals), CACHE_LINE_SIZE, NVIG*NGD2*sizeof(double));

        for (int i=0; i<NVIG; ++i) gamma_table[i] = gamma_table_vals + i*NGD2;

        double Av = (vg_max-vg_min) / double (NVIG);
        for (int i=0; i<NVIG; ++i) {
            double v = vg_min + i*Av;

            double vv[32];
            calcF(vv, L+7, v);

            //first value is the exponential
            gamma_table[i][0] = exp(-v);
            //rest of values are the derivatives up to order 6
            for (int j=0; j<7; ++j) gamma_table[i][j+1] = vv[L+j];
        }
    }


    GammaFunction IncompleteGamma;
    GammaFunction IncompleteGammas[4*LMAX+3];


    void InitIncompleteGammas() {
        IncompleteGamma.InitIncompleteGammaTable();

        for (int l=0; l<=4*LMAX+2; ++l) IncompleteGammas[l].InitIncompleteGammaTable(l);
    }



}

