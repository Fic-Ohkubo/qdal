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
#include "../defs.hpp"
#include "../math/angular.hpp"
#include "../math/affine.hpp"
using namespace std;


double AngularIntegral(int nx, int ny, int nz);

namespace LibAngular {
    int ExyzC[LMAX+2][LMAX+2][MMAX][MMAX];
    int ExyzT[LMAX+2][LMAX+2];

    SHTerm CartList[LMAX+2][MMAX];
    SphHar SHList  [LMAX+2][2*LMAX+1];

    double SHrot [LMAX+1][2*LMAX+1][2*LMAX+1];
    double SHrotI[LMAX+1][2*LMAX+1][2*LMAX+1];
}



void LibAngular::SphHar::SetNPS(int n) {
    nps = n;
    //T = new SHTerm[n];
}

void LibAngular::SphHar::SetT(int t, double c, int x, int y, int z) {
    T[t].cN = c;
    T[t].nx = x;
    T[t].ny = y;
    T[t].nz = z;

    //find the term
    //int l = x+y+z;
    int mm = nmC[l];

    for (int m=0; m<mm; ++m) {
        if ( (CartList[l][m].nx == x) && (CartList[l][m].ny == y) && (CartList[l][m].nz == z) ) {
            T[t].nc = m;
            return;
        }
    }
    //cout << "Error in angular: " << l << "  " << x << " " << y << " " << z << endl;
    throw 23;
}


//inicia las tablas de armonicos esfericos y las seminormaliza
void LibAngular::InitSHList() {

    for (int l=0; l<=LMAX; ++l) {

        for (int m=0; m<=2*l; ++m){
            SHList[l][m].l = l;
            SHList[l][m].m = m;
        }
    }

    //SP
    {
        for (int m=0; m<4; ++m){
            SHList[LSP][m].l = LSP;
            SHList[LSP][m].m = m;
        }
    }



    //inicia las tablas con datos; en adelante se usara un algoritmo para proyectar los
    //armonicos esfericos en coordenadas cartesianas

    //s
    {
        SHList[0][0].SetNPS(1);
        SHList[0][0].SetT(0, 1, 0,0,0);
    }

    //p
    if (LMAX > 0) {
        SHList[1][0].SetNPS(1); //y
        SHList[1][0].SetT(0, 1, 0,1,0);

        SHList[1][1].SetNPS(1); //z
        SHList[1][1].SetT(0, 1, 0,0,1);

        SHList[1][2].SetNPS(1); //x
        SHList[1][2].SetT(0, 1, 1,0,0);
    }

    //sp
    if (LMAX > 0) {
        SHList[LSP][0].SetNPS(1);
        SHList[LSP][0].SetT(0, 1, 0,0,0);

        SHList[LSP][1].SetNPS(1); //y
        SHList[LSP][1].SetT(0, 1, 0,1,0);

        SHList[LSP][2].SetNPS(1); //z
        SHList[LSP][2].SetT(0, 1, 0,0,1);

        SHList[LSP][3].SetNPS(1); //x
        SHList[LSP][3].SetT(0, 1, 1,0,0);
    }


    //d
    if (LMAX > 1) {
        SHList[2][0].SetNPS(1); //2i xy
        SHList[2][0].SetT(0, 1, 1,1,0);

        SHList[2][1].SetNPS(1); //yz
        SHList[2][1].SetT(0, 1, 0,1,1);

        SHList[2][2].SetNPS(3); //2z^2-x^2-y^2 = 3z^2 - r^2
        SHList[2][2].SetT(0,  2, 0,0,2);
        SHList[2][2].SetT(1, -1, 2,0,0);
        SHList[2][2].SetT(2, -1, 0,2,0);

        SHList[2][3].SetNPS(1); //xz
        SHList[2][3].SetT(0, 1, 1,0,1);

        SHList[2][4].SetNPS(2); //x^2-y^2
        SHList[2][4].SetT(0,  1, 2,0,0);
        SHList[2][4].SetT(1, -1, 0,2,0);
    }

    //f
    if (LMAX > 2) {
        // i (3x^2 - y^2) y
        SHList[3][0].SetNPS(2);
        SHList[3][0].SetT(0, -1, 0,3,0);
        SHList[3][0].SetT(1,  3, 2,1,0);

        //2i z xy
        SHList[3][1].SetNPS(1);
        SHList[3][1].SetT(0, 1, 1,1,1);

        // i (4z^2-x^2-y^2) y
        SHList[3][2].SetNPS(3);
        SHList[3][2].SetT(0,  4, 0,1,2);
        SHList[3][2].SetT(1, -1, 2,1,0);
        SHList[3][2].SetT(2, -1, 0,3,0);

        // (2z^2-3x^2-3y^2) z m=0
        SHList[3][3].SetNPS(3);
        SHList[3][3].SetT(0,  2, 0,0,3);
        SHList[3][3].SetT(1, -3, 2,0,1);
        SHList[3][3].SetT(2, -3, 0,2,1);

        // x(4z^2-x^2-y^2)
        SHList[3][4].SetNPS(3);
        SHList[3][4].SetT(0,  4, 1,0,2);
        SHList[3][4].SetT(1, -1, 3,0,0);
        SHList[3][4].SetT(2, -1, 1,2,0);

        // z(x^2-y^2)
        SHList[3][5].SetNPS(2);
        SHList[3][5].SetT(0,  1, 2,0,1);
        SHList[3][5].SetT(1, -1, 0,2,1);

        // x^3 - 3xy^2
        SHList[3][6].SetNPS(2);
        SHList[3][6].SetT(0,  1, 3,0,0);
        SHList[3][6].SetT(1, -3, 1,2,0);
    }

    //g
    if (LMAX > 3) {
        //xy (x^2 - y^2)
        SHList[4][0].SetNPS(2);
        SHList[4][0].SetT(0,  1, 3,1,0);
        SHList[4][0].SetT(1, -1, 1,3,0);

        //z y (3x^2 - y^2)
        SHList[4][1].SetNPS(2);
        SHList[4][1].SetT(0, 3, 2,1,1);
        SHList[4][1].SetT(1,-1, 0,3,1);

        //x y (7z^2 - r^2)
        SHList[4][2].SetNPS(3);
        SHList[4][2].SetT(0,  6, 1,1,2);
        SHList[4][2].SetT(1, -1, 3,1,0);
        SHList[4][2].SetT(2, -1, 1,3,0);

        //y z (7z^2 -3r^2)
        SHList[4][3].SetNPS(3);
        SHList[4][3].SetT(0,  4, 0,1,3);
        SHList[4][3].SetT(1, -3, 2,1,1);
        SHList[4][3].SetT(2, -3, 0,3,1);

        //3x^4 + 6x^2y^2 + 3y^4 - 24x^2z^2 - 24y^2z^2 + 8z^4 = 35z^4 - 30z^2 r^2 + 3r^4
        //= 35z^4 - 30z^2 r^2 + 3r^4
        //m=0
        SHList[4][4].SetNPS(6);
        SHList[4][4].SetT(0,   3, 4,0,0);
        SHList[4][4].SetT(1,   6, 2,2,0);
        SHList[4][4].SetT(2,   3, 0,4,0);
        SHList[4][4].SetT(3, -24, 2,0,2);
        SHList[4][4].SetT(4, -24, 0,2,2);
        SHList[4][4].SetT(5,   8, 0,0,4);

        //x z (7z^2 - 3r^2)
        SHList[4][5].SetNPS(3);
        SHList[4][5].SetT(0,  4, 1,0,3);
        SHList[4][5].SetT(1, -3, 3,0,1);
        SHList[4][5].SetT(2, -3, 1,2,1);

        //(x^2-y^2)*(7z^2-r^2) = -x^4 + y^4 - 6y^2 z^2 + 6x^2 z^2
        SHList[4][6].SetNPS(4);
        SHList[4][6].SetT(0, -1, 4,0,0);
        SHList[4][6].SetT(1,  1, 0,4,0);
        SHList[4][6].SetT(2, -6, 0,2,2);
        SHList[4][6].SetT(3,  6, 2,0,2);

        //x z (x^2 - 3y^2)
        SHList[4][7].SetNPS(2);
        SHList[4][7].SetT(0,  1, 3,0,1);
        SHList[4][7].SetT(1, -3, 1,2,1);

        //x^4 - 6x^2y^2 + y^4
        SHList[4][8].SetNPS(3);
        SHList[4][8].SetT(0,  1, 4,0,0);
        SHList[4][8].SetT(1,  1, 0,4,0);
        SHList[4][8].SetT(2, -6, 2,2,0);
    }

    //h
    if (LMAX > 4) {
        //i (5 x^4 y - 10 x^2 y^3 + y^5)
        SHList[5][0].SetNPS(3);
        SHList[5][0].SetT(0,   1, 0,5,0);
        SHList[5][0].SetT(1, -10, 2,3,0);
        SHList[5][0].SetT(2,   5, 4,1,0);

        //i z (x^3 y - x y^3)
        SHList[5][1].SetNPS(2);
        SHList[5][1].SetT(0,  1, 3,1,1);
        SHList[5][1].SetT(1, -1, 1,3,1);

        //i (y^5 - 2 y^3 x^2 - 3 y x^4 - 8 y^3 z^2 + 24 y x^2 z^2)
        SHList[5][2].SetNPS(5);
        SHList[5][2].SetT(0,  1, 0,5,0);
        SHList[5][2].SetT(1, -2, 2,3,0);
        SHList[5][2].SetT(2, -3, 4,1,0);
        SHList[5][2].SetT(3, -8, 0,3,2);
        SHList[5][2].SetT(4, 24, 2,1,2);

        //-x^3 yz - x y^3 z + 2xyz^3
        SHList[5][3].SetNPS(3);
        SHList[5][3].SetT(0, -1, 3,1,1);
        SHList[5][3].SetT(1, -1, 1,3,1);
        SHList[5][3].SetT(2,  2, 1,1,3);

        //i (x^4 y + 2 x^2 y^3 + y^5 - 12 x^2 y z^2 - 12 y^3 z^2 + 8 y z^4)
        SHList[5][4].SetNPS(6);
        SHList[5][4].SetT(0,   1, 0,5,0);
        SHList[5][4].SetT(1,   2, 2,3,0);
        SHList[5][4].SetT(2,   1, 4,1,0);
        SHList[5][4].SetT(3, -12, 0,3,2);
        SHList[5][4].SetT(4, -12, 2,1,2);
        SHList[5][4].SetT(5,   8, 0,1,4);

        //15 x^4 z + 30 x^2 y^2 z + 15 y^4 z - 40 x^2 z^3 - 40 y^2 z^3 + 8 z^5
        //63 z^5 - 70 z^3 r^2 + 15 z r^4
        // m=0
        SHList[5][5].SetNPS(6);
        SHList[5][5].SetT(0,   8, 0,0,5);
        SHList[5][5].SetT(1, -40, 2,0,3);
        SHList[5][5].SetT(2, -40, 0,2,3);
        SHList[5][5].SetT(3,  15, 4,0,1);
        SHList[5][5].SetT(4,  30, 2,2,1);
        SHList[5][5].SetT(5,  15, 0,4,1);

        //x^5 + 2 x^3 y^2 + x y^4 - 12 x^3 z^2 - 12 x y^2 z^2 + 8 x z^4
        SHList[5][6].SetNPS(6);
        SHList[5][6].SetT(0,   1, 5,0,0);
        SHList[5][6].SetT(1,   2, 3,2,0);
        SHList[5][6].SetT(2,   1, 1,4,0);
        SHList[5][6].SetT(3, -12, 3,0,2);
        SHList[5][6].SetT(4, -12, 1,2,2);
        SHList[5][6].SetT(5,   8, 1,0,4);

        //i (-x^4 z + y^4 z + 2 x^2 y^3 - 2 y^2 x^3)
        SHList[5][7].SetNPS(4);
        SHList[5][7].SetT(0, -1, 4,0,1);
        SHList[5][7].SetT(1,  1, 0,4,1);
        SHList[5][7].SetT(2,  2, 2,0,3);
        SHList[5][7].SetT(3, -2, 0,2,3);

        //3 x y^4 + 2 x^3 y^2 - x^5 - 24 x y^2 z^2 + 8 x^3 z^2
        SHList[5][8].SetNPS(5);
        SHList[5][8].SetT(0, -1, 5,0,0);
        SHList[5][8].SetT(1,  2, 3,2,0);
        SHList[5][8].SetT(2,  3, 1,4,0);
        SHList[5][8].SetT(3,  8, 3,0,2);
        SHList[5][8].SetT(4,-24, 1,2,2);

        //z (x^4 - 6 x^2 y^2 + y^4)
        SHList[5][9].SetNPS(3);
        SHList[5][9].SetT(0,  1, 4,0,1);
        SHList[5][9].SetT(1,  1, 0,4,1);
        SHList[5][9].SetT(2, -6, 2,2,1);

        //x^5 - 10 x^3 y^2 + 5 x y^4
        SHList[5][10].SetNPS(3);
        SHList[5][10].SetT(0,   1, 5,0,0);
        SHList[5][10].SetT(1, -10, 3,2,0);
        SHList[5][10].SetT(2,   5, 1,4,0);
    }


    //normalization
    for (int l=0; l<=LMAX+1; ++l) {
        for (int m=0; m<nmS[l]; ++m) {
            double sum = 0;

            for (int p=0; p<SHList[l][m].nps; ++p)
                for (int q=0; q<SHList[l][m].nps; ++q)
                    sum += SHList[l][m].T[p].cN * SHList[l][m].T[q].cN *
                    AngularIntegral( SHList[l][m].T[p].nx + SHList[l][m].T[q].nx ,
                                     SHList[l][m].T[p].ny + SHList[l][m].T[q].ny ,
                                     SHList[l][m].T[p].nz + SHList[l][m].T[q].nz);

            sum = 1/sqrt(sum);

            for (int p=0; p<SHList[l][m].nps; ++p)
                SHList[l][m].T[p].cN *= sum;

        }
    }


    // computes the axes permutation matrices for real spherical harmonics
    // *******************************************************************
    for (int l=0; l<=LMAX; ++l) {
        for (int m=0; m<2*l+1; ++m) {
            for (int n=0; n<2*l+1; ++n) {
                //(x,y,z) ->(y,z,x)
                //calculate the overlap between m and rotated n
                double sum = 0;

                for (int p=0; p<SHList[l][m].nps; ++p) {
                    for (int q=0; q<SHList[l][n].nps; ++q) {
                        SHTerm & Ta = SHList[l][m].T[p];
                        SHTerm & Tb = SHList[l][n].T[q];
                        sum += Ta.cN * Tb.cN * AngularIntegral( Ta.nx+Tb.ny , Ta.ny+Tb.nz , Ta.nz+Tb.nx);
                    }
                }
                SHrot[l][m][n] = sum;
            }
        }
    }

    for (int l=0; l<=LMAX; ++l) {

        for (int m=0; m<2*l+1; ++m) {
            for (int n=0; n<2*l+1; ++n) {
                double sum = 0;
                for (int o=0; o<2*l+1; ++o) {
                    sum += SHrot[l][m][o] * SHrot[l][o][n];
                }
                SHrotI[l][m][n] = sum;
            }
        }
    }

}

//inicia la lista de cartesianas
void LibAngular::InitCartList () {

    //CartList = new SHTerm*[LMAX+2];

    //añade cada termino a la lista
    for (int l=0; l<=LMAX; ++l) {
        //CartList[l] = new SHTerm[(l+1)*(l+2)/2];

        int t = 0;

        //cout << l << "  " << sqrt(0.5*hGamma[l+1]) << endl;

        for (int a=0; a<=l; ++a) {
            for (int b=0; b<=l-a; ++b) {
                CartList[l][t].nx = l-a-b;
                CartList[l][t].ny = b;
                CartList[l][t].nz = a;
                CartList[l][t].cN = 1./sqrt( AngularIntegral(2*(l-a-b), 2*b, 2*a)); // hGamma[i+1]/(2*hGamma[a] * hGamma[b] * hGamma[i-a-b]) );

                ++t;
            }
        }
    }


    //SP
    {
        //CartList[LSP] = new SHTerm[4];

        CartList[LSP][0].nx = 0;
        CartList[LSP][0].ny = 0;
        CartList[LSP][0].nz = 0;
        CartList[LSP][0].cN = 1./sqrt( AngularIntegral(0, 0, 0) ); // hGamma[i+1]/(2*hGamma[a] * hGamma[b] * hGamma[i-a-b]) );

        CartList[LSP][1].nx = 1;
        CartList[LSP][1].ny = 0;
        CartList[LSP][1].nz = 0;
        CartList[LSP][1].cN = 1./sqrt( AngularIntegral(2, 0, 0) ); // hGamma[i+1]/(2*hGamma[a] * hGamma[b] * hGamma[i-a-b]) );

        CartList[LSP][2].nx = 0;
        CartList[LSP][2].ny = 1;
        CartList[LSP][2].nz = 0;
        CartList[LSP][2].cN = 1./sqrt( AngularIntegral(0, 2, 0) ); // hGamma[i+1]/(2*hGamma[a] * hGamma[b] * hGamma[i-a-b]) );

        CartList[LSP][3].nx = 0;
        CartList[LSP][3].ny = 0;
        CartList[LSP][3].nz = 1;
        CartList[LSP][3].cN = 1./sqrt( AngularIntegral(0, 0, 2) ); // hGamma[i+1]/(2*hGamma[a] * hGamma[b] * hGamma[i-a-b]) );
    }


    //this is used for 1-electron integrals

    /*
    ExyzT = new int*[LMAX+2];
    for (int La=0; La<=LMAX+1; ++La)
        ExyzT[La] = new int[LMAX+2];
    */

    /*
    ExyzC = new int***[LMAX+2];
    for (int La=0; La<=LMAX+1; ++La) {
        ExyzC[La] = new int**[LMAX+2];
        for (int Lb=0; Lb<=LMAX+1; ++Lb) {

            ExyzC[La][Lb] = new int*[nmC[La]];
            for (int a=0; a<nmC[La]; ++a) {
                ExyzC[La][Lb][a] = new int[nmC[Lb]];
            }
        }
    }
    */

    //memoria que ocupa
    for (int La=0; La<=LMAX+1; ++La) {
        for (int Lb=0; Lb<=LMAX+1; ++Lb) {

            ExyzT[La][Lb]=0;

            for (int a=0; a<nmC[La]; ++a) {
                int ax = CartList[La][a].nx;
                int ay = CartList[La][a].ny;
                int az = CartList[La][a].nz;

                for (int b=0; b<nmC[Lb]; ++b) {
                    int bx = CartList[Lb][b].nx;
                    int by = CartList[Lb][b].ny;
                    int bz = CartList[Lb][b].nz;

                    ExyzC[La][Lb][a][b] = (ax+bx+1) * (ay+by+1) * (az+bz+1);
                    ExyzT[La][Lb] += ExyzC[La][Lb][a][b];
                }
            }

        }
    }

}


//convierte un momento angular (en numero) al caracter correspondiente
char LibAngular::L2S (int L) {
    //if (L == LSP) return 'L';

    if (L == 0) return 'S';
    if (L == 1) return 'P';
    if (L == 2) return 'D';
    if (L == 3) return 'F';
    if (L == 4) return 'G';
    if (L == 5) return 'H';
    if (L == 6) return 'I';
    if (L == 7) return 'K';
    return L + 68;
}

//convierte el momento angular de char a entero
int LibAngular::S2L (char S) {
    char shift = 255-32;
    char T = S & shift; //ucase(S);
    //if (T == 'L') return LSP;

    if (T == 'S') return 0;
    if (T == 'P') return 1;
    if (T == 'D') return 2;
    if (T == 'F') return 3;
    if (T == 'G') return 4;
    if (T == 'H') return 5;
    if (T == 'I') return 6;
    if (T == 'K') return 7;
    return T - 68;
}




double AngularIntegral(int nx, int ny, int nz) {
       if ((nx%2 == 1) || (ny%2 == 1) || (nz%2 == 1)) return 0;
       return 2*(LibAngular::hGamma[nx/2]*LibAngular::hGamma[ny/2]*LibAngular::hGamma[nz/2])/LibAngular::hGamma[(nx+ny+nz+2)/2];
}

