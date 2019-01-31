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


#include "newcontractions.hpp"
#include "../integrals/rotations.hpp"

static inline int Tpos (int M1, int M2, int M3, int M4, int v1, int v2, int v3, int v4) {
    return v4 * M3*M2*M1 + v3*M2*M1 + v2*M1 + v1;
}

static inline int Tpos (int M1, int M2, int v1, int v2) {
    return v2*M1 + v1;
}



template <int maxL1, int maxL2> inline void Mprod (cacheline64 * __restrict__ MR, const RotationMatrices64 * __restrict__ RM) {
    const RotationMatrices64 & RM8 = *RM;
    cacheline64 MM[maxL1*maxL2];

    //SP
    if      (maxL1==4) {
        //S
        for (int i=0; i<1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = MR[ij];
            }
        }
        //P
        for (int i=0; i<3; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<3; ++k) {
                    int ki = 3*k+i;
                    int kj = Tpos(maxL1,maxL2,k+1,j);
                    sum  += RM8.PP[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i+1,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==3) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 3*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.PP[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==5) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 5*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.DD[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==7) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 7*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.FF[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==9) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 9*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.GG[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    /*
    else if (maxL1==11) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 11*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.HH[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    */
    else {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = MR[ij];
            }
        }
    }

    if      (maxL2==4) {
        //S
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<1; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = MM[ij];
            }
        }
        //P
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<3; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<3; ++k) {
                    int kj = 3*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k+1);
                    sum  += RM8.PP[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j+1);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==3) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 3*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.PP[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==5) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 5*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.DD[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==7) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 7*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.FF[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==9) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 9*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.GG[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    /*
    else if (maxL2==11) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 11*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.HH[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    */
    else {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = MM[ij];
            }
        }
    }
}

template <int maxL1, int maxL2> inline void MprodT(cacheline64 * __restrict__ MR, const RotationMatrices64 * __restrict__ RM) {
    const RotationMatrices64 & RM8 = *RM;
    cacheline64 MM[maxL1*maxL2];

    if      (maxL1==4) {
        //S
        for (int i=0; i<1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = MR[ij];
            }
        }
        //P
        for (int i=0; i<3; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<3; ++k) {
                    int ik = 3*i+k;
                    int kj = Tpos(maxL1,maxL2,k+1,j);
                    sum  += RM8.PP[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i+1,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==3) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 3*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.PP[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==5) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 5*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.DD[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==7) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 7*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.FF[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==9) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 9*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.GG[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    /*
    else if (maxL1==11) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 11*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.HH[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    */
    else {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = MR[ij];
            }
        }
    }

    if      (maxL2==4) {
        //S
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<1; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = MM[ij];
            }
        }
        //P
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<3; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<3; ++k) {
                    int jk = 3*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k+1);
                    sum  += RM8.PP[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j+1);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==3) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 3*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.PP[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==5) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 5*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.DD[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==7) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 7*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.FF[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==9) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 9*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.GG[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    /*
    else if (maxL2==11) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline64 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 11*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.HH[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    */
    else {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = MM[ij];
            }
        }
    }
}


template <int maxLa, int maxLb, int maxLc, int maxLd> void jABCD(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mb=0; mb<maxLb; ++mb) {
            cacheline64 sum; sum = 0;
            for (int mc=0; mc<maxLc; ++mc) {
                for (int md=0; md<maxLd; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int cd = Tpos(maxLc,maxLd,mc,md);
                    sum += T[pos] * DD[cd];
                }
            }
            int ab = Tpos(maxLa,maxLb,ma,mb);
            FF[ab]+=sum*4;
        }
    }

}

template <int maxLa, int maxLb, int maxLc, int maxLd> void jCDAB(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int mc=0; mc<maxLc; ++mc) {
        for (int md=0; md<maxLd; ++md) {
            cacheline64 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mb=0; mb<maxLb; ++mb) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int ab = Tpos(maxLa,maxLb,ma,mb);
                    sum += T[pos] * DD[ab];
                }
            }
            int cd = Tpos(maxLc,maxLd,mc,md);
            FF[cd]+=sum*4;
        }
    }

}


template <int maxLa, int maxLb, int maxLc, int maxLd> void xACBD(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mc=0; mc<maxLc; ++mc) {
            cacheline64 sum; sum = 0;
            for (int mb=0; mb<maxLb; ++mb) {
                for (int md=0; md<maxLd; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int bd = Tpos(maxLb,maxLd,mb,md);
                    sum += T[pos] * DD[bd];
                }
            }
            int ac = Tpos(maxLa,maxLc,ma,mc);
            FF[ac]+=sum;
        }
    }

}

template <int maxLa, int maxLb, int maxLc, int maxLd> void xADBC(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int md=0; md<maxLd; ++md) {
            cacheline64 sum; sum = 0;
            for (int mb=0; mb<maxLb; ++mb) {
                for (int mc=0; mc<maxLc; ++mc) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int bc = Tpos(maxLb,maxLc,mb,mc);
                    sum+=T[pos] * DD[bc];
                }
            }
            int ad = Tpos(maxLa,maxLd,ma,md);
            FF[ad]+=sum;
        }
    }

}

template <int maxLa, int maxLb, int maxLc, int maxLd> void xBDAC(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int mb=0; mb<maxLb; ++mb) {
        for (int md=0; md<maxLd; ++md) {
            cacheline64 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mc=0; mc<maxLc; ++mc) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int ac = Tpos(maxLa,maxLc,ma,mc);
                    sum+=T[pos] * DD[ac];
                }
            }
            int bd = Tpos(maxLb,maxLd,mb,md);
            FF[bd]+=sum;
        }
    }

}

template <int maxLa, int maxLb, int maxLc, int maxLd> void xBCAD(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int mb=0; mb<maxLb; ++mb) {
        for (int mc=0; mc<maxLc; ++mc) {
            cacheline64 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int md=0; md<maxLd; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int ad = Tpos(maxLa,maxLd,ma,md);
                    sum+=T[pos] * DD[ad];
                }
            }
            int bc = Tpos(maxLb,maxLc,mb,mc);
            FF[bc]+=sum;
        }
    }

}


template <int maxLa, int maxLc, int maxLd> void jAACD(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mb=0; mb<maxLa; ++mb) {
            cacheline64 sum; sum = 0;
            for (int mc=0; mc<maxLc; ++mc) {
                for (int md=0; md<maxLd; ++md) {
                    int pos = Tpos (maxLa, maxLa, maxLc, maxLd, ma, mb, mc, md);
                    int cd = Tpos(maxLc,maxLd,mc,md);
                    sum += T[pos] * DD[cd];
                }
            }
            int ab = Tpos(maxLa, maxLa, ma, mb);
            FF[ab]+=sum*2;
        }
    }

}

template <int maxLa, int maxLc, int maxLd> void jCDAA(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int mc=0; mc<maxLc; ++mc) {
        for (int md=0; md<maxLd; ++md) {
            cacheline64 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mb=0; mb<maxLa; ++mb) {
                    int pos = Tpos (maxLa, maxLa, maxLc, maxLd, ma, mb, mc, md);
                    int ab = Tpos(maxLa,maxLa,ma,mb);
                    sum += T[pos] * DD[ab];
                }
            }
            int cd = Tpos(maxLc,maxLd,mc,md);
            FF[cd]+=sum*2;
        }
    }

}


template <int maxLa, int maxLb, int maxLc> void jABCC(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mb=0; mb<maxLb; ++mb) {
            cacheline64 sum; sum = 0;
            for (int mc=0; mc<maxLc; ++mc) {
                for (int md=0; md<maxLc; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLc, ma, mb, mc, md);
                    int cd = Tpos(maxLc,maxLc,mc,md);
                    sum += T[pos] * DD[cd];
                }
            }
            int ab = Tpos(maxLa,maxLb,ma,mb);
            FF[ab]+=sum*2;
        }
    }

}

template <int maxLa, int maxLb, int maxLc> void jCCAB(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int mc=0; mc<maxLc; ++mc) {
        for (int md=0; md<maxLc; ++md) {
            cacheline64 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mb=0; mb<maxLb; ++mb) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLc, ma, mb, mc, md);
                    int ab = Tpos(maxLa,maxLb,ma,mb);
                    sum += T[pos] * DD[ab];
                }
            }
            int cd = Tpos(maxLc,maxLc,mc,md);
            FF[cd]+=sum*2;
        }
    }

}


template <int maxLa, int maxLc> void jAACC(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mb=0; mb<maxLa; ++mb) {
            cacheline64 sum; sum = 0;
            for (int mc=0; mc<maxLc; ++mc) {
                for (int md=0; md<maxLc; ++md) {
                    int pos = Tpos (maxLa, maxLa, maxLc, maxLc, ma, mb, mc, md);
                    int cd = Tpos(maxLc,maxLc,mc,md);
                    sum += T[pos] * DD[cd];
                }
            }
            int ab = Tpos(maxLa,maxLa,ma,mb);
            FF[ab]+=sum;
        }
    }

}

template <int maxLa, int maxLc> void jCCAA(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int mc=0; mc<maxLc; ++mc) {
        for (int md=0; md<maxLc; ++md) {
            cacheline64 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mb=0; mb<maxLa; ++mb) {
                    int pos = Tpos (maxLa, maxLa, maxLc, maxLc, ma, mb, mc, md);
                    int ab = Tpos(maxLa,maxLa,ma,mb);
                    sum += T[pos] * DD[ab];
                }
            }
            int cd = Tpos(maxLc,maxLc,mc,md);
            FF[cd]+=sum;
        }
    }

}


template <int maxLa, int maxLb> void xAABB(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mc=0; mc<maxLa; ++mc) {
            cacheline64 sum; sum = 0;
            for (int mb=0; mb<maxLb; ++mb) {
                for (int md=0; md<maxLb; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLa, maxLb, ma, mb, mc, md);
                    int bd = Tpos(maxLb,maxLb,mb,md);
                    sum+=T[pos] * DD[bd];
                }
            }
            int ac = Tpos(maxLa, maxLa,ma, mc);
            FF[ac]+=sum*0.5;
        }
    }

}

template <int maxLa, int maxLb> void xBBAA(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int mb=0; mb<maxLb; ++mb) {
        for (int md=0; md<maxLb; ++md) {
            cacheline64 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mc=0; mc<maxLa; ++mc) {
                    int pos = Tpos (maxLa, maxLb, maxLa, maxLb, ma, mb, mc, md);
                    int ac = Tpos(maxLa,maxLa,ma,mc);
                    sum+=T[pos] * DD[ac];
                }
            }
            int bd = Tpos(maxLb, maxLb, mb, md);
            FF[bd]+=sum*0.5;
        }
    }

}

template <int maxLa, int maxLb> void xABAB(cacheline64 * __restrict__ FF, const cacheline64 * __restrict__ T, const cacheline64 * __restrict__ DD) {

    for (int mc=0; mc<maxLa; ++mc) {
        for (int mb=0; mb<maxLb; ++mb) {
            cacheline64 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int md=0; md<maxLb; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLa, maxLb, ma, mb, mc, md);
                    int ad = Tpos(maxLa,maxLb,ma,md);
                    sum += T[pos] * DD[ad];
                }
            }
            int cb = Tpos(maxLa,maxLb,mc,mb);
            FF[cb] += sum;
        }
    }

}

// FLOAT PACK CONTRACTIONS


template <int maxL1, int maxL2> inline void Mprod (cacheline32 * __restrict__ MR, const RotationMatrices32 * __restrict__ RM) {
    const RotationMatrices32 & RM8 = *RM;
    cacheline32 MM[maxL1*maxL2];

    //SP
    if      (maxL1==4) {
        //S
        for (int i=0; i<1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = MR[ij];
            }
        }
        //P
        for (int i=0; i<3; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<3; ++k) {
                    int ki = 3*k+i;
                    int kj = Tpos(maxL1,maxL2,k+1,j);
                    sum  += RM8.PP[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i+1,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==3) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 3*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.PP[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==5) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 5*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.DD[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==7) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 7*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.FF[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==9) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 9*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.GG[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    /*
    else if (maxL1==11) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ki = 11*k+i;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.HH[ki] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    */
    else {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = MR[ij];
            }
        }
    }

    if      (maxL2==4) {
        //S
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<1; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = MM[ij];
            }
        }
        //P
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<3; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<3; ++k) {
                    int kj = 3*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k+1);
                    sum  += RM8.PP[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j+1);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==3) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 3*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.PP[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==5) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 5*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.DD[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==7) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 7*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.FF[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==9) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 9*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.GG[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    /*
    else if (maxL2==11) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int kj = 11*k+j;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.HH[kj] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    */
    else {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = MM[ij];
            }
        }
    }
}

template <int maxL1, int maxL2> inline void MprodT(cacheline32 * __restrict__ MR, const RotationMatrices32 * __restrict__ RM) {
    const RotationMatrices32 & RM8 = *RM;
    cacheline32 MM[maxL1*maxL2];

    if      (maxL1==4) {
        //S
        for (int i=0; i<1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = MR[ij];
            }
        }
        //P
        for (int i=0; i<3; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<3; ++k) {
                    int ik = 3*i+k;
                    int kj = Tpos(maxL1,maxL2,k+1,j);
                    sum  += RM8.PP[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i+1,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==3) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 3*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.PP[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==5) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 5*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.DD[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==7) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 7*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.FF[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    else if (maxL1==9) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 9*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.GG[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    /*
    else if (maxL1==11) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL1; ++k) {
                    int ik = 11*i+k;
                    int kj = Tpos(maxL1,maxL2,k,j);
                    sum  += RM8.HH[ik] * MR[kj];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = sum;
            }
        }
    }
    */
    else {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MM[ij] = MR[ij];
            }
        }
    }

    if      (maxL2==4) {
        //S
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<1; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = MM[ij];
            }
        }
        //P
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<3; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<3; ++k) {
                    int jk = 3*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k+1);
                    sum  += RM8.PP[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j+1);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==3) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 3*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.PP[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==5) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 5*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.DD[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==7) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 7*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.FF[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    else if (maxL2==9) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 9*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.GG[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    /*
    else if (maxL2==11) {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                cacheline32 sum; sum = 0;
                for (int k=0; k<maxL2; ++k) {
                    int jk = 11*j+k;
                    int ik = Tpos(maxL1,maxL2,i,k);
                    sum  += RM8.HH[jk] * MM[ik];
                }
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = sum;
            }
        }
    }
    */
    else {
        for (int i=0; i<maxL1; ++i) {
            for (int j=0; j<maxL2; ++j) {
                int ij = Tpos(maxL1,maxL2,i,j);
                MR[ij] = MM[ij];
            }
        }
    }
}


template <int maxLa, int maxLb, int maxLc, int maxLd> void jABCD(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mb=0; mb<maxLb; ++mb) {
            cacheline32 sum; sum = 0;
            for (int mc=0; mc<maxLc; ++mc) {
                for (int md=0; md<maxLd; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int cd = Tpos(maxLc,maxLd,mc,md);
                    sum += T[pos] * DD[cd];
                }
            }
            int ab = Tpos(maxLa,maxLb,ma,mb);
            FF[ab]+=sum*4;
        }
    }

}

template <int maxLa, int maxLb, int maxLc, int maxLd> void jCDAB(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int mc=0; mc<maxLc; ++mc) {
        for (int md=0; md<maxLd; ++md) {
            cacheline32 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mb=0; mb<maxLb; ++mb) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int ab = Tpos(maxLa,maxLb,ma,mb);
                    sum += T[pos] * DD[ab];
                }
            }
            int cd = Tpos(maxLc,maxLd,mc,md);
            FF[cd]+=sum*4;
        }
    }

}


template <int maxLa, int maxLb, int maxLc, int maxLd> void xACBD(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mc=0; mc<maxLc; ++mc) {
            cacheline32 sum; sum = 0;
            for (int mb=0; mb<maxLb; ++mb) {
                for (int md=0; md<maxLd; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int bd = Tpos(maxLb,maxLd,mb,md);
                    sum += T[pos] * DD[bd];
                }
            }
            int ac = Tpos(maxLa,maxLc,ma,mc);
            FF[ac]+=sum;
        }
    }

}

template <int maxLa, int maxLb, int maxLc, int maxLd> void xADBC(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int md=0; md<maxLd; ++md) {
            cacheline32 sum; sum = 0;
            for (int mb=0; mb<maxLb; ++mb) {
                for (int mc=0; mc<maxLc; ++mc) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int bc = Tpos(maxLb,maxLc,mb,mc);
                    sum+=T[pos] * DD[bc];
                }
            }
            int ad = Tpos(maxLa,maxLd,ma,md);
            FF[ad]+=sum;
        }
    }

}

template <int maxLa, int maxLb, int maxLc, int maxLd> void xBDAC(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int mb=0; mb<maxLb; ++mb) {
        for (int md=0; md<maxLd; ++md) {
            cacheline32 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mc=0; mc<maxLc; ++mc) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int ac = Tpos(maxLa,maxLc,ma,mc);
                    sum+=T[pos] * DD[ac];
                }
            }
            int bd = Tpos(maxLb,maxLd,mb,md);
            FF[bd]+=sum;
        }
    }

}

template <int maxLa, int maxLb, int maxLc, int maxLd> void xBCAD(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int mb=0; mb<maxLb; ++mb) {
        for (int mc=0; mc<maxLc; ++mc) {
            cacheline32 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int md=0; md<maxLd; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLd, ma, mb, mc, md);
                    int ad = Tpos(maxLa,maxLd,ma,md);
                    sum+=T[pos] * DD[ad];
                }
            }
            int bc = Tpos(maxLb,maxLc,mb,mc);
            FF[bc]+=sum;
        }
    }

}


template <int maxLa, int maxLc, int maxLd> void jAACD(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mb=0; mb<maxLa; ++mb) {
            cacheline32 sum; sum = 0;
            for (int mc=0; mc<maxLc; ++mc) {
                for (int md=0; md<maxLd; ++md) {
                    int pos = Tpos (maxLa, maxLa, maxLc, maxLd, ma, mb, mc, md);
                    int cd = Tpos(maxLc,maxLd,mc,md);
                    sum += T[pos] * DD[cd];
                }
            }
            int ab = Tpos(maxLa, maxLa, ma, mb);
            FF[ab]+=sum*2;
        }
    }

}

template <int maxLa, int maxLc, int maxLd> void jCDAA(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int mc=0; mc<maxLc; ++mc) {
        for (int md=0; md<maxLd; ++md) {
            cacheline32 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mb=0; mb<maxLa; ++mb) {
                    int pos = Tpos (maxLa, maxLa, maxLc, maxLd, ma, mb, mc, md);
                    int ab = Tpos(maxLa,maxLa,ma,mb);
                    sum += T[pos] * DD[ab];
                }
            }
            int cd = Tpos(maxLc,maxLd,mc,md);
            FF[cd]+=sum*2;
        }
    }

}


template <int maxLa, int maxLb, int maxLc> void jABCC(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mb=0; mb<maxLb; ++mb) {
            cacheline32 sum; sum = 0;
            for (int mc=0; mc<maxLc; ++mc) {
                for (int md=0; md<maxLc; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLc, ma, mb, mc, md);
                    int cd = Tpos(maxLc,maxLc,mc,md);
                    sum += T[pos] * DD[cd];
                }
            }
            int ab = Tpos(maxLa,maxLb,ma,mb);
            FF[ab]+=sum*2;
        }
    }

}

template <int maxLa, int maxLb, int maxLc> void jCCAB(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int mc=0; mc<maxLc; ++mc) {
        for (int md=0; md<maxLc; ++md) {
            cacheline32 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mb=0; mb<maxLb; ++mb) {
                    int pos = Tpos (maxLa, maxLb, maxLc, maxLc, ma, mb, mc, md);
                    int ab = Tpos(maxLa,maxLb,ma,mb);
                    sum += T[pos] * DD[ab];
                }
            }
            int cd = Tpos(maxLc,maxLc,mc,md);
            FF[cd]+=sum*2;
        }
    }

}


template <int maxLa, int maxLc> void jAACC(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mb=0; mb<maxLa; ++mb) {
            cacheline32 sum; sum = 0;
            for (int mc=0; mc<maxLc; ++mc) {
                for (int md=0; md<maxLc; ++md) {
                    int pos = Tpos (maxLa, maxLa, maxLc, maxLc, ma, mb, mc, md);
                    int cd = Tpos(maxLc,maxLc,mc,md);
                    sum += T[pos] * DD[cd];
                }
            }
            int ab = Tpos(maxLa,maxLa,ma,mb);
            FF[ab]+=sum;
        }
    }

}

template <int maxLa, int maxLc> void jCCAA(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int mc=0; mc<maxLc; ++mc) {
        for (int md=0; md<maxLc; ++md) {
            cacheline32 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mb=0; mb<maxLa; ++mb) {
                    int pos = Tpos (maxLa, maxLa, maxLc, maxLc, ma, mb, mc, md);
                    int ab = Tpos(maxLa,maxLa,ma,mb);
                    sum += T[pos] * DD[ab];
                }
            }
            int cd = Tpos(maxLc,maxLc,mc,md);
            FF[cd]+=sum;
        }
    }

}


template <int maxLa, int maxLb> void xAABB(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int ma=0; ma<maxLa; ++ma) {
        for (int mc=0; mc<maxLa; ++mc) {
            cacheline32 sum; sum = 0;
            for (int mb=0; mb<maxLb; ++mb) {
                for (int md=0; md<maxLb; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLa, maxLb, ma, mb, mc, md);
                    int bd = Tpos(maxLb,maxLb,mb,md);
                    sum+=T[pos] * DD[bd];
                }
            }
            int ac = Tpos(maxLa, maxLa,ma, mc);
            FF[ac]+=sum*0.5;
        }
    }

}

template <int maxLa, int maxLb> void xBBAA(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int mb=0; mb<maxLb; ++mb) {
        for (int md=0; md<maxLb; ++md) {
            cacheline32 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int mc=0; mc<maxLa; ++mc) {
                    int pos = Tpos (maxLa, maxLb, maxLa, maxLb, ma, mb, mc, md);
                    int ac = Tpos(maxLa,maxLa,ma,mc);
                    sum+=T[pos] * DD[ac];
                }
            }
            int bd = Tpos(maxLb, maxLb, mb, md);
            FF[bd]+=sum*0.5;
        }
    }

}

template <int maxLa, int maxLb> void xABAB(cacheline32 * __restrict__ FF, const cacheline32 * __restrict__ T, const cacheline32 * __restrict__ DD) {

    for (int mc=0; mc<maxLa; ++mc) {
        for (int mb=0; mb<maxLb; ++mb) {
            cacheline32 sum; sum = 0;
            for (int ma=0; ma<maxLa; ++ma) {
                for (int md=0; md<maxLb; ++md) {
                    int pos = Tpos (maxLa, maxLb, maxLa, maxLb, ma, mb, mc, md);
                    int ad = Tpos(maxLa,maxLb,ma,md);
                    sum += T[pos] * DD[ad];
                }
            }
            int cb = Tpos(maxLa,maxLb,mc,mb);
            FF[cb] += sum;
        }
    }

}





DMD::DMD() {

    // 64

    //JXC NN
    {
        J64abcd [0][0][0][0] = jABCD  <1,1,1,1>;
        J64cdab [0][0][0][0] = jCDAB  <1,1,1,1>;
        X64acbd [0][0][0][0] = xACBD  <1,1,1,1>;
        X64adbc [0][0][0][0] = xADBC  <1,1,1,1>;
        X64bdac [0][0][0][0] = xBDAC  <1,1,1,1>;
        X64bcad [0][0][0][0] = xBCAD  <1,1,1,1>;
        J64abcd [0][0][1][0] = jABCD  <1,1,3,1>;
        J64cdab [0][0][1][0] = jCDAB  <1,1,3,1>;
        X64acbd [0][0][1][0] = xACBD  <1,1,3,1>;
        X64adbc [0][0][1][0] = xADBC  <1,1,3,1>;
        X64bdac [0][0][1][0] = xBDAC  <1,1,3,1>;
        X64bcad [0][0][1][0] = xBCAD  <1,1,3,1>;
        J64abcd [0][0][1][1] = jABCD  <1,1,3,3>;
        J64cdab [0][0][1][1] = jCDAB  <1,1,3,3>;
        X64acbd [0][0][1][1] = xACBD  <1,1,3,3>;
        X64adbc [0][0][1][1] = xADBC  <1,1,3,3>;
        X64bdac [0][0][1][1] = xBDAC  <1,1,3,3>;
        X64bcad [0][0][1][1] = xBCAD  <1,1,3,3>;
        J64abcd [0][0][2][0] = jABCD  <1,1,5,1>;
        J64cdab [0][0][2][0] = jCDAB  <1,1,5,1>;
        X64acbd [0][0][2][0] = xACBD  <1,1,5,1>;
        X64adbc [0][0][2][0] = xADBC  <1,1,5,1>;
        X64bdac [0][0][2][0] = xBDAC  <1,1,5,1>;
        X64bcad [0][0][2][0] = xBCAD  <1,1,5,1>;
        J64abcd [0][0][2][1] = jABCD  <1,1,5,3>;
        J64cdab [0][0][2][1] = jCDAB  <1,1,5,3>;
        X64acbd [0][0][2][1] = xACBD  <1,1,5,3>;
        X64adbc [0][0][2][1] = xADBC  <1,1,5,3>;
        X64bdac [0][0][2][1] = xBDAC  <1,1,5,3>;
        X64bcad [0][0][2][1] = xBCAD  <1,1,5,3>;
        J64abcd [0][0][2][2] = jABCD  <1,1,5,5>;
        J64cdab [0][0][2][2] = jCDAB  <1,1,5,5>;
        X64acbd [0][0][2][2] = xACBD  <1,1,5,5>;
        X64adbc [0][0][2][2] = xADBC  <1,1,5,5>;
        X64bdac [0][0][2][2] = xBDAC  <1,1,5,5>;
        X64bcad [0][0][2][2] = xBCAD  <1,1,5,5>;
        J64abcd [0][0][3][0] = jABCD  <1,1,7,1>;
        J64cdab [0][0][3][0] = jCDAB  <1,1,7,1>;
        X64acbd [0][0][3][0] = xACBD  <1,1,7,1>;
        X64adbc [0][0][3][0] = xADBC  <1,1,7,1>;
        X64bdac [0][0][3][0] = xBDAC  <1,1,7,1>;
        X64bcad [0][0][3][0] = xBCAD  <1,1,7,1>;
        J64abcd [0][0][3][1] = jABCD  <1,1,7,3>;
        J64cdab [0][0][3][1] = jCDAB  <1,1,7,3>;
        X64acbd [0][0][3][1] = xACBD  <1,1,7,3>;
        X64adbc [0][0][3][1] = xADBC  <1,1,7,3>;
        X64bdac [0][0][3][1] = xBDAC  <1,1,7,3>;
        X64bcad [0][0][3][1] = xBCAD  <1,1,7,3>;
        J64abcd [0][0][3][2] = jABCD  <1,1,7,5>;
        J64cdab [0][0][3][2] = jCDAB  <1,1,7,5>;
        X64acbd [0][0][3][2] = xACBD  <1,1,7,5>;
        X64adbc [0][0][3][2] = xADBC  <1,1,7,5>;
        X64bdac [0][0][3][2] = xBDAC  <1,1,7,5>;
        X64bcad [0][0][3][2] = xBCAD  <1,1,7,5>;
        J64abcd [0][0][3][3] = jABCD  <1,1,7,7>;
        J64cdab [0][0][3][3] = jCDAB  <1,1,7,7>;
        X64acbd [0][0][3][3] = xACBD  <1,1,7,7>;
        X64adbc [0][0][3][3] = xADBC  <1,1,7,7>;
        X64bdac [0][0][3][3] = xBDAC  <1,1,7,7>;
        X64bcad [0][0][3][3] = xBCAD  <1,1,7,7>;
        J64abcd [0][0][4][0] = jABCD  <1,1,9,1>;
        J64cdab [0][0][4][0] = jCDAB  <1,1,9,1>;
        X64acbd [0][0][4][0] = xACBD  <1,1,9,1>;
        X64adbc [0][0][4][0] = xADBC  <1,1,9,1>;
        X64bdac [0][0][4][0] = xBDAC  <1,1,9,1>;
        X64bcad [0][0][4][0] = xBCAD  <1,1,9,1>;
        J64abcd [0][0][4][1] = jABCD  <1,1,9,3>;
        J64cdab [0][0][4][1] = jCDAB  <1,1,9,3>;
        X64acbd [0][0][4][1] = xACBD  <1,1,9,3>;
        X64adbc [0][0][4][1] = xADBC  <1,1,9,3>;
        X64bdac [0][0][4][1] = xBDAC  <1,1,9,3>;
        X64bcad [0][0][4][1] = xBCAD  <1,1,9,3>;
        J64abcd [0][0][4][2] = jABCD  <1,1,9,5>;
        J64cdab [0][0][4][2] = jCDAB  <1,1,9,5>;
        X64acbd [0][0][4][2] = xACBD  <1,1,9,5>;
        X64adbc [0][0][4][2] = xADBC  <1,1,9,5>;
        X64bdac [0][0][4][2] = xBDAC  <1,1,9,5>;
        X64bcad [0][0][4][2] = xBCAD  <1,1,9,5>;
        J64abcd [0][0][4][3] = jABCD  <1,1,9,7>;
        J64cdab [0][0][4][3] = jCDAB  <1,1,9,7>;
        X64acbd [0][0][4][3] = xACBD  <1,1,9,7>;
        X64adbc [0][0][4][3] = xADBC  <1,1,9,7>;
        X64bdac [0][0][4][3] = xBDAC  <1,1,9,7>;
        X64bcad [0][0][4][3] = xBCAD  <1,1,9,7>;
        J64abcd [0][0][4][4] = jABCD  <1,1,9,9>;
        J64cdab [0][0][4][4] = jCDAB  <1,1,9,9>;
        X64acbd [0][0][4][4] = xACBD  <1,1,9,9>;
        X64adbc [0][0][4][4] = xADBC  <1,1,9,9>;
        X64bdac [0][0][4][4] = xBDAC  <1,1,9,9>;
        X64bcad [0][0][4][4] = xBCAD  <1,1,9,9>;

        J64abcd [1][0][0][0] = jABCD  <3,1,1,1>;
        J64cdab [1][0][0][0] = jCDAB  <3,1,1,1>;
        X64acbd [1][0][0][0] = xACBD  <3,1,1,1>;
        X64adbc [1][0][0][0] = xADBC  <3,1,1,1>;
        X64bdac [1][0][0][0] = xBDAC  <3,1,1,1>;
        X64bcad [1][0][0][0] = xBCAD  <3,1,1,1>;
        J64abcd [1][0][1][0] = jABCD  <3,1,3,1>;
        J64cdab [1][0][1][0] = jCDAB  <3,1,3,1>;
        X64acbd [1][0][1][0] = xACBD  <3,1,3,1>;
        X64adbc [1][0][1][0] = xADBC  <3,1,3,1>;
        X64bdac [1][0][1][0] = xBDAC  <3,1,3,1>;
        X64bcad [1][0][1][0] = xBCAD  <3,1,3,1>;
        J64abcd [1][0][1][1] = jABCD  <3,1,3,3>;
        J64cdab [1][0][1][1] = jCDAB  <3,1,3,3>;
        X64acbd [1][0][1][1] = xACBD  <3,1,3,3>;
        X64adbc [1][0][1][1] = xADBC  <3,1,3,3>;
        X64bdac [1][0][1][1] = xBDAC  <3,1,3,3>;
        X64bcad [1][0][1][1] = xBCAD  <3,1,3,3>;
        J64abcd [1][0][2][0] = jABCD  <3,1,5,1>;
        J64cdab [1][0][2][0] = jCDAB  <3,1,5,1>;
        X64acbd [1][0][2][0] = xACBD  <3,1,5,1>;
        X64adbc [1][0][2][0] = xADBC  <3,1,5,1>;
        X64bdac [1][0][2][0] = xBDAC  <3,1,5,1>;
        X64bcad [1][0][2][0] = xBCAD  <3,1,5,1>;
        J64abcd [1][0][2][1] = jABCD  <3,1,5,3>;
        J64cdab [1][0][2][1] = jCDAB  <3,1,5,3>;
        X64acbd [1][0][2][1] = xACBD  <3,1,5,3>;
        X64adbc [1][0][2][1] = xADBC  <3,1,5,3>;
        X64bdac [1][0][2][1] = xBDAC  <3,1,5,3>;
        X64bcad [1][0][2][1] = xBCAD  <3,1,5,3>;
        J64abcd [1][0][2][2] = jABCD  <3,1,5,5>;
        J64cdab [1][0][2][2] = jCDAB  <3,1,5,5>;
        X64acbd [1][0][2][2] = xACBD  <3,1,5,5>;
        X64adbc [1][0][2][2] = xADBC  <3,1,5,5>;
        X64bdac [1][0][2][2] = xBDAC  <3,1,5,5>;
        X64bcad [1][0][2][2] = xBCAD  <3,1,5,5>;
        J64abcd [1][0][3][0] = jABCD  <3,1,7,1>;
        J64cdab [1][0][3][0] = jCDAB  <3,1,7,1>;
        X64acbd [1][0][3][0] = xACBD  <3,1,7,1>;
        X64adbc [1][0][3][0] = xADBC  <3,1,7,1>;
        X64bdac [1][0][3][0] = xBDAC  <3,1,7,1>;
        X64bcad [1][0][3][0] = xBCAD  <3,1,7,1>;
        J64abcd [1][0][3][1] = jABCD  <3,1,7,3>;
        J64cdab [1][0][3][1] = jCDAB  <3,1,7,3>;
        X64acbd [1][0][3][1] = xACBD  <3,1,7,3>;
        X64adbc [1][0][3][1] = xADBC  <3,1,7,3>;
        X64bdac [1][0][3][1] = xBDAC  <3,1,7,3>;
        X64bcad [1][0][3][1] = xBCAD  <3,1,7,3>;
        J64abcd [1][0][3][2] = jABCD  <3,1,7,5>;
        J64cdab [1][0][3][2] = jCDAB  <3,1,7,5>;
        X64acbd [1][0][3][2] = xACBD  <3,1,7,5>;
        X64adbc [1][0][3][2] = xADBC  <3,1,7,5>;
        X64bdac [1][0][3][2] = xBDAC  <3,1,7,5>;
        X64bcad [1][0][3][2] = xBCAD  <3,1,7,5>;
        J64abcd [1][0][3][3] = jABCD  <3,1,7,7>;
        J64cdab [1][0][3][3] = jCDAB  <3,1,7,7>;
        X64acbd [1][0][3][3] = xACBD  <3,1,7,7>;
        X64adbc [1][0][3][3] = xADBC  <3,1,7,7>;
        X64bdac [1][0][3][3] = xBDAC  <3,1,7,7>;
        X64bcad [1][0][3][3] = xBCAD  <3,1,7,7>;
        J64abcd [1][0][4][0] = jABCD  <3,1,9,1>;
        J64cdab [1][0][4][0] = jCDAB  <3,1,9,1>;
        X64acbd [1][0][4][0] = xACBD  <3,1,9,1>;
        X64adbc [1][0][4][0] = xADBC  <3,1,9,1>;
        X64bdac [1][0][4][0] = xBDAC  <3,1,9,1>;
        X64bcad [1][0][4][0] = xBCAD  <3,1,9,1>;
        J64abcd [1][0][4][1] = jABCD  <3,1,9,3>;
        J64cdab [1][0][4][1] = jCDAB  <3,1,9,3>;
        X64acbd [1][0][4][1] = xACBD  <3,1,9,3>;
        X64adbc [1][0][4][1] = xADBC  <3,1,9,3>;
        X64bdac [1][0][4][1] = xBDAC  <3,1,9,3>;
        X64bcad [1][0][4][1] = xBCAD  <3,1,9,3>;
        J64abcd [1][0][4][2] = jABCD  <3,1,9,5>;
        J64cdab [1][0][4][2] = jCDAB  <3,1,9,5>;
        X64acbd [1][0][4][2] = xACBD  <3,1,9,5>;
        X64adbc [1][0][4][2] = xADBC  <3,1,9,5>;
        X64bdac [1][0][4][2] = xBDAC  <3,1,9,5>;
        X64bcad [1][0][4][2] = xBCAD  <3,1,9,5>;
        J64abcd [1][0][4][3] = jABCD  <3,1,9,7>;
        J64cdab [1][0][4][3] = jCDAB  <3,1,9,7>;
        X64acbd [1][0][4][3] = xACBD  <3,1,9,7>;
        X64adbc [1][0][4][3] = xADBC  <3,1,9,7>;
        X64bdac [1][0][4][3] = xBDAC  <3,1,9,7>;
        X64bcad [1][0][4][3] = xBCAD  <3,1,9,7>;
        J64abcd [1][0][4][4] = jABCD  <3,1,9,9>;
        J64cdab [1][0][4][4] = jCDAB  <3,1,9,9>;
        X64acbd [1][0][4][4] = xACBD  <3,1,9,9>;
        X64adbc [1][0][4][4] = xADBC  <3,1,9,9>;
        X64bdac [1][0][4][4] = xBDAC  <3,1,9,9>;
        X64bcad [1][0][4][4] = xBCAD  <3,1,9,9>;

        J64abcd [1][1][0][0] = jABCD  <3,3,1,1>;
        J64cdab [1][1][0][0] = jCDAB  <3,3,1,1>;
        X64acbd [1][1][0][0] = xACBD  <3,3,1,1>;
        X64adbc [1][1][0][0] = xADBC  <3,3,1,1>;
        X64bdac [1][1][0][0] = xBDAC  <3,3,1,1>;
        X64bcad [1][1][0][0] = xBCAD  <3,3,1,1>;
        J64abcd [1][1][1][0] = jABCD  <3,3,3,1>;
        J64cdab [1][1][1][0] = jCDAB  <3,3,3,1>;
        X64acbd [1][1][1][0] = xACBD  <3,3,3,1>;
        X64adbc [1][1][1][0] = xADBC  <3,3,3,1>;
        X64bdac [1][1][1][0] = xBDAC  <3,3,3,1>;
        X64bcad [1][1][1][0] = xBCAD  <3,3,3,1>;
        J64abcd [1][1][1][1] = jABCD  <3,3,3,3>;
        J64cdab [1][1][1][1] = jCDAB  <3,3,3,3>;
        X64acbd [1][1][1][1] = xACBD  <3,3,3,3>;
        X64adbc [1][1][1][1] = xADBC  <3,3,3,3>;
        X64bdac [1][1][1][1] = xBDAC  <3,3,3,3>;
        X64bcad [1][1][1][1] = xBCAD  <3,3,3,3>;
        J64abcd [1][1][2][0] = jABCD  <3,3,5,1>;
        J64cdab [1][1][2][0] = jCDAB  <3,3,5,1>;
        X64acbd [1][1][2][0] = xACBD  <3,3,5,1>;
        X64adbc [1][1][2][0] = xADBC  <3,3,5,1>;
        X64bdac [1][1][2][0] = xBDAC  <3,3,5,1>;
        X64bcad [1][1][2][0] = xBCAD  <3,3,5,1>;
        J64abcd [1][1][2][1] = jABCD  <3,3,5,3>;
        J64cdab [1][1][2][1] = jCDAB  <3,3,5,3>;
        X64acbd [1][1][2][1] = xACBD  <3,3,5,3>;
        X64adbc [1][1][2][1] = xADBC  <3,3,5,3>;
        X64bdac [1][1][2][1] = xBDAC  <3,3,5,3>;
        X64bcad [1][1][2][1] = xBCAD  <3,3,5,3>;
        J64abcd [1][1][2][2] = jABCD  <3,3,5,5>;
        J64cdab [1][1][2][2] = jCDAB  <3,3,5,5>;
        X64acbd [1][1][2][2] = xACBD  <3,3,5,5>;
        X64adbc [1][1][2][2] = xADBC  <3,3,5,5>;
        X64bdac [1][1][2][2] = xBDAC  <3,3,5,5>;
        X64bcad [1][1][2][2] = xBCAD  <3,3,5,5>;
        J64abcd [1][1][3][0] = jABCD  <3,3,7,1>;
        J64cdab [1][1][3][0] = jCDAB  <3,3,7,1>;
        X64acbd [1][1][3][0] = xACBD  <3,3,7,1>;
        X64adbc [1][1][3][0] = xADBC  <3,3,7,1>;
        X64bdac [1][1][3][0] = xBDAC  <3,3,7,1>;
        X64bcad [1][1][3][0] = xBCAD  <3,3,7,1>;
        J64abcd [1][1][3][1] = jABCD  <3,3,7,3>;
        J64cdab [1][1][3][1] = jCDAB  <3,3,7,3>;
        X64acbd [1][1][3][1] = xACBD  <3,3,7,3>;
        X64adbc [1][1][3][1] = xADBC  <3,3,7,3>;
        X64bdac [1][1][3][1] = xBDAC  <3,3,7,3>;
        X64bcad [1][1][3][1] = xBCAD  <3,3,7,3>;
        J64abcd [1][1][3][2] = jABCD  <3,3,7,5>;
        J64cdab [1][1][3][2] = jCDAB  <3,3,7,5>;
        X64acbd [1][1][3][2] = xACBD  <3,3,7,5>;
        X64adbc [1][1][3][2] = xADBC  <3,3,7,5>;
        X64bdac [1][1][3][2] = xBDAC  <3,3,7,5>;
        X64bcad [1][1][3][2] = xBCAD  <3,3,7,5>;
        J64abcd [1][1][3][3] = jABCD  <3,3,7,7>;
        J64cdab [1][1][3][3] = jCDAB  <3,3,7,7>;
        X64acbd [1][1][3][3] = xACBD  <3,3,7,7>;
        X64adbc [1][1][3][3] = xADBC  <3,3,7,7>;
        X64bdac [1][1][3][3] = xBDAC  <3,3,7,7>;
        X64bcad [1][1][3][3] = xBCAD  <3,3,7,7>;
        J64abcd [1][1][4][0] = jABCD  <3,3,9,1>;
        J64cdab [1][1][4][0] = jCDAB  <3,3,9,1>;
        X64acbd [1][1][4][0] = xACBD  <3,3,9,1>;
        X64adbc [1][1][4][0] = xADBC  <3,3,9,1>;
        X64bdac [1][1][4][0] = xBDAC  <3,3,9,1>;
        X64bcad [1][1][4][0] = xBCAD  <3,3,9,1>;
        J64abcd [1][1][4][1] = jABCD  <3,3,9,3>;
        J64cdab [1][1][4][1] = jCDAB  <3,3,9,3>;
        X64acbd [1][1][4][1] = xACBD  <3,3,9,3>;
        X64adbc [1][1][4][1] = xADBC  <3,3,9,3>;
        X64bdac [1][1][4][1] = xBDAC  <3,3,9,3>;
        X64bcad [1][1][4][1] = xBCAD  <3,3,9,3>;
        J64abcd [1][1][4][2] = jABCD  <3,3,9,5>;
        J64cdab [1][1][4][2] = jCDAB  <3,3,9,5>;
        X64acbd [1][1][4][2] = xACBD  <3,3,9,5>;
        X64adbc [1][1][4][2] = xADBC  <3,3,9,5>;
        X64bdac [1][1][4][2] = xBDAC  <3,3,9,5>;
        X64bcad [1][1][4][2] = xBCAD  <3,3,9,5>;
        J64abcd [1][1][4][3] = jABCD  <3,3,9,7>;
        J64cdab [1][1][4][3] = jCDAB  <3,3,9,7>;
        X64acbd [1][1][4][3] = xACBD  <3,3,9,7>;
        X64adbc [1][1][4][3] = xADBC  <3,3,9,7>;
        X64bdac [1][1][4][3] = xBDAC  <3,3,9,7>;
        X64bcad [1][1][4][3] = xBCAD  <3,3,9,7>;
        J64abcd [1][1][4][4] = jABCD  <3,3,9,9>;
        J64cdab [1][1][4][4] = jCDAB  <3,3,9,9>;
        X64acbd [1][1][4][4] = xACBD  <3,3,9,9>;
        X64adbc [1][1][4][4] = xADBC  <3,3,9,9>;
        X64bdac [1][1][4][4] = xBDAC  <3,3,9,9>;
        X64bcad [1][1][4][4] = xBCAD  <3,3,9,9>;

        J64abcd [2][0][0][0] = jABCD  <5,1,1,1>;
        J64cdab [2][0][0][0] = jCDAB  <5,1,1,1>;
        X64acbd [2][0][0][0] = xACBD  <5,1,1,1>;
        X64adbc [2][0][0][0] = xADBC  <5,1,1,1>;
        X64bdac [2][0][0][0] = xBDAC  <5,1,1,1>;
        X64bcad [2][0][0][0] = xBCAD  <5,1,1,1>;
        J64abcd [2][0][1][0] = jABCD  <5,1,3,1>;
        J64cdab [2][0][1][0] = jCDAB  <5,1,3,1>;
        X64acbd [2][0][1][0] = xACBD  <5,1,3,1>;
        X64adbc [2][0][1][0] = xADBC  <5,1,3,1>;
        X64bdac [2][0][1][0] = xBDAC  <5,1,3,1>;
        X64bcad [2][0][1][0] = xBCAD  <5,1,3,1>;
        J64abcd [2][0][1][1] = jABCD  <5,1,3,3>;
        J64cdab [2][0][1][1] = jCDAB  <5,1,3,3>;
        X64acbd [2][0][1][1] = xACBD  <5,1,3,3>;
        X64adbc [2][0][1][1] = xADBC  <5,1,3,3>;
        X64bdac [2][0][1][1] = xBDAC  <5,1,3,3>;
        X64bcad [2][0][1][1] = xBCAD  <5,1,3,3>;
        J64abcd [2][0][2][0] = jABCD  <5,1,5,1>;
        J64cdab [2][0][2][0] = jCDAB  <5,1,5,1>;
        X64acbd [2][0][2][0] = xACBD  <5,1,5,1>;
        X64adbc [2][0][2][0] = xADBC  <5,1,5,1>;
        X64bdac [2][0][2][0] = xBDAC  <5,1,5,1>;
        X64bcad [2][0][2][0] = xBCAD  <5,1,5,1>;
        J64abcd [2][0][2][1] = jABCD  <5,1,5,3>;
        J64cdab [2][0][2][1] = jCDAB  <5,1,5,3>;
        X64acbd [2][0][2][1] = xACBD  <5,1,5,3>;
        X64adbc [2][0][2][1] = xADBC  <5,1,5,3>;
        X64bdac [2][0][2][1] = xBDAC  <5,1,5,3>;
        X64bcad [2][0][2][1] = xBCAD  <5,1,5,3>;
        J64abcd [2][0][2][2] = jABCD  <5,1,5,5>;
        J64cdab [2][0][2][2] = jCDAB  <5,1,5,5>;
        X64acbd [2][0][2][2] = xACBD  <5,1,5,5>;
        X64adbc [2][0][2][2] = xADBC  <5,1,5,5>;
        X64bdac [2][0][2][2] = xBDAC  <5,1,5,5>;
        X64bcad [2][0][2][2] = xBCAD  <5,1,5,5>;
        J64abcd [2][0][3][0] = jABCD  <5,1,7,1>;
        J64cdab [2][0][3][0] = jCDAB  <5,1,7,1>;
        X64acbd [2][0][3][0] = xACBD  <5,1,7,1>;
        X64adbc [2][0][3][0] = xADBC  <5,1,7,1>;
        X64bdac [2][0][3][0] = xBDAC  <5,1,7,1>;
        X64bcad [2][0][3][0] = xBCAD  <5,1,7,1>;
        J64abcd [2][0][3][1] = jABCD  <5,1,7,3>;
        J64cdab [2][0][3][1] = jCDAB  <5,1,7,3>;
        X64acbd [2][0][3][1] = xACBD  <5,1,7,3>;
        X64adbc [2][0][3][1] = xADBC  <5,1,7,3>;
        X64bdac [2][0][3][1] = xBDAC  <5,1,7,3>;
        X64bcad [2][0][3][1] = xBCAD  <5,1,7,3>;
        J64abcd [2][0][3][2] = jABCD  <5,1,7,5>;
        J64cdab [2][0][3][2] = jCDAB  <5,1,7,5>;
        X64acbd [2][0][3][2] = xACBD  <5,1,7,5>;
        X64adbc [2][0][3][2] = xADBC  <5,1,7,5>;
        X64bdac [2][0][3][2] = xBDAC  <5,1,7,5>;
        X64bcad [2][0][3][2] = xBCAD  <5,1,7,5>;
        J64abcd [2][0][3][3] = jABCD  <5,1,7,7>;
        J64cdab [2][0][3][3] = jCDAB  <5,1,7,7>;
        X64acbd [2][0][3][3] = xACBD  <5,1,7,7>;
        X64adbc [2][0][3][3] = xADBC  <5,1,7,7>;
        X64bdac [2][0][3][3] = xBDAC  <5,1,7,7>;
        X64bcad [2][0][3][3] = xBCAD  <5,1,7,7>;
        J64abcd [2][0][4][0] = jABCD  <5,1,9,1>;
        J64cdab [2][0][4][0] = jCDAB  <5,1,9,1>;
        X64acbd [2][0][4][0] = xACBD  <5,1,9,1>;
        X64adbc [2][0][4][0] = xADBC  <5,1,9,1>;
        X64bdac [2][0][4][0] = xBDAC  <5,1,9,1>;
        X64bcad [2][0][4][0] = xBCAD  <5,1,9,1>;
        J64abcd [2][0][4][1] = jABCD  <5,1,9,3>;
        J64cdab [2][0][4][1] = jCDAB  <5,1,9,3>;
        X64acbd [2][0][4][1] = xACBD  <5,1,9,3>;
        X64adbc [2][0][4][1] = xADBC  <5,1,9,3>;
        X64bdac [2][0][4][1] = xBDAC  <5,1,9,3>;
        X64bcad [2][0][4][1] = xBCAD  <5,1,9,3>;
        J64abcd [2][0][4][2] = jABCD  <5,1,9,5>;
        J64cdab [2][0][4][2] = jCDAB  <5,1,9,5>;
        X64acbd [2][0][4][2] = xACBD  <5,1,9,5>;
        X64adbc [2][0][4][2] = xADBC  <5,1,9,5>;
        X64bdac [2][0][4][2] = xBDAC  <5,1,9,5>;
        X64bcad [2][0][4][2] = xBCAD  <5,1,9,5>;
        J64abcd [2][0][4][3] = jABCD  <5,1,9,7>;
        J64cdab [2][0][4][3] = jCDAB  <5,1,9,7>;
        X64acbd [2][0][4][3] = xACBD  <5,1,9,7>;
        X64adbc [2][0][4][3] = xADBC  <5,1,9,7>;
        X64bdac [2][0][4][3] = xBDAC  <5,1,9,7>;
        X64bcad [2][0][4][3] = xBCAD  <5,1,9,7>;
        J64abcd [2][0][4][4] = jABCD  <5,1,9,9>;
        J64cdab [2][0][4][4] = jCDAB  <5,1,9,9>;
        X64acbd [2][0][4][4] = xACBD  <5,1,9,9>;
        X64adbc [2][0][4][4] = xADBC  <5,1,9,9>;
        X64bdac [2][0][4][4] = xBDAC  <5,1,9,9>;
        X64bcad [2][0][4][4] = xBCAD  <5,1,9,9>;

        J64abcd [2][1][0][0] = jABCD  <5,3,1,1>;
        J64cdab [2][1][0][0] = jCDAB  <5,3,1,1>;
        X64acbd [2][1][0][0] = xACBD  <5,3,1,1>;
        X64adbc [2][1][0][0] = xADBC  <5,3,1,1>;
        X64bdac [2][1][0][0] = xBDAC  <5,3,1,1>;
        X64bcad [2][1][0][0] = xBCAD  <5,3,1,1>;
        J64abcd [2][1][1][0] = jABCD  <5,3,3,1>;
        J64cdab [2][1][1][0] = jCDAB  <5,3,3,1>;
        X64acbd [2][1][1][0] = xACBD  <5,3,3,1>;
        X64adbc [2][1][1][0] = xADBC  <5,3,3,1>;
        X64bdac [2][1][1][0] = xBDAC  <5,3,3,1>;
        X64bcad [2][1][1][0] = xBCAD  <5,3,3,1>;
        J64abcd [2][1][1][1] = jABCD  <5,3,3,3>;
        J64cdab [2][1][1][1] = jCDAB  <5,3,3,3>;
        X64acbd [2][1][1][1] = xACBD  <5,3,3,3>;
        X64adbc [2][1][1][1] = xADBC  <5,3,3,3>;
        X64bdac [2][1][1][1] = xBDAC  <5,3,3,3>;
        X64bcad [2][1][1][1] = xBCAD  <5,3,3,3>;
        J64abcd [2][1][2][0] = jABCD  <5,3,5,1>;
        J64cdab [2][1][2][0] = jCDAB  <5,3,5,1>;
        X64acbd [2][1][2][0] = xACBD  <5,3,5,1>;
        X64adbc [2][1][2][0] = xADBC  <5,3,5,1>;
        X64bdac [2][1][2][0] = xBDAC  <5,3,5,1>;
        X64bcad [2][1][2][0] = xBCAD  <5,3,5,1>;
        J64abcd [2][1][2][1] = jABCD  <5,3,5,3>;
        J64cdab [2][1][2][1] = jCDAB  <5,3,5,3>;
        X64acbd [2][1][2][1] = xACBD  <5,3,5,3>;
        X64adbc [2][1][2][1] = xADBC  <5,3,5,3>;
        X64bdac [2][1][2][1] = xBDAC  <5,3,5,3>;
        X64bcad [2][1][2][1] = xBCAD  <5,3,5,3>;
        J64abcd [2][1][2][2] = jABCD  <5,3,5,5>;
        J64cdab [2][1][2][2] = jCDAB  <5,3,5,5>;
        X64acbd [2][1][2][2] = xACBD  <5,3,5,5>;
        X64adbc [2][1][2][2] = xADBC  <5,3,5,5>;
        X64bdac [2][1][2][2] = xBDAC  <5,3,5,5>;
        X64bcad [2][1][2][2] = xBCAD  <5,3,5,5>;
        J64abcd [2][1][3][0] = jABCD  <5,3,7,1>;
        J64cdab [2][1][3][0] = jCDAB  <5,3,7,1>;
        X64acbd [2][1][3][0] = xACBD  <5,3,7,1>;
        X64adbc [2][1][3][0] = xADBC  <5,3,7,1>;
        X64bdac [2][1][3][0] = xBDAC  <5,3,7,1>;
        X64bcad [2][1][3][0] = xBCAD  <5,3,7,1>;
        J64abcd [2][1][3][1] = jABCD  <5,3,7,3>;
        J64cdab [2][1][3][1] = jCDAB  <5,3,7,3>;
        X64acbd [2][1][3][1] = xACBD  <5,3,7,3>;
        X64adbc [2][1][3][1] = xADBC  <5,3,7,3>;
        X64bdac [2][1][3][1] = xBDAC  <5,3,7,3>;
        X64bcad [2][1][3][1] = xBCAD  <5,3,7,3>;
        J64abcd [2][1][3][2] = jABCD  <5,3,7,5>;
        J64cdab [2][1][3][2] = jCDAB  <5,3,7,5>;
        X64acbd [2][1][3][2] = xACBD  <5,3,7,5>;
        X64adbc [2][1][3][2] = xADBC  <5,3,7,5>;
        X64bdac [2][1][3][2] = xBDAC  <5,3,7,5>;
        X64bcad [2][1][3][2] = xBCAD  <5,3,7,5>;
        J64abcd [2][1][3][3] = jABCD  <5,3,7,7>;
        J64cdab [2][1][3][3] = jCDAB  <5,3,7,7>;
        X64acbd [2][1][3][3] = xACBD  <5,3,7,7>;
        X64adbc [2][1][3][3] = xADBC  <5,3,7,7>;
        X64bdac [2][1][3][3] = xBDAC  <5,3,7,7>;
        X64bcad [2][1][3][3] = xBCAD  <5,3,7,7>;
        J64abcd [2][1][4][0] = jABCD  <5,3,9,1>;
        J64cdab [2][1][4][0] = jCDAB  <5,3,9,1>;
        X64acbd [2][1][4][0] = xACBD  <5,3,9,1>;
        X64adbc [2][1][4][0] = xADBC  <5,3,9,1>;
        X64bdac [2][1][4][0] = xBDAC  <5,3,9,1>;
        X64bcad [2][1][4][0] = xBCAD  <5,3,9,1>;
        J64abcd [2][1][4][1] = jABCD  <5,3,9,3>;
        J64cdab [2][1][4][1] = jCDAB  <5,3,9,3>;
        X64acbd [2][1][4][1] = xACBD  <5,3,9,3>;
        X64adbc [2][1][4][1] = xADBC  <5,3,9,3>;
        X64bdac [2][1][4][1] = xBDAC  <5,3,9,3>;
        X64bcad [2][1][4][1] = xBCAD  <5,3,9,3>;
        J64abcd [2][1][4][2] = jABCD  <5,3,9,5>;
        J64cdab [2][1][4][2] = jCDAB  <5,3,9,5>;
        X64acbd [2][1][4][2] = xACBD  <5,3,9,5>;
        X64adbc [2][1][4][2] = xADBC  <5,3,9,5>;
        X64bdac [2][1][4][2] = xBDAC  <5,3,9,5>;
        X64bcad [2][1][4][2] = xBCAD  <5,3,9,5>;
        J64abcd [2][1][4][3] = jABCD  <5,3,9,7>;
        J64cdab [2][1][4][3] = jCDAB  <5,3,9,7>;
        X64acbd [2][1][4][3] = xACBD  <5,3,9,7>;
        X64adbc [2][1][4][3] = xADBC  <5,3,9,7>;
        X64bdac [2][1][4][3] = xBDAC  <5,3,9,7>;
        X64bcad [2][1][4][3] = xBCAD  <5,3,9,7>;
        J64abcd [2][1][4][4] = jABCD  <5,3,9,9>;
        J64cdab [2][1][4][4] = jCDAB  <5,3,9,9>;
        X64acbd [2][1][4][4] = xACBD  <5,3,9,9>;
        X64adbc [2][1][4][4] = xADBC  <5,3,9,9>;
        X64bdac [2][1][4][4] = xBDAC  <5,3,9,9>;
        X64bcad [2][1][4][4] = xBCAD  <5,3,9,9>;

        J64abcd [2][2][0][0] = jABCD  <5,5,1,1>;
        J64cdab [2][2][0][0] = jCDAB  <5,5,1,1>;
        X64acbd [2][2][0][0] = xACBD  <5,5,1,1>;
        X64adbc [2][2][0][0] = xADBC  <5,5,1,1>;
        X64bdac [2][2][0][0] = xBDAC  <5,5,1,1>;
        X64bcad [2][2][0][0] = xBCAD  <5,5,1,1>;
        J64abcd [2][2][1][0] = jABCD  <5,5,3,1>;
        J64cdab [2][2][1][0] = jCDAB  <5,5,3,1>;
        X64acbd [2][2][1][0] = xACBD  <5,5,3,1>;
        X64adbc [2][2][1][0] = xADBC  <5,5,3,1>;
        X64bdac [2][2][1][0] = xBDAC  <5,5,3,1>;
        X64bcad [2][2][1][0] = xBCAD  <5,5,3,1>;
        J64abcd [2][2][1][1] = jABCD  <5,5,3,3>;
        J64cdab [2][2][1][1] = jCDAB  <5,5,3,3>;
        X64acbd [2][2][1][1] = xACBD  <5,5,3,3>;
        X64adbc [2][2][1][1] = xADBC  <5,5,3,3>;
        X64bdac [2][2][1][1] = xBDAC  <5,5,3,3>;
        X64bcad [2][2][1][1] = xBCAD  <5,5,3,3>;
        J64abcd [2][2][2][0] = jABCD  <5,5,5,1>;
        J64cdab [2][2][2][0] = jCDAB  <5,5,5,1>;
        X64acbd [2][2][2][0] = xACBD  <5,5,5,1>;
        X64adbc [2][2][2][0] = xADBC  <5,5,5,1>;
        X64bdac [2][2][2][0] = xBDAC  <5,5,5,1>;
        X64bcad [2][2][2][0] = xBCAD  <5,5,5,1>;
        J64abcd [2][2][2][1] = jABCD  <5,5,5,3>;
        J64cdab [2][2][2][1] = jCDAB  <5,5,5,3>;
        X64acbd [2][2][2][1] = xACBD  <5,5,5,3>;
        X64adbc [2][2][2][1] = xADBC  <5,5,5,3>;
        X64bdac [2][2][2][1] = xBDAC  <5,5,5,3>;
        X64bcad [2][2][2][1] = xBCAD  <5,5,5,3>;
        J64abcd [2][2][2][2] = jABCD  <5,5,5,5>;
        J64cdab [2][2][2][2] = jCDAB  <5,5,5,5>;
        X64acbd [2][2][2][2] = xACBD  <5,5,5,5>;
        X64adbc [2][2][2][2] = xADBC  <5,5,5,5>;
        X64bdac [2][2][2][2] = xBDAC  <5,5,5,5>;
        X64bcad [2][2][2][2] = xBCAD  <5,5,5,5>;
        J64abcd [2][2][3][0] = jABCD  <5,5,7,1>;
        J64cdab [2][2][3][0] = jCDAB  <5,5,7,1>;
        X64acbd [2][2][3][0] = xACBD  <5,5,7,1>;
        X64adbc [2][2][3][0] = xADBC  <5,5,7,1>;
        X64bdac [2][2][3][0] = xBDAC  <5,5,7,1>;
        X64bcad [2][2][3][0] = xBCAD  <5,5,7,1>;
        J64abcd [2][2][3][1] = jABCD  <5,5,7,3>;
        J64cdab [2][2][3][1] = jCDAB  <5,5,7,3>;
        X64acbd [2][2][3][1] = xACBD  <5,5,7,3>;
        X64adbc [2][2][3][1] = xADBC  <5,5,7,3>;
        X64bdac [2][2][3][1] = xBDAC  <5,5,7,3>;
        X64bcad [2][2][3][1] = xBCAD  <5,5,7,3>;
        J64abcd [2][2][3][2] = jABCD  <5,5,7,5>;
        J64cdab [2][2][3][2] = jCDAB  <5,5,7,5>;
        X64acbd [2][2][3][2] = xACBD  <5,5,7,5>;
        X64adbc [2][2][3][2] = xADBC  <5,5,7,5>;
        X64bdac [2][2][3][2] = xBDAC  <5,5,7,5>;
        X64bcad [2][2][3][2] = xBCAD  <5,5,7,5>;
        J64abcd [2][2][3][3] = jABCD  <5,5,7,7>;
        J64cdab [2][2][3][3] = jCDAB  <5,5,7,7>;
        X64acbd [2][2][3][3] = xACBD  <5,5,7,7>;
        X64adbc [2][2][3][3] = xADBC  <5,5,7,7>;
        X64bdac [2][2][3][3] = xBDAC  <5,5,7,7>;
        X64bcad [2][2][3][3] = xBCAD  <5,5,7,7>;
        J64abcd [2][2][4][0] = jABCD  <5,5,9,1>;
        J64cdab [2][2][4][0] = jCDAB  <5,5,9,1>;
        X64acbd [2][2][4][0] = xACBD  <5,5,9,1>;
        X64adbc [2][2][4][0] = xADBC  <5,5,9,1>;
        X64bdac [2][2][4][0] = xBDAC  <5,5,9,1>;
        X64bcad [2][2][4][0] = xBCAD  <5,5,9,1>;
        J64abcd [2][2][4][1] = jABCD  <5,5,9,3>;
        J64cdab [2][2][4][1] = jCDAB  <5,5,9,3>;
        X64acbd [2][2][4][1] = xACBD  <5,5,9,3>;
        X64adbc [2][2][4][1] = xADBC  <5,5,9,3>;
        X64bdac [2][2][4][1] = xBDAC  <5,5,9,3>;
        X64bcad [2][2][4][1] = xBCAD  <5,5,9,3>;
        J64abcd [2][2][4][2] = jABCD  <5,5,9,5>;
        J64cdab [2][2][4][2] = jCDAB  <5,5,9,5>;
        X64acbd [2][2][4][2] = xACBD  <5,5,9,5>;
        X64adbc [2][2][4][2] = xADBC  <5,5,9,5>;
        X64bdac [2][2][4][2] = xBDAC  <5,5,9,5>;
        X64bcad [2][2][4][2] = xBCAD  <5,5,9,5>;
        J64abcd [2][2][4][3] = jABCD  <5,5,9,7>;
        J64cdab [2][2][4][3] = jCDAB  <5,5,9,7>;
        X64acbd [2][2][4][3] = xACBD  <5,5,9,7>;
        X64adbc [2][2][4][3] = xADBC  <5,5,9,7>;
        X64bdac [2][2][4][3] = xBDAC  <5,5,9,7>;
        X64bcad [2][2][4][3] = xBCAD  <5,5,9,7>;
        J64abcd [2][2][4][4] = jABCD  <5,5,9,9>;
        J64cdab [2][2][4][4] = jCDAB  <5,5,9,9>;
        X64acbd [2][2][4][4] = xACBD  <5,5,9,9>;
        X64adbc [2][2][4][4] = xADBC  <5,5,9,9>;
        X64bdac [2][2][4][4] = xBDAC  <5,5,9,9>;
        X64bcad [2][2][4][4] = xBCAD  <5,5,9,9>;

        J64abcd [3][0][0][0] = jABCD  <7,1,1,1>;
        J64cdab [3][0][0][0] = jCDAB  <7,1,1,1>;
        X64acbd [3][0][0][0] = xACBD  <7,1,1,1>;
        X64adbc [3][0][0][0] = xADBC  <7,1,1,1>;
        X64bdac [3][0][0][0] = xBDAC  <7,1,1,1>;
        X64bcad [3][0][0][0] = xBCAD  <7,1,1,1>;
        J64abcd [3][0][1][0] = jABCD  <7,1,3,1>;
        J64cdab [3][0][1][0] = jCDAB  <7,1,3,1>;
        X64acbd [3][0][1][0] = xACBD  <7,1,3,1>;
        X64adbc [3][0][1][0] = xADBC  <7,1,3,1>;
        X64bdac [3][0][1][0] = xBDAC  <7,1,3,1>;
        X64bcad [3][0][1][0] = xBCAD  <7,1,3,1>;
        J64abcd [3][0][1][1] = jABCD  <7,1,3,3>;
        J64cdab [3][0][1][1] = jCDAB  <7,1,3,3>;
        X64acbd [3][0][1][1] = xACBD  <7,1,3,3>;
        X64adbc [3][0][1][1] = xADBC  <7,1,3,3>;
        X64bdac [3][0][1][1] = xBDAC  <7,1,3,3>;
        X64bcad [3][0][1][1] = xBCAD  <7,1,3,3>;
        J64abcd [3][0][2][0] = jABCD  <7,1,5,1>;
        J64cdab [3][0][2][0] = jCDAB  <7,1,5,1>;
        X64acbd [3][0][2][0] = xACBD  <7,1,5,1>;
        X64adbc [3][0][2][0] = xADBC  <7,1,5,1>;
        X64bdac [3][0][2][0] = xBDAC  <7,1,5,1>;
        X64bcad [3][0][2][0] = xBCAD  <7,1,5,1>;
        J64abcd [3][0][2][1] = jABCD  <7,1,5,3>;
        J64cdab [3][0][2][1] = jCDAB  <7,1,5,3>;
        X64acbd [3][0][2][1] = xACBD  <7,1,5,3>;
        X64adbc [3][0][2][1] = xADBC  <7,1,5,3>;
        X64bdac [3][0][2][1] = xBDAC  <7,1,5,3>;
        X64bcad [3][0][2][1] = xBCAD  <7,1,5,3>;
        J64abcd [3][0][2][2] = jABCD  <7,1,5,5>;
        J64cdab [3][0][2][2] = jCDAB  <7,1,5,5>;
        X64acbd [3][0][2][2] = xACBD  <7,1,5,5>;
        X64adbc [3][0][2][2] = xADBC  <7,1,5,5>;
        X64bdac [3][0][2][2] = xBDAC  <7,1,5,5>;
        X64bcad [3][0][2][2] = xBCAD  <7,1,5,5>;
        J64abcd [3][0][3][0] = jABCD  <7,1,7,1>;
        J64cdab [3][0][3][0] = jCDAB  <7,1,7,1>;
        X64acbd [3][0][3][0] = xACBD  <7,1,7,1>;
        X64adbc [3][0][3][0] = xADBC  <7,1,7,1>;
        X64bdac [3][0][3][0] = xBDAC  <7,1,7,1>;
        X64bcad [3][0][3][0] = xBCAD  <7,1,7,1>;
        J64abcd [3][0][3][1] = jABCD  <7,1,7,3>;
        J64cdab [3][0][3][1] = jCDAB  <7,1,7,3>;
        X64acbd [3][0][3][1] = xACBD  <7,1,7,3>;
        X64adbc [3][0][3][1] = xADBC  <7,1,7,3>;
        X64bdac [3][0][3][1] = xBDAC  <7,1,7,3>;
        X64bcad [3][0][3][1] = xBCAD  <7,1,7,3>;
        J64abcd [3][0][3][2] = jABCD  <7,1,7,5>;
        J64cdab [3][0][3][2] = jCDAB  <7,1,7,5>;
        X64acbd [3][0][3][2] = xACBD  <7,1,7,5>;
        X64adbc [3][0][3][2] = xADBC  <7,1,7,5>;
        X64bdac [3][0][3][2] = xBDAC  <7,1,7,5>;
        X64bcad [3][0][3][2] = xBCAD  <7,1,7,5>;
        J64abcd [3][0][3][3] = jABCD  <7,1,7,7>;
        J64cdab [3][0][3][3] = jCDAB  <7,1,7,7>;
        X64acbd [3][0][3][3] = xACBD  <7,1,7,7>;
        X64adbc [3][0][3][3] = xADBC  <7,1,7,7>;
        X64bdac [3][0][3][3] = xBDAC  <7,1,7,7>;
        X64bcad [3][0][3][3] = xBCAD  <7,1,7,7>;
        J64abcd [3][0][4][0] = jABCD  <7,1,9,1>;
        J64cdab [3][0][4][0] = jCDAB  <7,1,9,1>;
        X64acbd [3][0][4][0] = xACBD  <7,1,9,1>;
        X64adbc [3][0][4][0] = xADBC  <7,1,9,1>;
        X64bdac [3][0][4][0] = xBDAC  <7,1,9,1>;
        X64bcad [3][0][4][0] = xBCAD  <7,1,9,1>;
        J64abcd [3][0][4][1] = jABCD  <7,1,9,3>;
        J64cdab [3][0][4][1] = jCDAB  <7,1,9,3>;
        X64acbd [3][0][4][1] = xACBD  <7,1,9,3>;
        X64adbc [3][0][4][1] = xADBC  <7,1,9,3>;
        X64bdac [3][0][4][1] = xBDAC  <7,1,9,3>;
        X64bcad [3][0][4][1] = xBCAD  <7,1,9,3>;
        J64abcd [3][0][4][2] = jABCD  <7,1,9,5>;
        J64cdab [3][0][4][2] = jCDAB  <7,1,9,5>;
        X64acbd [3][0][4][2] = xACBD  <7,1,9,5>;
        X64adbc [3][0][4][2] = xADBC  <7,1,9,5>;
        X64bdac [3][0][4][2] = xBDAC  <7,1,9,5>;
        X64bcad [3][0][4][2] = xBCAD  <7,1,9,5>;
        J64abcd [3][0][4][3] = jABCD  <7,1,9,7>;
        J64cdab [3][0][4][3] = jCDAB  <7,1,9,7>;
        X64acbd [3][0][4][3] = xACBD  <7,1,9,7>;
        X64adbc [3][0][4][3] = xADBC  <7,1,9,7>;
        X64bdac [3][0][4][3] = xBDAC  <7,1,9,7>;
        X64bcad [3][0][4][3] = xBCAD  <7,1,9,7>;
        J64abcd [3][0][4][4] = jABCD  <7,1,9,9>;
        J64cdab [3][0][4][4] = jCDAB  <7,1,9,9>;
        X64acbd [3][0][4][4] = xACBD  <7,1,9,9>;
        X64adbc [3][0][4][4] = xADBC  <7,1,9,9>;
        X64bdac [3][0][4][4] = xBDAC  <7,1,9,9>;
        X64bcad [3][0][4][4] = xBCAD  <7,1,9,9>;

        J64abcd [3][1][0][0] = jABCD  <7,3,1,1>;
        J64cdab [3][1][0][0] = jCDAB  <7,3,1,1>;
        X64acbd [3][1][0][0] = xACBD  <7,3,1,1>;
        X64adbc [3][1][0][0] = xADBC  <7,3,1,1>;
        X64bdac [3][1][0][0] = xBDAC  <7,3,1,1>;
        X64bcad [3][1][0][0] = xBCAD  <7,3,1,1>;
        J64abcd [3][1][1][0] = jABCD  <7,3,3,1>;
        J64cdab [3][1][1][0] = jCDAB  <7,3,3,1>;
        X64acbd [3][1][1][0] = xACBD  <7,3,3,1>;
        X64adbc [3][1][1][0] = xADBC  <7,3,3,1>;
        X64bdac [3][1][1][0] = xBDAC  <7,3,3,1>;
        X64bcad [3][1][1][0] = xBCAD  <7,3,3,1>;
        J64abcd [3][1][1][1] = jABCD  <7,3,3,3>;
        J64cdab [3][1][1][1] = jCDAB  <7,3,3,3>;
        X64acbd [3][1][1][1] = xACBD  <7,3,3,3>;
        X64adbc [3][1][1][1] = xADBC  <7,3,3,3>;
        X64bdac [3][1][1][1] = xBDAC  <7,3,3,3>;
        X64bcad [3][1][1][1] = xBCAD  <7,3,3,3>;
        J64abcd [3][1][2][0] = jABCD  <7,3,5,1>;
        J64cdab [3][1][2][0] = jCDAB  <7,3,5,1>;
        X64acbd [3][1][2][0] = xACBD  <7,3,5,1>;
        X64adbc [3][1][2][0] = xADBC  <7,3,5,1>;
        X64bdac [3][1][2][0] = xBDAC  <7,3,5,1>;
        X64bcad [3][1][2][0] = xBCAD  <7,3,5,1>;
        J64abcd [3][1][2][1] = jABCD  <7,3,5,3>;
        J64cdab [3][1][2][1] = jCDAB  <7,3,5,3>;
        X64acbd [3][1][2][1] = xACBD  <7,3,5,3>;
        X64adbc [3][1][2][1] = xADBC  <7,3,5,3>;
        X64bdac [3][1][2][1] = xBDAC  <7,3,5,3>;
        X64bcad [3][1][2][1] = xBCAD  <7,3,5,3>;
        J64abcd [3][1][2][2] = jABCD  <7,3,5,5>;
        J64cdab [3][1][2][2] = jCDAB  <7,3,5,5>;
        X64acbd [3][1][2][2] = xACBD  <7,3,5,5>;
        X64adbc [3][1][2][2] = xADBC  <7,3,5,5>;
        X64bdac [3][1][2][2] = xBDAC  <7,3,5,5>;
        X64bcad [3][1][2][2] = xBCAD  <7,3,5,5>;
        J64abcd [3][1][3][0] = jABCD  <7,3,7,1>;
        J64cdab [3][1][3][0] = jCDAB  <7,3,7,1>;
        X64acbd [3][1][3][0] = xACBD  <7,3,7,1>;
        X64adbc [3][1][3][0] = xADBC  <7,3,7,1>;
        X64bdac [3][1][3][0] = xBDAC  <7,3,7,1>;
        X64bcad [3][1][3][0] = xBCAD  <7,3,7,1>;
        J64abcd [3][1][3][1] = jABCD  <7,3,7,3>;
        J64cdab [3][1][3][1] = jCDAB  <7,3,7,3>;
        X64acbd [3][1][3][1] = xACBD  <7,3,7,3>;
        X64adbc [3][1][3][1] = xADBC  <7,3,7,3>;
        X64bdac [3][1][3][1] = xBDAC  <7,3,7,3>;
        X64bcad [3][1][3][1] = xBCAD  <7,3,7,3>;
        J64abcd [3][1][3][2] = jABCD  <7,3,7,5>;
        J64cdab [3][1][3][2] = jCDAB  <7,3,7,5>;
        X64acbd [3][1][3][2] = xACBD  <7,3,7,5>;
        X64adbc [3][1][3][2] = xADBC  <7,3,7,5>;
        X64bdac [3][1][3][2] = xBDAC  <7,3,7,5>;
        X64bcad [3][1][3][2] = xBCAD  <7,3,7,5>;
        J64abcd [3][1][3][3] = jABCD  <7,3,7,7>;
        J64cdab [3][1][3][3] = jCDAB  <7,3,7,7>;
        X64acbd [3][1][3][3] = xACBD  <7,3,7,7>;
        X64adbc [3][1][3][3] = xADBC  <7,3,7,7>;
        X64bdac [3][1][3][3] = xBDAC  <7,3,7,7>;
        X64bcad [3][1][3][3] = xBCAD  <7,3,7,7>;
        J64abcd [3][1][4][0] = jABCD  <7,3,9,1>;
        J64cdab [3][1][4][0] = jCDAB  <7,3,9,1>;
        X64acbd [3][1][4][0] = xACBD  <7,3,9,1>;
        X64adbc [3][1][4][0] = xADBC  <7,3,9,1>;
        X64bdac [3][1][4][0] = xBDAC  <7,3,9,1>;
        X64bcad [3][1][4][0] = xBCAD  <7,3,9,1>;
        J64abcd [3][1][4][1] = jABCD  <7,3,9,3>;
        J64cdab [3][1][4][1] = jCDAB  <7,3,9,3>;
        X64acbd [3][1][4][1] = xACBD  <7,3,9,3>;
        X64adbc [3][1][4][1] = xADBC  <7,3,9,3>;
        X64bdac [3][1][4][1] = xBDAC  <7,3,9,3>;
        X64bcad [3][1][4][1] = xBCAD  <7,3,9,3>;
        J64abcd [3][1][4][2] = jABCD  <7,3,9,5>;
        J64cdab [3][1][4][2] = jCDAB  <7,3,9,5>;
        X64acbd [3][1][4][2] = xACBD  <7,3,9,5>;
        X64adbc [3][1][4][2] = xADBC  <7,3,9,5>;
        X64bdac [3][1][4][2] = xBDAC  <7,3,9,5>;
        X64bcad [3][1][4][2] = xBCAD  <7,3,9,5>;
        J64abcd [3][1][4][3] = jABCD  <7,3,9,7>;
        J64cdab [3][1][4][3] = jCDAB  <7,3,9,7>;
        X64acbd [3][1][4][3] = xACBD  <7,3,9,7>;
        X64adbc [3][1][4][3] = xADBC  <7,3,9,7>;
        X64bdac [3][1][4][3] = xBDAC  <7,3,9,7>;
        X64bcad [3][1][4][3] = xBCAD  <7,3,9,7>;
        J64abcd [3][1][4][4] = jABCD  <7,3,9,9>;
        J64cdab [3][1][4][4] = jCDAB  <7,3,9,9>;
        X64acbd [3][1][4][4] = xACBD  <7,3,9,9>;
        X64adbc [3][1][4][4] = xADBC  <7,3,9,9>;
        X64bdac [3][1][4][4] = xBDAC  <7,3,9,9>;
        X64bcad [3][1][4][4] = xBCAD  <7,3,9,9>;

        J64abcd [3][2][0][0] = jABCD  <7,5,1,1>;
        J64cdab [3][2][0][0] = jCDAB  <7,5,1,1>;
        X64acbd [3][2][0][0] = xACBD  <7,5,1,1>;
        X64adbc [3][2][0][0] = xADBC  <7,5,1,1>;
        X64bdac [3][2][0][0] = xBDAC  <7,5,1,1>;
        X64bcad [3][2][0][0] = xBCAD  <7,5,1,1>;
        J64abcd [3][2][1][0] = jABCD  <7,5,3,1>;
        J64cdab [3][2][1][0] = jCDAB  <7,5,3,1>;
        X64acbd [3][2][1][0] = xACBD  <7,5,3,1>;
        X64adbc [3][2][1][0] = xADBC  <7,5,3,1>;
        X64bdac [3][2][1][0] = xBDAC  <7,5,3,1>;
        X64bcad [3][2][1][0] = xBCAD  <7,5,3,1>;
        J64abcd [3][2][1][1] = jABCD  <7,5,3,3>;
        J64cdab [3][2][1][1] = jCDAB  <7,5,3,3>;
        X64acbd [3][2][1][1] = xACBD  <7,5,3,3>;
        X64adbc [3][2][1][1] = xADBC  <7,5,3,3>;
        X64bdac [3][2][1][1] = xBDAC  <7,5,3,3>;
        X64bcad [3][2][1][1] = xBCAD  <7,5,3,3>;
        J64abcd [3][2][2][0] = jABCD  <7,5,5,1>;
        J64cdab [3][2][2][0] = jCDAB  <7,5,5,1>;
        X64acbd [3][2][2][0] = xACBD  <7,5,5,1>;
        X64adbc [3][2][2][0] = xADBC  <7,5,5,1>;
        X64bdac [3][2][2][0] = xBDAC  <7,5,5,1>;
        X64bcad [3][2][2][0] = xBCAD  <7,5,5,1>;
        J64abcd [3][2][2][1] = jABCD  <7,5,5,3>;
        J64cdab [3][2][2][1] = jCDAB  <7,5,5,3>;
        X64acbd [3][2][2][1] = xACBD  <7,5,5,3>;
        X64adbc [3][2][2][1] = xADBC  <7,5,5,3>;
        X64bdac [3][2][2][1] = xBDAC  <7,5,5,3>;
        X64bcad [3][2][2][1] = xBCAD  <7,5,5,3>;
        J64abcd [3][2][2][2] = jABCD  <7,5,5,5>;
        J64cdab [3][2][2][2] = jCDAB  <7,5,5,5>;
        X64acbd [3][2][2][2] = xACBD  <7,5,5,5>;
        X64adbc [3][2][2][2] = xADBC  <7,5,5,5>;
        X64bdac [3][2][2][2] = xBDAC  <7,5,5,5>;
        X64bcad [3][2][2][2] = xBCAD  <7,5,5,5>;
        J64abcd [3][2][3][0] = jABCD  <7,5,7,1>;
        J64cdab [3][2][3][0] = jCDAB  <7,5,7,1>;
        X64acbd [3][2][3][0] = xACBD  <7,5,7,1>;
        X64adbc [3][2][3][0] = xADBC  <7,5,7,1>;
        X64bdac [3][2][3][0] = xBDAC  <7,5,7,1>;
        X64bcad [3][2][3][0] = xBCAD  <7,5,7,1>;
        J64abcd [3][2][3][1] = jABCD  <7,5,7,3>;
        J64cdab [3][2][3][1] = jCDAB  <7,5,7,3>;
        X64acbd [3][2][3][1] = xACBD  <7,5,7,3>;
        X64adbc [3][2][3][1] = xADBC  <7,5,7,3>;
        X64bdac [3][2][3][1] = xBDAC  <7,5,7,3>;
        X64bcad [3][2][3][1] = xBCAD  <7,5,7,3>;
        J64abcd [3][2][3][2] = jABCD  <7,5,7,5>;
        J64cdab [3][2][3][2] = jCDAB  <7,5,7,5>;
        X64acbd [3][2][3][2] = xACBD  <7,5,7,5>;
        X64adbc [3][2][3][2] = xADBC  <7,5,7,5>;
        X64bdac [3][2][3][2] = xBDAC  <7,5,7,5>;
        X64bcad [3][2][3][2] = xBCAD  <7,5,7,5>;
        J64abcd [3][2][3][3] = jABCD  <7,5,7,7>;
        J64cdab [3][2][3][3] = jCDAB  <7,5,7,7>;
        X64acbd [3][2][3][3] = xACBD  <7,5,7,7>;
        X64adbc [3][2][3][3] = xADBC  <7,5,7,7>;
        X64bdac [3][2][3][3] = xBDAC  <7,5,7,7>;
        X64bcad [3][2][3][3] = xBCAD  <7,5,7,7>;
        J64abcd [3][2][4][0] = jABCD  <7,5,9,1>;
        J64cdab [3][2][4][0] = jCDAB  <7,5,9,1>;
        X64acbd [3][2][4][0] = xACBD  <7,5,9,1>;
        X64adbc [3][2][4][0] = xADBC  <7,5,9,1>;
        X64bdac [3][2][4][0] = xBDAC  <7,5,9,1>;
        X64bcad [3][2][4][0] = xBCAD  <7,5,9,1>;
        J64abcd [3][2][4][1] = jABCD  <7,5,9,3>;
        J64cdab [3][2][4][1] = jCDAB  <7,5,9,3>;
        X64acbd [3][2][4][1] = xACBD  <7,5,9,3>;
        X64adbc [3][2][4][1] = xADBC  <7,5,9,3>;
        X64bdac [3][2][4][1] = xBDAC  <7,5,9,3>;
        X64bcad [3][2][4][1] = xBCAD  <7,5,9,3>;
        J64abcd [3][2][4][2] = jABCD  <7,5,9,5>;
        J64cdab [3][2][4][2] = jCDAB  <7,5,9,5>;
        X64acbd [3][2][4][2] = xACBD  <7,5,9,5>;
        X64adbc [3][2][4][2] = xADBC  <7,5,9,5>;
        X64bdac [3][2][4][2] = xBDAC  <7,5,9,5>;
        X64bcad [3][2][4][2] = xBCAD  <7,5,9,5>;
        J64abcd [3][2][4][3] = jABCD  <7,5,9,7>;
        J64cdab [3][2][4][3] = jCDAB  <7,5,9,7>;
        X64acbd [3][2][4][3] = xACBD  <7,5,9,7>;
        X64adbc [3][2][4][3] = xADBC  <7,5,9,7>;
        X64bdac [3][2][4][3] = xBDAC  <7,5,9,7>;
        X64bcad [3][2][4][3] = xBCAD  <7,5,9,7>;
        J64abcd [3][2][4][4] = jABCD  <7,5,9,9>;
        J64cdab [3][2][4][4] = jCDAB  <7,5,9,9>;
        X64acbd [3][2][4][4] = xACBD  <7,5,9,9>;
        X64adbc [3][2][4][4] = xADBC  <7,5,9,9>;
        X64bdac [3][2][4][4] = xBDAC  <7,5,9,9>;
        X64bcad [3][2][4][4] = xBCAD  <7,5,9,9>;

        J64abcd [3][3][0][0] = jABCD  <7,7,1,1>;
        J64cdab [3][3][0][0] = jCDAB  <7,7,1,1>;
        X64acbd [3][3][0][0] = xACBD  <7,7,1,1>;
        X64adbc [3][3][0][0] = xADBC  <7,7,1,1>;
        X64bdac [3][3][0][0] = xBDAC  <7,7,1,1>;
        X64bcad [3][3][0][0] = xBCAD  <7,7,1,1>;
        J64abcd [3][3][1][0] = jABCD  <7,7,3,1>;
        J64cdab [3][3][1][0] = jCDAB  <7,7,3,1>;
        X64acbd [3][3][1][0] = xACBD  <7,7,3,1>;
        X64adbc [3][3][1][0] = xADBC  <7,7,3,1>;
        X64bdac [3][3][1][0] = xBDAC  <7,7,3,1>;
        X64bcad [3][3][1][0] = xBCAD  <7,7,3,1>;
        J64abcd [3][3][1][1] = jABCD  <7,7,3,3>;
        J64cdab [3][3][1][1] = jCDAB  <7,7,3,3>;
        X64acbd [3][3][1][1] = xACBD  <7,7,3,3>;
        X64adbc [3][3][1][1] = xADBC  <7,7,3,3>;
        X64bdac [3][3][1][1] = xBDAC  <7,7,3,3>;
        X64bcad [3][3][1][1] = xBCAD  <7,7,3,3>;
        J64abcd [3][3][2][0] = jABCD  <7,7,5,1>;
        J64cdab [3][3][2][0] = jCDAB  <7,7,5,1>;
        X64acbd [3][3][2][0] = xACBD  <7,7,5,1>;
        X64adbc [3][3][2][0] = xADBC  <7,7,5,1>;
        X64bdac [3][3][2][0] = xBDAC  <7,7,5,1>;
        X64bcad [3][3][2][0] = xBCAD  <7,7,5,1>;
        J64abcd [3][3][2][1] = jABCD  <7,7,5,3>;
        J64cdab [3][3][2][1] = jCDAB  <7,7,5,3>;
        X64acbd [3][3][2][1] = xACBD  <7,7,5,3>;
        X64adbc [3][3][2][1] = xADBC  <7,7,5,3>;
        X64bdac [3][3][2][1] = xBDAC  <7,7,5,3>;
        X64bcad [3][3][2][1] = xBCAD  <7,7,5,3>;
        J64abcd [3][3][2][2] = jABCD  <7,7,5,5>;
        J64cdab [3][3][2][2] = jCDAB  <7,7,5,5>;
        X64acbd [3][3][2][2] = xACBD  <7,7,5,5>;
        X64adbc [3][3][2][2] = xADBC  <7,7,5,5>;
        X64bdac [3][3][2][2] = xBDAC  <7,7,5,5>;
        X64bcad [3][3][2][2] = xBCAD  <7,7,5,5>;
        J64abcd [3][3][3][0] = jABCD  <7,7,7,1>;
        J64cdab [3][3][3][0] = jCDAB  <7,7,7,1>;
        X64acbd [3][3][3][0] = xACBD  <7,7,7,1>;
        X64adbc [3][3][3][0] = xADBC  <7,7,7,1>;
        X64bdac [3][3][3][0] = xBDAC  <7,7,7,1>;
        X64bcad [3][3][3][0] = xBCAD  <7,7,7,1>;
        J64abcd [3][3][3][1] = jABCD  <7,7,7,3>;
        J64cdab [3][3][3][1] = jCDAB  <7,7,7,3>;
        X64acbd [3][3][3][1] = xACBD  <7,7,7,3>;
        X64adbc [3][3][3][1] = xADBC  <7,7,7,3>;
        X64bdac [3][3][3][1] = xBDAC  <7,7,7,3>;
        X64bcad [3][3][3][1] = xBCAD  <7,7,7,3>;
        J64abcd [3][3][3][2] = jABCD  <7,7,7,5>;
        J64cdab [3][3][3][2] = jCDAB  <7,7,7,5>;
        X64acbd [3][3][3][2] = xACBD  <7,7,7,5>;
        X64adbc [3][3][3][2] = xADBC  <7,7,7,5>;
        X64bdac [3][3][3][2] = xBDAC  <7,7,7,5>;
        X64bcad [3][3][3][2] = xBCAD  <7,7,7,5>;
        J64abcd [3][3][3][3] = jABCD  <7,7,7,7>;
        J64cdab [3][3][3][3] = jCDAB  <7,7,7,7>;
        X64acbd [3][3][3][3] = xACBD  <7,7,7,7>;
        X64adbc [3][3][3][3] = xADBC  <7,7,7,7>;
        X64bdac [3][3][3][3] = xBDAC  <7,7,7,7>;
        X64bcad [3][3][3][3] = xBCAD  <7,7,7,7>;
        J64abcd [3][3][4][0] = jABCD  <7,7,9,1>;
        J64cdab [3][3][4][0] = jCDAB  <7,7,9,1>;
        X64acbd [3][3][4][0] = xACBD  <7,7,9,1>;
        X64adbc [3][3][4][0] = xADBC  <7,7,9,1>;
        X64bdac [3][3][4][0] = xBDAC  <7,7,9,1>;
        X64bcad [3][3][4][0] = xBCAD  <7,7,9,1>;
        J64abcd [3][3][4][1] = jABCD  <7,7,9,3>;
        J64cdab [3][3][4][1] = jCDAB  <7,7,9,3>;
        X64acbd [3][3][4][1] = xACBD  <7,7,9,3>;
        X64adbc [3][3][4][1] = xADBC  <7,7,9,3>;
        X64bdac [3][3][4][1] = xBDAC  <7,7,9,3>;
        X64bcad [3][3][4][1] = xBCAD  <7,7,9,3>;
        J64abcd [3][3][4][2] = jABCD  <7,7,9,5>;
        J64cdab [3][3][4][2] = jCDAB  <7,7,9,5>;
        X64acbd [3][3][4][2] = xACBD  <7,7,9,5>;
        X64adbc [3][3][4][2] = xADBC  <7,7,9,5>;
        X64bdac [3][3][4][2] = xBDAC  <7,7,9,5>;
        X64bcad [3][3][4][2] = xBCAD  <7,7,9,5>;
        J64abcd [3][3][4][3] = jABCD  <7,7,9,7>;
        J64cdab [3][3][4][3] = jCDAB  <7,7,9,7>;
        X64acbd [3][3][4][3] = xACBD  <7,7,9,7>;
        X64adbc [3][3][4][3] = xADBC  <7,7,9,7>;
        X64bdac [3][3][4][3] = xBDAC  <7,7,9,7>;
        X64bcad [3][3][4][3] = xBCAD  <7,7,9,7>;
        J64abcd [3][3][4][4] = jABCD  <7,7,9,9>;
        J64cdab [3][3][4][4] = jCDAB  <7,7,9,9>;
        X64acbd [3][3][4][4] = xACBD  <7,7,9,9>;
        X64adbc [3][3][4][4] = xADBC  <7,7,9,9>;
        X64bdac [3][3][4][4] = xBDAC  <7,7,9,9>;
        X64bcad [3][3][4][4] = xBCAD  <7,7,9,9>;

        J64abcd [4][0][0][0] = jABCD  <9,1,1,1>;
        J64cdab [4][0][0][0] = jCDAB  <9,1,1,1>;
        X64acbd [4][0][0][0] = xACBD  <9,1,1,1>;
        X64adbc [4][0][0][0] = xADBC  <9,1,1,1>;
        X64bdac [4][0][0][0] = xBDAC  <9,1,1,1>;
        X64bcad [4][0][0][0] = xBCAD  <9,1,1,1>;
        J64abcd [4][0][1][0] = jABCD  <9,1,3,1>;
        J64cdab [4][0][1][0] = jCDAB  <9,1,3,1>;
        X64acbd [4][0][1][0] = xACBD  <9,1,3,1>;
        X64adbc [4][0][1][0] = xADBC  <9,1,3,1>;
        X64bdac [4][0][1][0] = xBDAC  <9,1,3,1>;
        X64bcad [4][0][1][0] = xBCAD  <9,1,3,1>;
        J64abcd [4][0][1][1] = jABCD  <9,1,3,3>;
        J64cdab [4][0][1][1] = jCDAB  <9,1,3,3>;
        X64acbd [4][0][1][1] = xACBD  <9,1,3,3>;
        X64adbc [4][0][1][1] = xADBC  <9,1,3,3>;
        X64bdac [4][0][1][1] = xBDAC  <9,1,3,3>;
        X64bcad [4][0][1][1] = xBCAD  <9,1,3,3>;
        J64abcd [4][0][2][0] = jABCD  <9,1,5,1>;
        J64cdab [4][0][2][0] = jCDAB  <9,1,5,1>;
        X64acbd [4][0][2][0] = xACBD  <9,1,5,1>;
        X64adbc [4][0][2][0] = xADBC  <9,1,5,1>;
        X64bdac [4][0][2][0] = xBDAC  <9,1,5,1>;
        X64bcad [4][0][2][0] = xBCAD  <9,1,5,1>;
        J64abcd [4][0][2][1] = jABCD  <9,1,5,3>;
        J64cdab [4][0][2][1] = jCDAB  <9,1,5,3>;
        X64acbd [4][0][2][1] = xACBD  <9,1,5,3>;
        X64adbc [4][0][2][1] = xADBC  <9,1,5,3>;
        X64bdac [4][0][2][1] = xBDAC  <9,1,5,3>;
        X64bcad [4][0][2][1] = xBCAD  <9,1,5,3>;
        J64abcd [4][0][2][2] = jABCD  <9,1,5,5>;
        J64cdab [4][0][2][2] = jCDAB  <9,1,5,5>;
        X64acbd [4][0][2][2] = xACBD  <9,1,5,5>;
        X64adbc [4][0][2][2] = xADBC  <9,1,5,5>;
        X64bdac [4][0][2][2] = xBDAC  <9,1,5,5>;
        X64bcad [4][0][2][2] = xBCAD  <9,1,5,5>;
        J64abcd [4][0][3][0] = jABCD  <9,1,7,1>;
        J64cdab [4][0][3][0] = jCDAB  <9,1,7,1>;
        X64acbd [4][0][3][0] = xACBD  <9,1,7,1>;
        X64adbc [4][0][3][0] = xADBC  <9,1,7,1>;
        X64bdac [4][0][3][0] = xBDAC  <9,1,7,1>;
        X64bcad [4][0][3][0] = xBCAD  <9,1,7,1>;
        J64abcd [4][0][3][1] = jABCD  <9,1,7,3>;
        J64cdab [4][0][3][1] = jCDAB  <9,1,7,3>;
        X64acbd [4][0][3][1] = xACBD  <9,1,7,3>;
        X64adbc [4][0][3][1] = xADBC  <9,1,7,3>;
        X64bdac [4][0][3][1] = xBDAC  <9,1,7,3>;
        X64bcad [4][0][3][1] = xBCAD  <9,1,7,3>;
        J64abcd [4][0][3][2] = jABCD  <9,1,7,5>;
        J64cdab [4][0][3][2] = jCDAB  <9,1,7,5>;
        X64acbd [4][0][3][2] = xACBD  <9,1,7,5>;
        X64adbc [4][0][3][2] = xADBC  <9,1,7,5>;
        X64bdac [4][0][3][2] = xBDAC  <9,1,7,5>;
        X64bcad [4][0][3][2] = xBCAD  <9,1,7,5>;
        J64abcd [4][0][3][3] = jABCD  <9,1,7,7>;
        J64cdab [4][0][3][3] = jCDAB  <9,1,7,7>;
        X64acbd [4][0][3][3] = xACBD  <9,1,7,7>;
        X64adbc [4][0][3][3] = xADBC  <9,1,7,7>;
        X64bdac [4][0][3][3] = xBDAC  <9,1,7,7>;
        X64bcad [4][0][3][3] = xBCAD  <9,1,7,7>;
        J64abcd [4][0][4][0] = jABCD  <9,1,9,1>;
        J64cdab [4][0][4][0] = jCDAB  <9,1,9,1>;
        X64acbd [4][0][4][0] = xACBD  <9,1,9,1>;
        X64adbc [4][0][4][0] = xADBC  <9,1,9,1>;
        X64bdac [4][0][4][0] = xBDAC  <9,1,9,1>;
        X64bcad [4][0][4][0] = xBCAD  <9,1,9,1>;
        J64abcd [4][0][4][1] = jABCD  <9,1,9,3>;
        J64cdab [4][0][4][1] = jCDAB  <9,1,9,3>;
        X64acbd [4][0][4][1] = xACBD  <9,1,9,3>;
        X64adbc [4][0][4][1] = xADBC  <9,1,9,3>;
        X64bdac [4][0][4][1] = xBDAC  <9,1,9,3>;
        X64bcad [4][0][4][1] = xBCAD  <9,1,9,3>;
        J64abcd [4][0][4][2] = jABCD  <9,1,9,5>;
        J64cdab [4][0][4][2] = jCDAB  <9,1,9,5>;
        X64acbd [4][0][4][2] = xACBD  <9,1,9,5>;
        X64adbc [4][0][4][2] = xADBC  <9,1,9,5>;
        X64bdac [4][0][4][2] = xBDAC  <9,1,9,5>;
        X64bcad [4][0][4][2] = xBCAD  <9,1,9,5>;
        J64abcd [4][0][4][3] = jABCD  <9,1,9,7>;
        J64cdab [4][0][4][3] = jCDAB  <9,1,9,7>;
        X64acbd [4][0][4][3] = xACBD  <9,1,9,7>;
        X64adbc [4][0][4][3] = xADBC  <9,1,9,7>;
        X64bdac [4][0][4][3] = xBDAC  <9,1,9,7>;
        X64bcad [4][0][4][3] = xBCAD  <9,1,9,7>;
        J64abcd [4][0][4][4] = jABCD  <9,1,9,9>;
        J64cdab [4][0][4][4] = jCDAB  <9,1,9,9>;
        X64acbd [4][0][4][4] = xACBD  <9,1,9,9>;
        X64adbc [4][0][4][4] = xADBC  <9,1,9,9>;
        X64bdac [4][0][4][4] = xBDAC  <9,1,9,9>;
        X64bcad [4][0][4][4] = xBCAD  <9,1,9,9>;

        J64abcd [4][1][0][0] = jABCD  <9,3,1,1>;
        J64cdab [4][1][0][0] = jCDAB  <9,3,1,1>;
        X64acbd [4][1][0][0] = xACBD  <9,3,1,1>;
        X64adbc [4][1][0][0] = xADBC  <9,3,1,1>;
        X64bdac [4][1][0][0] = xBDAC  <9,3,1,1>;
        X64bcad [4][1][0][0] = xBCAD  <9,3,1,1>;
        J64abcd [4][1][1][0] = jABCD  <9,3,3,1>;
        J64cdab [4][1][1][0] = jCDAB  <9,3,3,1>;
        X64acbd [4][1][1][0] = xACBD  <9,3,3,1>;
        X64adbc [4][1][1][0] = xADBC  <9,3,3,1>;
        X64bdac [4][1][1][0] = xBDAC  <9,3,3,1>;
        X64bcad [4][1][1][0] = xBCAD  <9,3,3,1>;
        J64abcd [4][1][1][1] = jABCD  <9,3,3,3>;
        J64cdab [4][1][1][1] = jCDAB  <9,3,3,3>;
        X64acbd [4][1][1][1] = xACBD  <9,3,3,3>;
        X64adbc [4][1][1][1] = xADBC  <9,3,3,3>;
        X64bdac [4][1][1][1] = xBDAC  <9,3,3,3>;
        X64bcad [4][1][1][1] = xBCAD  <9,3,3,3>;
        J64abcd [4][1][2][0] = jABCD  <9,3,5,1>;
        J64cdab [4][1][2][0] = jCDAB  <9,3,5,1>;
        X64acbd [4][1][2][0] = xACBD  <9,3,5,1>;
        X64adbc [4][1][2][0] = xADBC  <9,3,5,1>;
        X64bdac [4][1][2][0] = xBDAC  <9,3,5,1>;
        X64bcad [4][1][2][0] = xBCAD  <9,3,5,1>;
        J64abcd [4][1][2][1] = jABCD  <9,3,5,3>;
        J64cdab [4][1][2][1] = jCDAB  <9,3,5,3>;
        X64acbd [4][1][2][1] = xACBD  <9,3,5,3>;
        X64adbc [4][1][2][1] = xADBC  <9,3,5,3>;
        X64bdac [4][1][2][1] = xBDAC  <9,3,5,3>;
        X64bcad [4][1][2][1] = xBCAD  <9,3,5,3>;
        J64abcd [4][1][2][2] = jABCD  <9,3,5,5>;
        J64cdab [4][1][2][2] = jCDAB  <9,3,5,5>;
        X64acbd [4][1][2][2] = xACBD  <9,3,5,5>;
        X64adbc [4][1][2][2] = xADBC  <9,3,5,5>;
        X64bdac [4][1][2][2] = xBDAC  <9,3,5,5>;
        X64bcad [4][1][2][2] = xBCAD  <9,3,5,5>;
        J64abcd [4][1][3][0] = jABCD  <9,3,7,1>;
        J64cdab [4][1][3][0] = jCDAB  <9,3,7,1>;
        X64acbd [4][1][3][0] = xACBD  <9,3,7,1>;
        X64adbc [4][1][3][0] = xADBC  <9,3,7,1>;
        X64bdac [4][1][3][0] = xBDAC  <9,3,7,1>;
        X64bcad [4][1][3][0] = xBCAD  <9,3,7,1>;
        J64abcd [4][1][3][1] = jABCD  <9,3,7,3>;
        J64cdab [4][1][3][1] = jCDAB  <9,3,7,3>;
        X64acbd [4][1][3][1] = xACBD  <9,3,7,3>;
        X64adbc [4][1][3][1] = xADBC  <9,3,7,3>;
        X64bdac [4][1][3][1] = xBDAC  <9,3,7,3>;
        X64bcad [4][1][3][1] = xBCAD  <9,3,7,3>;
        J64abcd [4][1][3][2] = jABCD  <9,3,7,5>;
        J64cdab [4][1][3][2] = jCDAB  <9,3,7,5>;
        X64acbd [4][1][3][2] = xACBD  <9,3,7,5>;
        X64adbc [4][1][3][2] = xADBC  <9,3,7,5>;
        X64bdac [4][1][3][2] = xBDAC  <9,3,7,5>;
        X64bcad [4][1][3][2] = xBCAD  <9,3,7,5>;
        J64abcd [4][1][3][3] = jABCD  <9,3,7,7>;
        J64cdab [4][1][3][3] = jCDAB  <9,3,7,7>;
        X64acbd [4][1][3][3] = xACBD  <9,3,7,7>;
        X64adbc [4][1][3][3] = xADBC  <9,3,7,7>;
        X64bdac [4][1][3][3] = xBDAC  <9,3,7,7>;
        X64bcad [4][1][3][3] = xBCAD  <9,3,7,7>;
        J64abcd [4][1][4][0] = jABCD  <9,3,9,1>;
        J64cdab [4][1][4][0] = jCDAB  <9,3,9,1>;
        X64acbd [4][1][4][0] = xACBD  <9,3,9,1>;
        X64adbc [4][1][4][0] = xADBC  <9,3,9,1>;
        X64bdac [4][1][4][0] = xBDAC  <9,3,9,1>;
        X64bcad [4][1][4][0] = xBCAD  <9,3,9,1>;
        J64abcd [4][1][4][1] = jABCD  <9,3,9,3>;
        J64cdab [4][1][4][1] = jCDAB  <9,3,9,3>;
        X64acbd [4][1][4][1] = xACBD  <9,3,9,3>;
        X64adbc [4][1][4][1] = xADBC  <9,3,9,3>;
        X64bdac [4][1][4][1] = xBDAC  <9,3,9,3>;
        X64bcad [4][1][4][1] = xBCAD  <9,3,9,3>;
        J64abcd [4][1][4][2] = jABCD  <9,3,9,5>;
        J64cdab [4][1][4][2] = jCDAB  <9,3,9,5>;
        X64acbd [4][1][4][2] = xACBD  <9,3,9,5>;
        X64adbc [4][1][4][2] = xADBC  <9,3,9,5>;
        X64bdac [4][1][4][2] = xBDAC  <9,3,9,5>;
        X64bcad [4][1][4][2] = xBCAD  <9,3,9,5>;
        J64abcd [4][1][4][3] = jABCD  <9,3,9,7>;
        J64cdab [4][1][4][3] = jCDAB  <9,3,9,7>;
        X64acbd [4][1][4][3] = xACBD  <9,3,9,7>;
        X64adbc [4][1][4][3] = xADBC  <9,3,9,7>;
        X64bdac [4][1][4][3] = xBDAC  <9,3,9,7>;
        X64bcad [4][1][4][3] = xBCAD  <9,3,9,7>;
        J64abcd [4][1][4][4] = jABCD  <9,3,9,9>;
        J64cdab [4][1][4][4] = jCDAB  <9,3,9,9>;
        X64acbd [4][1][4][4] = xACBD  <9,3,9,9>;
        X64adbc [4][1][4][4] = xADBC  <9,3,9,9>;
        X64bdac [4][1][4][4] = xBDAC  <9,3,9,9>;
        X64bcad [4][1][4][4] = xBCAD  <9,3,9,9>;

        J64abcd [4][2][0][0] = jABCD  <9,5,1,1>;
        J64cdab [4][2][0][0] = jCDAB  <9,5,1,1>;
        X64acbd [4][2][0][0] = xACBD  <9,5,1,1>;
        X64adbc [4][2][0][0] = xADBC  <9,5,1,1>;
        X64bdac [4][2][0][0] = xBDAC  <9,5,1,1>;
        X64bcad [4][2][0][0] = xBCAD  <9,5,1,1>;
        J64abcd [4][2][1][0] = jABCD  <9,5,3,1>;
        J64cdab [4][2][1][0] = jCDAB  <9,5,3,1>;
        X64acbd [4][2][1][0] = xACBD  <9,5,3,1>;
        X64adbc [4][2][1][0] = xADBC  <9,5,3,1>;
        X64bdac [4][2][1][0] = xBDAC  <9,5,3,1>;
        X64bcad [4][2][1][0] = xBCAD  <9,5,3,1>;
        J64abcd [4][2][1][1] = jABCD  <9,5,3,3>;
        J64cdab [4][2][1][1] = jCDAB  <9,5,3,3>;
        X64acbd [4][2][1][1] = xACBD  <9,5,3,3>;
        X64adbc [4][2][1][1] = xADBC  <9,5,3,3>;
        X64bdac [4][2][1][1] = xBDAC  <9,5,3,3>;
        X64bcad [4][2][1][1] = xBCAD  <9,5,3,3>;
        J64abcd [4][2][2][0] = jABCD  <9,5,5,1>;
        J64cdab [4][2][2][0] = jCDAB  <9,5,5,1>;
        X64acbd [4][2][2][0] = xACBD  <9,5,5,1>;
        X64adbc [4][2][2][0] = xADBC  <9,5,5,1>;
        X64bdac [4][2][2][0] = xBDAC  <9,5,5,1>;
        X64bcad [4][2][2][0] = xBCAD  <9,5,5,1>;
        J64abcd [4][2][2][1] = jABCD  <9,5,5,3>;
        J64cdab [4][2][2][1] = jCDAB  <9,5,5,3>;
        X64acbd [4][2][2][1] = xACBD  <9,5,5,3>;
        X64adbc [4][2][2][1] = xADBC  <9,5,5,3>;
        X64bdac [4][2][2][1] = xBDAC  <9,5,5,3>;
        X64bcad [4][2][2][1] = xBCAD  <9,5,5,3>;
        J64abcd [4][2][2][2] = jABCD  <9,5,5,5>;
        J64cdab [4][2][2][2] = jCDAB  <9,5,5,5>;
        X64acbd [4][2][2][2] = xACBD  <9,5,5,5>;
        X64adbc [4][2][2][2] = xADBC  <9,5,5,5>;
        X64bdac [4][2][2][2] = xBDAC  <9,5,5,5>;
        X64bcad [4][2][2][2] = xBCAD  <9,5,5,5>;
        J64abcd [4][2][3][0] = jABCD  <9,5,7,1>;
        J64cdab [4][2][3][0] = jCDAB  <9,5,7,1>;
        X64acbd [4][2][3][0] = xACBD  <9,5,7,1>;
        X64adbc [4][2][3][0] = xADBC  <9,5,7,1>;
        X64bdac [4][2][3][0] = xBDAC  <9,5,7,1>;
        X64bcad [4][2][3][0] = xBCAD  <9,5,7,1>;
        J64abcd [4][2][3][1] = jABCD  <9,5,7,3>;
        J64cdab [4][2][3][1] = jCDAB  <9,5,7,3>;
        X64acbd [4][2][3][1] = xACBD  <9,5,7,3>;
        X64adbc [4][2][3][1] = xADBC  <9,5,7,3>;
        X64bdac [4][2][3][1] = xBDAC  <9,5,7,3>;
        X64bcad [4][2][3][1] = xBCAD  <9,5,7,3>;
        J64abcd [4][2][3][2] = jABCD  <9,5,7,5>;
        J64cdab [4][2][3][2] = jCDAB  <9,5,7,5>;
        X64acbd [4][2][3][2] = xACBD  <9,5,7,5>;
        X64adbc [4][2][3][2] = xADBC  <9,5,7,5>;
        X64bdac [4][2][3][2] = xBDAC  <9,5,7,5>;
        X64bcad [4][2][3][2] = xBCAD  <9,5,7,5>;
        J64abcd [4][2][3][3] = jABCD  <9,5,7,7>;
        J64cdab [4][2][3][3] = jCDAB  <9,5,7,7>;
        X64acbd [4][2][3][3] = xACBD  <9,5,7,7>;
        X64adbc [4][2][3][3] = xADBC  <9,5,7,7>;
        X64bdac [4][2][3][3] = xBDAC  <9,5,7,7>;
        X64bcad [4][2][3][3] = xBCAD  <9,5,7,7>;
        J64abcd [4][2][4][0] = jABCD  <9,5,9,1>;
        J64cdab [4][2][4][0] = jCDAB  <9,5,9,1>;
        X64acbd [4][2][4][0] = xACBD  <9,5,9,1>;
        X64adbc [4][2][4][0] = xADBC  <9,5,9,1>;
        X64bdac [4][2][4][0] = xBDAC  <9,5,9,1>;
        X64bcad [4][2][4][0] = xBCAD  <9,5,9,1>;
        J64abcd [4][2][4][1] = jABCD  <9,5,9,3>;
        J64cdab [4][2][4][1] = jCDAB  <9,5,9,3>;
        X64acbd [4][2][4][1] = xACBD  <9,5,9,3>;
        X64adbc [4][2][4][1] = xADBC  <9,5,9,3>;
        X64bdac [4][2][4][1] = xBDAC  <9,5,9,3>;
        X64bcad [4][2][4][1] = xBCAD  <9,5,9,3>;
        J64abcd [4][2][4][2] = jABCD  <9,5,9,5>;
        J64cdab [4][2][4][2] = jCDAB  <9,5,9,5>;
        X64acbd [4][2][4][2] = xACBD  <9,5,9,5>;
        X64adbc [4][2][4][2] = xADBC  <9,5,9,5>;
        X64bdac [4][2][4][2] = xBDAC  <9,5,9,5>;
        X64bcad [4][2][4][2] = xBCAD  <9,5,9,5>;
        J64abcd [4][2][4][3] = jABCD  <9,5,9,7>;
        J64cdab [4][2][4][3] = jCDAB  <9,5,9,7>;
        X64acbd [4][2][4][3] = xACBD  <9,5,9,7>;
        X64adbc [4][2][4][3] = xADBC  <9,5,9,7>;
        X64bdac [4][2][4][3] = xBDAC  <9,5,9,7>;
        X64bcad [4][2][4][3] = xBCAD  <9,5,9,7>;
        J64abcd [4][2][4][4] = jABCD  <9,5,9,9>;
        J64cdab [4][2][4][4] = jCDAB  <9,5,9,9>;
        X64acbd [4][2][4][4] = xACBD  <9,5,9,9>;
        X64adbc [4][2][4][4] = xADBC  <9,5,9,9>;
        X64bdac [4][2][4][4] = xBDAC  <9,5,9,9>;
        X64bcad [4][2][4][4] = xBCAD  <9,5,9,9>;

        J64abcd [4][3][0][0] = jABCD  <9,7,1,1>;
        J64cdab [4][3][0][0] = jCDAB  <9,7,1,1>;
        X64acbd [4][3][0][0] = xACBD  <9,7,1,1>;
        X64adbc [4][3][0][0] = xADBC  <9,7,1,1>;
        X64bdac [4][3][0][0] = xBDAC  <9,7,1,1>;
        X64bcad [4][3][0][0] = xBCAD  <9,7,1,1>;
        J64abcd [4][3][1][0] = jABCD  <9,7,3,1>;
        J64cdab [4][3][1][0] = jCDAB  <9,7,3,1>;
        X64acbd [4][3][1][0] = xACBD  <9,7,3,1>;
        X64adbc [4][3][1][0] = xADBC  <9,7,3,1>;
        X64bdac [4][3][1][0] = xBDAC  <9,7,3,1>;
        X64bcad [4][3][1][0] = xBCAD  <9,7,3,1>;
        J64abcd [4][3][1][1] = jABCD  <9,7,3,3>;
        J64cdab [4][3][1][1] = jCDAB  <9,7,3,3>;
        X64acbd [4][3][1][1] = xACBD  <9,7,3,3>;
        X64adbc [4][3][1][1] = xADBC  <9,7,3,3>;
        X64bdac [4][3][1][1] = xBDAC  <9,7,3,3>;
        X64bcad [4][3][1][1] = xBCAD  <9,7,3,3>;
        J64abcd [4][3][2][0] = jABCD  <9,7,5,1>;
        J64cdab [4][3][2][0] = jCDAB  <9,7,5,1>;
        X64acbd [4][3][2][0] = xACBD  <9,7,5,1>;
        X64adbc [4][3][2][0] = xADBC  <9,7,5,1>;
        X64bdac [4][3][2][0] = xBDAC  <9,7,5,1>;
        X64bcad [4][3][2][0] = xBCAD  <9,7,5,1>;
        J64abcd [4][3][2][1] = jABCD  <9,7,5,3>;
        J64cdab [4][3][2][1] = jCDAB  <9,7,5,3>;
        X64acbd [4][3][2][1] = xACBD  <9,7,5,3>;
        X64adbc [4][3][2][1] = xADBC  <9,7,5,3>;
        X64bdac [4][3][2][1] = xBDAC  <9,7,5,3>;
        X64bcad [4][3][2][1] = xBCAD  <9,7,5,3>;
        J64abcd [4][3][2][2] = jABCD  <9,7,5,5>;
        J64cdab [4][3][2][2] = jCDAB  <9,7,5,5>;
        X64acbd [4][3][2][2] = xACBD  <9,7,5,5>;
        X64adbc [4][3][2][2] = xADBC  <9,7,5,5>;
        X64bdac [4][3][2][2] = xBDAC  <9,7,5,5>;
        X64bcad [4][3][2][2] = xBCAD  <9,7,5,5>;
        J64abcd [4][3][3][0] = jABCD  <9,7,7,1>;
        J64cdab [4][3][3][0] = jCDAB  <9,7,7,1>;
        X64acbd [4][3][3][0] = xACBD  <9,7,7,1>;
        X64adbc [4][3][3][0] = xADBC  <9,7,7,1>;
        X64bdac [4][3][3][0] = xBDAC  <9,7,7,1>;
        X64bcad [4][3][3][0] = xBCAD  <9,7,7,1>;
        J64abcd [4][3][3][1] = jABCD  <9,7,7,3>;
        J64cdab [4][3][3][1] = jCDAB  <9,7,7,3>;
        X64acbd [4][3][3][1] = xACBD  <9,7,7,3>;
        X64adbc [4][3][3][1] = xADBC  <9,7,7,3>;
        X64bdac [4][3][3][1] = xBDAC  <9,7,7,3>;
        X64bcad [4][3][3][1] = xBCAD  <9,7,7,3>;
        J64abcd [4][3][3][2] = jABCD  <9,7,7,5>;
        J64cdab [4][3][3][2] = jCDAB  <9,7,7,5>;
        X64acbd [4][3][3][2] = xACBD  <9,7,7,5>;
        X64adbc [4][3][3][2] = xADBC  <9,7,7,5>;
        X64bdac [4][3][3][2] = xBDAC  <9,7,7,5>;
        X64bcad [4][3][3][2] = xBCAD  <9,7,7,5>;
        J64abcd [4][3][3][3] = jABCD  <9,7,7,7>;
        J64cdab [4][3][3][3] = jCDAB  <9,7,7,7>;
        X64acbd [4][3][3][3] = xACBD  <9,7,7,7>;
        X64adbc [4][3][3][3] = xADBC  <9,7,7,7>;
        X64bdac [4][3][3][3] = xBDAC  <9,7,7,7>;
        X64bcad [4][3][3][3] = xBCAD  <9,7,7,7>;
        J64abcd [4][3][4][0] = jABCD  <9,7,9,1>;
        J64cdab [4][3][4][0] = jCDAB  <9,7,9,1>;
        X64acbd [4][3][4][0] = xACBD  <9,7,9,1>;
        X64adbc [4][3][4][0] = xADBC  <9,7,9,1>;
        X64bdac [4][3][4][0] = xBDAC  <9,7,9,1>;
        X64bcad [4][3][4][0] = xBCAD  <9,7,9,1>;
        J64abcd [4][3][4][1] = jABCD  <9,7,9,3>;
        J64cdab [4][3][4][1] = jCDAB  <9,7,9,3>;
        X64acbd [4][3][4][1] = xACBD  <9,7,9,3>;
        X64adbc [4][3][4][1] = xADBC  <9,7,9,3>;
        X64bdac [4][3][4][1] = xBDAC  <9,7,9,3>;
        X64bcad [4][3][4][1] = xBCAD  <9,7,9,3>;
        J64abcd [4][3][4][2] = jABCD  <9,7,9,5>;
        J64cdab [4][3][4][2] = jCDAB  <9,7,9,5>;
        X64acbd [4][3][4][2] = xACBD  <9,7,9,5>;
        X64adbc [4][3][4][2] = xADBC  <9,7,9,5>;
        X64bdac [4][3][4][2] = xBDAC  <9,7,9,5>;
        X64bcad [4][3][4][2] = xBCAD  <9,7,9,5>;
        J64abcd [4][3][4][3] = jABCD  <9,7,9,7>;
        J64cdab [4][3][4][3] = jCDAB  <9,7,9,7>;
        X64acbd [4][3][4][3] = xACBD  <9,7,9,7>;
        X64adbc [4][3][4][3] = xADBC  <9,7,9,7>;
        X64bdac [4][3][4][3] = xBDAC  <9,7,9,7>;
        X64bcad [4][3][4][3] = xBCAD  <9,7,9,7>;
        J64abcd [4][3][4][4] = jABCD  <9,7,9,9>;
        J64cdab [4][3][4][4] = jCDAB  <9,7,9,9>;
        X64acbd [4][3][4][4] = xACBD  <9,7,9,9>;
        X64adbc [4][3][4][4] = xADBC  <9,7,9,9>;
        X64bdac [4][3][4][4] = xBDAC  <9,7,9,9>;
        X64bcad [4][3][4][4] = xBCAD  <9,7,9,9>;

        J64abcd [4][4][0][0] = jABCD  <9,9,1,1>;
        J64cdab [4][4][0][0] = jCDAB  <9,9,1,1>;
        X64acbd [4][4][0][0] = xACBD  <9,9,1,1>;
        X64adbc [4][4][0][0] = xADBC  <9,9,1,1>;
        X64bdac [4][4][0][0] = xBDAC  <9,9,1,1>;
        X64bcad [4][4][0][0] = xBCAD  <9,9,1,1>;
        J64abcd [4][4][1][0] = jABCD  <9,9,3,1>;
        J64cdab [4][4][1][0] = jCDAB  <9,9,3,1>;
        X64acbd [4][4][1][0] = xACBD  <9,9,3,1>;
        X64adbc [4][4][1][0] = xADBC  <9,9,3,1>;
        X64bdac [4][4][1][0] = xBDAC  <9,9,3,1>;
        X64bcad [4][4][1][0] = xBCAD  <9,9,3,1>;
        J64abcd [4][4][1][1] = jABCD  <9,9,3,3>;
        J64cdab [4][4][1][1] = jCDAB  <9,9,3,3>;
        X64acbd [4][4][1][1] = xACBD  <9,9,3,3>;
        X64adbc [4][4][1][1] = xADBC  <9,9,3,3>;
        X64bdac [4][4][1][1] = xBDAC  <9,9,3,3>;
        X64bcad [4][4][1][1] = xBCAD  <9,9,3,3>;
        J64abcd [4][4][2][0] = jABCD  <9,9,5,1>;
        J64cdab [4][4][2][0] = jCDAB  <9,9,5,1>;
        X64acbd [4][4][2][0] = xACBD  <9,9,5,1>;
        X64adbc [4][4][2][0] = xADBC  <9,9,5,1>;
        X64bdac [4][4][2][0] = xBDAC  <9,9,5,1>;
        X64bcad [4][4][2][0] = xBCAD  <9,9,5,1>;
        J64abcd [4][4][2][1] = jABCD  <9,9,5,3>;
        J64cdab [4][4][2][1] = jCDAB  <9,9,5,3>;
        X64acbd [4][4][2][1] = xACBD  <9,9,5,3>;
        X64adbc [4][4][2][1] = xADBC  <9,9,5,3>;
        X64bdac [4][4][2][1] = xBDAC  <9,9,5,3>;
        X64bcad [4][4][2][1] = xBCAD  <9,9,5,3>;
        J64abcd [4][4][2][2] = jABCD  <9,9,5,5>;
        J64cdab [4][4][2][2] = jCDAB  <9,9,5,5>;
        X64acbd [4][4][2][2] = xACBD  <9,9,5,5>;
        X64adbc [4][4][2][2] = xADBC  <9,9,5,5>;
        X64bdac [4][4][2][2] = xBDAC  <9,9,5,5>;
        X64bcad [4][4][2][2] = xBCAD  <9,9,5,5>;
        J64abcd [4][4][3][0] = jABCD  <9,9,7,1>;
        J64cdab [4][4][3][0] = jCDAB  <9,9,7,1>;
        X64acbd [4][4][3][0] = xACBD  <9,9,7,1>;
        X64adbc [4][4][3][0] = xADBC  <9,9,7,1>;
        X64bdac [4][4][3][0] = xBDAC  <9,9,7,1>;
        X64bcad [4][4][3][0] = xBCAD  <9,9,7,1>;
        J64abcd [4][4][3][1] = jABCD  <9,9,7,3>;
        J64cdab [4][4][3][1] = jCDAB  <9,9,7,3>;
        X64acbd [4][4][3][1] = xACBD  <9,9,7,3>;
        X64adbc [4][4][3][1] = xADBC  <9,9,7,3>;
        X64bdac [4][4][3][1] = xBDAC  <9,9,7,3>;
        X64bcad [4][4][3][1] = xBCAD  <9,9,7,3>;
        J64abcd [4][4][3][2] = jABCD  <9,9,7,5>;
        J64cdab [4][4][3][2] = jCDAB  <9,9,7,5>;
        X64acbd [4][4][3][2] = xACBD  <9,9,7,5>;
        X64adbc [4][4][3][2] = xADBC  <9,9,7,5>;
        X64bdac [4][4][3][2] = xBDAC  <9,9,7,5>;
        X64bcad [4][4][3][2] = xBCAD  <9,9,7,5>;
        J64abcd [4][4][3][3] = jABCD  <9,9,7,7>;
        J64cdab [4][4][3][3] = jCDAB  <9,9,7,7>;
        X64acbd [4][4][3][3] = xACBD  <9,9,7,7>;
        X64adbc [4][4][3][3] = xADBC  <9,9,7,7>;
        X64bdac [4][4][3][3] = xBDAC  <9,9,7,7>;
        X64bcad [4][4][3][3] = xBCAD  <9,9,7,7>;
        J64abcd [4][4][4][0] = jABCD  <9,9,9,1>;
        J64cdab [4][4][4][0] = jCDAB  <9,9,9,1>;
        X64acbd [4][4][4][0] = xACBD  <9,9,9,1>;
        X64adbc [4][4][4][0] = xADBC  <9,9,9,1>;
        X64bdac [4][4][4][0] = xBDAC  <9,9,9,1>;
        X64bcad [4][4][4][0] = xBCAD  <9,9,9,1>;
        J64abcd [4][4][4][1] = jABCD  <9,9,9,3>;
        J64cdab [4][4][4][1] = jCDAB  <9,9,9,3>;
        X64acbd [4][4][4][1] = xACBD  <9,9,9,3>;
        X64adbc [4][4][4][1] = xADBC  <9,9,9,3>;
        X64bdac [4][4][4][1] = xBDAC  <9,9,9,3>;
        X64bcad [4][4][4][1] = xBCAD  <9,9,9,3>;
        J64abcd [4][4][4][2] = jABCD  <9,9,9,5>;
        J64cdab [4][4][4][2] = jCDAB  <9,9,9,5>;
        X64acbd [4][4][4][2] = xACBD  <9,9,9,5>;
        X64adbc [4][4][4][2] = xADBC  <9,9,9,5>;
        X64bdac [4][4][4][2] = xBDAC  <9,9,9,5>;
        X64bcad [4][4][4][2] = xBCAD  <9,9,9,5>;
        J64abcd [4][4][4][3] = jABCD  <9,9,9,7>;
        J64cdab [4][4][4][3] = jCDAB  <9,9,9,7>;
        X64acbd [4][4][4][3] = xACBD  <9,9,9,7>;
        X64adbc [4][4][4][3] = xADBC  <9,9,9,7>;
        X64bdac [4][4][4][3] = xBDAC  <9,9,9,7>;
        X64bcad [4][4][4][3] = xBCAD  <9,9,9,7>;
        J64abcd [4][4][4][4] = jABCD  <9,9,9,9>;
        J64cdab [4][4][4][4] = jCDAB  <9,9,9,9>;
        X64acbd [4][4][4][4] = xACBD  <9,9,9,9>;
        X64adbc [4][4][4][4] = xADBC  <9,9,9,9>;
        X64bdac [4][4][4][4] = xBDAC  <9,9,9,9>;
        X64bcad [4][4][4][4] = xBCAD  <9,9,9,9>;
    }

    //JXC NS
    {
        J64abcc [0][0][0] = jABCC  <1,1,1>;
        J64ccab [0][0][0] = jCCAB  <1,1,1>;
        X64acbc [0][0][0] = xACBD  <1,1,1,1>;
        X64bcac [0][0][0] = xBDAC  <1,1,1,1>;
        J64abcc [0][0][1] = jABCC  <1,1,3>;
        J64ccab [0][0][1] = jCCAB  <1,1,3>;
        X64acbc [0][0][1] = xACBD  <1,1,3,3>;
        X64bcac [0][0][1] = xBDAC  <1,1,3,3>;
        J64abcc [0][0][2] = jABCC  <1,1,5>;
        J64ccab [0][0][2] = jCCAB  <1,1,5>;
        X64acbc [0][0][2] = xACBD  <1,1,5,5>;
        X64bcac [0][0][2] = xBDAC  <1,1,5,5>;
        J64abcc [0][0][3] = jABCC  <1,1,7>;
        J64ccab [0][0][3] = jCCAB  <1,1,7>;
        X64acbc [0][0][3] = xACBD  <1,1,7,7>;
        X64bcac [0][0][3] = xBDAC  <1,1,7,7>;
        J64abcc [0][0][4] = jABCC  <1,1,9>;
        J64ccab [0][0][4] = jCCAB  <1,1,9>;
        X64acbc [0][0][4] = xACBD  <1,1,9,9>;
        X64bcac [0][0][4] = xBDAC  <1,1,9,9>;

        J64abcc [1][0][0] = jABCC  <3,1,1>;
        J64ccab [1][0][0] = jCCAB  <3,1,1>;
        X64acbc [1][0][0] = xACBD  <3,1,1,1>;
        X64bcac [1][0][0] = xBDAC  <3,1,1,1>;
        J64abcc [1][0][1] = jABCC  <3,1,3>;
        J64ccab [1][0][1] = jCCAB  <3,1,3>;
        X64acbc [1][0][1] = xACBD  <3,1,3,3>;
        X64bcac [1][0][1] = xBDAC  <3,1,3,3>;
        J64abcc [1][0][2] = jABCC  <3,1,5>;
        J64ccab [1][0][2] = jCCAB  <3,1,5>;
        X64acbc [1][0][2] = xACBD  <3,1,5,5>;
        X64bcac [1][0][2] = xBDAC  <3,1,5,5>;
        J64abcc [1][0][3] = jABCC  <3,1,7>;
        J64ccab [1][0][3] = jCCAB  <3,1,7>;
        X64acbc [1][0][3] = xACBD  <3,1,7,7>;
        X64bcac [1][0][3] = xBDAC  <3,1,7,7>;
        J64abcc [1][0][4] = jABCC  <3,1,9>;
        J64ccab [1][0][4] = jCCAB  <3,1,9>;
        X64acbc [1][0][4] = xACBD  <3,1,9,9>;
        X64bcac [1][0][4] = xBDAC  <3,1,9,9>;

        J64abcc [1][1][0] = jABCC  <3,3,1>;
        J64ccab [1][1][0] = jCCAB  <3,3,1>;
        X64acbc [1][1][0] = xACBD  <3,3,1,1>;
        X64bcac [1][1][0] = xBDAC  <3,3,1,1>;
        J64abcc [1][1][1] = jABCC  <3,3,3>;
        J64ccab [1][1][1] = jCCAB  <3,3,3>;
        X64acbc [1][1][1] = xACBD  <3,3,3,3>;
        X64bcac [1][1][1] = xBDAC  <3,3,3,3>;
        J64abcc [1][1][2] = jABCC  <3,3,5>;
        J64ccab [1][1][2] = jCCAB  <3,3,5>;
        X64acbc [1][1][2] = xACBD  <3,3,5,5>;
        X64bcac [1][1][2] = xBDAC  <3,3,5,5>;
        J64abcc [1][1][3] = jABCC  <3,3,7>;
        J64ccab [1][1][3] = jCCAB  <3,3,7>;
        X64acbc [1][1][3] = xACBD  <3,3,7,7>;
        X64bcac [1][1][3] = xBDAC  <3,3,7,7>;
        J64abcc [1][1][4] = jABCC  <3,3,9>;
        J64ccab [1][1][4] = jCCAB  <3,3,9>;
        X64acbc [1][1][4] = xACBD  <3,3,9,9>;
        X64bcac [1][1][4] = xBDAC  <3,3,9,9>;

        J64abcc [2][0][0] = jABCC  <5,1,1>;
        J64ccab [2][0][0] = jCCAB  <5,1,1>;
        X64acbc [2][0][0] = xACBD  <5,1,1,1>;
        X64bcac [2][0][0] = xBDAC  <5,1,1,1>;
        J64abcc [2][0][1] = jABCC  <5,1,3>;
        J64ccab [2][0][1] = jCCAB  <5,1,3>;
        X64acbc [2][0][1] = xACBD  <5,1,3,3>;
        X64bcac [2][0][1] = xBDAC  <5,1,3,3>;
        J64abcc [2][0][2] = jABCC  <5,1,5>;
        J64ccab [2][0][2] = jCCAB  <5,1,5>;
        X64acbc [2][0][2] = xACBD  <5,1,5,5>;
        X64bcac [2][0][2] = xBDAC  <5,1,5,5>;
        J64abcc [2][0][3] = jABCC  <5,1,7>;
        J64ccab [2][0][3] = jCCAB  <5,1,7>;
        X64acbc [2][0][3] = xACBD  <5,1,7,7>;
        X64bcac [2][0][3] = xBDAC  <5,1,7,7>;
        J64abcc [2][0][4] = jABCC  <5,1,9>;
        J64ccab [2][0][4] = jCCAB  <5,1,9>;
        X64acbc [2][0][4] = xACBD  <5,1,9,9>;
        X64bcac [2][0][4] = xBDAC  <5,1,9,9>;

        J64abcc [2][1][0] = jABCC  <5,3,1>;
        J64ccab [2][1][0] = jCCAB  <5,3,1>;
        X64acbc [2][1][0] = xACBD  <5,3,1,1>;
        X64bcac [2][1][0] = xBDAC  <5,3,1,1>;
        J64abcc [2][1][1] = jABCC  <5,3,3>;
        J64ccab [2][1][1] = jCCAB  <5,3,3>;
        X64acbc [2][1][1] = xACBD  <5,3,3,3>;
        X64bcac [2][1][1] = xBDAC  <5,3,3,3>;
        J64abcc [2][1][2] = jABCC  <5,3,5>;
        J64ccab [2][1][2] = jCCAB  <5,3,5>;
        X64acbc [2][1][2] = xACBD  <5,3,5,5>;
        X64bcac [2][1][2] = xBDAC  <5,3,5,5>;
        J64abcc [2][1][3] = jABCC  <5,3,7>;
        J64ccab [2][1][3] = jCCAB  <5,3,7>;
        X64acbc [2][1][3] = xACBD  <5,3,7,7>;
        X64bcac [2][1][3] = xBDAC  <5,3,7,7>;
        J64abcc [2][1][4] = jABCC  <5,3,9>;
        J64ccab [2][1][4] = jCCAB  <5,3,9>;
        X64acbc [2][1][4] = xACBD  <5,3,9,9>;
        X64bcac [2][1][4] = xBDAC  <5,3,9,9>;

        J64abcc [2][2][0] = jABCC  <5,5,1>;
        J64ccab [2][2][0] = jCCAB  <5,5,1>;
        X64acbc [2][2][0] = xACBD  <5,5,1,1>;
        X64bcac [2][2][0] = xBDAC  <5,5,1,1>;
        J64abcc [2][2][1] = jABCC  <5,5,3>;
        J64ccab [2][2][1] = jCCAB  <5,5,3>;
        X64acbc [2][2][1] = xACBD  <5,5,3,3>;
        X64bcac [2][2][1] = xBDAC  <5,5,3,3>;
        J64abcc [2][2][2] = jABCC  <5,5,5>;
        J64ccab [2][2][2] = jCCAB  <5,5,5>;
        X64acbc [2][2][2] = xACBD  <5,5,5,5>;
        X64bcac [2][2][2] = xBDAC  <5,5,5,5>;
        J64abcc [2][2][3] = jABCC  <5,5,7>;
        J64ccab [2][2][3] = jCCAB  <5,5,7>;
        X64acbc [2][2][3] = xACBD  <5,5,7,7>;
        X64bcac [2][2][3] = xBDAC  <5,5,7,7>;
        J64abcc [2][2][4] = jABCC  <5,5,9>;
        J64ccab [2][2][4] = jCCAB  <5,5,9>;
        X64acbc [2][2][4] = xACBD  <5,5,9,9>;
        X64bcac [2][2][4] = xBDAC  <5,5,9,9>;

        J64abcc [3][0][0] = jABCC  <7,1,1>;
        J64ccab [3][0][0] = jCCAB  <7,1,1>;
        X64acbc [3][0][0] = xACBD  <7,1,1,1>;
        X64bcac [3][0][0] = xBDAC  <7,1,1,1>;
        J64abcc [3][0][1] = jABCC  <7,1,3>;
        J64ccab [3][0][1] = jCCAB  <7,1,3>;
        X64acbc [3][0][1] = xACBD  <7,1,3,3>;
        X64bcac [3][0][1] = xBDAC  <7,1,3,3>;
        J64abcc [3][0][2] = jABCC  <7,1,5>;
        J64ccab [3][0][2] = jCCAB  <7,1,5>;
        X64acbc [3][0][2] = xACBD  <7,1,5,5>;
        X64bcac [3][0][2] = xBDAC  <7,1,5,5>;
        J64abcc [3][0][3] = jABCC  <7,1,7>;
        J64ccab [3][0][3] = jCCAB  <7,1,7>;
        X64acbc [3][0][3] = xACBD  <7,1,7,7>;
        X64bcac [3][0][3] = xBDAC  <7,1,7,7>;
        J64abcc [3][0][4] = jABCC  <7,1,9>;
        J64ccab [3][0][4] = jCCAB  <7,1,9>;
        X64acbc [3][0][4] = xACBD  <7,1,9,9>;
        X64bcac [3][0][4] = xBDAC  <7,1,9,9>;

        J64abcc [3][1][0] = jABCC  <7,3,1>;
        J64ccab [3][1][0] = jCCAB  <7,3,1>;
        X64acbc [3][1][0] = xACBD  <7,3,1,1>;
        X64bcac [3][1][0] = xBDAC  <7,3,1,1>;
        J64abcc [3][1][1] = jABCC  <7,3,3>;
        J64ccab [3][1][1] = jCCAB  <7,3,3>;
        X64acbc [3][1][1] = xACBD  <7,3,3,3>;
        X64bcac [3][1][1] = xBDAC  <7,3,3,3>;
        J64abcc [3][1][2] = jABCC  <7,3,5>;
        J64ccab [3][1][2] = jCCAB  <7,3,5>;
        X64acbc [3][1][2] = xACBD  <7,3,5,5>;
        X64bcac [3][1][2] = xBDAC  <7,3,5,5>;
        J64abcc [3][1][3] = jABCC  <7,3,7>;
        J64ccab [3][1][3] = jCCAB  <7,3,7>;
        X64acbc [3][1][3] = xACBD  <7,3,7,7>;
        X64bcac [3][1][3] = xBDAC  <7,3,7,7>;
        J64abcc [3][1][4] = jABCC  <7,3,9>;
        J64ccab [3][1][4] = jCCAB  <7,3,9>;
        X64acbc [3][1][4] = xACBD  <7,3,9,9>;
        X64bcac [3][1][4] = xBDAC  <7,3,9,9>;

        J64abcc [3][2][0] = jABCC  <7,5,1>;
        J64ccab [3][2][0] = jCCAB  <7,5,1>;
        X64acbc [3][2][0] = xACBD  <7,5,1,1>;
        X64bcac [3][2][0] = xBDAC  <7,5,1,1>;
        J64abcc [3][2][1] = jABCC  <7,5,3>;
        J64ccab [3][2][1] = jCCAB  <7,5,3>;
        X64acbc [3][2][1] = xACBD  <7,5,3,3>;
        X64bcac [3][2][1] = xBDAC  <7,5,3,3>;
        J64abcc [3][2][2] = jABCC  <7,5,5>;
        J64ccab [3][2][2] = jCCAB  <7,5,5>;
        X64acbc [3][2][2] = xACBD  <7,5,5,5>;
        X64bcac [3][2][2] = xBDAC  <7,5,5,5>;
        J64abcc [3][2][3] = jABCC  <7,5,7>;
        J64ccab [3][2][3] = jCCAB  <7,5,7>;
        X64acbc [3][2][3] = xACBD  <7,5,7,7>;
        X64bcac [3][2][3] = xBDAC  <7,5,7,7>;
        J64abcc [3][2][4] = jABCC  <7,5,9>;
        J64ccab [3][2][4] = jCCAB  <7,5,9>;
        X64acbc [3][2][4] = xACBD  <7,5,9,9>;
        X64bcac [3][2][4] = xBDAC  <7,5,9,9>;

        J64abcc [3][3][0] = jABCC  <7,7,1>;
        J64ccab [3][3][0] = jCCAB  <7,7,1>;
        X64acbc [3][3][0] = xACBD  <7,7,1,1>;
        X64bcac [3][3][0] = xBDAC  <7,7,1,1>;
        J64abcc [3][3][1] = jABCC  <7,7,3>;
        J64ccab [3][3][1] = jCCAB  <7,7,3>;
        X64acbc [3][3][1] = xACBD  <7,7,3,3>;
        X64bcac [3][3][1] = xBDAC  <7,7,3,3>;
        J64abcc [3][3][2] = jABCC  <7,7,5>;
        J64ccab [3][3][2] = jCCAB  <7,7,5>;
        X64acbc [3][3][2] = xACBD  <7,7,5,5>;
        X64bcac [3][3][2] = xBDAC  <7,7,5,5>;
        J64abcc [3][3][3] = jABCC  <7,7,7>;
        J64ccab [3][3][3] = jCCAB  <7,7,7>;
        X64acbc [3][3][3] = xACBD  <7,7,7,7>;
        X64bcac [3][3][3] = xBDAC  <7,7,7,7>;
        J64abcc [3][3][4] = jABCC  <7,7,9>;
        J64ccab [3][3][4] = jCCAB  <7,7,9>;
        X64acbc [3][3][4] = xACBD  <7,7,9,9>;
        X64bcac [3][3][4] = xBDAC  <7,7,9,9>;

        J64abcc [4][0][0] = jABCC  <9,1,1>;
        J64ccab [4][0][0] = jCCAB  <9,1,1>;
        X64acbc [4][0][0] = xACBD  <9,1,1,1>;
        X64bcac [4][0][0] = xBDAC  <9,1,1,1>;
        J64abcc [4][0][1] = jABCC  <9,1,3>;
        J64ccab [4][0][1] = jCCAB  <9,1,3>;
        X64acbc [4][0][1] = xACBD  <9,1,3,3>;
        X64bcac [4][0][1] = xBDAC  <9,1,3,3>;
        J64abcc [4][0][2] = jABCC  <9,1,5>;
        J64ccab [4][0][2] = jCCAB  <9,1,5>;
        X64acbc [4][0][2] = xACBD  <9,1,5,5>;
        X64bcac [4][0][2] = xBDAC  <9,1,5,5>;
        J64abcc [4][0][3] = jABCC  <9,1,7>;
        J64ccab [4][0][3] = jCCAB  <9,1,7>;
        X64acbc [4][0][3] = xACBD  <9,1,7,7>;
        X64bcac [4][0][3] = xBDAC  <9,1,7,7>;
        J64abcc [4][0][4] = jABCC  <9,1,9>;
        J64ccab [4][0][4] = jCCAB  <9,1,9>;
        X64acbc [4][0][4] = xACBD  <9,1,9,9>;
        X64bcac [4][0][4] = xBDAC  <9,1,9,9>;

        J64abcc [4][1][0] = jABCC  <9,3,1>;
        J64ccab [4][1][0] = jCCAB  <9,3,1>;
        X64acbc [4][1][0] = xACBD  <9,3,1,1>;
        X64bcac [4][1][0] = xBDAC  <9,3,1,1>;
        J64abcc [4][1][1] = jABCC  <9,3,3>;
        J64ccab [4][1][1] = jCCAB  <9,3,3>;
        X64acbc [4][1][1] = xACBD  <9,3,3,3>;
        X64bcac [4][1][1] = xBDAC  <9,3,3,3>;
        J64abcc [4][1][2] = jABCC  <9,3,5>;
        J64ccab [4][1][2] = jCCAB  <9,3,5>;
        X64acbc [4][1][2] = xACBD  <9,3,5,5>;
        X64bcac [4][1][2] = xBDAC  <9,3,5,5>;
        J64abcc [4][1][3] = jABCC  <9,3,7>;
        J64ccab [4][1][3] = jCCAB  <9,3,7>;
        X64acbc [4][1][3] = xACBD  <9,3,7,7>;
        X64bcac [4][1][3] = xBDAC  <9,3,7,7>;
        J64abcc [4][1][4] = jABCC  <9,3,9>;
        J64ccab [4][1][4] = jCCAB  <9,3,9>;
        X64acbc [4][1][4] = xACBD  <9,3,9,9>;
        X64bcac [4][1][4] = xBDAC  <9,3,9,9>;

        J64abcc [4][2][0] = jABCC  <9,5,1>;
        J64ccab [4][2][0] = jCCAB  <9,5,1>;
        X64acbc [4][2][0] = xACBD  <9,5,1,1>;
        X64bcac [4][2][0] = xBDAC  <9,5,1,1>;
        J64abcc [4][2][1] = jABCC  <9,5,3>;
        J64ccab [4][2][1] = jCCAB  <9,5,3>;
        X64acbc [4][2][1] = xACBD  <9,5,3,3>;
        X64bcac [4][2][1] = xBDAC  <9,5,3,3>;
        J64abcc [4][2][2] = jABCC  <9,5,5>;
        J64ccab [4][2][2] = jCCAB  <9,5,5>;
        X64acbc [4][2][2] = xACBD  <9,5,5,5>;
        X64bcac [4][2][2] = xBDAC  <9,5,5,5>;
        J64abcc [4][2][3] = jABCC  <9,5,7>;
        J64ccab [4][2][3] = jCCAB  <9,5,7>;
        X64acbc [4][2][3] = xACBD  <9,5,7,7>;
        X64bcac [4][2][3] = xBDAC  <9,5,7,7>;
        J64abcc [4][2][4] = jABCC  <9,5,9>;
        J64ccab [4][2][4] = jCCAB  <9,5,9>;
        X64acbc [4][2][4] = xACBD  <9,5,9,9>;
        X64bcac [4][2][4] = xBDAC  <9,5,9,9>;

        J64abcc [4][3][0] = jABCC  <9,7,1>;
        J64ccab [4][3][0] = jCCAB  <9,7,1>;
        X64acbc [4][3][0] = xACBD  <9,7,1,1>;
        X64bcac [4][3][0] = xBDAC  <9,7,1,1>;
        J64abcc [4][3][1] = jABCC  <9,7,3>;
        J64ccab [4][3][1] = jCCAB  <9,7,3>;
        X64acbc [4][3][1] = xACBD  <9,7,3,3>;
        X64bcac [4][3][1] = xBDAC  <9,7,3,3>;
        J64abcc [4][3][2] = jABCC  <9,7,5>;
        J64ccab [4][3][2] = jCCAB  <9,7,5>;
        X64acbc [4][3][2] = xACBD  <9,7,5,5>;
        X64bcac [4][3][2] = xBDAC  <9,7,5,5>;
        J64abcc [4][3][3] = jABCC  <9,7,7>;
        J64ccab [4][3][3] = jCCAB  <9,7,7>;
        X64acbc [4][3][3] = xACBD  <9,7,7,7>;
        X64bcac [4][3][3] = xBDAC  <9,7,7,7>;
        J64abcc [4][3][4] = jABCC  <9,7,9>;
        J64ccab [4][3][4] = jCCAB  <9,7,9>;
        X64acbc [4][3][4] = xACBD  <9,7,9,9>;
        X64bcac [4][3][4] = xBDAC  <9,7,9,9>;
        J64abcc [4][4][0] = jABCC  <9,9,1>;
        J64ccab [4][4][0] = jCCAB  <9,9,1>;
        X64acbc [4][4][0] = xACBD  <9,9,1,1>;
        X64bcac [4][4][0] = xBDAC  <9,9,1,1>;
        J64abcc [4][4][1] = jABCC  <9,9,3>;
        J64ccab [4][4][1] = jCCAB  <9,9,3>;
        X64acbc [4][4][1] = xACBD  <9,9,3,3>;
        X64bcac [4][4][1] = xBDAC  <9,9,3,3>;
        J64abcc [4][4][2] = jABCC  <9,9,5>;
        J64ccab [4][4][2] = jCCAB  <9,9,5>;
        X64acbc [4][4][2] = xACBD  <9,9,5,5>;
        X64bcac [4][4][2] = xBDAC  <9,9,5,5>;
        J64abcc [4][4][3] = jABCC  <9,9,7>;
        J64ccab [4][4][3] = jCCAB  <9,9,7>;
        X64acbc [4][4][3] = xACBD  <9,9,7,7>;
        X64bcac [4][4][3] = xBDAC  <9,9,7,7>;
        J64abcc [4][4][4] = jABCC  <9,9,9>;
        J64ccab [4][4][4] = jCCAB  <9,9,9>;
        X64acbc [4][4][4] = xACBD  <9,9,9,9>;
        X64bcac [4][4][4] = xBDAC  <9,9,9,9>;
    }

    //JXC SN
    {
        J64aacd [0][0][0] = jAACD  <1,1,1>;
        J64cdaa [0][0][0] = jCDAA  <1,1,1>;
        X64acad [0][0][0] = xACBD  <1,1,1,1>;
        X64adac [0][0][0] = xADBC  <1,1,1,1>;
        J64aacd [0][1][0] = jAACD  <1,3,1>;
        J64cdaa [0][1][0] = jCDAA  <1,3,1>;
        X64acad [0][1][0] = xACBD  <1,1,3,1>;
        X64adac [0][1][0] = xADBC  <1,1,3,1>;
        J64aacd [0][1][1] = jAACD  <1,3,3>;
        J64cdaa [0][1][1] = jCDAA  <1,3,3>;
        X64acad [0][1][1] = xACBD  <1,1,3,3>;
        X64adac [0][1][1] = xADBC  <1,1,3,3>;
        J64aacd [0][2][0] = jAACD  <1,5,1>;
        J64cdaa [0][2][0] = jCDAA  <1,5,1>;
        X64acad [0][2][0] = xACBD  <1,1,5,1>;
        X64adac [0][2][0] = xADBC  <1,1,5,1>;
        J64aacd [0][2][1] = jAACD  <1,5,3>;
        J64cdaa [0][2][1] = jCDAA  <1,5,3>;
        X64acad [0][2][1] = xACBD  <1,1,5,3>;
        X64adac [0][2][1] = xADBC  <1,1,5,3>;
        J64aacd [0][2][2] = jAACD  <1,5,5>;
        J64cdaa [0][2][2] = jCDAA  <1,5,5>;
        X64acad [0][2][2] = xACBD  <1,1,5,5>;
        X64adac [0][2][2] = xADBC  <1,1,5,5>;
        J64aacd [0][3][0] = jAACD  <1,7,1>;
        J64cdaa [0][3][0] = jCDAA  <1,7,1>;
        X64acad [0][3][0] = xACBD  <1,1,7,1>;
        X64adac [0][3][0] = xADBC  <1,1,7,1>;
        J64aacd [0][3][1] = jAACD  <1,7,3>;
        J64cdaa [0][3][1] = jCDAA  <1,7,3>;
        X64acad [0][3][1] = xACBD  <1,1,7,3>;
        X64adac [0][3][1] = xADBC  <1,1,7,3>;
        J64aacd [0][3][2] = jAACD  <1,7,5>;
        J64cdaa [0][3][2] = jCDAA  <1,7,5>;
        X64acad [0][3][2] = xACBD  <1,1,7,5>;
        X64adac [0][3][2] = xADBC  <1,1,7,5>;
        J64aacd [0][3][3] = jAACD  <1,7,7>;
        J64cdaa [0][3][3] = jCDAA  <1,7,7>;
        X64acad [0][3][3] = xACBD  <1,1,7,7>;
        X64adac [0][3][3] = xADBC  <1,1,7,7>;
        J64aacd [0][4][0] = jAACD  <1,9,1>;
        J64cdaa [0][4][0] = jCDAA  <1,9,1>;
        X64acad [0][4][0] = xACBD  <1,1,9,1>;
        X64adac [0][4][0] = xADBC  <1,1,9,1>;
        J64aacd [0][4][1] = jAACD  <1,9,3>;
        J64cdaa [0][4][1] = jCDAA  <1,9,3>;
        X64acad [0][4][1] = xACBD  <1,1,9,3>;
        X64adac [0][4][1] = xADBC  <1,1,9,3>;
        J64aacd [0][4][2] = jAACD  <1,9,5>;
        J64cdaa [0][4][2] = jCDAA  <1,9,5>;
        X64acad [0][4][2] = xACBD  <1,1,9,5>;
        X64adac [0][4][2] = xADBC  <1,1,9,5>;
        J64aacd [0][4][3] = jAACD  <1,9,7>;
        J64cdaa [0][4][3] = jCDAA  <1,9,7>;
        X64acad [0][4][3] = xACBD  <1,1,9,7>;
        X64adac [0][4][3] = xADBC  <1,1,9,7>;
        J64aacd [0][4][4] = jAACD  <1,9,9>;
        J64cdaa [0][4][4] = jCDAA  <1,9,9>;
        X64acad [0][4][4] = xACBD  <1,1,9,9>;
        X64adac [0][4][4] = xADBC  <1,1,9,9>;
        J64aacd [1][0][0] = jAACD  <3,1,1>;
        J64cdaa [1][0][0] = jCDAA  <3,1,1>;
        X64acad [1][0][0] = xACBD  <3,3,1,1>;
        X64adac [1][0][0] = xADBC  <3,3,1,1>;
        J64aacd [1][1][0] = jAACD  <3,3,1>;
        J64cdaa [1][1][0] = jCDAA  <3,3,1>;
        X64acad [1][1][0] = xACBD  <3,3,3,1>;
        X64adac [1][1][0] = xADBC  <3,3,3,1>;
        J64aacd [1][1][1] = jAACD  <3,3,3>;
        J64cdaa [1][1][1] = jCDAA  <3,3,3>;
        X64acad [1][1][1] = xACBD  <3,3,3,3>;
        X64adac [1][1][1] = xADBC  <3,3,3,3>;
        J64aacd [1][2][0] = jAACD  <3,5,1>;
        J64cdaa [1][2][0] = jCDAA  <3,5,1>;
        X64acad [1][2][0] = xACBD  <3,3,5,1>;
        X64adac [1][2][0] = xADBC  <3,3,5,1>;
        J64aacd [1][2][1] = jAACD  <3,5,3>;
        J64cdaa [1][2][1] = jCDAA  <3,5,3>;
        X64acad [1][2][1] = xACBD  <3,3,5,3>;
        X64adac [1][2][1] = xADBC  <3,3,5,3>;
        J64aacd [1][2][2] = jAACD  <3,5,5>;
        J64cdaa [1][2][2] = jCDAA  <3,5,5>;
        X64acad [1][2][2] = xACBD  <3,3,5,5>;
        X64adac [1][2][2] = xADBC  <3,3,5,5>;
        J64aacd [1][3][0] = jAACD  <3,7,1>;
        J64cdaa [1][3][0] = jCDAA  <3,7,1>;
        X64acad [1][3][0] = xACBD  <3,3,7,1>;
        X64adac [1][3][0] = xADBC  <3,3,7,1>;
        J64aacd [1][3][1] = jAACD  <3,7,3>;
        J64cdaa [1][3][1] = jCDAA  <3,7,3>;
        X64acad [1][3][1] = xACBD  <3,3,7,3>;
        X64adac [1][3][1] = xADBC  <3,3,7,3>;
        J64aacd [1][3][2] = jAACD  <3,7,5>;
        J64cdaa [1][3][2] = jCDAA  <3,7,5>;
        X64acad [1][3][2] = xACBD  <3,3,7,5>;
        X64adac [1][3][2] = xADBC  <3,3,7,5>;
        J64aacd [1][3][3] = jAACD  <3,7,7>;
        J64cdaa [1][3][3] = jCDAA  <3,7,7>;
        X64acad [1][3][3] = xACBD  <3,3,7,7>;
        X64adac [1][3][3] = xADBC  <3,3,7,7>;
        J64aacd [1][4][0] = jAACD  <3,9,1>;
        J64cdaa [1][4][0] = jCDAA  <3,9,1>;
        X64acad [1][4][0] = xACBD  <3,3,9,1>;
        X64adac [1][4][0] = xADBC  <3,3,9,1>;
        J64aacd [1][4][1] = jAACD  <3,9,3>;
        J64cdaa [1][4][1] = jCDAA  <3,9,3>;
        X64acad [1][4][1] = xACBD  <3,3,9,3>;
        X64adac [1][4][1] = xADBC  <3,3,9,3>;
        J64aacd [1][4][2] = jAACD  <3,9,5>;
        J64cdaa [1][4][2] = jCDAA  <3,9,5>;
        X64acad [1][4][2] = xACBD  <3,3,9,5>;
        X64adac [1][4][2] = xADBC  <3,3,9,5>;
        J64aacd [1][4][3] = jAACD  <3,9,7>;
        J64cdaa [1][4][3] = jCDAA  <3,9,7>;
        X64acad [1][4][3] = xACBD  <3,3,9,7>;
        X64adac [1][4][3] = xADBC  <3,3,9,7>;
        J64aacd [1][4][4] = jAACD  <3,9,9>;
        J64cdaa [1][4][4] = jCDAA  <3,9,9>;
        X64acad [1][4][4] = xACBD  <3,3,9,9>;
        X64adac [1][4][4] = xADBC  <3,3,9,9>;
        J64aacd [2][0][0] = jAACD  <5,1,1>;
        J64cdaa [2][0][0] = jCDAA  <5,1,1>;
        X64acad [2][0][0] = xACBD  <5,5,1,1>;
        X64adac [2][0][0] = xADBC  <5,5,1,1>;
        J64aacd [2][1][0] = jAACD  <5,3,1>;
        J64cdaa [2][1][0] = jCDAA  <5,3,1>;
        X64acad [2][1][0] = xACBD  <5,5,3,1>;
        X64adac [2][1][0] = xADBC  <5,5,3,1>;
        J64aacd [2][1][1] = jAACD  <5,3,3>;
        J64cdaa [2][1][1] = jCDAA  <5,3,3>;
        X64acad [2][1][1] = xACBD  <5,5,3,3>;
        X64adac [2][1][1] = xADBC  <5,5,3,3>;
        J64aacd [2][2][0] = jAACD  <5,5,1>;
        J64cdaa [2][2][0] = jCDAA  <5,5,1>;
        X64acad [2][2][0] = xACBD  <5,5,5,1>;
        X64adac [2][2][0] = xADBC  <5,5,5,1>;
        J64aacd [2][2][1] = jAACD  <5,5,3>;
        J64cdaa [2][2][1] = jCDAA  <5,5,3>;
        X64acad [2][2][1] = xACBD  <5,5,5,3>;
        X64adac [2][2][1] = xADBC  <5,5,5,3>;
        J64aacd [2][2][2] = jAACD  <5,5,5>;
        J64cdaa [2][2][2] = jCDAA  <5,5,5>;
        X64acad [2][2][2] = xACBD  <5,5,5,5>;
        X64adac [2][2][2] = xADBC  <5,5,5,5>;
        J64aacd [2][3][0] = jAACD  <5,7,1>;
        J64cdaa [2][3][0] = jCDAA  <5,7,1>;
        X64acad [2][3][0] = xACBD  <5,5,7,1>;
        X64adac [2][3][0] = xADBC  <5,5,7,1>;
        J64aacd [2][3][1] = jAACD  <5,7,3>;
        J64cdaa [2][3][1] = jCDAA  <5,7,3>;
        X64acad [2][3][1] = xACBD  <5,5,7,3>;
        X64adac [2][3][1] = xADBC  <5,5,7,3>;
        J64aacd [2][3][2] = jAACD  <5,7,5>;
        J64cdaa [2][3][2] = jCDAA  <5,7,5>;
        X64acad [2][3][2] = xACBD  <5,5,7,5>;
        X64adac [2][3][2] = xADBC  <5,5,7,5>;
        J64aacd [2][3][3] = jAACD  <5,7,7>;
        J64cdaa [2][3][3] = jCDAA  <5,7,7>;
        X64acad [2][3][3] = xACBD  <5,5,7,7>;
        X64adac [2][3][3] = xADBC  <5,5,7,7>;
        J64aacd [2][4][0] = jAACD  <5,9,1>;
        J64cdaa [2][4][0] = jCDAA  <5,9,1>;
        X64acad [2][4][0] = xACBD  <5,5,9,1>;
        X64adac [2][4][0] = xADBC  <5,5,9,1>;
        J64aacd [2][4][1] = jAACD  <5,9,3>;
        J64cdaa [2][4][1] = jCDAA  <5,9,3>;
        X64acad [2][4][1] = xACBD  <5,5,9,3>;
        X64adac [2][4][1] = xADBC  <5,5,9,3>;
        J64aacd [2][4][2] = jAACD  <5,9,5>;
        J64cdaa [2][4][2] = jCDAA  <5,9,5>;
        X64acad [2][4][2] = xACBD  <5,5,9,5>;
        X64adac [2][4][2] = xADBC  <5,5,9,5>;
        J64aacd [2][4][3] = jAACD  <5,9,7>;
        J64cdaa [2][4][3] = jCDAA  <5,9,7>;
        X64acad [2][4][3] = xACBD  <5,5,9,7>;
        X64adac [2][4][3] = xADBC  <5,5,9,7>;
        J64aacd [2][4][4] = jAACD  <5,9,9>;
        J64cdaa [2][4][4] = jCDAA  <5,9,9>;
        X64acad [2][4][4] = xACBD  <5,5,9,9>;
        X64adac [2][4][4] = xADBC  <5,5,9,9>;
        J64aacd [3][0][0] = jAACD  <7,1,1>;
        J64cdaa [3][0][0] = jCDAA  <7,1,1>;
        X64acad [3][0][0] = xACBD  <7,7,1,1>;
        X64adac [3][0][0] = xADBC  <7,7,1,1>;
        J64aacd [3][1][0] = jAACD  <7,3,1>;
        J64cdaa [3][1][0] = jCDAA  <7,3,1>;
        X64acad [3][1][0] = xACBD  <7,7,3,1>;
        X64adac [3][1][0] = xADBC  <7,7,3,1>;
        J64aacd [3][1][1] = jAACD  <7,3,3>;
        J64cdaa [3][1][1] = jCDAA  <7,3,3>;
        X64acad [3][1][1] = xACBD  <7,7,3,3>;
        X64adac [3][1][1] = xADBC  <7,7,3,3>;
        J64aacd [3][2][0] = jAACD  <7,5,1>;
        J64cdaa [3][2][0] = jCDAA  <7,5,1>;
        X64acad [3][2][0] = xACBD  <7,7,5,1>;
        X64adac [3][2][0] = xADBC  <7,7,5,1>;
        J64aacd [3][2][1] = jAACD  <7,5,3>;
        J64cdaa [3][2][1] = jCDAA  <7,5,3>;
        X64acad [3][2][1] = xACBD  <7,7,5,3>;
        X64adac [3][2][1] = xADBC  <7,7,5,3>;
        J64aacd [3][2][2] = jAACD  <7,5,5>;
        J64cdaa [3][2][2] = jCDAA  <7,5,5>;
        X64acad [3][2][2] = xACBD  <7,7,5,5>;
        X64adac [3][2][2] = xADBC  <7,7,5,5>;
        J64aacd [3][3][0] = jAACD  <7,7,1>;
        J64cdaa [3][3][0] = jCDAA  <7,7,1>;
        X64acad [3][3][0] = xACBD  <7,7,7,1>;
        X64adac [3][3][0] = xADBC  <7,7,7,1>;
        J64aacd [3][3][1] = jAACD  <7,7,3>;
        J64cdaa [3][3][1] = jCDAA  <7,7,3>;
        X64acad [3][3][1] = xACBD  <7,7,7,3>;
        X64adac [3][3][1] = xADBC  <7,7,7,3>;
        J64aacd [3][3][2] = jAACD  <7,7,5>;
        J64cdaa [3][3][2] = jCDAA  <7,7,5>;
        X64acad [3][3][2] = xACBD  <7,7,7,5>;
        X64adac [3][3][2] = xADBC  <7,7,7,5>;
        J64aacd [3][3][3] = jAACD  <7,7,7>;
        J64cdaa [3][3][3] = jCDAA  <7,7,7>;
        X64acad [3][3][3] = xACBD  <7,7,7,7>;
        X64adac [3][3][3] = xADBC  <7,7,7,7>;
        J64aacd [3][4][0] = jAACD  <7,9,1>;
        J64cdaa [3][4][0] = jCDAA  <7,9,1>;
        X64acad [3][4][0] = xACBD  <7,7,9,1>;
        X64adac [3][4][0] = xADBC  <7,7,9,1>;
        J64aacd [3][4][1] = jAACD  <7,9,3>;
        J64cdaa [3][4][1] = jCDAA  <7,9,3>;
        X64acad [3][4][1] = xACBD  <7,7,9,3>;
        X64adac [3][4][1] = xADBC  <7,7,9,3>;
        J64aacd [3][4][2] = jAACD  <7,9,5>;
        J64cdaa [3][4][2] = jCDAA  <7,9,5>;
        X64acad [3][4][2] = xACBD  <7,7,9,5>;
        X64adac [3][4][2] = xADBC  <7,7,9,5>;
        J64aacd [3][4][3] = jAACD  <7,9,7>;
        J64cdaa [3][4][3] = jCDAA  <7,9,7>;
        X64acad [3][4][3] = xACBD  <7,7,9,7>;
        X64adac [3][4][3] = xADBC  <7,7,9,7>;
        J64aacd [3][4][4] = jAACD  <7,9,9>;
        J64cdaa [3][4][4] = jCDAA  <7,9,9>;
        X64acad [3][4][4] = xACBD  <7,7,9,9>;
        X64adac [3][4][4] = xADBC  <7,7,9,9>;
        J64aacd [4][0][0] = jAACD  <9,1,1>;
        J64cdaa [4][0][0] = jCDAA  <9,1,1>;
        X64acad [4][0][0] = xACBD  <9,9,1,1>;
        X64adac [4][0][0] = xADBC  <9,9,1,1>;
        J64aacd [4][1][0] = jAACD  <9,3,1>;
        J64cdaa [4][1][0] = jCDAA  <9,3,1>;
        X64acad [4][1][0] = xACBD  <9,9,3,1>;
        X64adac [4][1][0] = xADBC  <9,9,3,1>;
        J64aacd [4][1][1] = jAACD  <9,3,3>;
        J64cdaa [4][1][1] = jCDAA  <9,3,3>;
        X64acad [4][1][1] = xACBD  <9,9,3,3>;
        X64adac [4][1][1] = xADBC  <9,9,3,3>;
        J64aacd [4][2][0] = jAACD  <9,5,1>;
        J64cdaa [4][2][0] = jCDAA  <9,5,1>;
        X64acad [4][2][0] = xACBD  <9,9,5,1>;
        X64adac [4][2][0] = xADBC  <9,9,5,1>;
        J64aacd [4][2][1] = jAACD  <9,5,3>;
        J64cdaa [4][2][1] = jCDAA  <9,5,3>;
        X64acad [4][2][1] = xACBD  <9,9,5,3>;
        X64adac [4][2][1] = xADBC  <9,9,5,3>;
        J64aacd [4][2][2] = jAACD  <9,5,5>;
        J64cdaa [4][2][2] = jCDAA  <9,5,5>;
        X64acad [4][2][2] = xACBD  <9,9,5,5>;
        X64adac [4][2][2] = xADBC  <9,9,5,5>;
        J64aacd [4][3][0] = jAACD  <9,7,1>;
        J64cdaa [4][3][0] = jCDAA  <9,7,1>;
        X64acad [4][3][0] = xACBD  <9,9,7,1>;
        X64adac [4][3][0] = xADBC  <9,9,7,1>;
        J64aacd [4][3][1] = jAACD  <9,7,3>;
        J64cdaa [4][3][1] = jCDAA  <9,7,3>;
        X64acad [4][3][1] = xACBD  <9,9,7,3>;
        X64adac [4][3][1] = xADBC  <9,9,7,3>;
        J64aacd [4][3][2] = jAACD  <9,7,5>;
        J64cdaa [4][3][2] = jCDAA  <9,7,5>;
        X64acad [4][3][2] = xACBD  <9,9,7,5>;
        X64adac [4][3][2] = xADBC  <9,9,7,5>;
        J64aacd [4][3][3] = jAACD  <9,7,7>;
        J64cdaa [4][3][3] = jCDAA  <9,7,7>;
        X64acad [4][3][3] = xACBD  <9,9,7,7>;
        X64adac [4][3][3] = xADBC  <9,9,7,7>;
        J64aacd [4][4][0] = jAACD  <9,9,1>;
        J64cdaa [4][4][0] = jCDAA  <9,9,1>;
        X64acad [4][4][0] = xACBD  <9,9,9,1>;
        X64adac [4][4][0] = xADBC  <9,9,9,1>;
        J64aacd [4][4][1] = jAACD  <9,9,3>;
        J64cdaa [4][4][1] = jCDAA  <9,9,3>;
        X64acad [4][4][1] = xACBD  <9,9,9,3>;
        X64adac [4][4][1] = xADBC  <9,9,9,3>;
        J64aacd [4][4][2] = jAACD  <9,9,5>;
        J64cdaa [4][4][2] = jCDAA  <9,9,5>;
        X64acad [4][4][2] = xACBD  <9,9,9,5>;
        X64adac [4][4][2] = xADBC  <9,9,9,5>;
        J64aacd [4][4][3] = jAACD  <9,9,7>;
        J64cdaa [4][4][3] = jCDAA  <9,9,7>;
        X64acad [4][4][3] = xACBD  <9,9,9,7>;
        X64adac [4][4][3] = xADBC  <9,9,9,7>;
        J64aacd [4][4][4] = jAACD  <9,9,9>;
        J64cdaa [4][4][4] = jCDAA  <9,9,9>;
        X64acad [4][4][4] = xACBD  <9,9,9,9>;
        X64adac [4][4][4] = xADBC  <9,9,9,9>;
    }

    //JXC SS
    {
        J64aacc [0][0] = jAACC  <1,1>;
        J64ccaa [0][0] = jCCAA  <1,1>;
        X64acac [0][0] = xACBD  <1,1,1,1>;
        J64aacc [0][1] = jAACC  <1,3>;
        J64ccaa [0][1] = jCCAA  <1,3>;
        X64acac [0][1] = xACBD  <1,1,3,3>;
        J64aacc [0][2] = jAACC  <1,5>;
        J64ccaa [0][2] = jCCAA  <1,5>;
        X64acac [0][2] = xACBD  <1,1,5,5>;
        J64aacc [0][3] = jAACC  <1,7>;
        J64ccaa [0][3] = jCCAA  <1,7>;
        X64acac [0][3] = xACBD  <1,1,7,7>;
        J64aacc [0][4] = jAACC  <1,9>;
        J64ccaa [0][4] = jCCAA  <1,9>;
        X64acac [0][4] = xACBD  <1,1,9,9>;
        J64aacc [1][0] = jAACC  <3,1>;
        J64ccaa [1][0] = jCCAA  <3,1>;
        X64acac [1][0] = xACBD  <3,3,1,1>;
        J64aacc [1][1] = jAACC  <3,3>;
        J64ccaa [1][1] = jCCAA  <3,3>;
        X64acac [1][1] = xACBD  <3,3,3,3>;
        J64aacc [1][2] = jAACC  <3,5>;
        J64ccaa [1][2] = jCCAA  <3,5>;
        X64acac [1][2] = xACBD  <3,3,5,5>;
        J64aacc [1][3] = jAACC  <3,7>;
        J64ccaa [1][3] = jCCAA  <3,7>;
        X64acac [1][3] = xACBD  <3,3,7,7>;
        J64aacc [1][4] = jAACC  <3,9>;
        J64ccaa [1][4] = jCCAA  <3,9>;
        X64acac [1][4] = xACBD  <3,3,9,9>;
        J64aacc [2][0] = jAACC  <5,1>;
        J64ccaa [2][0] = jCCAA  <5,1>;
        X64acac [2][0] = xACBD  <5,5,1,1>;
        J64aacc [2][1] = jAACC  <5,3>;
        J64ccaa [2][1] = jCCAA  <5,3>;
        X64acac [2][1] = xACBD  <5,5,3,3>;
        J64aacc [2][2] = jAACC  <5,5>;
        J64ccaa [2][2] = jCCAA  <5,5>;
        X64acac [2][2] = xACBD  <5,5,5,5>;
        J64aacc [2][3] = jAACC  <5,7>;
        J64ccaa [2][3] = jCCAA  <5,7>;
        X64acac [2][3] = xACBD  <5,5,7,7>;
        J64aacc [2][4] = jAACC  <5,9>;
        J64ccaa [2][4] = jCCAA  <5,9>;
        X64acac [2][4] = xACBD  <5,5,9,9>;
        J64aacc [3][0] = jAACC  <7,1>;
        J64ccaa [3][0] = jCCAA  <7,1>;
        X64acac [3][0] = xACBD  <7,7,1,1>;
        J64aacc [3][1] = jAACC  <7,3>;
        J64ccaa [3][1] = jCCAA  <7,3>;
        X64acac [3][1] = xACBD  <7,7,3,3>;
        J64aacc [3][2] = jAACC  <7,5>;
        J64ccaa [3][2] = jCCAA  <7,5>;
        X64acac [3][2] = xACBD  <7,7,5,5>;
        J64aacc [3][3] = jAACC  <7,7>;
        J64ccaa [3][3] = jCCAA  <7,7>;
        X64acac [3][3] = xACBD  <7,7,7,7>;
        J64aacc [3][4] = jAACC  <7,9>;
        J64ccaa [3][4] = jCCAA  <7,9>;
        X64acac [3][4] = xACBD  <7,7,9,9>;
        J64aacc [4][0] = jAACC  <9,1>;
        J64ccaa [4][0] = jCCAA  <9,1>;
        X64acac [4][0] = xACBD  <9,9,1,1>;
        J64aacc [4][1] = jAACC  <9,3>;
        J64ccaa [4][1] = jCCAA  <9,3>;
        X64acac [4][1] = xACBD  <9,9,3,3>;
        J64aacc [4][2] = jAACC  <9,5>;
        J64ccaa [4][2] = jCCAA  <9,5>;
        X64acac [4][2] = xACBD  <9,9,5,5>;
        J64aacc [4][3] = jAACC  <9,7>;
        J64ccaa [4][3] = jCCAA  <9,7>;
        X64acac [4][3] = xACBD  <9,9,7,7>;
        J64aacc [4][4] = jAACC  <9,9>;
        J64ccaa [4][4] = jCCAA  <9,9>;
        X64acac [4][4] = xACBD  <9,9,9,9>;
    }

    //JXC S
    {
        J64abab [0][0] = jABCD  <1,1,1,1>;
        X64aabb [0][0] = xAABB  <1,1>;
        X64bbaa [0][0] = xBBAA  <1,1>;
        X64abba [0][0] = xABAB  <1,1>;
        J64abab [1][0] = jABCD  <3,1,3,1>;
        X64aabb [1][0] = xAABB  <3,1>;
        X64bbaa [1][0] = xBBAA  <3,1>;
        X64abba [1][0] = xABAB  <3,1>;
        J64abab [1][1] = jABCD  <3,3,3,3>;
        X64aabb [1][1] = xAABB  <3,3>;
        X64bbaa [1][1] = xBBAA  <3,3>;
        X64abba [1][1] = xABAB  <3,3>;
        J64abab [2][0] = jABCD  <5,1,5,1>;
        X64aabb [2][0] = xAABB  <5,1>;
        X64bbaa [2][0] = xBBAA  <5,1>;
        X64abba [2][0] = xABAB  <5,1>;
        J64abab [2][1] = jABCD  <5,3,5,3>;
        X64aabb [2][1] = xAABB  <5,3>;
        X64bbaa [2][1] = xBBAA  <5,3>;
        X64abba [2][1] = xABAB  <5,3>;
        J64abab [2][2] = jABCD  <5,5,5,5>;
        X64aabb [2][2] = xAABB  <5,5>;
        X64bbaa [2][2] = xBBAA  <5,5>;
        X64abba [2][2] = xABAB  <5,5>;
        J64abab [3][0] = jABCD  <7,1,7,1>;
        X64aabb [3][0] = xAABB  <7,1>;
        X64bbaa [3][0] = xBBAA  <7,1>;
        X64abba [3][0] = xABAB  <7,1>;
        J64abab [3][1] = jABCD  <7,3,7,3>;
        X64aabb [3][1] = xAABB  <7,3>;
        X64bbaa [3][1] = xBBAA  <7,3>;
        X64abba [3][1] = xABAB  <7,3>;
        J64abab [3][2] = jABCD  <7,5,7,5>;
        X64aabb [3][2] = xAABB  <7,5>;
        X64bbaa [3][2] = xBBAA  <7,5>;
        X64abba [3][2] = xABAB  <7,5>;
        J64abab [3][3] = jABCD  <7,7,7,7>;
        X64aabb [3][3] = xAABB  <7,7>;
        X64bbaa [3][3] = xBBAA  <7,7>;
        X64abba [3][3] = xABAB  <7,7>;
        J64abab [4][0] = jABCD  <9,1,9,1>;
        X64aabb [4][0] = xAABB  <9,1>;
        X64bbaa [4][0] = xBBAA  <9,1>;
        X64abba [4][0] = xABAB  <9,1>;
        J64abab [4][1] = jABCD  <9,3,9,3>;
        X64aabb [4][1] = xAABB  <9,3>;
        X64bbaa [4][1] = xBBAA  <9,3>;
        X64abba [4][1] = xABAB  <9,3>;
        J64abab [4][2] = jABCD  <9,5,9,5>;
        X64aabb [4][2] = xAABB  <9,5>;
        X64bbaa [4][2] = xBBAA  <9,5>;
        X64abba [4][2] = xABAB  <9,5>;
        J64abab [4][3] = jABCD  <9,7,9,7>;
        X64aabb [4][3] = xAABB  <9,7>;
        X64bbaa [4][3] = xBBAA  <9,7>;
        X64abba [4][3] = xABAB  <9,7>;
        J64abab [4][4] = jABCD  <9,9,9,9>;
        X64aabb [4][4] = xAABB  <9,9>;
        X64bbaa [4][4] = xBBAA  <9,9>;
        X64abba [4][4] = xABAB  <9,9>;
    }


    // 32

    //JXC NN
    {
        J32abcd [0][0][0][0] = jABCD  <1,1,1,1>;
        J32cdab [0][0][0][0] = jCDAB  <1,1,1,1>;
        X32acbd [0][0][0][0] = xACBD  <1,1,1,1>;
        X32adbc [0][0][0][0] = xADBC  <1,1,1,1>;
        X32bdac [0][0][0][0] = xBDAC  <1,1,1,1>;
        X32bcad [0][0][0][0] = xBCAD  <1,1,1,1>;
        J32abcd [0][0][1][0] = jABCD  <1,1,3,1>;
        J32cdab [0][0][1][0] = jCDAB  <1,1,3,1>;
        X32acbd [0][0][1][0] = xACBD  <1,1,3,1>;
        X32adbc [0][0][1][0] = xADBC  <1,1,3,1>;
        X32bdac [0][0][1][0] = xBDAC  <1,1,3,1>;
        X32bcad [0][0][1][0] = xBCAD  <1,1,3,1>;
        J32abcd [0][0][1][1] = jABCD  <1,1,3,3>;
        J32cdab [0][0][1][1] = jCDAB  <1,1,3,3>;
        X32acbd [0][0][1][1] = xACBD  <1,1,3,3>;
        X32adbc [0][0][1][1] = xADBC  <1,1,3,3>;
        X32bdac [0][0][1][1] = xBDAC  <1,1,3,3>;
        X32bcad [0][0][1][1] = xBCAD  <1,1,3,3>;
        J32abcd [0][0][2][0] = jABCD  <1,1,5,1>;
        J32cdab [0][0][2][0] = jCDAB  <1,1,5,1>;
        X32acbd [0][0][2][0] = xACBD  <1,1,5,1>;
        X32adbc [0][0][2][0] = xADBC  <1,1,5,1>;
        X32bdac [0][0][2][0] = xBDAC  <1,1,5,1>;
        X32bcad [0][0][2][0] = xBCAD  <1,1,5,1>;
        J32abcd [0][0][2][1] = jABCD  <1,1,5,3>;
        J32cdab [0][0][2][1] = jCDAB  <1,1,5,3>;
        X32acbd [0][0][2][1] = xACBD  <1,1,5,3>;
        X32adbc [0][0][2][1] = xADBC  <1,1,5,3>;
        X32bdac [0][0][2][1] = xBDAC  <1,1,5,3>;
        X32bcad [0][0][2][1] = xBCAD  <1,1,5,3>;
        J32abcd [0][0][2][2] = jABCD  <1,1,5,5>;
        J32cdab [0][0][2][2] = jCDAB  <1,1,5,5>;
        X32acbd [0][0][2][2] = xACBD  <1,1,5,5>;
        X32adbc [0][0][2][2] = xADBC  <1,1,5,5>;
        X32bdac [0][0][2][2] = xBDAC  <1,1,5,5>;
        X32bcad [0][0][2][2] = xBCAD  <1,1,5,5>;
        J32abcd [0][0][3][0] = jABCD  <1,1,7,1>;
        J32cdab [0][0][3][0] = jCDAB  <1,1,7,1>;
        X32acbd [0][0][3][0] = xACBD  <1,1,7,1>;
        X32adbc [0][0][3][0] = xADBC  <1,1,7,1>;
        X32bdac [0][0][3][0] = xBDAC  <1,1,7,1>;
        X32bcad [0][0][3][0] = xBCAD  <1,1,7,1>;
        J32abcd [0][0][3][1] = jABCD  <1,1,7,3>;
        J32cdab [0][0][3][1] = jCDAB  <1,1,7,3>;
        X32acbd [0][0][3][1] = xACBD  <1,1,7,3>;
        X32adbc [0][0][3][1] = xADBC  <1,1,7,3>;
        X32bdac [0][0][3][1] = xBDAC  <1,1,7,3>;
        X32bcad [0][0][3][1] = xBCAD  <1,1,7,3>;
        J32abcd [0][0][3][2] = jABCD  <1,1,7,5>;
        J32cdab [0][0][3][2] = jCDAB  <1,1,7,5>;
        X32acbd [0][0][3][2] = xACBD  <1,1,7,5>;
        X32adbc [0][0][3][2] = xADBC  <1,1,7,5>;
        X32bdac [0][0][3][2] = xBDAC  <1,1,7,5>;
        X32bcad [0][0][3][2] = xBCAD  <1,1,7,5>;
        J32abcd [0][0][3][3] = jABCD  <1,1,7,7>;
        J32cdab [0][0][3][3] = jCDAB  <1,1,7,7>;
        X32acbd [0][0][3][3] = xACBD  <1,1,7,7>;
        X32adbc [0][0][3][3] = xADBC  <1,1,7,7>;
        X32bdac [0][0][3][3] = xBDAC  <1,1,7,7>;
        X32bcad [0][0][3][3] = xBCAD  <1,1,7,7>;
        J32abcd [0][0][4][0] = jABCD  <1,1,9,1>;
        J32cdab [0][0][4][0] = jCDAB  <1,1,9,1>;
        X32acbd [0][0][4][0] = xACBD  <1,1,9,1>;
        X32adbc [0][0][4][0] = xADBC  <1,1,9,1>;
        X32bdac [0][0][4][0] = xBDAC  <1,1,9,1>;
        X32bcad [0][0][4][0] = xBCAD  <1,1,9,1>;
        J32abcd [0][0][4][1] = jABCD  <1,1,9,3>;
        J32cdab [0][0][4][1] = jCDAB  <1,1,9,3>;
        X32acbd [0][0][4][1] = xACBD  <1,1,9,3>;
        X32adbc [0][0][4][1] = xADBC  <1,1,9,3>;
        X32bdac [0][0][4][1] = xBDAC  <1,1,9,3>;
        X32bcad [0][0][4][1] = xBCAD  <1,1,9,3>;
        J32abcd [0][0][4][2] = jABCD  <1,1,9,5>;
        J32cdab [0][0][4][2] = jCDAB  <1,1,9,5>;
        X32acbd [0][0][4][2] = xACBD  <1,1,9,5>;
        X32adbc [0][0][4][2] = xADBC  <1,1,9,5>;
        X32bdac [0][0][4][2] = xBDAC  <1,1,9,5>;
        X32bcad [0][0][4][2] = xBCAD  <1,1,9,5>;
        J32abcd [0][0][4][3] = jABCD  <1,1,9,7>;
        J32cdab [0][0][4][3] = jCDAB  <1,1,9,7>;
        X32acbd [0][0][4][3] = xACBD  <1,1,9,7>;
        X32adbc [0][0][4][3] = xADBC  <1,1,9,7>;
        X32bdac [0][0][4][3] = xBDAC  <1,1,9,7>;
        X32bcad [0][0][4][3] = xBCAD  <1,1,9,7>;
        J32abcd [0][0][4][4] = jABCD  <1,1,9,9>;
        J32cdab [0][0][4][4] = jCDAB  <1,1,9,9>;
        X32acbd [0][0][4][4] = xACBD  <1,1,9,9>;
        X32adbc [0][0][4][4] = xADBC  <1,1,9,9>;
        X32bdac [0][0][4][4] = xBDAC  <1,1,9,9>;
        X32bcad [0][0][4][4] = xBCAD  <1,1,9,9>;

        J32abcd [1][0][0][0] = jABCD  <3,1,1,1>;
        J32cdab [1][0][0][0] = jCDAB  <3,1,1,1>;
        X32acbd [1][0][0][0] = xACBD  <3,1,1,1>;
        X32adbc [1][0][0][0] = xADBC  <3,1,1,1>;
        X32bdac [1][0][0][0] = xBDAC  <3,1,1,1>;
        X32bcad [1][0][0][0] = xBCAD  <3,1,1,1>;
        J32abcd [1][0][1][0] = jABCD  <3,1,3,1>;
        J32cdab [1][0][1][0] = jCDAB  <3,1,3,1>;
        X32acbd [1][0][1][0] = xACBD  <3,1,3,1>;
        X32adbc [1][0][1][0] = xADBC  <3,1,3,1>;
        X32bdac [1][0][1][0] = xBDAC  <3,1,3,1>;
        X32bcad [1][0][1][0] = xBCAD  <3,1,3,1>;
        J32abcd [1][0][1][1] = jABCD  <3,1,3,3>;
        J32cdab [1][0][1][1] = jCDAB  <3,1,3,3>;
        X32acbd [1][0][1][1] = xACBD  <3,1,3,3>;
        X32adbc [1][0][1][1] = xADBC  <3,1,3,3>;
        X32bdac [1][0][1][1] = xBDAC  <3,1,3,3>;
        X32bcad [1][0][1][1] = xBCAD  <3,1,3,3>;
        J32abcd [1][0][2][0] = jABCD  <3,1,5,1>;
        J32cdab [1][0][2][0] = jCDAB  <3,1,5,1>;
        X32acbd [1][0][2][0] = xACBD  <3,1,5,1>;
        X32adbc [1][0][2][0] = xADBC  <3,1,5,1>;
        X32bdac [1][0][2][0] = xBDAC  <3,1,5,1>;
        X32bcad [1][0][2][0] = xBCAD  <3,1,5,1>;
        J32abcd [1][0][2][1] = jABCD  <3,1,5,3>;
        J32cdab [1][0][2][1] = jCDAB  <3,1,5,3>;
        X32acbd [1][0][2][1] = xACBD  <3,1,5,3>;
        X32adbc [1][0][2][1] = xADBC  <3,1,5,3>;
        X32bdac [1][0][2][1] = xBDAC  <3,1,5,3>;
        X32bcad [1][0][2][1] = xBCAD  <3,1,5,3>;
        J32abcd [1][0][2][2] = jABCD  <3,1,5,5>;
        J32cdab [1][0][2][2] = jCDAB  <3,1,5,5>;
        X32acbd [1][0][2][2] = xACBD  <3,1,5,5>;
        X32adbc [1][0][2][2] = xADBC  <3,1,5,5>;
        X32bdac [1][0][2][2] = xBDAC  <3,1,5,5>;
        X32bcad [1][0][2][2] = xBCAD  <3,1,5,5>;
        J32abcd [1][0][3][0] = jABCD  <3,1,7,1>;
        J32cdab [1][0][3][0] = jCDAB  <3,1,7,1>;
        X32acbd [1][0][3][0] = xACBD  <3,1,7,1>;
        X32adbc [1][0][3][0] = xADBC  <3,1,7,1>;
        X32bdac [1][0][3][0] = xBDAC  <3,1,7,1>;
        X32bcad [1][0][3][0] = xBCAD  <3,1,7,1>;
        J32abcd [1][0][3][1] = jABCD  <3,1,7,3>;
        J32cdab [1][0][3][1] = jCDAB  <3,1,7,3>;
        X32acbd [1][0][3][1] = xACBD  <3,1,7,3>;
        X32adbc [1][0][3][1] = xADBC  <3,1,7,3>;
        X32bdac [1][0][3][1] = xBDAC  <3,1,7,3>;
        X32bcad [1][0][3][1] = xBCAD  <3,1,7,3>;
        J32abcd [1][0][3][2] = jABCD  <3,1,7,5>;
        J32cdab [1][0][3][2] = jCDAB  <3,1,7,5>;
        X32acbd [1][0][3][2] = xACBD  <3,1,7,5>;
        X32adbc [1][0][3][2] = xADBC  <3,1,7,5>;
        X32bdac [1][0][3][2] = xBDAC  <3,1,7,5>;
        X32bcad [1][0][3][2] = xBCAD  <3,1,7,5>;
        J32abcd [1][0][3][3] = jABCD  <3,1,7,7>;
        J32cdab [1][0][3][3] = jCDAB  <3,1,7,7>;
        X32acbd [1][0][3][3] = xACBD  <3,1,7,7>;
        X32adbc [1][0][3][3] = xADBC  <3,1,7,7>;
        X32bdac [1][0][3][3] = xBDAC  <3,1,7,7>;
        X32bcad [1][0][3][3] = xBCAD  <3,1,7,7>;
        J32abcd [1][0][4][0] = jABCD  <3,1,9,1>;
        J32cdab [1][0][4][0] = jCDAB  <3,1,9,1>;
        X32acbd [1][0][4][0] = xACBD  <3,1,9,1>;
        X32adbc [1][0][4][0] = xADBC  <3,1,9,1>;
        X32bdac [1][0][4][0] = xBDAC  <3,1,9,1>;
        X32bcad [1][0][4][0] = xBCAD  <3,1,9,1>;
        J32abcd [1][0][4][1] = jABCD  <3,1,9,3>;
        J32cdab [1][0][4][1] = jCDAB  <3,1,9,3>;
        X32acbd [1][0][4][1] = xACBD  <3,1,9,3>;
        X32adbc [1][0][4][1] = xADBC  <3,1,9,3>;
        X32bdac [1][0][4][1] = xBDAC  <3,1,9,3>;
        X32bcad [1][0][4][1] = xBCAD  <3,1,9,3>;
        J32abcd [1][0][4][2] = jABCD  <3,1,9,5>;
        J32cdab [1][0][4][2] = jCDAB  <3,1,9,5>;
        X32acbd [1][0][4][2] = xACBD  <3,1,9,5>;
        X32adbc [1][0][4][2] = xADBC  <3,1,9,5>;
        X32bdac [1][0][4][2] = xBDAC  <3,1,9,5>;
        X32bcad [1][0][4][2] = xBCAD  <3,1,9,5>;
        J32abcd [1][0][4][3] = jABCD  <3,1,9,7>;
        J32cdab [1][0][4][3] = jCDAB  <3,1,9,7>;
        X32acbd [1][0][4][3] = xACBD  <3,1,9,7>;
        X32adbc [1][0][4][3] = xADBC  <3,1,9,7>;
        X32bdac [1][0][4][3] = xBDAC  <3,1,9,7>;
        X32bcad [1][0][4][3] = xBCAD  <3,1,9,7>;
        J32abcd [1][0][4][4] = jABCD  <3,1,9,9>;
        J32cdab [1][0][4][4] = jCDAB  <3,1,9,9>;
        X32acbd [1][0][4][4] = xACBD  <3,1,9,9>;
        X32adbc [1][0][4][4] = xADBC  <3,1,9,9>;
        X32bdac [1][0][4][4] = xBDAC  <3,1,9,9>;
        X32bcad [1][0][4][4] = xBCAD  <3,1,9,9>;

        J32abcd [1][1][0][0] = jABCD  <3,3,1,1>;
        J32cdab [1][1][0][0] = jCDAB  <3,3,1,1>;
        X32acbd [1][1][0][0] = xACBD  <3,3,1,1>;
        X32adbc [1][1][0][0] = xADBC  <3,3,1,1>;
        X32bdac [1][1][0][0] = xBDAC  <3,3,1,1>;
        X32bcad [1][1][0][0] = xBCAD  <3,3,1,1>;
        J32abcd [1][1][1][0] = jABCD  <3,3,3,1>;
        J32cdab [1][1][1][0] = jCDAB  <3,3,3,1>;
        X32acbd [1][1][1][0] = xACBD  <3,3,3,1>;
        X32adbc [1][1][1][0] = xADBC  <3,3,3,1>;
        X32bdac [1][1][1][0] = xBDAC  <3,3,3,1>;
        X32bcad [1][1][1][0] = xBCAD  <3,3,3,1>;
        J32abcd [1][1][1][1] = jABCD  <3,3,3,3>;
        J32cdab [1][1][1][1] = jCDAB  <3,3,3,3>;
        X32acbd [1][1][1][1] = xACBD  <3,3,3,3>;
        X32adbc [1][1][1][1] = xADBC  <3,3,3,3>;
        X32bdac [1][1][1][1] = xBDAC  <3,3,3,3>;
        X32bcad [1][1][1][1] = xBCAD  <3,3,3,3>;
        J32abcd [1][1][2][0] = jABCD  <3,3,5,1>;
        J32cdab [1][1][2][0] = jCDAB  <3,3,5,1>;
        X32acbd [1][1][2][0] = xACBD  <3,3,5,1>;
        X32adbc [1][1][2][0] = xADBC  <3,3,5,1>;
        X32bdac [1][1][2][0] = xBDAC  <3,3,5,1>;
        X32bcad [1][1][2][0] = xBCAD  <3,3,5,1>;
        J32abcd [1][1][2][1] = jABCD  <3,3,5,3>;
        J32cdab [1][1][2][1] = jCDAB  <3,3,5,3>;
        X32acbd [1][1][2][1] = xACBD  <3,3,5,3>;
        X32adbc [1][1][2][1] = xADBC  <3,3,5,3>;
        X32bdac [1][1][2][1] = xBDAC  <3,3,5,3>;
        X32bcad [1][1][2][1] = xBCAD  <3,3,5,3>;
        J32abcd [1][1][2][2] = jABCD  <3,3,5,5>;
        J32cdab [1][1][2][2] = jCDAB  <3,3,5,5>;
        X32acbd [1][1][2][2] = xACBD  <3,3,5,5>;
        X32adbc [1][1][2][2] = xADBC  <3,3,5,5>;
        X32bdac [1][1][2][2] = xBDAC  <3,3,5,5>;
        X32bcad [1][1][2][2] = xBCAD  <3,3,5,5>;
        J32abcd [1][1][3][0] = jABCD  <3,3,7,1>;
        J32cdab [1][1][3][0] = jCDAB  <3,3,7,1>;
        X32acbd [1][1][3][0] = xACBD  <3,3,7,1>;
        X32adbc [1][1][3][0] = xADBC  <3,3,7,1>;
        X32bdac [1][1][3][0] = xBDAC  <3,3,7,1>;
        X32bcad [1][1][3][0] = xBCAD  <3,3,7,1>;
        J32abcd [1][1][3][1] = jABCD  <3,3,7,3>;
        J32cdab [1][1][3][1] = jCDAB  <3,3,7,3>;
        X32acbd [1][1][3][1] = xACBD  <3,3,7,3>;
        X32adbc [1][1][3][1] = xADBC  <3,3,7,3>;
        X32bdac [1][1][3][1] = xBDAC  <3,3,7,3>;
        X32bcad [1][1][3][1] = xBCAD  <3,3,7,3>;
        J32abcd [1][1][3][2] = jABCD  <3,3,7,5>;
        J32cdab [1][1][3][2] = jCDAB  <3,3,7,5>;
        X32acbd [1][1][3][2] = xACBD  <3,3,7,5>;
        X32adbc [1][1][3][2] = xADBC  <3,3,7,5>;
        X32bdac [1][1][3][2] = xBDAC  <3,3,7,5>;
        X32bcad [1][1][3][2] = xBCAD  <3,3,7,5>;
        J32abcd [1][1][3][3] = jABCD  <3,3,7,7>;
        J32cdab [1][1][3][3] = jCDAB  <3,3,7,7>;
        X32acbd [1][1][3][3] = xACBD  <3,3,7,7>;
        X32adbc [1][1][3][3] = xADBC  <3,3,7,7>;
        X32bdac [1][1][3][3] = xBDAC  <3,3,7,7>;
        X32bcad [1][1][3][3] = xBCAD  <3,3,7,7>;
        J32abcd [1][1][4][0] = jABCD  <3,3,9,1>;
        J32cdab [1][1][4][0] = jCDAB  <3,3,9,1>;
        X32acbd [1][1][4][0] = xACBD  <3,3,9,1>;
        X32adbc [1][1][4][0] = xADBC  <3,3,9,1>;
        X32bdac [1][1][4][0] = xBDAC  <3,3,9,1>;
        X32bcad [1][1][4][0] = xBCAD  <3,3,9,1>;
        J32abcd [1][1][4][1] = jABCD  <3,3,9,3>;
        J32cdab [1][1][4][1] = jCDAB  <3,3,9,3>;
        X32acbd [1][1][4][1] = xACBD  <3,3,9,3>;
        X32adbc [1][1][4][1] = xADBC  <3,3,9,3>;
        X32bdac [1][1][4][1] = xBDAC  <3,3,9,3>;
        X32bcad [1][1][4][1] = xBCAD  <3,3,9,3>;
        J32abcd [1][1][4][2] = jABCD  <3,3,9,5>;
        J32cdab [1][1][4][2] = jCDAB  <3,3,9,5>;
        X32acbd [1][1][4][2] = xACBD  <3,3,9,5>;
        X32adbc [1][1][4][2] = xADBC  <3,3,9,5>;
        X32bdac [1][1][4][2] = xBDAC  <3,3,9,5>;
        X32bcad [1][1][4][2] = xBCAD  <3,3,9,5>;
        J32abcd [1][1][4][3] = jABCD  <3,3,9,7>;
        J32cdab [1][1][4][3] = jCDAB  <3,3,9,7>;
        X32acbd [1][1][4][3] = xACBD  <3,3,9,7>;
        X32adbc [1][1][4][3] = xADBC  <3,3,9,7>;
        X32bdac [1][1][4][3] = xBDAC  <3,3,9,7>;
        X32bcad [1][1][4][3] = xBCAD  <3,3,9,7>;
        J32abcd [1][1][4][4] = jABCD  <3,3,9,9>;
        J32cdab [1][1][4][4] = jCDAB  <3,3,9,9>;
        X32acbd [1][1][4][4] = xACBD  <3,3,9,9>;
        X32adbc [1][1][4][4] = xADBC  <3,3,9,9>;
        X32bdac [1][1][4][4] = xBDAC  <3,3,9,9>;
        X32bcad [1][1][4][4] = xBCAD  <3,3,9,9>;

        J32abcd [2][0][0][0] = jABCD  <5,1,1,1>;
        J32cdab [2][0][0][0] = jCDAB  <5,1,1,1>;
        X32acbd [2][0][0][0] = xACBD  <5,1,1,1>;
        X32adbc [2][0][0][0] = xADBC  <5,1,1,1>;
        X32bdac [2][0][0][0] = xBDAC  <5,1,1,1>;
        X32bcad [2][0][0][0] = xBCAD  <5,1,1,1>;
        J32abcd [2][0][1][0] = jABCD  <5,1,3,1>;
        J32cdab [2][0][1][0] = jCDAB  <5,1,3,1>;
        X32acbd [2][0][1][0] = xACBD  <5,1,3,1>;
        X32adbc [2][0][1][0] = xADBC  <5,1,3,1>;
        X32bdac [2][0][1][0] = xBDAC  <5,1,3,1>;
        X32bcad [2][0][1][0] = xBCAD  <5,1,3,1>;
        J32abcd [2][0][1][1] = jABCD  <5,1,3,3>;
        J32cdab [2][0][1][1] = jCDAB  <5,1,3,3>;
        X32acbd [2][0][1][1] = xACBD  <5,1,3,3>;
        X32adbc [2][0][1][1] = xADBC  <5,1,3,3>;
        X32bdac [2][0][1][1] = xBDAC  <5,1,3,3>;
        X32bcad [2][0][1][1] = xBCAD  <5,1,3,3>;
        J32abcd [2][0][2][0] = jABCD  <5,1,5,1>;
        J32cdab [2][0][2][0] = jCDAB  <5,1,5,1>;
        X32acbd [2][0][2][0] = xACBD  <5,1,5,1>;
        X32adbc [2][0][2][0] = xADBC  <5,1,5,1>;
        X32bdac [2][0][2][0] = xBDAC  <5,1,5,1>;
        X32bcad [2][0][2][0] = xBCAD  <5,1,5,1>;
        J32abcd [2][0][2][1] = jABCD  <5,1,5,3>;
        J32cdab [2][0][2][1] = jCDAB  <5,1,5,3>;
        X32acbd [2][0][2][1] = xACBD  <5,1,5,3>;
        X32adbc [2][0][2][1] = xADBC  <5,1,5,3>;
        X32bdac [2][0][2][1] = xBDAC  <5,1,5,3>;
        X32bcad [2][0][2][1] = xBCAD  <5,1,5,3>;
        J32abcd [2][0][2][2] = jABCD  <5,1,5,5>;
        J32cdab [2][0][2][2] = jCDAB  <5,1,5,5>;
        X32acbd [2][0][2][2] = xACBD  <5,1,5,5>;
        X32adbc [2][0][2][2] = xADBC  <5,1,5,5>;
        X32bdac [2][0][2][2] = xBDAC  <5,1,5,5>;
        X32bcad [2][0][2][2] = xBCAD  <5,1,5,5>;
        J32abcd [2][0][3][0] = jABCD  <5,1,7,1>;
        J32cdab [2][0][3][0] = jCDAB  <5,1,7,1>;
        X32acbd [2][0][3][0] = xACBD  <5,1,7,1>;
        X32adbc [2][0][3][0] = xADBC  <5,1,7,1>;
        X32bdac [2][0][3][0] = xBDAC  <5,1,7,1>;
        X32bcad [2][0][3][0] = xBCAD  <5,1,7,1>;
        J32abcd [2][0][3][1] = jABCD  <5,1,7,3>;
        J32cdab [2][0][3][1] = jCDAB  <5,1,7,3>;
        X32acbd [2][0][3][1] = xACBD  <5,1,7,3>;
        X32adbc [2][0][3][1] = xADBC  <5,1,7,3>;
        X32bdac [2][0][3][1] = xBDAC  <5,1,7,3>;
        X32bcad [2][0][3][1] = xBCAD  <5,1,7,3>;
        J32abcd [2][0][3][2] = jABCD  <5,1,7,5>;
        J32cdab [2][0][3][2] = jCDAB  <5,1,7,5>;
        X32acbd [2][0][3][2] = xACBD  <5,1,7,5>;
        X32adbc [2][0][3][2] = xADBC  <5,1,7,5>;
        X32bdac [2][0][3][2] = xBDAC  <5,1,7,5>;
        X32bcad [2][0][3][2] = xBCAD  <5,1,7,5>;
        J32abcd [2][0][3][3] = jABCD  <5,1,7,7>;
        J32cdab [2][0][3][3] = jCDAB  <5,1,7,7>;
        X32acbd [2][0][3][3] = xACBD  <5,1,7,7>;
        X32adbc [2][0][3][3] = xADBC  <5,1,7,7>;
        X32bdac [2][0][3][3] = xBDAC  <5,1,7,7>;
        X32bcad [2][0][3][3] = xBCAD  <5,1,7,7>;
        J32abcd [2][0][4][0] = jABCD  <5,1,9,1>;
        J32cdab [2][0][4][0] = jCDAB  <5,1,9,1>;
        X32acbd [2][0][4][0] = xACBD  <5,1,9,1>;
        X32adbc [2][0][4][0] = xADBC  <5,1,9,1>;
        X32bdac [2][0][4][0] = xBDAC  <5,1,9,1>;
        X32bcad [2][0][4][0] = xBCAD  <5,1,9,1>;
        J32abcd [2][0][4][1] = jABCD  <5,1,9,3>;
        J32cdab [2][0][4][1] = jCDAB  <5,1,9,3>;
        X32acbd [2][0][4][1] = xACBD  <5,1,9,3>;
        X32adbc [2][0][4][1] = xADBC  <5,1,9,3>;
        X32bdac [2][0][4][1] = xBDAC  <5,1,9,3>;
        X32bcad [2][0][4][1] = xBCAD  <5,1,9,3>;
        J32abcd [2][0][4][2] = jABCD  <5,1,9,5>;
        J32cdab [2][0][4][2] = jCDAB  <5,1,9,5>;
        X32acbd [2][0][4][2] = xACBD  <5,1,9,5>;
        X32adbc [2][0][4][2] = xADBC  <5,1,9,5>;
        X32bdac [2][0][4][2] = xBDAC  <5,1,9,5>;
        X32bcad [2][0][4][2] = xBCAD  <5,1,9,5>;
        J32abcd [2][0][4][3] = jABCD  <5,1,9,7>;
        J32cdab [2][0][4][3] = jCDAB  <5,1,9,7>;
        X32acbd [2][0][4][3] = xACBD  <5,1,9,7>;
        X32adbc [2][0][4][3] = xADBC  <5,1,9,7>;
        X32bdac [2][0][4][3] = xBDAC  <5,1,9,7>;
        X32bcad [2][0][4][3] = xBCAD  <5,1,9,7>;
        J32abcd [2][0][4][4] = jABCD  <5,1,9,9>;
        J32cdab [2][0][4][4] = jCDAB  <5,1,9,9>;
        X32acbd [2][0][4][4] = xACBD  <5,1,9,9>;
        X32adbc [2][0][4][4] = xADBC  <5,1,9,9>;
        X32bdac [2][0][4][4] = xBDAC  <5,1,9,9>;
        X32bcad [2][0][4][4] = xBCAD  <5,1,9,9>;

        J32abcd [2][1][0][0] = jABCD  <5,3,1,1>;
        J32cdab [2][1][0][0] = jCDAB  <5,3,1,1>;
        X32acbd [2][1][0][0] = xACBD  <5,3,1,1>;
        X32adbc [2][1][0][0] = xADBC  <5,3,1,1>;
        X32bdac [2][1][0][0] = xBDAC  <5,3,1,1>;
        X32bcad [2][1][0][0] = xBCAD  <5,3,1,1>;
        J32abcd [2][1][1][0] = jABCD  <5,3,3,1>;
        J32cdab [2][1][1][0] = jCDAB  <5,3,3,1>;
        X32acbd [2][1][1][0] = xACBD  <5,3,3,1>;
        X32adbc [2][1][1][0] = xADBC  <5,3,3,1>;
        X32bdac [2][1][1][0] = xBDAC  <5,3,3,1>;
        X32bcad [2][1][1][0] = xBCAD  <5,3,3,1>;
        J32abcd [2][1][1][1] = jABCD  <5,3,3,3>;
        J32cdab [2][1][1][1] = jCDAB  <5,3,3,3>;
        X32acbd [2][1][1][1] = xACBD  <5,3,3,3>;
        X32adbc [2][1][1][1] = xADBC  <5,3,3,3>;
        X32bdac [2][1][1][1] = xBDAC  <5,3,3,3>;
        X32bcad [2][1][1][1] = xBCAD  <5,3,3,3>;
        J32abcd [2][1][2][0] = jABCD  <5,3,5,1>;
        J32cdab [2][1][2][0] = jCDAB  <5,3,5,1>;
        X32acbd [2][1][2][0] = xACBD  <5,3,5,1>;
        X32adbc [2][1][2][0] = xADBC  <5,3,5,1>;
        X32bdac [2][1][2][0] = xBDAC  <5,3,5,1>;
        X32bcad [2][1][2][0] = xBCAD  <5,3,5,1>;
        J32abcd [2][1][2][1] = jABCD  <5,3,5,3>;
        J32cdab [2][1][2][1] = jCDAB  <5,3,5,3>;
        X32acbd [2][1][2][1] = xACBD  <5,3,5,3>;
        X32adbc [2][1][2][1] = xADBC  <5,3,5,3>;
        X32bdac [2][1][2][1] = xBDAC  <5,3,5,3>;
        X32bcad [2][1][2][1] = xBCAD  <5,3,5,3>;
        J32abcd [2][1][2][2] = jABCD  <5,3,5,5>;
        J32cdab [2][1][2][2] = jCDAB  <5,3,5,5>;
        X32acbd [2][1][2][2] = xACBD  <5,3,5,5>;
        X32adbc [2][1][2][2] = xADBC  <5,3,5,5>;
        X32bdac [2][1][2][2] = xBDAC  <5,3,5,5>;
        X32bcad [2][1][2][2] = xBCAD  <5,3,5,5>;
        J32abcd [2][1][3][0] = jABCD  <5,3,7,1>;
        J32cdab [2][1][3][0] = jCDAB  <5,3,7,1>;
        X32acbd [2][1][3][0] = xACBD  <5,3,7,1>;
        X32adbc [2][1][3][0] = xADBC  <5,3,7,1>;
        X32bdac [2][1][3][0] = xBDAC  <5,3,7,1>;
        X32bcad [2][1][3][0] = xBCAD  <5,3,7,1>;
        J32abcd [2][1][3][1] = jABCD  <5,3,7,3>;
        J32cdab [2][1][3][1] = jCDAB  <5,3,7,3>;
        X32acbd [2][1][3][1] = xACBD  <5,3,7,3>;
        X32adbc [2][1][3][1] = xADBC  <5,3,7,3>;
        X32bdac [2][1][3][1] = xBDAC  <5,3,7,3>;
        X32bcad [2][1][3][1] = xBCAD  <5,3,7,3>;
        J32abcd [2][1][3][2] = jABCD  <5,3,7,5>;
        J32cdab [2][1][3][2] = jCDAB  <5,3,7,5>;
        X32acbd [2][1][3][2] = xACBD  <5,3,7,5>;
        X32adbc [2][1][3][2] = xADBC  <5,3,7,5>;
        X32bdac [2][1][3][2] = xBDAC  <5,3,7,5>;
        X32bcad [2][1][3][2] = xBCAD  <5,3,7,5>;
        J32abcd [2][1][3][3] = jABCD  <5,3,7,7>;
        J32cdab [2][1][3][3] = jCDAB  <5,3,7,7>;
        X32acbd [2][1][3][3] = xACBD  <5,3,7,7>;
        X32adbc [2][1][3][3] = xADBC  <5,3,7,7>;
        X32bdac [2][1][3][3] = xBDAC  <5,3,7,7>;
        X32bcad [2][1][3][3] = xBCAD  <5,3,7,7>;
        J32abcd [2][1][4][0] = jABCD  <5,3,9,1>;
        J32cdab [2][1][4][0] = jCDAB  <5,3,9,1>;
        X32acbd [2][1][4][0] = xACBD  <5,3,9,1>;
        X32adbc [2][1][4][0] = xADBC  <5,3,9,1>;
        X32bdac [2][1][4][0] = xBDAC  <5,3,9,1>;
        X32bcad [2][1][4][0] = xBCAD  <5,3,9,1>;
        J32abcd [2][1][4][1] = jABCD  <5,3,9,3>;
        J32cdab [2][1][4][1] = jCDAB  <5,3,9,3>;
        X32acbd [2][1][4][1] = xACBD  <5,3,9,3>;
        X32adbc [2][1][4][1] = xADBC  <5,3,9,3>;
        X32bdac [2][1][4][1] = xBDAC  <5,3,9,3>;
        X32bcad [2][1][4][1] = xBCAD  <5,3,9,3>;
        J32abcd [2][1][4][2] = jABCD  <5,3,9,5>;
        J32cdab [2][1][4][2] = jCDAB  <5,3,9,5>;
        X32acbd [2][1][4][2] = xACBD  <5,3,9,5>;
        X32adbc [2][1][4][2] = xADBC  <5,3,9,5>;
        X32bdac [2][1][4][2] = xBDAC  <5,3,9,5>;
        X32bcad [2][1][4][2] = xBCAD  <5,3,9,5>;
        J32abcd [2][1][4][3] = jABCD  <5,3,9,7>;
        J32cdab [2][1][4][3] = jCDAB  <5,3,9,7>;
        X32acbd [2][1][4][3] = xACBD  <5,3,9,7>;
        X32adbc [2][1][4][3] = xADBC  <5,3,9,7>;
        X32bdac [2][1][4][3] = xBDAC  <5,3,9,7>;
        X32bcad [2][1][4][3] = xBCAD  <5,3,9,7>;
        J32abcd [2][1][4][4] = jABCD  <5,3,9,9>;
        J32cdab [2][1][4][4] = jCDAB  <5,3,9,9>;
        X32acbd [2][1][4][4] = xACBD  <5,3,9,9>;
        X32adbc [2][1][4][4] = xADBC  <5,3,9,9>;
        X32bdac [2][1][4][4] = xBDAC  <5,3,9,9>;
        X32bcad [2][1][4][4] = xBCAD  <5,3,9,9>;

        J32abcd [2][2][0][0] = jABCD  <5,5,1,1>;
        J32cdab [2][2][0][0] = jCDAB  <5,5,1,1>;
        X32acbd [2][2][0][0] = xACBD  <5,5,1,1>;
        X32adbc [2][2][0][0] = xADBC  <5,5,1,1>;
        X32bdac [2][2][0][0] = xBDAC  <5,5,1,1>;
        X32bcad [2][2][0][0] = xBCAD  <5,5,1,1>;
        J32abcd [2][2][1][0] = jABCD  <5,5,3,1>;
        J32cdab [2][2][1][0] = jCDAB  <5,5,3,1>;
        X32acbd [2][2][1][0] = xACBD  <5,5,3,1>;
        X32adbc [2][2][1][0] = xADBC  <5,5,3,1>;
        X32bdac [2][2][1][0] = xBDAC  <5,5,3,1>;
        X32bcad [2][2][1][0] = xBCAD  <5,5,3,1>;
        J32abcd [2][2][1][1] = jABCD  <5,5,3,3>;
        J32cdab [2][2][1][1] = jCDAB  <5,5,3,3>;
        X32acbd [2][2][1][1] = xACBD  <5,5,3,3>;
        X32adbc [2][2][1][1] = xADBC  <5,5,3,3>;
        X32bdac [2][2][1][1] = xBDAC  <5,5,3,3>;
        X32bcad [2][2][1][1] = xBCAD  <5,5,3,3>;
        J32abcd [2][2][2][0] = jABCD  <5,5,5,1>;
        J32cdab [2][2][2][0] = jCDAB  <5,5,5,1>;
        X32acbd [2][2][2][0] = xACBD  <5,5,5,1>;
        X32adbc [2][2][2][0] = xADBC  <5,5,5,1>;
        X32bdac [2][2][2][0] = xBDAC  <5,5,5,1>;
        X32bcad [2][2][2][0] = xBCAD  <5,5,5,1>;
        J32abcd [2][2][2][1] = jABCD  <5,5,5,3>;
        J32cdab [2][2][2][1] = jCDAB  <5,5,5,3>;
        X32acbd [2][2][2][1] = xACBD  <5,5,5,3>;
        X32adbc [2][2][2][1] = xADBC  <5,5,5,3>;
        X32bdac [2][2][2][1] = xBDAC  <5,5,5,3>;
        X32bcad [2][2][2][1] = xBCAD  <5,5,5,3>;
        J32abcd [2][2][2][2] = jABCD  <5,5,5,5>;
        J32cdab [2][2][2][2] = jCDAB  <5,5,5,5>;
        X32acbd [2][2][2][2] = xACBD  <5,5,5,5>;
        X32adbc [2][2][2][2] = xADBC  <5,5,5,5>;
        X32bdac [2][2][2][2] = xBDAC  <5,5,5,5>;
        X32bcad [2][2][2][2] = xBCAD  <5,5,5,5>;
        J32abcd [2][2][3][0] = jABCD  <5,5,7,1>;
        J32cdab [2][2][3][0] = jCDAB  <5,5,7,1>;
        X32acbd [2][2][3][0] = xACBD  <5,5,7,1>;
        X32adbc [2][2][3][0] = xADBC  <5,5,7,1>;
        X32bdac [2][2][3][0] = xBDAC  <5,5,7,1>;
        X32bcad [2][2][3][0] = xBCAD  <5,5,7,1>;
        J32abcd [2][2][3][1] = jABCD  <5,5,7,3>;
        J32cdab [2][2][3][1] = jCDAB  <5,5,7,3>;
        X32acbd [2][2][3][1] = xACBD  <5,5,7,3>;
        X32adbc [2][2][3][1] = xADBC  <5,5,7,3>;
        X32bdac [2][2][3][1] = xBDAC  <5,5,7,3>;
        X32bcad [2][2][3][1] = xBCAD  <5,5,7,3>;
        J32abcd [2][2][3][2] = jABCD  <5,5,7,5>;
        J32cdab [2][2][3][2] = jCDAB  <5,5,7,5>;
        X32acbd [2][2][3][2] = xACBD  <5,5,7,5>;
        X32adbc [2][2][3][2] = xADBC  <5,5,7,5>;
        X32bdac [2][2][3][2] = xBDAC  <5,5,7,5>;
        X32bcad [2][2][3][2] = xBCAD  <5,5,7,5>;
        J32abcd [2][2][3][3] = jABCD  <5,5,7,7>;
        J32cdab [2][2][3][3] = jCDAB  <5,5,7,7>;
        X32acbd [2][2][3][3] = xACBD  <5,5,7,7>;
        X32adbc [2][2][3][3] = xADBC  <5,5,7,7>;
        X32bdac [2][2][3][3] = xBDAC  <5,5,7,7>;
        X32bcad [2][2][3][3] = xBCAD  <5,5,7,7>;
        J32abcd [2][2][4][0] = jABCD  <5,5,9,1>;
        J32cdab [2][2][4][0] = jCDAB  <5,5,9,1>;
        X32acbd [2][2][4][0] = xACBD  <5,5,9,1>;
        X32adbc [2][2][4][0] = xADBC  <5,5,9,1>;
        X32bdac [2][2][4][0] = xBDAC  <5,5,9,1>;
        X32bcad [2][2][4][0] = xBCAD  <5,5,9,1>;
        J32abcd [2][2][4][1] = jABCD  <5,5,9,3>;
        J32cdab [2][2][4][1] = jCDAB  <5,5,9,3>;
        X32acbd [2][2][4][1] = xACBD  <5,5,9,3>;
        X32adbc [2][2][4][1] = xADBC  <5,5,9,3>;
        X32bdac [2][2][4][1] = xBDAC  <5,5,9,3>;
        X32bcad [2][2][4][1] = xBCAD  <5,5,9,3>;
        J32abcd [2][2][4][2] = jABCD  <5,5,9,5>;
        J32cdab [2][2][4][2] = jCDAB  <5,5,9,5>;
        X32acbd [2][2][4][2] = xACBD  <5,5,9,5>;
        X32adbc [2][2][4][2] = xADBC  <5,5,9,5>;
        X32bdac [2][2][4][2] = xBDAC  <5,5,9,5>;
        X32bcad [2][2][4][2] = xBCAD  <5,5,9,5>;
        J32abcd [2][2][4][3] = jABCD  <5,5,9,7>;
        J32cdab [2][2][4][3] = jCDAB  <5,5,9,7>;
        X32acbd [2][2][4][3] = xACBD  <5,5,9,7>;
        X32adbc [2][2][4][3] = xADBC  <5,5,9,7>;
        X32bdac [2][2][4][3] = xBDAC  <5,5,9,7>;
        X32bcad [2][2][4][3] = xBCAD  <5,5,9,7>;
        J32abcd [2][2][4][4] = jABCD  <5,5,9,9>;
        J32cdab [2][2][4][4] = jCDAB  <5,5,9,9>;
        X32acbd [2][2][4][4] = xACBD  <5,5,9,9>;
        X32adbc [2][2][4][4] = xADBC  <5,5,9,9>;
        X32bdac [2][2][4][4] = xBDAC  <5,5,9,9>;
        X32bcad [2][2][4][4] = xBCAD  <5,5,9,9>;

        J32abcd [3][0][0][0] = jABCD  <7,1,1,1>;
        J32cdab [3][0][0][0] = jCDAB  <7,1,1,1>;
        X32acbd [3][0][0][0] = xACBD  <7,1,1,1>;
        X32adbc [3][0][0][0] = xADBC  <7,1,1,1>;
        X32bdac [3][0][0][0] = xBDAC  <7,1,1,1>;
        X32bcad [3][0][0][0] = xBCAD  <7,1,1,1>;
        J32abcd [3][0][1][0] = jABCD  <7,1,3,1>;
        J32cdab [3][0][1][0] = jCDAB  <7,1,3,1>;
        X32acbd [3][0][1][0] = xACBD  <7,1,3,1>;
        X32adbc [3][0][1][0] = xADBC  <7,1,3,1>;
        X32bdac [3][0][1][0] = xBDAC  <7,1,3,1>;
        X32bcad [3][0][1][0] = xBCAD  <7,1,3,1>;
        J32abcd [3][0][1][1] = jABCD  <7,1,3,3>;
        J32cdab [3][0][1][1] = jCDAB  <7,1,3,3>;
        X32acbd [3][0][1][1] = xACBD  <7,1,3,3>;
        X32adbc [3][0][1][1] = xADBC  <7,1,3,3>;
        X32bdac [3][0][1][1] = xBDAC  <7,1,3,3>;
        X32bcad [3][0][1][1] = xBCAD  <7,1,3,3>;
        J32abcd [3][0][2][0] = jABCD  <7,1,5,1>;
        J32cdab [3][0][2][0] = jCDAB  <7,1,5,1>;
        X32acbd [3][0][2][0] = xACBD  <7,1,5,1>;
        X32adbc [3][0][2][0] = xADBC  <7,1,5,1>;
        X32bdac [3][0][2][0] = xBDAC  <7,1,5,1>;
        X32bcad [3][0][2][0] = xBCAD  <7,1,5,1>;
        J32abcd [3][0][2][1] = jABCD  <7,1,5,3>;
        J32cdab [3][0][2][1] = jCDAB  <7,1,5,3>;
        X32acbd [3][0][2][1] = xACBD  <7,1,5,3>;
        X32adbc [3][0][2][1] = xADBC  <7,1,5,3>;
        X32bdac [3][0][2][1] = xBDAC  <7,1,5,3>;
        X32bcad [3][0][2][1] = xBCAD  <7,1,5,3>;
        J32abcd [3][0][2][2] = jABCD  <7,1,5,5>;
        J32cdab [3][0][2][2] = jCDAB  <7,1,5,5>;
        X32acbd [3][0][2][2] = xACBD  <7,1,5,5>;
        X32adbc [3][0][2][2] = xADBC  <7,1,5,5>;
        X32bdac [3][0][2][2] = xBDAC  <7,1,5,5>;
        X32bcad [3][0][2][2] = xBCAD  <7,1,5,5>;
        J32abcd [3][0][3][0] = jABCD  <7,1,7,1>;
        J32cdab [3][0][3][0] = jCDAB  <7,1,7,1>;
        X32acbd [3][0][3][0] = xACBD  <7,1,7,1>;
        X32adbc [3][0][3][0] = xADBC  <7,1,7,1>;
        X32bdac [3][0][3][0] = xBDAC  <7,1,7,1>;
        X32bcad [3][0][3][0] = xBCAD  <7,1,7,1>;
        J32abcd [3][0][3][1] = jABCD  <7,1,7,3>;
        J32cdab [3][0][3][1] = jCDAB  <7,1,7,3>;
        X32acbd [3][0][3][1] = xACBD  <7,1,7,3>;
        X32adbc [3][0][3][1] = xADBC  <7,1,7,3>;
        X32bdac [3][0][3][1] = xBDAC  <7,1,7,3>;
        X32bcad [3][0][3][1] = xBCAD  <7,1,7,3>;
        J32abcd [3][0][3][2] = jABCD  <7,1,7,5>;
        J32cdab [3][0][3][2] = jCDAB  <7,1,7,5>;
        X32acbd [3][0][3][2] = xACBD  <7,1,7,5>;
        X32adbc [3][0][3][2] = xADBC  <7,1,7,5>;
        X32bdac [3][0][3][2] = xBDAC  <7,1,7,5>;
        X32bcad [3][0][3][2] = xBCAD  <7,1,7,5>;
        J32abcd [3][0][3][3] = jABCD  <7,1,7,7>;
        J32cdab [3][0][3][3] = jCDAB  <7,1,7,7>;
        X32acbd [3][0][3][3] = xACBD  <7,1,7,7>;
        X32adbc [3][0][3][3] = xADBC  <7,1,7,7>;
        X32bdac [3][0][3][3] = xBDAC  <7,1,7,7>;
        X32bcad [3][0][3][3] = xBCAD  <7,1,7,7>;
        J32abcd [3][0][4][0] = jABCD  <7,1,9,1>;
        J32cdab [3][0][4][0] = jCDAB  <7,1,9,1>;
        X32acbd [3][0][4][0] = xACBD  <7,1,9,1>;
        X32adbc [3][0][4][0] = xADBC  <7,1,9,1>;
        X32bdac [3][0][4][0] = xBDAC  <7,1,9,1>;
        X32bcad [3][0][4][0] = xBCAD  <7,1,9,1>;
        J32abcd [3][0][4][1] = jABCD  <7,1,9,3>;
        J32cdab [3][0][4][1] = jCDAB  <7,1,9,3>;
        X32acbd [3][0][4][1] = xACBD  <7,1,9,3>;
        X32adbc [3][0][4][1] = xADBC  <7,1,9,3>;
        X32bdac [3][0][4][1] = xBDAC  <7,1,9,3>;
        X32bcad [3][0][4][1] = xBCAD  <7,1,9,3>;
        J32abcd [3][0][4][2] = jABCD  <7,1,9,5>;
        J32cdab [3][0][4][2] = jCDAB  <7,1,9,5>;
        X32acbd [3][0][4][2] = xACBD  <7,1,9,5>;
        X32adbc [3][0][4][2] = xADBC  <7,1,9,5>;
        X32bdac [3][0][4][2] = xBDAC  <7,1,9,5>;
        X32bcad [3][0][4][2] = xBCAD  <7,1,9,5>;
        J32abcd [3][0][4][3] = jABCD  <7,1,9,7>;
        J32cdab [3][0][4][3] = jCDAB  <7,1,9,7>;
        X32acbd [3][0][4][3] = xACBD  <7,1,9,7>;
        X32adbc [3][0][4][3] = xADBC  <7,1,9,7>;
        X32bdac [3][0][4][3] = xBDAC  <7,1,9,7>;
        X32bcad [3][0][4][3] = xBCAD  <7,1,9,7>;
        J32abcd [3][0][4][4] = jABCD  <7,1,9,9>;
        J32cdab [3][0][4][4] = jCDAB  <7,1,9,9>;
        X32acbd [3][0][4][4] = xACBD  <7,1,9,9>;
        X32adbc [3][0][4][4] = xADBC  <7,1,9,9>;
        X32bdac [3][0][4][4] = xBDAC  <7,1,9,9>;
        X32bcad [3][0][4][4] = xBCAD  <7,1,9,9>;

        J32abcd [3][1][0][0] = jABCD  <7,3,1,1>;
        J32cdab [3][1][0][0] = jCDAB  <7,3,1,1>;
        X32acbd [3][1][0][0] = xACBD  <7,3,1,1>;
        X32adbc [3][1][0][0] = xADBC  <7,3,1,1>;
        X32bdac [3][1][0][0] = xBDAC  <7,3,1,1>;
        X32bcad [3][1][0][0] = xBCAD  <7,3,1,1>;
        J32abcd [3][1][1][0] = jABCD  <7,3,3,1>;
        J32cdab [3][1][1][0] = jCDAB  <7,3,3,1>;
        X32acbd [3][1][1][0] = xACBD  <7,3,3,1>;
        X32adbc [3][1][1][0] = xADBC  <7,3,3,1>;
        X32bdac [3][1][1][0] = xBDAC  <7,3,3,1>;
        X32bcad [3][1][1][0] = xBCAD  <7,3,3,1>;
        J32abcd [3][1][1][1] = jABCD  <7,3,3,3>;
        J32cdab [3][1][1][1] = jCDAB  <7,3,3,3>;
        X32acbd [3][1][1][1] = xACBD  <7,3,3,3>;
        X32adbc [3][1][1][1] = xADBC  <7,3,3,3>;
        X32bdac [3][1][1][1] = xBDAC  <7,3,3,3>;
        X32bcad [3][1][1][1] = xBCAD  <7,3,3,3>;
        J32abcd [3][1][2][0] = jABCD  <7,3,5,1>;
        J32cdab [3][1][2][0] = jCDAB  <7,3,5,1>;
        X32acbd [3][1][2][0] = xACBD  <7,3,5,1>;
        X32adbc [3][1][2][0] = xADBC  <7,3,5,1>;
        X32bdac [3][1][2][0] = xBDAC  <7,3,5,1>;
        X32bcad [3][1][2][0] = xBCAD  <7,3,5,1>;
        J32abcd [3][1][2][1] = jABCD  <7,3,5,3>;
        J32cdab [3][1][2][1] = jCDAB  <7,3,5,3>;
        X32acbd [3][1][2][1] = xACBD  <7,3,5,3>;
        X32adbc [3][1][2][1] = xADBC  <7,3,5,3>;
        X32bdac [3][1][2][1] = xBDAC  <7,3,5,3>;
        X32bcad [3][1][2][1] = xBCAD  <7,3,5,3>;
        J32abcd [3][1][2][2] = jABCD  <7,3,5,5>;
        J32cdab [3][1][2][2] = jCDAB  <7,3,5,5>;
        X32acbd [3][1][2][2] = xACBD  <7,3,5,5>;
        X32adbc [3][1][2][2] = xADBC  <7,3,5,5>;
        X32bdac [3][1][2][2] = xBDAC  <7,3,5,5>;
        X32bcad [3][1][2][2] = xBCAD  <7,3,5,5>;
        J32abcd [3][1][3][0] = jABCD  <7,3,7,1>;
        J32cdab [3][1][3][0] = jCDAB  <7,3,7,1>;
        X32acbd [3][1][3][0] = xACBD  <7,3,7,1>;
        X32adbc [3][1][3][0] = xADBC  <7,3,7,1>;
        X32bdac [3][1][3][0] = xBDAC  <7,3,7,1>;
        X32bcad [3][1][3][0] = xBCAD  <7,3,7,1>;
        J32abcd [3][1][3][1] = jABCD  <7,3,7,3>;
        J32cdab [3][1][3][1] = jCDAB  <7,3,7,3>;
        X32acbd [3][1][3][1] = xACBD  <7,3,7,3>;
        X32adbc [3][1][3][1] = xADBC  <7,3,7,3>;
        X32bdac [3][1][3][1] = xBDAC  <7,3,7,3>;
        X32bcad [3][1][3][1] = xBCAD  <7,3,7,3>;
        J32abcd [3][1][3][2] = jABCD  <7,3,7,5>;
        J32cdab [3][1][3][2] = jCDAB  <7,3,7,5>;
        X32acbd [3][1][3][2] = xACBD  <7,3,7,5>;
        X32adbc [3][1][3][2] = xADBC  <7,3,7,5>;
        X32bdac [3][1][3][2] = xBDAC  <7,3,7,5>;
        X32bcad [3][1][3][2] = xBCAD  <7,3,7,5>;
        J32abcd [3][1][3][3] = jABCD  <7,3,7,7>;
        J32cdab [3][1][3][3] = jCDAB  <7,3,7,7>;
        X32acbd [3][1][3][3] = xACBD  <7,3,7,7>;
        X32adbc [3][1][3][3] = xADBC  <7,3,7,7>;
        X32bdac [3][1][3][3] = xBDAC  <7,3,7,7>;
        X32bcad [3][1][3][3] = xBCAD  <7,3,7,7>;
        J32abcd [3][1][4][0] = jABCD  <7,3,9,1>;
        J32cdab [3][1][4][0] = jCDAB  <7,3,9,1>;
        X32acbd [3][1][4][0] = xACBD  <7,3,9,1>;
        X32adbc [3][1][4][0] = xADBC  <7,3,9,1>;
        X32bdac [3][1][4][0] = xBDAC  <7,3,9,1>;
        X32bcad [3][1][4][0] = xBCAD  <7,3,9,1>;
        J32abcd [3][1][4][1] = jABCD  <7,3,9,3>;
        J32cdab [3][1][4][1] = jCDAB  <7,3,9,3>;
        X32acbd [3][1][4][1] = xACBD  <7,3,9,3>;
        X32adbc [3][1][4][1] = xADBC  <7,3,9,3>;
        X32bdac [3][1][4][1] = xBDAC  <7,3,9,3>;
        X32bcad [3][1][4][1] = xBCAD  <7,3,9,3>;
        J32abcd [3][1][4][2] = jABCD  <7,3,9,5>;
        J32cdab [3][1][4][2] = jCDAB  <7,3,9,5>;
        X32acbd [3][1][4][2] = xACBD  <7,3,9,5>;
        X32adbc [3][1][4][2] = xADBC  <7,3,9,5>;
        X32bdac [3][1][4][2] = xBDAC  <7,3,9,5>;
        X32bcad [3][1][4][2] = xBCAD  <7,3,9,5>;
        J32abcd [3][1][4][3] = jABCD  <7,3,9,7>;
        J32cdab [3][1][4][3] = jCDAB  <7,3,9,7>;
        X32acbd [3][1][4][3] = xACBD  <7,3,9,7>;
        X32adbc [3][1][4][3] = xADBC  <7,3,9,7>;
        X32bdac [3][1][4][3] = xBDAC  <7,3,9,7>;
        X32bcad [3][1][4][3] = xBCAD  <7,3,9,7>;
        J32abcd [3][1][4][4] = jABCD  <7,3,9,9>;
        J32cdab [3][1][4][4] = jCDAB  <7,3,9,9>;
        X32acbd [3][1][4][4] = xACBD  <7,3,9,9>;
        X32adbc [3][1][4][4] = xADBC  <7,3,9,9>;
        X32bdac [3][1][4][4] = xBDAC  <7,3,9,9>;
        X32bcad [3][1][4][4] = xBCAD  <7,3,9,9>;

        J32abcd [3][2][0][0] = jABCD  <7,5,1,1>;
        J32cdab [3][2][0][0] = jCDAB  <7,5,1,1>;
        X32acbd [3][2][0][0] = xACBD  <7,5,1,1>;
        X32adbc [3][2][0][0] = xADBC  <7,5,1,1>;
        X32bdac [3][2][0][0] = xBDAC  <7,5,1,1>;
        X32bcad [3][2][0][0] = xBCAD  <7,5,1,1>;
        J32abcd [3][2][1][0] = jABCD  <7,5,3,1>;
        J32cdab [3][2][1][0] = jCDAB  <7,5,3,1>;
        X32acbd [3][2][1][0] = xACBD  <7,5,3,1>;
        X32adbc [3][2][1][0] = xADBC  <7,5,3,1>;
        X32bdac [3][2][1][0] = xBDAC  <7,5,3,1>;
        X32bcad [3][2][1][0] = xBCAD  <7,5,3,1>;
        J32abcd [3][2][1][1] = jABCD  <7,5,3,3>;
        J32cdab [3][2][1][1] = jCDAB  <7,5,3,3>;
        X32acbd [3][2][1][1] = xACBD  <7,5,3,3>;
        X32adbc [3][2][1][1] = xADBC  <7,5,3,3>;
        X32bdac [3][2][1][1] = xBDAC  <7,5,3,3>;
        X32bcad [3][2][1][1] = xBCAD  <7,5,3,3>;
        J32abcd [3][2][2][0] = jABCD  <7,5,5,1>;
        J32cdab [3][2][2][0] = jCDAB  <7,5,5,1>;
        X32acbd [3][2][2][0] = xACBD  <7,5,5,1>;
        X32adbc [3][2][2][0] = xADBC  <7,5,5,1>;
        X32bdac [3][2][2][0] = xBDAC  <7,5,5,1>;
        X32bcad [3][2][2][0] = xBCAD  <7,5,5,1>;
        J32abcd [3][2][2][1] = jABCD  <7,5,5,3>;
        J32cdab [3][2][2][1] = jCDAB  <7,5,5,3>;
        X32acbd [3][2][2][1] = xACBD  <7,5,5,3>;
        X32adbc [3][2][2][1] = xADBC  <7,5,5,3>;
        X32bdac [3][2][2][1] = xBDAC  <7,5,5,3>;
        X32bcad [3][2][2][1] = xBCAD  <7,5,5,3>;
        J32abcd [3][2][2][2] = jABCD  <7,5,5,5>;
        J32cdab [3][2][2][2] = jCDAB  <7,5,5,5>;
        X32acbd [3][2][2][2] = xACBD  <7,5,5,5>;
        X32adbc [3][2][2][2] = xADBC  <7,5,5,5>;
        X32bdac [3][2][2][2] = xBDAC  <7,5,5,5>;
        X32bcad [3][2][2][2] = xBCAD  <7,5,5,5>;
        J32abcd [3][2][3][0] = jABCD  <7,5,7,1>;
        J32cdab [3][2][3][0] = jCDAB  <7,5,7,1>;
        X32acbd [3][2][3][0] = xACBD  <7,5,7,1>;
        X32adbc [3][2][3][0] = xADBC  <7,5,7,1>;
        X32bdac [3][2][3][0] = xBDAC  <7,5,7,1>;
        X32bcad [3][2][3][0] = xBCAD  <7,5,7,1>;
        J32abcd [3][2][3][1] = jABCD  <7,5,7,3>;
        J32cdab [3][2][3][1] = jCDAB  <7,5,7,3>;
        X32acbd [3][2][3][1] = xACBD  <7,5,7,3>;
        X32adbc [3][2][3][1] = xADBC  <7,5,7,3>;
        X32bdac [3][2][3][1] = xBDAC  <7,5,7,3>;
        X32bcad [3][2][3][1] = xBCAD  <7,5,7,3>;
        J32abcd [3][2][3][2] = jABCD  <7,5,7,5>;
        J32cdab [3][2][3][2] = jCDAB  <7,5,7,5>;
        X32acbd [3][2][3][2] = xACBD  <7,5,7,5>;
        X32adbc [3][2][3][2] = xADBC  <7,5,7,5>;
        X32bdac [3][2][3][2] = xBDAC  <7,5,7,5>;
        X32bcad [3][2][3][2] = xBCAD  <7,5,7,5>;
        J32abcd [3][2][3][3] = jABCD  <7,5,7,7>;
        J32cdab [3][2][3][3] = jCDAB  <7,5,7,7>;
        X32acbd [3][2][3][3] = xACBD  <7,5,7,7>;
        X32adbc [3][2][3][3] = xADBC  <7,5,7,7>;
        X32bdac [3][2][3][3] = xBDAC  <7,5,7,7>;
        X32bcad [3][2][3][3] = xBCAD  <7,5,7,7>;
        J32abcd [3][2][4][0] = jABCD  <7,5,9,1>;
        J32cdab [3][2][4][0] = jCDAB  <7,5,9,1>;
        X32acbd [3][2][4][0] = xACBD  <7,5,9,1>;
        X32adbc [3][2][4][0] = xADBC  <7,5,9,1>;
        X32bdac [3][2][4][0] = xBDAC  <7,5,9,1>;
        X32bcad [3][2][4][0] = xBCAD  <7,5,9,1>;
        J32abcd [3][2][4][1] = jABCD  <7,5,9,3>;
        J32cdab [3][2][4][1] = jCDAB  <7,5,9,3>;
        X32acbd [3][2][4][1] = xACBD  <7,5,9,3>;
        X32adbc [3][2][4][1] = xADBC  <7,5,9,3>;
        X32bdac [3][2][4][1] = xBDAC  <7,5,9,3>;
        X32bcad [3][2][4][1] = xBCAD  <7,5,9,3>;
        J32abcd [3][2][4][2] = jABCD  <7,5,9,5>;
        J32cdab [3][2][4][2] = jCDAB  <7,5,9,5>;
        X32acbd [3][2][4][2] = xACBD  <7,5,9,5>;
        X32adbc [3][2][4][2] = xADBC  <7,5,9,5>;
        X32bdac [3][2][4][2] = xBDAC  <7,5,9,5>;
        X32bcad [3][2][4][2] = xBCAD  <7,5,9,5>;
        J32abcd [3][2][4][3] = jABCD  <7,5,9,7>;
        J32cdab [3][2][4][3] = jCDAB  <7,5,9,7>;
        X32acbd [3][2][4][3] = xACBD  <7,5,9,7>;
        X32adbc [3][2][4][3] = xADBC  <7,5,9,7>;
        X32bdac [3][2][4][3] = xBDAC  <7,5,9,7>;
        X32bcad [3][2][4][3] = xBCAD  <7,5,9,7>;
        J32abcd [3][2][4][4] = jABCD  <7,5,9,9>;
        J32cdab [3][2][4][4] = jCDAB  <7,5,9,9>;
        X32acbd [3][2][4][4] = xACBD  <7,5,9,9>;
        X32adbc [3][2][4][4] = xADBC  <7,5,9,9>;
        X32bdac [3][2][4][4] = xBDAC  <7,5,9,9>;
        X32bcad [3][2][4][4] = xBCAD  <7,5,9,9>;

        J32abcd [3][3][0][0] = jABCD  <7,7,1,1>;
        J32cdab [3][3][0][0] = jCDAB  <7,7,1,1>;
        X32acbd [3][3][0][0] = xACBD  <7,7,1,1>;
        X32adbc [3][3][0][0] = xADBC  <7,7,1,1>;
        X32bdac [3][3][0][0] = xBDAC  <7,7,1,1>;
        X32bcad [3][3][0][0] = xBCAD  <7,7,1,1>;
        J32abcd [3][3][1][0] = jABCD  <7,7,3,1>;
        J32cdab [3][3][1][0] = jCDAB  <7,7,3,1>;
        X32acbd [3][3][1][0] = xACBD  <7,7,3,1>;
        X32adbc [3][3][1][0] = xADBC  <7,7,3,1>;
        X32bdac [3][3][1][0] = xBDAC  <7,7,3,1>;
        X32bcad [3][3][1][0] = xBCAD  <7,7,3,1>;
        J32abcd [3][3][1][1] = jABCD  <7,7,3,3>;
        J32cdab [3][3][1][1] = jCDAB  <7,7,3,3>;
        X32acbd [3][3][1][1] = xACBD  <7,7,3,3>;
        X32adbc [3][3][1][1] = xADBC  <7,7,3,3>;
        X32bdac [3][3][1][1] = xBDAC  <7,7,3,3>;
        X32bcad [3][3][1][1] = xBCAD  <7,7,3,3>;
        J32abcd [3][3][2][0] = jABCD  <7,7,5,1>;
        J32cdab [3][3][2][0] = jCDAB  <7,7,5,1>;
        X32acbd [3][3][2][0] = xACBD  <7,7,5,1>;
        X32adbc [3][3][2][0] = xADBC  <7,7,5,1>;
        X32bdac [3][3][2][0] = xBDAC  <7,7,5,1>;
        X32bcad [3][3][2][0] = xBCAD  <7,7,5,1>;
        J32abcd [3][3][2][1] = jABCD  <7,7,5,3>;
        J32cdab [3][3][2][1] = jCDAB  <7,7,5,3>;
        X32acbd [3][3][2][1] = xACBD  <7,7,5,3>;
        X32adbc [3][3][2][1] = xADBC  <7,7,5,3>;
        X32bdac [3][3][2][1] = xBDAC  <7,7,5,3>;
        X32bcad [3][3][2][1] = xBCAD  <7,7,5,3>;
        J32abcd [3][3][2][2] = jABCD  <7,7,5,5>;
        J32cdab [3][3][2][2] = jCDAB  <7,7,5,5>;
        X32acbd [3][3][2][2] = xACBD  <7,7,5,5>;
        X32adbc [3][3][2][2] = xADBC  <7,7,5,5>;
        X32bdac [3][3][2][2] = xBDAC  <7,7,5,5>;
        X32bcad [3][3][2][2] = xBCAD  <7,7,5,5>;
        J32abcd [3][3][3][0] = jABCD  <7,7,7,1>;
        J32cdab [3][3][3][0] = jCDAB  <7,7,7,1>;
        X32acbd [3][3][3][0] = xACBD  <7,7,7,1>;
        X32adbc [3][3][3][0] = xADBC  <7,7,7,1>;
        X32bdac [3][3][3][0] = xBDAC  <7,7,7,1>;
        X32bcad [3][3][3][0] = xBCAD  <7,7,7,1>;
        J32abcd [3][3][3][1] = jABCD  <7,7,7,3>;
        J32cdab [3][3][3][1] = jCDAB  <7,7,7,3>;
        X32acbd [3][3][3][1] = xACBD  <7,7,7,3>;
        X32adbc [3][3][3][1] = xADBC  <7,7,7,3>;
        X32bdac [3][3][3][1] = xBDAC  <7,7,7,3>;
        X32bcad [3][3][3][1] = xBCAD  <7,7,7,3>;
        J32abcd [3][3][3][2] = jABCD  <7,7,7,5>;
        J32cdab [3][3][3][2] = jCDAB  <7,7,7,5>;
        X32acbd [3][3][3][2] = xACBD  <7,7,7,5>;
        X32adbc [3][3][3][2] = xADBC  <7,7,7,5>;
        X32bdac [3][3][3][2] = xBDAC  <7,7,7,5>;
        X32bcad [3][3][3][2] = xBCAD  <7,7,7,5>;
        J32abcd [3][3][3][3] = jABCD  <7,7,7,7>;
        J32cdab [3][3][3][3] = jCDAB  <7,7,7,7>;
        X32acbd [3][3][3][3] = xACBD  <7,7,7,7>;
        X32adbc [3][3][3][3] = xADBC  <7,7,7,7>;
        X32bdac [3][3][3][3] = xBDAC  <7,7,7,7>;
        X32bcad [3][3][3][3] = xBCAD  <7,7,7,7>;
        J32abcd [3][3][4][0] = jABCD  <7,7,9,1>;
        J32cdab [3][3][4][0] = jCDAB  <7,7,9,1>;
        X32acbd [3][3][4][0] = xACBD  <7,7,9,1>;
        X32adbc [3][3][4][0] = xADBC  <7,7,9,1>;
        X32bdac [3][3][4][0] = xBDAC  <7,7,9,1>;
        X32bcad [3][3][4][0] = xBCAD  <7,7,9,1>;
        J32abcd [3][3][4][1] = jABCD  <7,7,9,3>;
        J32cdab [3][3][4][1] = jCDAB  <7,7,9,3>;
        X32acbd [3][3][4][1] = xACBD  <7,7,9,3>;
        X32adbc [3][3][4][1] = xADBC  <7,7,9,3>;
        X32bdac [3][3][4][1] = xBDAC  <7,7,9,3>;
        X32bcad [3][3][4][1] = xBCAD  <7,7,9,3>;
        J32abcd [3][3][4][2] = jABCD  <7,7,9,5>;
        J32cdab [3][3][4][2] = jCDAB  <7,7,9,5>;
        X32acbd [3][3][4][2] = xACBD  <7,7,9,5>;
        X32adbc [3][3][4][2] = xADBC  <7,7,9,5>;
        X32bdac [3][3][4][2] = xBDAC  <7,7,9,5>;
        X32bcad [3][3][4][2] = xBCAD  <7,7,9,5>;
        J32abcd [3][3][4][3] = jABCD  <7,7,9,7>;
        J32cdab [3][3][4][3] = jCDAB  <7,7,9,7>;
        X32acbd [3][3][4][3] = xACBD  <7,7,9,7>;
        X32adbc [3][3][4][3] = xADBC  <7,7,9,7>;
        X32bdac [3][3][4][3] = xBDAC  <7,7,9,7>;
        X32bcad [3][3][4][3] = xBCAD  <7,7,9,7>;
        J32abcd [3][3][4][4] = jABCD  <7,7,9,9>;
        J32cdab [3][3][4][4] = jCDAB  <7,7,9,9>;
        X32acbd [3][3][4][4] = xACBD  <7,7,9,9>;
        X32adbc [3][3][4][4] = xADBC  <7,7,9,9>;
        X32bdac [3][3][4][4] = xBDAC  <7,7,9,9>;
        X32bcad [3][3][4][4] = xBCAD  <7,7,9,9>;

        J32abcd [4][0][0][0] = jABCD  <9,1,1,1>;
        J32cdab [4][0][0][0] = jCDAB  <9,1,1,1>;
        X32acbd [4][0][0][0] = xACBD  <9,1,1,1>;
        X32adbc [4][0][0][0] = xADBC  <9,1,1,1>;
        X32bdac [4][0][0][0] = xBDAC  <9,1,1,1>;
        X32bcad [4][0][0][0] = xBCAD  <9,1,1,1>;
        J32abcd [4][0][1][0] = jABCD  <9,1,3,1>;
        J32cdab [4][0][1][0] = jCDAB  <9,1,3,1>;
        X32acbd [4][0][1][0] = xACBD  <9,1,3,1>;
        X32adbc [4][0][1][0] = xADBC  <9,1,3,1>;
        X32bdac [4][0][1][0] = xBDAC  <9,1,3,1>;
        X32bcad [4][0][1][0] = xBCAD  <9,1,3,1>;
        J32abcd [4][0][1][1] = jABCD  <9,1,3,3>;
        J32cdab [4][0][1][1] = jCDAB  <9,1,3,3>;
        X32acbd [4][0][1][1] = xACBD  <9,1,3,3>;
        X32adbc [4][0][1][1] = xADBC  <9,1,3,3>;
        X32bdac [4][0][1][1] = xBDAC  <9,1,3,3>;
        X32bcad [4][0][1][1] = xBCAD  <9,1,3,3>;
        J32abcd [4][0][2][0] = jABCD  <9,1,5,1>;
        J32cdab [4][0][2][0] = jCDAB  <9,1,5,1>;
        X32acbd [4][0][2][0] = xACBD  <9,1,5,1>;
        X32adbc [4][0][2][0] = xADBC  <9,1,5,1>;
        X32bdac [4][0][2][0] = xBDAC  <9,1,5,1>;
        X32bcad [4][0][2][0] = xBCAD  <9,1,5,1>;
        J32abcd [4][0][2][1] = jABCD  <9,1,5,3>;
        J32cdab [4][0][2][1] = jCDAB  <9,1,5,3>;
        X32acbd [4][0][2][1] = xACBD  <9,1,5,3>;
        X32adbc [4][0][2][1] = xADBC  <9,1,5,3>;
        X32bdac [4][0][2][1] = xBDAC  <9,1,5,3>;
        X32bcad [4][0][2][1] = xBCAD  <9,1,5,3>;
        J32abcd [4][0][2][2] = jABCD  <9,1,5,5>;
        J32cdab [4][0][2][2] = jCDAB  <9,1,5,5>;
        X32acbd [4][0][2][2] = xACBD  <9,1,5,5>;
        X32adbc [4][0][2][2] = xADBC  <9,1,5,5>;
        X32bdac [4][0][2][2] = xBDAC  <9,1,5,5>;
        X32bcad [4][0][2][2] = xBCAD  <9,1,5,5>;
        J32abcd [4][0][3][0] = jABCD  <9,1,7,1>;
        J32cdab [4][0][3][0] = jCDAB  <9,1,7,1>;
        X32acbd [4][0][3][0] = xACBD  <9,1,7,1>;
        X32adbc [4][0][3][0] = xADBC  <9,1,7,1>;
        X32bdac [4][0][3][0] = xBDAC  <9,1,7,1>;
        X32bcad [4][0][3][0] = xBCAD  <9,1,7,1>;
        J32abcd [4][0][3][1] = jABCD  <9,1,7,3>;
        J32cdab [4][0][3][1] = jCDAB  <9,1,7,3>;
        X32acbd [4][0][3][1] = xACBD  <9,1,7,3>;
        X32adbc [4][0][3][1] = xADBC  <9,1,7,3>;
        X32bdac [4][0][3][1] = xBDAC  <9,1,7,3>;
        X32bcad [4][0][3][1] = xBCAD  <9,1,7,3>;
        J32abcd [4][0][3][2] = jABCD  <9,1,7,5>;
        J32cdab [4][0][3][2] = jCDAB  <9,1,7,5>;
        X32acbd [4][0][3][2] = xACBD  <9,1,7,5>;
        X32adbc [4][0][3][2] = xADBC  <9,1,7,5>;
        X32bdac [4][0][3][2] = xBDAC  <9,1,7,5>;
        X32bcad [4][0][3][2] = xBCAD  <9,1,7,5>;
        J32abcd [4][0][3][3] = jABCD  <9,1,7,7>;
        J32cdab [4][0][3][3] = jCDAB  <9,1,7,7>;
        X32acbd [4][0][3][3] = xACBD  <9,1,7,7>;
        X32adbc [4][0][3][3] = xADBC  <9,1,7,7>;
        X32bdac [4][0][3][3] = xBDAC  <9,1,7,7>;
        X32bcad [4][0][3][3] = xBCAD  <9,1,7,7>;
        J32abcd [4][0][4][0] = jABCD  <9,1,9,1>;
        J32cdab [4][0][4][0] = jCDAB  <9,1,9,1>;
        X32acbd [4][0][4][0] = xACBD  <9,1,9,1>;
        X32adbc [4][0][4][0] = xADBC  <9,1,9,1>;
        X32bdac [4][0][4][0] = xBDAC  <9,1,9,1>;
        X32bcad [4][0][4][0] = xBCAD  <9,1,9,1>;
        J32abcd [4][0][4][1] = jABCD  <9,1,9,3>;
        J32cdab [4][0][4][1] = jCDAB  <9,1,9,3>;
        X32acbd [4][0][4][1] = xACBD  <9,1,9,3>;
        X32adbc [4][0][4][1] = xADBC  <9,1,9,3>;
        X32bdac [4][0][4][1] = xBDAC  <9,1,9,3>;
        X32bcad [4][0][4][1] = xBCAD  <9,1,9,3>;
        J32abcd [4][0][4][2] = jABCD  <9,1,9,5>;
        J32cdab [4][0][4][2] = jCDAB  <9,1,9,5>;
        X32acbd [4][0][4][2] = xACBD  <9,1,9,5>;
        X32adbc [4][0][4][2] = xADBC  <9,1,9,5>;
        X32bdac [4][0][4][2] = xBDAC  <9,1,9,5>;
        X32bcad [4][0][4][2] = xBCAD  <9,1,9,5>;
        J32abcd [4][0][4][3] = jABCD  <9,1,9,7>;
        J32cdab [4][0][4][3] = jCDAB  <9,1,9,7>;
        X32acbd [4][0][4][3] = xACBD  <9,1,9,7>;
        X32adbc [4][0][4][3] = xADBC  <9,1,9,7>;
        X32bdac [4][0][4][3] = xBDAC  <9,1,9,7>;
        X32bcad [4][0][4][3] = xBCAD  <9,1,9,7>;
        J32abcd [4][0][4][4] = jABCD  <9,1,9,9>;
        J32cdab [4][0][4][4] = jCDAB  <9,1,9,9>;
        X32acbd [4][0][4][4] = xACBD  <9,1,9,9>;
        X32adbc [4][0][4][4] = xADBC  <9,1,9,9>;
        X32bdac [4][0][4][4] = xBDAC  <9,1,9,9>;
        X32bcad [4][0][4][4] = xBCAD  <9,1,9,9>;

        J32abcd [4][1][0][0] = jABCD  <9,3,1,1>;
        J32cdab [4][1][0][0] = jCDAB  <9,3,1,1>;
        X32acbd [4][1][0][0] = xACBD  <9,3,1,1>;
        X32adbc [4][1][0][0] = xADBC  <9,3,1,1>;
        X32bdac [4][1][0][0] = xBDAC  <9,3,1,1>;
        X32bcad [4][1][0][0] = xBCAD  <9,3,1,1>;
        J32abcd [4][1][1][0] = jABCD  <9,3,3,1>;
        J32cdab [4][1][1][0] = jCDAB  <9,3,3,1>;
        X32acbd [4][1][1][0] = xACBD  <9,3,3,1>;
        X32adbc [4][1][1][0] = xADBC  <9,3,3,1>;
        X32bdac [4][1][1][0] = xBDAC  <9,3,3,1>;
        X32bcad [4][1][1][0] = xBCAD  <9,3,3,1>;
        J32abcd [4][1][1][1] = jABCD  <9,3,3,3>;
        J32cdab [4][1][1][1] = jCDAB  <9,3,3,3>;
        X32acbd [4][1][1][1] = xACBD  <9,3,3,3>;
        X32adbc [4][1][1][1] = xADBC  <9,3,3,3>;
        X32bdac [4][1][1][1] = xBDAC  <9,3,3,3>;
        X32bcad [4][1][1][1] = xBCAD  <9,3,3,3>;
        J32abcd [4][1][2][0] = jABCD  <9,3,5,1>;
        J32cdab [4][1][2][0] = jCDAB  <9,3,5,1>;
        X32acbd [4][1][2][0] = xACBD  <9,3,5,1>;
        X32adbc [4][1][2][0] = xADBC  <9,3,5,1>;
        X32bdac [4][1][2][0] = xBDAC  <9,3,5,1>;
        X32bcad [4][1][2][0] = xBCAD  <9,3,5,1>;
        J32abcd [4][1][2][1] = jABCD  <9,3,5,3>;
        J32cdab [4][1][2][1] = jCDAB  <9,3,5,3>;
        X32acbd [4][1][2][1] = xACBD  <9,3,5,3>;
        X32adbc [4][1][2][1] = xADBC  <9,3,5,3>;
        X32bdac [4][1][2][1] = xBDAC  <9,3,5,3>;
        X32bcad [4][1][2][1] = xBCAD  <9,3,5,3>;
        J32abcd [4][1][2][2] = jABCD  <9,3,5,5>;
        J32cdab [4][1][2][2] = jCDAB  <9,3,5,5>;
        X32acbd [4][1][2][2] = xACBD  <9,3,5,5>;
        X32adbc [4][1][2][2] = xADBC  <9,3,5,5>;
        X32bdac [4][1][2][2] = xBDAC  <9,3,5,5>;
        X32bcad [4][1][2][2] = xBCAD  <9,3,5,5>;
        J32abcd [4][1][3][0] = jABCD  <9,3,7,1>;
        J32cdab [4][1][3][0] = jCDAB  <9,3,7,1>;
        X32acbd [4][1][3][0] = xACBD  <9,3,7,1>;
        X32adbc [4][1][3][0] = xADBC  <9,3,7,1>;
        X32bdac [4][1][3][0] = xBDAC  <9,3,7,1>;
        X32bcad [4][1][3][0] = xBCAD  <9,3,7,1>;
        J32abcd [4][1][3][1] = jABCD  <9,3,7,3>;
        J32cdab [4][1][3][1] = jCDAB  <9,3,7,3>;
        X32acbd [4][1][3][1] = xACBD  <9,3,7,3>;
        X32adbc [4][1][3][1] = xADBC  <9,3,7,3>;
        X32bdac [4][1][3][1] = xBDAC  <9,3,7,3>;
        X32bcad [4][1][3][1] = xBCAD  <9,3,7,3>;
        J32abcd [4][1][3][2] = jABCD  <9,3,7,5>;
        J32cdab [4][1][3][2] = jCDAB  <9,3,7,5>;
        X32acbd [4][1][3][2] = xACBD  <9,3,7,5>;
        X32adbc [4][1][3][2] = xADBC  <9,3,7,5>;
        X32bdac [4][1][3][2] = xBDAC  <9,3,7,5>;
        X32bcad [4][1][3][2] = xBCAD  <9,3,7,5>;
        J32abcd [4][1][3][3] = jABCD  <9,3,7,7>;
        J32cdab [4][1][3][3] = jCDAB  <9,3,7,7>;
        X32acbd [4][1][3][3] = xACBD  <9,3,7,7>;
        X32adbc [4][1][3][3] = xADBC  <9,3,7,7>;
        X32bdac [4][1][3][3] = xBDAC  <9,3,7,7>;
        X32bcad [4][1][3][3] = xBCAD  <9,3,7,7>;
        J32abcd [4][1][4][0] = jABCD  <9,3,9,1>;
        J32cdab [4][1][4][0] = jCDAB  <9,3,9,1>;
        X32acbd [4][1][4][0] = xACBD  <9,3,9,1>;
        X32adbc [4][1][4][0] = xADBC  <9,3,9,1>;
        X32bdac [4][1][4][0] = xBDAC  <9,3,9,1>;
        X32bcad [4][1][4][0] = xBCAD  <9,3,9,1>;
        J32abcd [4][1][4][1] = jABCD  <9,3,9,3>;
        J32cdab [4][1][4][1] = jCDAB  <9,3,9,3>;
        X32acbd [4][1][4][1] = xACBD  <9,3,9,3>;
        X32adbc [4][1][4][1] = xADBC  <9,3,9,3>;
        X32bdac [4][1][4][1] = xBDAC  <9,3,9,3>;
        X32bcad [4][1][4][1] = xBCAD  <9,3,9,3>;
        J32abcd [4][1][4][2] = jABCD  <9,3,9,5>;
        J32cdab [4][1][4][2] = jCDAB  <9,3,9,5>;
        X32acbd [4][1][4][2] = xACBD  <9,3,9,5>;
        X32adbc [4][1][4][2] = xADBC  <9,3,9,5>;
        X32bdac [4][1][4][2] = xBDAC  <9,3,9,5>;
        X32bcad [4][1][4][2] = xBCAD  <9,3,9,5>;
        J32abcd [4][1][4][3] = jABCD  <9,3,9,7>;
        J32cdab [4][1][4][3] = jCDAB  <9,3,9,7>;
        X32acbd [4][1][4][3] = xACBD  <9,3,9,7>;
        X32adbc [4][1][4][3] = xADBC  <9,3,9,7>;
        X32bdac [4][1][4][3] = xBDAC  <9,3,9,7>;
        X32bcad [4][1][4][3] = xBCAD  <9,3,9,7>;
        J32abcd [4][1][4][4] = jABCD  <9,3,9,9>;
        J32cdab [4][1][4][4] = jCDAB  <9,3,9,9>;
        X32acbd [4][1][4][4] = xACBD  <9,3,9,9>;
        X32adbc [4][1][4][4] = xADBC  <9,3,9,9>;
        X32bdac [4][1][4][4] = xBDAC  <9,3,9,9>;
        X32bcad [4][1][4][4] = xBCAD  <9,3,9,9>;

        J32abcd [4][2][0][0] = jABCD  <9,5,1,1>;
        J32cdab [4][2][0][0] = jCDAB  <9,5,1,1>;
        X32acbd [4][2][0][0] = xACBD  <9,5,1,1>;
        X32adbc [4][2][0][0] = xADBC  <9,5,1,1>;
        X32bdac [4][2][0][0] = xBDAC  <9,5,1,1>;
        X32bcad [4][2][0][0] = xBCAD  <9,5,1,1>;
        J32abcd [4][2][1][0] = jABCD  <9,5,3,1>;
        J32cdab [4][2][1][0] = jCDAB  <9,5,3,1>;
        X32acbd [4][2][1][0] = xACBD  <9,5,3,1>;
        X32adbc [4][2][1][0] = xADBC  <9,5,3,1>;
        X32bdac [4][2][1][0] = xBDAC  <9,5,3,1>;
        X32bcad [4][2][1][0] = xBCAD  <9,5,3,1>;
        J32abcd [4][2][1][1] = jABCD  <9,5,3,3>;
        J32cdab [4][2][1][1] = jCDAB  <9,5,3,3>;
        X32acbd [4][2][1][1] = xACBD  <9,5,3,3>;
        X32adbc [4][2][1][1] = xADBC  <9,5,3,3>;
        X32bdac [4][2][1][1] = xBDAC  <9,5,3,3>;
        X32bcad [4][2][1][1] = xBCAD  <9,5,3,3>;
        J32abcd [4][2][2][0] = jABCD  <9,5,5,1>;
        J32cdab [4][2][2][0] = jCDAB  <9,5,5,1>;
        X32acbd [4][2][2][0] = xACBD  <9,5,5,1>;
        X32adbc [4][2][2][0] = xADBC  <9,5,5,1>;
        X32bdac [4][2][2][0] = xBDAC  <9,5,5,1>;
        X32bcad [4][2][2][0] = xBCAD  <9,5,5,1>;
        J32abcd [4][2][2][1] = jABCD  <9,5,5,3>;
        J32cdab [4][2][2][1] = jCDAB  <9,5,5,3>;
        X32acbd [4][2][2][1] = xACBD  <9,5,5,3>;
        X32adbc [4][2][2][1] = xADBC  <9,5,5,3>;
        X32bdac [4][2][2][1] = xBDAC  <9,5,5,3>;
        X32bcad [4][2][2][1] = xBCAD  <9,5,5,3>;
        J32abcd [4][2][2][2] = jABCD  <9,5,5,5>;
        J32cdab [4][2][2][2] = jCDAB  <9,5,5,5>;
        X32acbd [4][2][2][2] = xACBD  <9,5,5,5>;
        X32adbc [4][2][2][2] = xADBC  <9,5,5,5>;
        X32bdac [4][2][2][2] = xBDAC  <9,5,5,5>;
        X32bcad [4][2][2][2] = xBCAD  <9,5,5,5>;
        J32abcd [4][2][3][0] = jABCD  <9,5,7,1>;
        J32cdab [4][2][3][0] = jCDAB  <9,5,7,1>;
        X32acbd [4][2][3][0] = xACBD  <9,5,7,1>;
        X32adbc [4][2][3][0] = xADBC  <9,5,7,1>;
        X32bdac [4][2][3][0] = xBDAC  <9,5,7,1>;
        X32bcad [4][2][3][0] = xBCAD  <9,5,7,1>;
        J32abcd [4][2][3][1] = jABCD  <9,5,7,3>;
        J32cdab [4][2][3][1] = jCDAB  <9,5,7,3>;
        X32acbd [4][2][3][1] = xACBD  <9,5,7,3>;
        X32adbc [4][2][3][1] = xADBC  <9,5,7,3>;
        X32bdac [4][2][3][1] = xBDAC  <9,5,7,3>;
        X32bcad [4][2][3][1] = xBCAD  <9,5,7,3>;
        J32abcd [4][2][3][2] = jABCD  <9,5,7,5>;
        J32cdab [4][2][3][2] = jCDAB  <9,5,7,5>;
        X32acbd [4][2][3][2] = xACBD  <9,5,7,5>;
        X32adbc [4][2][3][2] = xADBC  <9,5,7,5>;
        X32bdac [4][2][3][2] = xBDAC  <9,5,7,5>;
        X32bcad [4][2][3][2] = xBCAD  <9,5,7,5>;
        J32abcd [4][2][3][3] = jABCD  <9,5,7,7>;
        J32cdab [4][2][3][3] = jCDAB  <9,5,7,7>;
        X32acbd [4][2][3][3] = xACBD  <9,5,7,7>;
        X32adbc [4][2][3][3] = xADBC  <9,5,7,7>;
        X32bdac [4][2][3][3] = xBDAC  <9,5,7,7>;
        X32bcad [4][2][3][3] = xBCAD  <9,5,7,7>;
        J32abcd [4][2][4][0] = jABCD  <9,5,9,1>;
        J32cdab [4][2][4][0] = jCDAB  <9,5,9,1>;
        X32acbd [4][2][4][0] = xACBD  <9,5,9,1>;
        X32adbc [4][2][4][0] = xADBC  <9,5,9,1>;
        X32bdac [4][2][4][0] = xBDAC  <9,5,9,1>;
        X32bcad [4][2][4][0] = xBCAD  <9,5,9,1>;
        J32abcd [4][2][4][1] = jABCD  <9,5,9,3>;
        J32cdab [4][2][4][1] = jCDAB  <9,5,9,3>;
        X32acbd [4][2][4][1] = xACBD  <9,5,9,3>;
        X32adbc [4][2][4][1] = xADBC  <9,5,9,3>;
        X32bdac [4][2][4][1] = xBDAC  <9,5,9,3>;
        X32bcad [4][2][4][1] = xBCAD  <9,5,9,3>;
        J32abcd [4][2][4][2] = jABCD  <9,5,9,5>;
        J32cdab [4][2][4][2] = jCDAB  <9,5,9,5>;
        X32acbd [4][2][4][2] = xACBD  <9,5,9,5>;
        X32adbc [4][2][4][2] = xADBC  <9,5,9,5>;
        X32bdac [4][2][4][2] = xBDAC  <9,5,9,5>;
        X32bcad [4][2][4][2] = xBCAD  <9,5,9,5>;
        J32abcd [4][2][4][3] = jABCD  <9,5,9,7>;
        J32cdab [4][2][4][3] = jCDAB  <9,5,9,7>;
        X32acbd [4][2][4][3] = xACBD  <9,5,9,7>;
        X32adbc [4][2][4][3] = xADBC  <9,5,9,7>;
        X32bdac [4][2][4][3] = xBDAC  <9,5,9,7>;
        X32bcad [4][2][4][3] = xBCAD  <9,5,9,7>;
        J32abcd [4][2][4][4] = jABCD  <9,5,9,9>;
        J32cdab [4][2][4][4] = jCDAB  <9,5,9,9>;
        X32acbd [4][2][4][4] = xACBD  <9,5,9,9>;
        X32adbc [4][2][4][4] = xADBC  <9,5,9,9>;
        X32bdac [4][2][4][4] = xBDAC  <9,5,9,9>;
        X32bcad [4][2][4][4] = xBCAD  <9,5,9,9>;

        J32abcd [4][3][0][0] = jABCD  <9,7,1,1>;
        J32cdab [4][3][0][0] = jCDAB  <9,7,1,1>;
        X32acbd [4][3][0][0] = xACBD  <9,7,1,1>;
        X32adbc [4][3][0][0] = xADBC  <9,7,1,1>;
        X32bdac [4][3][0][0] = xBDAC  <9,7,1,1>;
        X32bcad [4][3][0][0] = xBCAD  <9,7,1,1>;
        J32abcd [4][3][1][0] = jABCD  <9,7,3,1>;
        J32cdab [4][3][1][0] = jCDAB  <9,7,3,1>;
        X32acbd [4][3][1][0] = xACBD  <9,7,3,1>;
        X32adbc [4][3][1][0] = xADBC  <9,7,3,1>;
        X32bdac [4][3][1][0] = xBDAC  <9,7,3,1>;
        X32bcad [4][3][1][0] = xBCAD  <9,7,3,1>;
        J32abcd [4][3][1][1] = jABCD  <9,7,3,3>;
        J32cdab [4][3][1][1] = jCDAB  <9,7,3,3>;
        X32acbd [4][3][1][1] = xACBD  <9,7,3,3>;
        X32adbc [4][3][1][1] = xADBC  <9,7,3,3>;
        X32bdac [4][3][1][1] = xBDAC  <9,7,3,3>;
        X32bcad [4][3][1][1] = xBCAD  <9,7,3,3>;
        J32abcd [4][3][2][0] = jABCD  <9,7,5,1>;
        J32cdab [4][3][2][0] = jCDAB  <9,7,5,1>;
        X32acbd [4][3][2][0] = xACBD  <9,7,5,1>;
        X32adbc [4][3][2][0] = xADBC  <9,7,5,1>;
        X32bdac [4][3][2][0] = xBDAC  <9,7,5,1>;
        X32bcad [4][3][2][0] = xBCAD  <9,7,5,1>;
        J32abcd [4][3][2][1] = jABCD  <9,7,5,3>;
        J32cdab [4][3][2][1] = jCDAB  <9,7,5,3>;
        X32acbd [4][3][2][1] = xACBD  <9,7,5,3>;
        X32adbc [4][3][2][1] = xADBC  <9,7,5,3>;
        X32bdac [4][3][2][1] = xBDAC  <9,7,5,3>;
        X32bcad [4][3][2][1] = xBCAD  <9,7,5,3>;
        J32abcd [4][3][2][2] = jABCD  <9,7,5,5>;
        J32cdab [4][3][2][2] = jCDAB  <9,7,5,5>;
        X32acbd [4][3][2][2] = xACBD  <9,7,5,5>;
        X32adbc [4][3][2][2] = xADBC  <9,7,5,5>;
        X32bdac [4][3][2][2] = xBDAC  <9,7,5,5>;
        X32bcad [4][3][2][2] = xBCAD  <9,7,5,5>;
        J32abcd [4][3][3][0] = jABCD  <9,7,7,1>;
        J32cdab [4][3][3][0] = jCDAB  <9,7,7,1>;
        X32acbd [4][3][3][0] = xACBD  <9,7,7,1>;
        X32adbc [4][3][3][0] = xADBC  <9,7,7,1>;
        X32bdac [4][3][3][0] = xBDAC  <9,7,7,1>;
        X32bcad [4][3][3][0] = xBCAD  <9,7,7,1>;
        J32abcd [4][3][3][1] = jABCD  <9,7,7,3>;
        J32cdab [4][3][3][1] = jCDAB  <9,7,7,3>;
        X32acbd [4][3][3][1] = xACBD  <9,7,7,3>;
        X32adbc [4][3][3][1] = xADBC  <9,7,7,3>;
        X32bdac [4][3][3][1] = xBDAC  <9,7,7,3>;
        X32bcad [4][3][3][1] = xBCAD  <9,7,7,3>;
        J32abcd [4][3][3][2] = jABCD  <9,7,7,5>;
        J32cdab [4][3][3][2] = jCDAB  <9,7,7,5>;
        X32acbd [4][3][3][2] = xACBD  <9,7,7,5>;
        X32adbc [4][3][3][2] = xADBC  <9,7,7,5>;
        X32bdac [4][3][3][2] = xBDAC  <9,7,7,5>;
        X32bcad [4][3][3][2] = xBCAD  <9,7,7,5>;
        J32abcd [4][3][3][3] = jABCD  <9,7,7,7>;
        J32cdab [4][3][3][3] = jCDAB  <9,7,7,7>;
        X32acbd [4][3][3][3] = xACBD  <9,7,7,7>;
        X32adbc [4][3][3][3] = xADBC  <9,7,7,7>;
        X32bdac [4][3][3][3] = xBDAC  <9,7,7,7>;
        X32bcad [4][3][3][3] = xBCAD  <9,7,7,7>;
        J32abcd [4][3][4][0] = jABCD  <9,7,9,1>;
        J32cdab [4][3][4][0] = jCDAB  <9,7,9,1>;
        X32acbd [4][3][4][0] = xACBD  <9,7,9,1>;
        X32adbc [4][3][4][0] = xADBC  <9,7,9,1>;
        X32bdac [4][3][4][0] = xBDAC  <9,7,9,1>;
        X32bcad [4][3][4][0] = xBCAD  <9,7,9,1>;
        J32abcd [4][3][4][1] = jABCD  <9,7,9,3>;
        J32cdab [4][3][4][1] = jCDAB  <9,7,9,3>;
        X32acbd [4][3][4][1] = xACBD  <9,7,9,3>;
        X32adbc [4][3][4][1] = xADBC  <9,7,9,3>;
        X32bdac [4][3][4][1] = xBDAC  <9,7,9,3>;
        X32bcad [4][3][4][1] = xBCAD  <9,7,9,3>;
        J32abcd [4][3][4][2] = jABCD  <9,7,9,5>;
        J32cdab [4][3][4][2] = jCDAB  <9,7,9,5>;
        X32acbd [4][3][4][2] = xACBD  <9,7,9,5>;
        X32adbc [4][3][4][2] = xADBC  <9,7,9,5>;
        X32bdac [4][3][4][2] = xBDAC  <9,7,9,5>;
        X32bcad [4][3][4][2] = xBCAD  <9,7,9,5>;
        J32abcd [4][3][4][3] = jABCD  <9,7,9,7>;
        J32cdab [4][3][4][3] = jCDAB  <9,7,9,7>;
        X32acbd [4][3][4][3] = xACBD  <9,7,9,7>;
        X32adbc [4][3][4][3] = xADBC  <9,7,9,7>;
        X32bdac [4][3][4][3] = xBDAC  <9,7,9,7>;
        X32bcad [4][3][4][3] = xBCAD  <9,7,9,7>;
        J32abcd [4][3][4][4] = jABCD  <9,7,9,9>;
        J32cdab [4][3][4][4] = jCDAB  <9,7,9,9>;
        X32acbd [4][3][4][4] = xACBD  <9,7,9,9>;
        X32adbc [4][3][4][4] = xADBC  <9,7,9,9>;
        X32bdac [4][3][4][4] = xBDAC  <9,7,9,9>;
        X32bcad [4][3][4][4] = xBCAD  <9,7,9,9>;

        J32abcd [4][4][0][0] = jABCD  <9,9,1,1>;
        J32cdab [4][4][0][0] = jCDAB  <9,9,1,1>;
        X32acbd [4][4][0][0] = xACBD  <9,9,1,1>;
        X32adbc [4][4][0][0] = xADBC  <9,9,1,1>;
        X32bdac [4][4][0][0] = xBDAC  <9,9,1,1>;
        X32bcad [4][4][0][0] = xBCAD  <9,9,1,1>;
        J32abcd [4][4][1][0] = jABCD  <9,9,3,1>;
        J32cdab [4][4][1][0] = jCDAB  <9,9,3,1>;
        X32acbd [4][4][1][0] = xACBD  <9,9,3,1>;
        X32adbc [4][4][1][0] = xADBC  <9,9,3,1>;
        X32bdac [4][4][1][0] = xBDAC  <9,9,3,1>;
        X32bcad [4][4][1][0] = xBCAD  <9,9,3,1>;
        J32abcd [4][4][1][1] = jABCD  <9,9,3,3>;
        J32cdab [4][4][1][1] = jCDAB  <9,9,3,3>;
        X32acbd [4][4][1][1] = xACBD  <9,9,3,3>;
        X32adbc [4][4][1][1] = xADBC  <9,9,3,3>;
        X32bdac [4][4][1][1] = xBDAC  <9,9,3,3>;
        X32bcad [4][4][1][1] = xBCAD  <9,9,3,3>;
        J32abcd [4][4][2][0] = jABCD  <9,9,5,1>;
        J32cdab [4][4][2][0] = jCDAB  <9,9,5,1>;
        X32acbd [4][4][2][0] = xACBD  <9,9,5,1>;
        X32adbc [4][4][2][0] = xADBC  <9,9,5,1>;
        X32bdac [4][4][2][0] = xBDAC  <9,9,5,1>;
        X32bcad [4][4][2][0] = xBCAD  <9,9,5,1>;
        J32abcd [4][4][2][1] = jABCD  <9,9,5,3>;
        J32cdab [4][4][2][1] = jCDAB  <9,9,5,3>;
        X32acbd [4][4][2][1] = xACBD  <9,9,5,3>;
        X32adbc [4][4][2][1] = xADBC  <9,9,5,3>;
        X32bdac [4][4][2][1] = xBDAC  <9,9,5,3>;
        X32bcad [4][4][2][1] = xBCAD  <9,9,5,3>;
        J32abcd [4][4][2][2] = jABCD  <9,9,5,5>;
        J32cdab [4][4][2][2] = jCDAB  <9,9,5,5>;
        X32acbd [4][4][2][2] = xACBD  <9,9,5,5>;
        X32adbc [4][4][2][2] = xADBC  <9,9,5,5>;
        X32bdac [4][4][2][2] = xBDAC  <9,9,5,5>;
        X32bcad [4][4][2][2] = xBCAD  <9,9,5,5>;
        J32abcd [4][4][3][0] = jABCD  <9,9,7,1>;
        J32cdab [4][4][3][0] = jCDAB  <9,9,7,1>;
        X32acbd [4][4][3][0] = xACBD  <9,9,7,1>;
        X32adbc [4][4][3][0] = xADBC  <9,9,7,1>;
        X32bdac [4][4][3][0] = xBDAC  <9,9,7,1>;
        X32bcad [4][4][3][0] = xBCAD  <9,9,7,1>;
        J32abcd [4][4][3][1] = jABCD  <9,9,7,3>;
        J32cdab [4][4][3][1] = jCDAB  <9,9,7,3>;
        X32acbd [4][4][3][1] = xACBD  <9,9,7,3>;
        X32adbc [4][4][3][1] = xADBC  <9,9,7,3>;
        X32bdac [4][4][3][1] = xBDAC  <9,9,7,3>;
        X32bcad [4][4][3][1] = xBCAD  <9,9,7,3>;
        J32abcd [4][4][3][2] = jABCD  <9,9,7,5>;
        J32cdab [4][4][3][2] = jCDAB  <9,9,7,5>;
        X32acbd [4][4][3][2] = xACBD  <9,9,7,5>;
        X32adbc [4][4][3][2] = xADBC  <9,9,7,5>;
        X32bdac [4][4][3][2] = xBDAC  <9,9,7,5>;
        X32bcad [4][4][3][2] = xBCAD  <9,9,7,5>;
        J32abcd [4][4][3][3] = jABCD  <9,9,7,7>;
        J32cdab [4][4][3][3] = jCDAB  <9,9,7,7>;
        X32acbd [4][4][3][3] = xACBD  <9,9,7,7>;
        X32adbc [4][4][3][3] = xADBC  <9,9,7,7>;
        X32bdac [4][4][3][3] = xBDAC  <9,9,7,7>;
        X32bcad [4][4][3][3] = xBCAD  <9,9,7,7>;
        J32abcd [4][4][4][0] = jABCD  <9,9,9,1>;
        J32cdab [4][4][4][0] = jCDAB  <9,9,9,1>;
        X32acbd [4][4][4][0] = xACBD  <9,9,9,1>;
        X32adbc [4][4][4][0] = xADBC  <9,9,9,1>;
        X32bdac [4][4][4][0] = xBDAC  <9,9,9,1>;
        X32bcad [4][4][4][0] = xBCAD  <9,9,9,1>;
        J32abcd [4][4][4][1] = jABCD  <9,9,9,3>;
        J32cdab [4][4][4][1] = jCDAB  <9,9,9,3>;
        X32acbd [4][4][4][1] = xACBD  <9,9,9,3>;
        X32adbc [4][4][4][1] = xADBC  <9,9,9,3>;
        X32bdac [4][4][4][1] = xBDAC  <9,9,9,3>;
        X32bcad [4][4][4][1] = xBCAD  <9,9,9,3>;
        J32abcd [4][4][4][2] = jABCD  <9,9,9,5>;
        J32cdab [4][4][4][2] = jCDAB  <9,9,9,5>;
        X32acbd [4][4][4][2] = xACBD  <9,9,9,5>;
        X32adbc [4][4][4][2] = xADBC  <9,9,9,5>;
        X32bdac [4][4][4][2] = xBDAC  <9,9,9,5>;
        X32bcad [4][4][4][2] = xBCAD  <9,9,9,5>;
        J32abcd [4][4][4][3] = jABCD  <9,9,9,7>;
        J32cdab [4][4][4][3] = jCDAB  <9,9,9,7>;
        X32acbd [4][4][4][3] = xACBD  <9,9,9,7>;
        X32adbc [4][4][4][3] = xADBC  <9,9,9,7>;
        X32bdac [4][4][4][3] = xBDAC  <9,9,9,7>;
        X32bcad [4][4][4][3] = xBCAD  <9,9,9,7>;
        J32abcd [4][4][4][4] = jABCD  <9,9,9,9>;
        J32cdab [4][4][4][4] = jCDAB  <9,9,9,9>;
        X32acbd [4][4][4][4] = xACBD  <9,9,9,9>;
        X32adbc [4][4][4][4] = xADBC  <9,9,9,9>;
        X32bdac [4][4][4][4] = xBDAC  <9,9,9,9>;
        X32bcad [4][4][4][4] = xBCAD  <9,9,9,9>;
    }

    //JXC NS
    {
        J32abcc [0][0][0] = jABCC  <1,1,1>;
        J32ccab [0][0][0] = jCCAB  <1,1,1>;
        X32acbc [0][0][0] = xACBD  <1,1,1,1>;
        X32bcac [0][0][0] = xBDAC  <1,1,1,1>;
        J32abcc [0][0][1] = jABCC  <1,1,3>;
        J32ccab [0][0][1] = jCCAB  <1,1,3>;
        X32acbc [0][0][1] = xACBD  <1,1,3,3>;
        X32bcac [0][0][1] = xBDAC  <1,1,3,3>;
        J32abcc [0][0][2] = jABCC  <1,1,5>;
        J32ccab [0][0][2] = jCCAB  <1,1,5>;
        X32acbc [0][0][2] = xACBD  <1,1,5,5>;
        X32bcac [0][0][2] = xBDAC  <1,1,5,5>;
        J32abcc [0][0][3] = jABCC  <1,1,7>;
        J32ccab [0][0][3] = jCCAB  <1,1,7>;
        X32acbc [0][0][3] = xACBD  <1,1,7,7>;
        X32bcac [0][0][3] = xBDAC  <1,1,7,7>;
        J32abcc [0][0][4] = jABCC  <1,1,9>;
        J32ccab [0][0][4] = jCCAB  <1,1,9>;
        X32acbc [0][0][4] = xACBD  <1,1,9,9>;
        X32bcac [0][0][4] = xBDAC  <1,1,9,9>;

        J32abcc [1][0][0] = jABCC  <3,1,1>;
        J32ccab [1][0][0] = jCCAB  <3,1,1>;
        X32acbc [1][0][0] = xACBD  <3,1,1,1>;
        X32bcac [1][0][0] = xBDAC  <3,1,1,1>;
        J32abcc [1][0][1] = jABCC  <3,1,3>;
        J32ccab [1][0][1] = jCCAB  <3,1,3>;
        X32acbc [1][0][1] = xACBD  <3,1,3,3>;
        X32bcac [1][0][1] = xBDAC  <3,1,3,3>;
        J32abcc [1][0][2] = jABCC  <3,1,5>;
        J32ccab [1][0][2] = jCCAB  <3,1,5>;
        X32acbc [1][0][2] = xACBD  <3,1,5,5>;
        X32bcac [1][0][2] = xBDAC  <3,1,5,5>;
        J32abcc [1][0][3] = jABCC  <3,1,7>;
        J32ccab [1][0][3] = jCCAB  <3,1,7>;
        X32acbc [1][0][3] = xACBD  <3,1,7,7>;
        X32bcac [1][0][3] = xBDAC  <3,1,7,7>;
        J32abcc [1][0][4] = jABCC  <3,1,9>;
        J32ccab [1][0][4] = jCCAB  <3,1,9>;
        X32acbc [1][0][4] = xACBD  <3,1,9,9>;
        X32bcac [1][0][4] = xBDAC  <3,1,9,9>;

        J32abcc [1][1][0] = jABCC  <3,3,1>;
        J32ccab [1][1][0] = jCCAB  <3,3,1>;
        X32acbc [1][1][0] = xACBD  <3,3,1,1>;
        X32bcac [1][1][0] = xBDAC  <3,3,1,1>;
        J32abcc [1][1][1] = jABCC  <3,3,3>;
        J32ccab [1][1][1] = jCCAB  <3,3,3>;
        X32acbc [1][1][1] = xACBD  <3,3,3,3>;
        X32bcac [1][1][1] = xBDAC  <3,3,3,3>;
        J32abcc [1][1][2] = jABCC  <3,3,5>;
        J32ccab [1][1][2] = jCCAB  <3,3,5>;
        X32acbc [1][1][2] = xACBD  <3,3,5,5>;
        X32bcac [1][1][2] = xBDAC  <3,3,5,5>;
        J32abcc [1][1][3] = jABCC  <3,3,7>;
        J32ccab [1][1][3] = jCCAB  <3,3,7>;
        X32acbc [1][1][3] = xACBD  <3,3,7,7>;
        X32bcac [1][1][3] = xBDAC  <3,3,7,7>;
        J32abcc [1][1][4] = jABCC  <3,3,9>;
        J32ccab [1][1][4] = jCCAB  <3,3,9>;
        X32acbc [1][1][4] = xACBD  <3,3,9,9>;
        X32bcac [1][1][4] = xBDAC  <3,3,9,9>;

        J32abcc [2][0][0] = jABCC  <5,1,1>;
        J32ccab [2][0][0] = jCCAB  <5,1,1>;
        X32acbc [2][0][0] = xACBD  <5,1,1,1>;
        X32bcac [2][0][0] = xBDAC  <5,1,1,1>;
        J32abcc [2][0][1] = jABCC  <5,1,3>;
        J32ccab [2][0][1] = jCCAB  <5,1,3>;
        X32acbc [2][0][1] = xACBD  <5,1,3,3>;
        X32bcac [2][0][1] = xBDAC  <5,1,3,3>;
        J32abcc [2][0][2] = jABCC  <5,1,5>;
        J32ccab [2][0][2] = jCCAB  <5,1,5>;
        X32acbc [2][0][2] = xACBD  <5,1,5,5>;
        X32bcac [2][0][2] = xBDAC  <5,1,5,5>;
        J32abcc [2][0][3] = jABCC  <5,1,7>;
        J32ccab [2][0][3] = jCCAB  <5,1,7>;
        X32acbc [2][0][3] = xACBD  <5,1,7,7>;
        X32bcac [2][0][3] = xBDAC  <5,1,7,7>;
        J32abcc [2][0][4] = jABCC  <5,1,9>;
        J32ccab [2][0][4] = jCCAB  <5,1,9>;
        X32acbc [2][0][4] = xACBD  <5,1,9,9>;
        X32bcac [2][0][4] = xBDAC  <5,1,9,9>;

        J32abcc [2][1][0] = jABCC  <5,3,1>;
        J32ccab [2][1][0] = jCCAB  <5,3,1>;
        X32acbc [2][1][0] = xACBD  <5,3,1,1>;
        X32bcac [2][1][0] = xBDAC  <5,3,1,1>;
        J32abcc [2][1][1] = jABCC  <5,3,3>;
        J32ccab [2][1][1] = jCCAB  <5,3,3>;
        X32acbc [2][1][1] = xACBD  <5,3,3,3>;
        X32bcac [2][1][1] = xBDAC  <5,3,3,3>;
        J32abcc [2][1][2] = jABCC  <5,3,5>;
        J32ccab [2][1][2] = jCCAB  <5,3,5>;
        X32acbc [2][1][2] = xACBD  <5,3,5,5>;
        X32bcac [2][1][2] = xBDAC  <5,3,5,5>;
        J32abcc [2][1][3] = jABCC  <5,3,7>;
        J32ccab [2][1][3] = jCCAB  <5,3,7>;
        X32acbc [2][1][3] = xACBD  <5,3,7,7>;
        X32bcac [2][1][3] = xBDAC  <5,3,7,7>;
        J32abcc [2][1][4] = jABCC  <5,3,9>;
        J32ccab [2][1][4] = jCCAB  <5,3,9>;
        X32acbc [2][1][4] = xACBD  <5,3,9,9>;
        X32bcac [2][1][4] = xBDAC  <5,3,9,9>;

        J32abcc [2][2][0] = jABCC  <5,5,1>;
        J32ccab [2][2][0] = jCCAB  <5,5,1>;
        X32acbc [2][2][0] = xACBD  <5,5,1,1>;
        X32bcac [2][2][0] = xBDAC  <5,5,1,1>;
        J32abcc [2][2][1] = jABCC  <5,5,3>;
        J32ccab [2][2][1] = jCCAB  <5,5,3>;
        X32acbc [2][2][1] = xACBD  <5,5,3,3>;
        X32bcac [2][2][1] = xBDAC  <5,5,3,3>;
        J32abcc [2][2][2] = jABCC  <5,5,5>;
        J32ccab [2][2][2] = jCCAB  <5,5,5>;
        X32acbc [2][2][2] = xACBD  <5,5,5,5>;
        X32bcac [2][2][2] = xBDAC  <5,5,5,5>;
        J32abcc [2][2][3] = jABCC  <5,5,7>;
        J32ccab [2][2][3] = jCCAB  <5,5,7>;
        X32acbc [2][2][3] = xACBD  <5,5,7,7>;
        X32bcac [2][2][3] = xBDAC  <5,5,7,7>;
        J32abcc [2][2][4] = jABCC  <5,5,9>;
        J32ccab [2][2][4] = jCCAB  <5,5,9>;
        X32acbc [2][2][4] = xACBD  <5,5,9,9>;
        X32bcac [2][2][4] = xBDAC  <5,5,9,9>;

        J32abcc [3][0][0] = jABCC  <7,1,1>;
        J32ccab [3][0][0] = jCCAB  <7,1,1>;
        X32acbc [3][0][0] = xACBD  <7,1,1,1>;
        X32bcac [3][0][0] = xBDAC  <7,1,1,1>;
        J32abcc [3][0][1] = jABCC  <7,1,3>;
        J32ccab [3][0][1] = jCCAB  <7,1,3>;
        X32acbc [3][0][1] = xACBD  <7,1,3,3>;
        X32bcac [3][0][1] = xBDAC  <7,1,3,3>;
        J32abcc [3][0][2] = jABCC  <7,1,5>;
        J32ccab [3][0][2] = jCCAB  <7,1,5>;
        X32acbc [3][0][2] = xACBD  <7,1,5,5>;
        X32bcac [3][0][2] = xBDAC  <7,1,5,5>;
        J32abcc [3][0][3] = jABCC  <7,1,7>;
        J32ccab [3][0][3] = jCCAB  <7,1,7>;
        X32acbc [3][0][3] = xACBD  <7,1,7,7>;
        X32bcac [3][0][3] = xBDAC  <7,1,7,7>;
        J32abcc [3][0][4] = jABCC  <7,1,9>;
        J32ccab [3][0][4] = jCCAB  <7,1,9>;
        X32acbc [3][0][4] = xACBD  <7,1,9,9>;
        X32bcac [3][0][4] = xBDAC  <7,1,9,9>;

        J32abcc [3][1][0] = jABCC  <7,3,1>;
        J32ccab [3][1][0] = jCCAB  <7,3,1>;
        X32acbc [3][1][0] = xACBD  <7,3,1,1>;
        X32bcac [3][1][0] = xBDAC  <7,3,1,1>;
        J32abcc [3][1][1] = jABCC  <7,3,3>;
        J32ccab [3][1][1] = jCCAB  <7,3,3>;
        X32acbc [3][1][1] = xACBD  <7,3,3,3>;
        X32bcac [3][1][1] = xBDAC  <7,3,3,3>;
        J32abcc [3][1][2] = jABCC  <7,3,5>;
        J32ccab [3][1][2] = jCCAB  <7,3,5>;
        X32acbc [3][1][2] = xACBD  <7,3,5,5>;
        X32bcac [3][1][2] = xBDAC  <7,3,5,5>;
        J32abcc [3][1][3] = jABCC  <7,3,7>;
        J32ccab [3][1][3] = jCCAB  <7,3,7>;
        X32acbc [3][1][3] = xACBD  <7,3,7,7>;
        X32bcac [3][1][3] = xBDAC  <7,3,7,7>;
        J32abcc [3][1][4] = jABCC  <7,3,9>;
        J32ccab [3][1][4] = jCCAB  <7,3,9>;
        X32acbc [3][1][4] = xACBD  <7,3,9,9>;
        X32bcac [3][1][4] = xBDAC  <7,3,9,9>;

        J32abcc [3][2][0] = jABCC  <7,5,1>;
        J32ccab [3][2][0] = jCCAB  <7,5,1>;
        X32acbc [3][2][0] = xACBD  <7,5,1,1>;
        X32bcac [3][2][0] = xBDAC  <7,5,1,1>;
        J32abcc [3][2][1] = jABCC  <7,5,3>;
        J32ccab [3][2][1] = jCCAB  <7,5,3>;
        X32acbc [3][2][1] = xACBD  <7,5,3,3>;
        X32bcac [3][2][1] = xBDAC  <7,5,3,3>;
        J32abcc [3][2][2] = jABCC  <7,5,5>;
        J32ccab [3][2][2] = jCCAB  <7,5,5>;
        X32acbc [3][2][2] = xACBD  <7,5,5,5>;
        X32bcac [3][2][2] = xBDAC  <7,5,5,5>;
        J32abcc [3][2][3] = jABCC  <7,5,7>;
        J32ccab [3][2][3] = jCCAB  <7,5,7>;
        X32acbc [3][2][3] = xACBD  <7,5,7,7>;
        X32bcac [3][2][3] = xBDAC  <7,5,7,7>;
        J32abcc [3][2][4] = jABCC  <7,5,9>;
        J32ccab [3][2][4] = jCCAB  <7,5,9>;
        X32acbc [3][2][4] = xACBD  <7,5,9,9>;
        X32bcac [3][2][4] = xBDAC  <7,5,9,9>;

        J32abcc [3][3][0] = jABCC  <7,7,1>;
        J32ccab [3][3][0] = jCCAB  <7,7,1>;
        X32acbc [3][3][0] = xACBD  <7,7,1,1>;
        X32bcac [3][3][0] = xBDAC  <7,7,1,1>;
        J32abcc [3][3][1] = jABCC  <7,7,3>;
        J32ccab [3][3][1] = jCCAB  <7,7,3>;
        X32acbc [3][3][1] = xACBD  <7,7,3,3>;
        X32bcac [3][3][1] = xBDAC  <7,7,3,3>;
        J32abcc [3][3][2] = jABCC  <7,7,5>;
        J32ccab [3][3][2] = jCCAB  <7,7,5>;
        X32acbc [3][3][2] = xACBD  <7,7,5,5>;
        X32bcac [3][3][2] = xBDAC  <7,7,5,5>;
        J32abcc [3][3][3] = jABCC  <7,7,7>;
        J32ccab [3][3][3] = jCCAB  <7,7,7>;
        X32acbc [3][3][3] = xACBD  <7,7,7,7>;
        X32bcac [3][3][3] = xBDAC  <7,7,7,7>;
        J32abcc [3][3][4] = jABCC  <7,7,9>;
        J32ccab [3][3][4] = jCCAB  <7,7,9>;
        X32acbc [3][3][4] = xACBD  <7,7,9,9>;
        X32bcac [3][3][4] = xBDAC  <7,7,9,9>;

        J32abcc [4][0][0] = jABCC  <9,1,1>;
        J32ccab [4][0][0] = jCCAB  <9,1,1>;
        X32acbc [4][0][0] = xACBD  <9,1,1,1>;
        X32bcac [4][0][0] = xBDAC  <9,1,1,1>;
        J32abcc [4][0][1] = jABCC  <9,1,3>;
        J32ccab [4][0][1] = jCCAB  <9,1,3>;
        X32acbc [4][0][1] = xACBD  <9,1,3,3>;
        X32bcac [4][0][1] = xBDAC  <9,1,3,3>;
        J32abcc [4][0][2] = jABCC  <9,1,5>;
        J32ccab [4][0][2] = jCCAB  <9,1,5>;
        X32acbc [4][0][2] = xACBD  <9,1,5,5>;
        X32bcac [4][0][2] = xBDAC  <9,1,5,5>;
        J32abcc [4][0][3] = jABCC  <9,1,7>;
        J32ccab [4][0][3] = jCCAB  <9,1,7>;
        X32acbc [4][0][3] = xACBD  <9,1,7,7>;
        X32bcac [4][0][3] = xBDAC  <9,1,7,7>;
        J32abcc [4][0][4] = jABCC  <9,1,9>;
        J32ccab [4][0][4] = jCCAB  <9,1,9>;
        X32acbc [4][0][4] = xACBD  <9,1,9,9>;
        X32bcac [4][0][4] = xBDAC  <9,1,9,9>;

        J32abcc [4][1][0] = jABCC  <9,3,1>;
        J32ccab [4][1][0] = jCCAB  <9,3,1>;
        X32acbc [4][1][0] = xACBD  <9,3,1,1>;
        X32bcac [4][1][0] = xBDAC  <9,3,1,1>;
        J32abcc [4][1][1] = jABCC  <9,3,3>;
        J32ccab [4][1][1] = jCCAB  <9,3,3>;
        X32acbc [4][1][1] = xACBD  <9,3,3,3>;
        X32bcac [4][1][1] = xBDAC  <9,3,3,3>;
        J32abcc [4][1][2] = jABCC  <9,3,5>;
        J32ccab [4][1][2] = jCCAB  <9,3,5>;
        X32acbc [4][1][2] = xACBD  <9,3,5,5>;
        X32bcac [4][1][2] = xBDAC  <9,3,5,5>;
        J32abcc [4][1][3] = jABCC  <9,3,7>;
        J32ccab [4][1][3] = jCCAB  <9,3,7>;
        X32acbc [4][1][3] = xACBD  <9,3,7,7>;
        X32bcac [4][1][3] = xBDAC  <9,3,7,7>;
        J32abcc [4][1][4] = jABCC  <9,3,9>;
        J32ccab [4][1][4] = jCCAB  <9,3,9>;
        X32acbc [4][1][4] = xACBD  <9,3,9,9>;
        X32bcac [4][1][4] = xBDAC  <9,3,9,9>;

        J32abcc [4][2][0] = jABCC  <9,5,1>;
        J32ccab [4][2][0] = jCCAB  <9,5,1>;
        X32acbc [4][2][0] = xACBD  <9,5,1,1>;
        X32bcac [4][2][0] = xBDAC  <9,5,1,1>;
        J32abcc [4][2][1] = jABCC  <9,5,3>;
        J32ccab [4][2][1] = jCCAB  <9,5,3>;
        X32acbc [4][2][1] = xACBD  <9,5,3,3>;
        X32bcac [4][2][1] = xBDAC  <9,5,3,3>;
        J32abcc [4][2][2] = jABCC  <9,5,5>;
        J32ccab [4][2][2] = jCCAB  <9,5,5>;
        X32acbc [4][2][2] = xACBD  <9,5,5,5>;
        X32bcac [4][2][2] = xBDAC  <9,5,5,5>;
        J32abcc [4][2][3] = jABCC  <9,5,7>;
        J32ccab [4][2][3] = jCCAB  <9,5,7>;
        X32acbc [4][2][3] = xACBD  <9,5,7,7>;
        X32bcac [4][2][3] = xBDAC  <9,5,7,7>;
        J32abcc [4][2][4] = jABCC  <9,5,9>;
        J32ccab [4][2][4] = jCCAB  <9,5,9>;
        X32acbc [4][2][4] = xACBD  <9,5,9,9>;
        X32bcac [4][2][4] = xBDAC  <9,5,9,9>;

        J32abcc [4][3][0] = jABCC  <9,7,1>;
        J32ccab [4][3][0] = jCCAB  <9,7,1>;
        X32acbc [4][3][0] = xACBD  <9,7,1,1>;
        X32bcac [4][3][0] = xBDAC  <9,7,1,1>;
        J32abcc [4][3][1] = jABCC  <9,7,3>;
        J32ccab [4][3][1] = jCCAB  <9,7,3>;
        X32acbc [4][3][1] = xACBD  <9,7,3,3>;
        X32bcac [4][3][1] = xBDAC  <9,7,3,3>;
        J32abcc [4][3][2] = jABCC  <9,7,5>;
        J32ccab [4][3][2] = jCCAB  <9,7,5>;
        X32acbc [4][3][2] = xACBD  <9,7,5,5>;
        X32bcac [4][3][2] = xBDAC  <9,7,5,5>;
        J32abcc [4][3][3] = jABCC  <9,7,7>;
        J32ccab [4][3][3] = jCCAB  <9,7,7>;
        X32acbc [4][3][3] = xACBD  <9,7,7,7>;
        X32bcac [4][3][3] = xBDAC  <9,7,7,7>;
        J32abcc [4][3][4] = jABCC  <9,7,9>;
        J32ccab [4][3][4] = jCCAB  <9,7,9>;
        X32acbc [4][3][4] = xACBD  <9,7,9,9>;
        X32bcac [4][3][4] = xBDAC  <9,7,9,9>;
        J32abcc [4][4][0] = jABCC  <9,9,1>;
        J32ccab [4][4][0] = jCCAB  <9,9,1>;
        X32acbc [4][4][0] = xACBD  <9,9,1,1>;
        X32bcac [4][4][0] = xBDAC  <9,9,1,1>;
        J32abcc [4][4][1] = jABCC  <9,9,3>;
        J32ccab [4][4][1] = jCCAB  <9,9,3>;
        X32acbc [4][4][1] = xACBD  <9,9,3,3>;
        X32bcac [4][4][1] = xBDAC  <9,9,3,3>;
        J32abcc [4][4][2] = jABCC  <9,9,5>;
        J32ccab [4][4][2] = jCCAB  <9,9,5>;
        X32acbc [4][4][2] = xACBD  <9,9,5,5>;
        X32bcac [4][4][2] = xBDAC  <9,9,5,5>;
        J32abcc [4][4][3] = jABCC  <9,9,7>;
        J32ccab [4][4][3] = jCCAB  <9,9,7>;
        X32acbc [4][4][3] = xACBD  <9,9,7,7>;
        X32bcac [4][4][3] = xBDAC  <9,9,7,7>;
        J32abcc [4][4][4] = jABCC  <9,9,9>;
        J32ccab [4][4][4] = jCCAB  <9,9,9>;
        X32acbc [4][4][4] = xACBD  <9,9,9,9>;
        X32bcac [4][4][4] = xBDAC  <9,9,9,9>;
    }

    //JXC SN
    {
        J32aacd [0][0][0] = jAACD  <1,1,1>;
        J32cdaa [0][0][0] = jCDAA  <1,1,1>;
        X32acad [0][0][0] = xACBD  <1,1,1,1>;
        X32adac [0][0][0] = xADBC  <1,1,1,1>;
        J32aacd [0][1][0] = jAACD  <1,3,1>;
        J32cdaa [0][1][0] = jCDAA  <1,3,1>;
        X32acad [0][1][0] = xACBD  <1,1,3,1>;
        X32adac [0][1][0] = xADBC  <1,1,3,1>;
        J32aacd [0][1][1] = jAACD  <1,3,3>;
        J32cdaa [0][1][1] = jCDAA  <1,3,3>;
        X32acad [0][1][1] = xACBD  <1,1,3,3>;
        X32adac [0][1][1] = xADBC  <1,1,3,3>;
        J32aacd [0][2][0] = jAACD  <1,5,1>;
        J32cdaa [0][2][0] = jCDAA  <1,5,1>;
        X32acad [0][2][0] = xACBD  <1,1,5,1>;
        X32adac [0][2][0] = xADBC  <1,1,5,1>;
        J32aacd [0][2][1] = jAACD  <1,5,3>;
        J32cdaa [0][2][1] = jCDAA  <1,5,3>;
        X32acad [0][2][1] = xACBD  <1,1,5,3>;
        X32adac [0][2][1] = xADBC  <1,1,5,3>;
        J32aacd [0][2][2] = jAACD  <1,5,5>;
        J32cdaa [0][2][2] = jCDAA  <1,5,5>;
        X32acad [0][2][2] = xACBD  <1,1,5,5>;
        X32adac [0][2][2] = xADBC  <1,1,5,5>;
        J32aacd [0][3][0] = jAACD  <1,7,1>;
        J32cdaa [0][3][0] = jCDAA  <1,7,1>;
        X32acad [0][3][0] = xACBD  <1,1,7,1>;
        X32adac [0][3][0] = xADBC  <1,1,7,1>;
        J32aacd [0][3][1] = jAACD  <1,7,3>;
        J32cdaa [0][3][1] = jCDAA  <1,7,3>;
        X32acad [0][3][1] = xACBD  <1,1,7,3>;
        X32adac [0][3][1] = xADBC  <1,1,7,3>;
        J32aacd [0][3][2] = jAACD  <1,7,5>;
        J32cdaa [0][3][2] = jCDAA  <1,7,5>;
        X32acad [0][3][2] = xACBD  <1,1,7,5>;
        X32adac [0][3][2] = xADBC  <1,1,7,5>;
        J32aacd [0][3][3] = jAACD  <1,7,7>;
        J32cdaa [0][3][3] = jCDAA  <1,7,7>;
        X32acad [0][3][3] = xACBD  <1,1,7,7>;
        X32adac [0][3][3] = xADBC  <1,1,7,7>;
        J32aacd [0][4][0] = jAACD  <1,9,1>;
        J32cdaa [0][4][0] = jCDAA  <1,9,1>;
        X32acad [0][4][0] = xACBD  <1,1,9,1>;
        X32adac [0][4][0] = xADBC  <1,1,9,1>;
        J32aacd [0][4][1] = jAACD  <1,9,3>;
        J32cdaa [0][4][1] = jCDAA  <1,9,3>;
        X32acad [0][4][1] = xACBD  <1,1,9,3>;
        X32adac [0][4][1] = xADBC  <1,1,9,3>;
        J32aacd [0][4][2] = jAACD  <1,9,5>;
        J32cdaa [0][4][2] = jCDAA  <1,9,5>;
        X32acad [0][4][2] = xACBD  <1,1,9,5>;
        X32adac [0][4][2] = xADBC  <1,1,9,5>;
        J32aacd [0][4][3] = jAACD  <1,9,7>;
        J32cdaa [0][4][3] = jCDAA  <1,9,7>;
        X32acad [0][4][3] = xACBD  <1,1,9,7>;
        X32adac [0][4][3] = xADBC  <1,1,9,7>;
        J32aacd [0][4][4] = jAACD  <1,9,9>;
        J32cdaa [0][4][4] = jCDAA  <1,9,9>;
        X32acad [0][4][4] = xACBD  <1,1,9,9>;
        X32adac [0][4][4] = xADBC  <1,1,9,9>;
        J32aacd [1][0][0] = jAACD  <3,1,1>;
        J32cdaa [1][0][0] = jCDAA  <3,1,1>;
        X32acad [1][0][0] = xACBD  <3,3,1,1>;
        X32adac [1][0][0] = xADBC  <3,3,1,1>;
        J32aacd [1][1][0] = jAACD  <3,3,1>;
        J32cdaa [1][1][0] = jCDAA  <3,3,1>;
        X32acad [1][1][0] = xACBD  <3,3,3,1>;
        X32adac [1][1][0] = xADBC  <3,3,3,1>;
        J32aacd [1][1][1] = jAACD  <3,3,3>;
        J32cdaa [1][1][1] = jCDAA  <3,3,3>;
        X32acad [1][1][1] = xACBD  <3,3,3,3>;
        X32adac [1][1][1] = xADBC  <3,3,3,3>;
        J32aacd [1][2][0] = jAACD  <3,5,1>;
        J32cdaa [1][2][0] = jCDAA  <3,5,1>;
        X32acad [1][2][0] = xACBD  <3,3,5,1>;
        X32adac [1][2][0] = xADBC  <3,3,5,1>;
        J32aacd [1][2][1] = jAACD  <3,5,3>;
        J32cdaa [1][2][1] = jCDAA  <3,5,3>;
        X32acad [1][2][1] = xACBD  <3,3,5,3>;
        X32adac [1][2][1] = xADBC  <3,3,5,3>;
        J32aacd [1][2][2] = jAACD  <3,5,5>;
        J32cdaa [1][2][2] = jCDAA  <3,5,5>;
        X32acad [1][2][2] = xACBD  <3,3,5,5>;
        X32adac [1][2][2] = xADBC  <3,3,5,5>;
        J32aacd [1][3][0] = jAACD  <3,7,1>;
        J32cdaa [1][3][0] = jCDAA  <3,7,1>;
        X32acad [1][3][0] = xACBD  <3,3,7,1>;
        X32adac [1][3][0] = xADBC  <3,3,7,1>;
        J32aacd [1][3][1] = jAACD  <3,7,3>;
        J32cdaa [1][3][1] = jCDAA  <3,7,3>;
        X32acad [1][3][1] = xACBD  <3,3,7,3>;
        X32adac [1][3][1] = xADBC  <3,3,7,3>;
        J32aacd [1][3][2] = jAACD  <3,7,5>;
        J32cdaa [1][3][2] = jCDAA  <3,7,5>;
        X32acad [1][3][2] = xACBD  <3,3,7,5>;
        X32adac [1][3][2] = xADBC  <3,3,7,5>;
        J32aacd [1][3][3] = jAACD  <3,7,7>;
        J32cdaa [1][3][3] = jCDAA  <3,7,7>;
        X32acad [1][3][3] = xACBD  <3,3,7,7>;
        X32adac [1][3][3] = xADBC  <3,3,7,7>;
        J32aacd [1][4][0] = jAACD  <3,9,1>;
        J32cdaa [1][4][0] = jCDAA  <3,9,1>;
        X32acad [1][4][0] = xACBD  <3,3,9,1>;
        X32adac [1][4][0] = xADBC  <3,3,9,1>;
        J32aacd [1][4][1] = jAACD  <3,9,3>;
        J32cdaa [1][4][1] = jCDAA  <3,9,3>;
        X32acad [1][4][1] = xACBD  <3,3,9,3>;
        X32adac [1][4][1] = xADBC  <3,3,9,3>;
        J32aacd [1][4][2] = jAACD  <3,9,5>;
        J32cdaa [1][4][2] = jCDAA  <3,9,5>;
        X32acad [1][4][2] = xACBD  <3,3,9,5>;
        X32adac [1][4][2] = xADBC  <3,3,9,5>;
        J32aacd [1][4][3] = jAACD  <3,9,7>;
        J32cdaa [1][4][3] = jCDAA  <3,9,7>;
        X32acad [1][4][3] = xACBD  <3,3,9,7>;
        X32adac [1][4][3] = xADBC  <3,3,9,7>;
        J32aacd [1][4][4] = jAACD  <3,9,9>;
        J32cdaa [1][4][4] = jCDAA  <3,9,9>;
        X32acad [1][4][4] = xACBD  <3,3,9,9>;
        X32adac [1][4][4] = xADBC  <3,3,9,9>;
        J32aacd [2][0][0] = jAACD  <5,1,1>;
        J32cdaa [2][0][0] = jCDAA  <5,1,1>;
        X32acad [2][0][0] = xACBD  <5,5,1,1>;
        X32adac [2][0][0] = xADBC  <5,5,1,1>;
        J32aacd [2][1][0] = jAACD  <5,3,1>;
        J32cdaa [2][1][0] = jCDAA  <5,3,1>;
        X32acad [2][1][0] = xACBD  <5,5,3,1>;
        X32adac [2][1][0] = xADBC  <5,5,3,1>;
        J32aacd [2][1][1] = jAACD  <5,3,3>;
        J32cdaa [2][1][1] = jCDAA  <5,3,3>;
        X32acad [2][1][1] = xACBD  <5,5,3,3>;
        X32adac [2][1][1] = xADBC  <5,5,3,3>;
        J32aacd [2][2][0] = jAACD  <5,5,1>;
        J32cdaa [2][2][0] = jCDAA  <5,5,1>;
        X32acad [2][2][0] = xACBD  <5,5,5,1>;
        X32adac [2][2][0] = xADBC  <5,5,5,1>;
        J32aacd [2][2][1] = jAACD  <5,5,3>;
        J32cdaa [2][2][1] = jCDAA  <5,5,3>;
        X32acad [2][2][1] = xACBD  <5,5,5,3>;
        X32adac [2][2][1] = xADBC  <5,5,5,3>;
        J32aacd [2][2][2] = jAACD  <5,5,5>;
        J32cdaa [2][2][2] = jCDAA  <5,5,5>;
        X32acad [2][2][2] = xACBD  <5,5,5,5>;
        X32adac [2][2][2] = xADBC  <5,5,5,5>;
        J32aacd [2][3][0] = jAACD  <5,7,1>;
        J32cdaa [2][3][0] = jCDAA  <5,7,1>;
        X32acad [2][3][0] = xACBD  <5,5,7,1>;
        X32adac [2][3][0] = xADBC  <5,5,7,1>;
        J32aacd [2][3][1] = jAACD  <5,7,3>;
        J32cdaa [2][3][1] = jCDAA  <5,7,3>;
        X32acad [2][3][1] = xACBD  <5,5,7,3>;
        X32adac [2][3][1] = xADBC  <5,5,7,3>;
        J32aacd [2][3][2] = jAACD  <5,7,5>;
        J32cdaa [2][3][2] = jCDAA  <5,7,5>;
        X32acad [2][3][2] = xACBD  <5,5,7,5>;
        X32adac [2][3][2] = xADBC  <5,5,7,5>;
        J32aacd [2][3][3] = jAACD  <5,7,7>;
        J32cdaa [2][3][3] = jCDAA  <5,7,7>;
        X32acad [2][3][3] = xACBD  <5,5,7,7>;
        X32adac [2][3][3] = xADBC  <5,5,7,7>;
        J32aacd [2][4][0] = jAACD  <5,9,1>;
        J32cdaa [2][4][0] = jCDAA  <5,9,1>;
        X32acad [2][4][0] = xACBD  <5,5,9,1>;
        X32adac [2][4][0] = xADBC  <5,5,9,1>;
        J32aacd [2][4][1] = jAACD  <5,9,3>;
        J32cdaa [2][4][1] = jCDAA  <5,9,3>;
        X32acad [2][4][1] = xACBD  <5,5,9,3>;
        X32adac [2][4][1] = xADBC  <5,5,9,3>;
        J32aacd [2][4][2] = jAACD  <5,9,5>;
        J32cdaa [2][4][2] = jCDAA  <5,9,5>;
        X32acad [2][4][2] = xACBD  <5,5,9,5>;
        X32adac [2][4][2] = xADBC  <5,5,9,5>;
        J32aacd [2][4][3] = jAACD  <5,9,7>;
        J32cdaa [2][4][3] = jCDAA  <5,9,7>;
        X32acad [2][4][3] = xACBD  <5,5,9,7>;
        X32adac [2][4][3] = xADBC  <5,5,9,7>;
        J32aacd [2][4][4] = jAACD  <5,9,9>;
        J32cdaa [2][4][4] = jCDAA  <5,9,9>;
        X32acad [2][4][4] = xACBD  <5,5,9,9>;
        X32adac [2][4][4] = xADBC  <5,5,9,9>;
        J32aacd [3][0][0] = jAACD  <7,1,1>;
        J32cdaa [3][0][0] = jCDAA  <7,1,1>;
        X32acad [3][0][0] = xACBD  <7,7,1,1>;
        X32adac [3][0][0] = xADBC  <7,7,1,1>;
        J32aacd [3][1][0] = jAACD  <7,3,1>;
        J32cdaa [3][1][0] = jCDAA  <7,3,1>;
        X32acad [3][1][0] = xACBD  <7,7,3,1>;
        X32adac [3][1][0] = xADBC  <7,7,3,1>;
        J32aacd [3][1][1] = jAACD  <7,3,3>;
        J32cdaa [3][1][1] = jCDAA  <7,3,3>;
        X32acad [3][1][1] = xACBD  <7,7,3,3>;
        X32adac [3][1][1] = xADBC  <7,7,3,3>;
        J32aacd [3][2][0] = jAACD  <7,5,1>;
        J32cdaa [3][2][0] = jCDAA  <7,5,1>;
        X32acad [3][2][0] = xACBD  <7,7,5,1>;
        X32adac [3][2][0] = xADBC  <7,7,5,1>;
        J32aacd [3][2][1] = jAACD  <7,5,3>;
        J32cdaa [3][2][1] = jCDAA  <7,5,3>;
        X32acad [3][2][1] = xACBD  <7,7,5,3>;
        X32adac [3][2][1] = xADBC  <7,7,5,3>;
        J32aacd [3][2][2] = jAACD  <7,5,5>;
        J32cdaa [3][2][2] = jCDAA  <7,5,5>;
        X32acad [3][2][2] = xACBD  <7,7,5,5>;
        X32adac [3][2][2] = xADBC  <7,7,5,5>;
        J32aacd [3][3][0] = jAACD  <7,7,1>;
        J32cdaa [3][3][0] = jCDAA  <7,7,1>;
        X32acad [3][3][0] = xACBD  <7,7,7,1>;
        X32adac [3][3][0] = xADBC  <7,7,7,1>;
        J32aacd [3][3][1] = jAACD  <7,7,3>;
        J32cdaa [3][3][1] = jCDAA  <7,7,3>;
        X32acad [3][3][1] = xACBD  <7,7,7,3>;
        X32adac [3][3][1] = xADBC  <7,7,7,3>;
        J32aacd [3][3][2] = jAACD  <7,7,5>;
        J32cdaa [3][3][2] = jCDAA  <7,7,5>;
        X32acad [3][3][2] = xACBD  <7,7,7,5>;
        X32adac [3][3][2] = xADBC  <7,7,7,5>;
        J32aacd [3][3][3] = jAACD  <7,7,7>;
        J32cdaa [3][3][3] = jCDAA  <7,7,7>;
        X32acad [3][3][3] = xACBD  <7,7,7,7>;
        X32adac [3][3][3] = xADBC  <7,7,7,7>;
        J32aacd [3][4][0] = jAACD  <7,9,1>;
        J32cdaa [3][4][0] = jCDAA  <7,9,1>;
        X32acad [3][4][0] = xACBD  <7,7,9,1>;
        X32adac [3][4][0] = xADBC  <7,7,9,1>;
        J32aacd [3][4][1] = jAACD  <7,9,3>;
        J32cdaa [3][4][1] = jCDAA  <7,9,3>;
        X32acad [3][4][1] = xACBD  <7,7,9,3>;
        X32adac [3][4][1] = xADBC  <7,7,9,3>;
        J32aacd [3][4][2] = jAACD  <7,9,5>;
        J32cdaa [3][4][2] = jCDAA  <7,9,5>;
        X32acad [3][4][2] = xACBD  <7,7,9,5>;
        X32adac [3][4][2] = xADBC  <7,7,9,5>;
        J32aacd [3][4][3] = jAACD  <7,9,7>;
        J32cdaa [3][4][3] = jCDAA  <7,9,7>;
        X32acad [3][4][3] = xACBD  <7,7,9,7>;
        X32adac [3][4][3] = xADBC  <7,7,9,7>;
        J32aacd [3][4][4] = jAACD  <7,9,9>;
        J32cdaa [3][4][4] = jCDAA  <7,9,9>;
        X32acad [3][4][4] = xACBD  <7,7,9,9>;
        X32adac [3][4][4] = xADBC  <7,7,9,9>;
        J32aacd [4][0][0] = jAACD  <9,1,1>;
        J32cdaa [4][0][0] = jCDAA  <9,1,1>;
        X32acad [4][0][0] = xACBD  <9,9,1,1>;
        X32adac [4][0][0] = xADBC  <9,9,1,1>;
        J32aacd [4][1][0] = jAACD  <9,3,1>;
        J32cdaa [4][1][0] = jCDAA  <9,3,1>;
        X32acad [4][1][0] = xACBD  <9,9,3,1>;
        X32adac [4][1][0] = xADBC  <9,9,3,1>;
        J32aacd [4][1][1] = jAACD  <9,3,3>;
        J32cdaa [4][1][1] = jCDAA  <9,3,3>;
        X32acad [4][1][1] = xACBD  <9,9,3,3>;
        X32adac [4][1][1] = xADBC  <9,9,3,3>;
        J32aacd [4][2][0] = jAACD  <9,5,1>;
        J32cdaa [4][2][0] = jCDAA  <9,5,1>;
        X32acad [4][2][0] = xACBD  <9,9,5,1>;
        X32adac [4][2][0] = xADBC  <9,9,5,1>;
        J32aacd [4][2][1] = jAACD  <9,5,3>;
        J32cdaa [4][2][1] = jCDAA  <9,5,3>;
        X32acad [4][2][1] = xACBD  <9,9,5,3>;
        X32adac [4][2][1] = xADBC  <9,9,5,3>;
        J32aacd [4][2][2] = jAACD  <9,5,5>;
        J32cdaa [4][2][2] = jCDAA  <9,5,5>;
        X32acad [4][2][2] = xACBD  <9,9,5,5>;
        X32adac [4][2][2] = xADBC  <9,9,5,5>;
        J32aacd [4][3][0] = jAACD  <9,7,1>;
        J32cdaa [4][3][0] = jCDAA  <9,7,1>;
        X32acad [4][3][0] = xACBD  <9,9,7,1>;
        X32adac [4][3][0] = xADBC  <9,9,7,1>;
        J32aacd [4][3][1] = jAACD  <9,7,3>;
        J32cdaa [4][3][1] = jCDAA  <9,7,3>;
        X32acad [4][3][1] = xACBD  <9,9,7,3>;
        X32adac [4][3][1] = xADBC  <9,9,7,3>;
        J32aacd [4][3][2] = jAACD  <9,7,5>;
        J32cdaa [4][3][2] = jCDAA  <9,7,5>;
        X32acad [4][3][2] = xACBD  <9,9,7,5>;
        X32adac [4][3][2] = xADBC  <9,9,7,5>;
        J32aacd [4][3][3] = jAACD  <9,7,7>;
        J32cdaa [4][3][3] = jCDAA  <9,7,7>;
        X32acad [4][3][3] = xACBD  <9,9,7,7>;
        X32adac [4][3][3] = xADBC  <9,9,7,7>;
        J32aacd [4][4][0] = jAACD  <9,9,1>;
        J32cdaa [4][4][0] = jCDAA  <9,9,1>;
        X32acad [4][4][0] = xACBD  <9,9,9,1>;
        X32adac [4][4][0] = xADBC  <9,9,9,1>;
        J32aacd [4][4][1] = jAACD  <9,9,3>;
        J32cdaa [4][4][1] = jCDAA  <9,9,3>;
        X32acad [4][4][1] = xACBD  <9,9,9,3>;
        X32adac [4][4][1] = xADBC  <9,9,9,3>;
        J32aacd [4][4][2] = jAACD  <9,9,5>;
        J32cdaa [4][4][2] = jCDAA  <9,9,5>;
        X32acad [4][4][2] = xACBD  <9,9,9,5>;
        X32adac [4][4][2] = xADBC  <9,9,9,5>;
        J32aacd [4][4][3] = jAACD  <9,9,7>;
        J32cdaa [4][4][3] = jCDAA  <9,9,7>;
        X32acad [4][4][3] = xACBD  <9,9,9,7>;
        X32adac [4][4][3] = xADBC  <9,9,9,7>;
        J32aacd [4][4][4] = jAACD  <9,9,9>;
        J32cdaa [4][4][4] = jCDAA  <9,9,9>;
        X32acad [4][4][4] = xACBD  <9,9,9,9>;
        X32adac [4][4][4] = xADBC  <9,9,9,9>;
    }

    //JXC SS
    {
        J32aacc [0][0] = jAACC  <1,1>;
        J32ccaa [0][0] = jCCAA  <1,1>;
        X32acac [0][0] = xACBD  <1,1,1,1>;
        J32aacc [0][1] = jAACC  <1,3>;
        J32ccaa [0][1] = jCCAA  <1,3>;
        X32acac [0][1] = xACBD  <1,1,3,3>;
        J32aacc [0][2] = jAACC  <1,5>;
        J32ccaa [0][2] = jCCAA  <1,5>;
        X32acac [0][2] = xACBD  <1,1,5,5>;
        J32aacc [0][3] = jAACC  <1,7>;
        J32ccaa [0][3] = jCCAA  <1,7>;
        X32acac [0][3] = xACBD  <1,1,7,7>;
        J32aacc [0][4] = jAACC  <1,9>;
        J32ccaa [0][4] = jCCAA  <1,9>;
        X32acac [0][4] = xACBD  <1,1,9,9>;
        J32aacc [1][0] = jAACC  <3,1>;
        J32ccaa [1][0] = jCCAA  <3,1>;
        X32acac [1][0] = xACBD  <3,3,1,1>;
        J32aacc [1][1] = jAACC  <3,3>;
        J32ccaa [1][1] = jCCAA  <3,3>;
        X32acac [1][1] = xACBD  <3,3,3,3>;
        J32aacc [1][2] = jAACC  <3,5>;
        J32ccaa [1][2] = jCCAA  <3,5>;
        X32acac [1][2] = xACBD  <3,3,5,5>;
        J32aacc [1][3] = jAACC  <3,7>;
        J32ccaa [1][3] = jCCAA  <3,7>;
        X32acac [1][3] = xACBD  <3,3,7,7>;
        J32aacc [1][4] = jAACC  <3,9>;
        J32ccaa [1][4] = jCCAA  <3,9>;
        X32acac [1][4] = xACBD  <3,3,9,9>;
        J32aacc [2][0] = jAACC  <5,1>;
        J32ccaa [2][0] = jCCAA  <5,1>;
        X32acac [2][0] = xACBD  <5,5,1,1>;
        J32aacc [2][1] = jAACC  <5,3>;
        J32ccaa [2][1] = jCCAA  <5,3>;
        X32acac [2][1] = xACBD  <5,5,3,3>;
        J32aacc [2][2] = jAACC  <5,5>;
        J32ccaa [2][2] = jCCAA  <5,5>;
        X32acac [2][2] = xACBD  <5,5,5,5>;
        J32aacc [2][3] = jAACC  <5,7>;
        J32ccaa [2][3] = jCCAA  <5,7>;
        X32acac [2][3] = xACBD  <5,5,7,7>;
        J32aacc [2][4] = jAACC  <5,9>;
        J32ccaa [2][4] = jCCAA  <5,9>;
        X32acac [2][4] = xACBD  <5,5,9,9>;
        J32aacc [3][0] = jAACC  <7,1>;
        J32ccaa [3][0] = jCCAA  <7,1>;
        X32acac [3][0] = xACBD  <7,7,1,1>;
        J32aacc [3][1] = jAACC  <7,3>;
        J32ccaa [3][1] = jCCAA  <7,3>;
        X32acac [3][1] = xACBD  <7,7,3,3>;
        J32aacc [3][2] = jAACC  <7,5>;
        J32ccaa [3][2] = jCCAA  <7,5>;
        X32acac [3][2] = xACBD  <7,7,5,5>;
        J32aacc [3][3] = jAACC  <7,7>;
        J32ccaa [3][3] = jCCAA  <7,7>;
        X32acac [3][3] = xACBD  <7,7,7,7>;
        J32aacc [3][4] = jAACC  <7,9>;
        J32ccaa [3][4] = jCCAA  <7,9>;
        X32acac [3][4] = xACBD  <7,7,9,9>;
        J32aacc [4][0] = jAACC  <9,1>;
        J32ccaa [4][0] = jCCAA  <9,1>;
        X32acac [4][0] = xACBD  <9,9,1,1>;
        J32aacc [4][1] = jAACC  <9,3>;
        J32ccaa [4][1] = jCCAA  <9,3>;
        X32acac [4][1] = xACBD  <9,9,3,3>;
        J32aacc [4][2] = jAACC  <9,5>;
        J32ccaa [4][2] = jCCAA  <9,5>;
        X32acac [4][2] = xACBD  <9,9,5,5>;
        J32aacc [4][3] = jAACC  <9,7>;
        J32ccaa [4][3] = jCCAA  <9,7>;
        X32acac [4][3] = xACBD  <9,9,7,7>;
        J32aacc [4][4] = jAACC  <9,9>;
        J32ccaa [4][4] = jCCAA  <9,9>;
        X32acac [4][4] = xACBD  <9,9,9,9>;
    }

    //JXC S
    {
        J32abab [0][0] = jABCD  <1,1,1,1>;
        X32aabb [0][0] = xAABB  <1,1>;
        X32bbaa [0][0] = xBBAA  <1,1>;
        X32abba [0][0] = xABAB  <1,1>;
        J32abab [1][0] = jABCD  <3,1,3,1>;
        X32aabb [1][0] = xAABB  <3,1>;
        X32bbaa [1][0] = xBBAA  <3,1>;
        X32abba [1][0] = xABAB  <3,1>;
        J32abab [1][1] = jABCD  <3,3,3,3>;
        X32aabb [1][1] = xAABB  <3,3>;
        X32bbaa [1][1] = xBBAA  <3,3>;
        X32abba [1][1] = xABAB  <3,3>;
        J32abab [2][0] = jABCD  <5,1,5,1>;
        X32aabb [2][0] = xAABB  <5,1>;
        X32bbaa [2][0] = xBBAA  <5,1>;
        X32abba [2][0] = xABAB  <5,1>;
        J32abab [2][1] = jABCD  <5,3,5,3>;
        X32aabb [2][1] = xAABB  <5,3>;
        X32bbaa [2][1] = xBBAA  <5,3>;
        X32abba [2][1] = xABAB  <5,3>;
        J32abab [2][2] = jABCD  <5,5,5,5>;
        X32aabb [2][2] = xAABB  <5,5>;
        X32bbaa [2][2] = xBBAA  <5,5>;
        X32abba [2][2] = xABAB  <5,5>;
        J32abab [3][0] = jABCD  <7,1,7,1>;
        X32aabb [3][0] = xAABB  <7,1>;
        X32bbaa [3][0] = xBBAA  <7,1>;
        X32abba [3][0] = xABAB  <7,1>;
        J32abab [3][1] = jABCD  <7,3,7,3>;
        X32aabb [3][1] = xAABB  <7,3>;
        X32bbaa [3][1] = xBBAA  <7,3>;
        X32abba [3][1] = xABAB  <7,3>;
        J32abab [3][2] = jABCD  <7,5,7,5>;
        X32aabb [3][2] = xAABB  <7,5>;
        X32bbaa [3][2] = xBBAA  <7,5>;
        X32abba [3][2] = xABAB  <7,5>;
        J32abab [3][3] = jABCD  <7,7,7,7>;
        X32aabb [3][3] = xAABB  <7,7>;
        X32bbaa [3][3] = xBBAA  <7,7>;
        X32abba [3][3] = xABAB  <7,7>;
        J32abab [4][0] = jABCD  <9,1,9,1>;
        X32aabb [4][0] = xAABB  <9,1>;
        X32bbaa [4][0] = xBBAA  <9,1>;
        X32abba [4][0] = xABAB  <9,1>;
        J32abab [4][1] = jABCD  <9,3,9,3>;
        X32aabb [4][1] = xAABB  <9,3>;
        X32bbaa [4][1] = xBBAA  <9,3>;
        X32abba [4][1] = xABAB  <9,3>;
        J32abab [4][2] = jABCD  <9,5,9,5>;
        X32aabb [4][2] = xAABB  <9,5>;
        X32bbaa [4][2] = xBBAA  <9,5>;
        X32abba [4][2] = xABAB  <9,5>;
        J32abab [4][3] = jABCD  <9,7,9,7>;
        X32aabb [4][3] = xAABB  <9,7>;
        X32bbaa [4][3] = xBBAA  <9,7>;
        X32abba [4][3] = xABAB  <9,7>;
        J32abab [4][4] = jABCD  <9,9,9,9>;
        X32aabb [4][4] = xAABB  <9,9>;
        X32bbaa [4][4] = xBBAA  <9,9>;
        X32abba [4][4] = xABAB  <9,9>;
    }






    //rotations 64
    {
        R64 [0][0] = Mprod  <1,1>;
        R64 [0][1] = Mprod  <1,3>;
        R64 [0][2] = Mprod  <1,5>;
        R64 [0][3] = Mprod  <1,7>;
        R64 [0][4] = Mprod  <1,9>;

        R64 [1][0] = Mprod  <3,1>;
        R64 [1][1] = Mprod  <3,3>;
        R64 [1][2] = Mprod  <3,5>;
        R64 [1][3] = Mprod  <3,7>;
        R64 [1][4] = Mprod  <3,9>;

        R64 [2][0] = Mprod  <5,1>;
        R64 [2][1] = Mprod  <5,3>;
        R64 [2][2] = Mprod  <5,5>;
        R64 [2][3] = Mprod  <5,7>;
        R64 [2][4] = Mprod  <5,9>;

        R64 [3][0] = Mprod  <7,1>;
        R64 [3][1] = Mprod  <7,3>;
        R64 [3][2] = Mprod  <7,5>;
        R64 [3][3] = Mprod  <7,7>;
        R64 [3][4] = Mprod  <7,9>;

        R64 [4][0] = Mprod  <9,1>;
        R64 [4][1] = Mprod  <9,3>;
        R64 [4][2] = Mprod  <9,5>;
        R64 [4][3] = Mprod  <9,7>;
        R64 [4][4] = Mprod  <9,9>;



        RT64 [0][0] = MprodT  <1,1>;
        RT64 [0][1] = MprodT  <1,3>;
        RT64 [0][2] = MprodT  <1,5>;
        RT64 [0][3] = MprodT  <1,7>;
        RT64 [0][4] = MprodT  <1,9>;

        RT64 [1][0] = MprodT  <3,1>;
        RT64 [1][1] = MprodT  <3,3>;
        RT64 [1][2] = MprodT  <3,5>;
        RT64 [1][3] = MprodT  <3,7>;
        RT64 [1][4] = MprodT  <3,9>;

        RT64 [2][0] = MprodT  <5,1>;
        RT64 [2][1] = MprodT  <5,3>;
        RT64 [2][2] = MprodT  <5,5>;
        RT64 [2][3] = MprodT  <5,7>;
        RT64 [2][4] = MprodT  <5,9>;

        RT64 [3][0] = MprodT  <7,1>;
        RT64 [3][1] = MprodT  <7,3>;
        RT64 [3][2] = MprodT  <7,5>;
        RT64 [3][3] = MprodT  <7,7>;
        RT64 [3][4] = MprodT  <7,9>;

        RT64 [4][0] = MprodT  <9,1>;
        RT64 [4][1] = MprodT  <9,3>;
        RT64 [4][2] = MprodT  <9,5>;
        RT64 [4][3] = MprodT  <9,7>;
        RT64 [4][4] = MprodT  <9,9>;
    }

    //rotations 32
    {
        R32 [0][0] = Mprod  <1,1>;
        R32 [0][1] = Mprod  <1,3>;
        R32 [0][2] = Mprod  <1,5>;
        R32 [0][3] = Mprod  <1,7>;
        R32 [0][4] = Mprod  <1,9>;

        R32 [1][0] = Mprod  <3,1>;
        R32 [1][1] = Mprod  <3,3>;
        R32 [1][2] = Mprod  <3,5>;
        R32 [1][3] = Mprod  <3,7>;
        R32 [1][4] = Mprod  <3,9>;

        R32 [2][0] = Mprod  <5,1>;
        R32 [2][1] = Mprod  <5,3>;
        R32 [2][2] = Mprod  <5,5>;
        R32 [2][3] = Mprod  <5,7>;
        R32 [2][4] = Mprod  <5,9>;

        R32 [3][0] = Mprod  <7,1>;
        R32 [3][1] = Mprod  <7,3>;
        R32 [3][2] = Mprod  <7,5>;
        R32 [3][3] = Mprod  <7,7>;
        R32 [3][4] = Mprod  <7,9>;

        R32 [4][0] = Mprod  <9,1>;
        R32 [4][1] = Mprod  <9,3>;
        R32 [4][2] = Mprod  <9,5>;
        R32 [4][3] = Mprod  <9,7>;
        R32 [4][4] = Mprod  <9,9>;



        RT32 [0][0] = MprodT  <1,1>;
        RT32 [0][1] = MprodT  <1,3>;
        RT32 [0][2] = MprodT  <1,5>;
        RT32 [0][3] = MprodT  <1,7>;
        RT32 [0][4] = MprodT  <1,9>;

        RT32 [1][0] = MprodT  <3,1>;
        RT32 [1][1] = MprodT  <3,3>;
        RT32 [1][2] = MprodT  <3,5>;
        RT32 [1][3] = MprodT  <3,7>;
        RT32 [1][4] = MprodT  <3,9>;

        RT32 [2][0] = MprodT  <5,1>;
        RT32 [2][1] = MprodT  <5,3>;
        RT32 [2][2] = MprodT  <5,5>;
        RT32 [2][3] = MprodT  <5,7>;
        RT32 [2][4] = MprodT  <5,9>;

        RT32 [3][0] = MprodT  <7,1>;
        RT32 [3][1] = MprodT  <7,3>;
        RT32 [3][2] = MprodT  <7,5>;
        RT32 [3][3] = MprodT  <7,7>;
        RT32 [3][4] = MprodT  <7,9>;

        RT32 [4][0] = MprodT  <9,1>;
        RT32 [4][1] = MprodT  <9,3>;
        RT32 [4][2] = MprodT  <9,5>;
        RT32 [4][3] = MprodT  <9,7>;
        RT32 [4][4] = MprodT  <9,9>;
    }
}



