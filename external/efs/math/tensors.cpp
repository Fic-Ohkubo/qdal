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
#include "../math/eigen.hpp"
#include "../math/tensors.hpp"
#include "../math/newsparse.hpp"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif



using namespace std;


void tensor2::setsize(int p) {
    if (n==p && m==p) return;

    if (n!=0 || m!=0) {
        delete[] c2;
	    delete[] c;
    }

	n = p;
	m = p;

	c  = new double*[n];
    c2 = new double[n*n];

	//asigna los punteros
    //#pragma omp parallel for
	for (int i=0;i<n;++i)
		c[i] = c2 + i*n;
}

void tensor2::setsize(int p, int q) {
    if (n==p && m==q) return;

    if (n!=0 || m!=0) {
        delete[] c2;
	    delete[] c;
    }

	n = p;
	m = q;

	c  = new double*[n];
    c2 = new double[n*m];

	//asigna los punteros
    //#pragma omp parallel for
	for (int i=0;i<n;++i)
		c[i] = c2 + i*m;
}

void tensor2::zeroize() {
    //#pragma omp parallel for
	for (int i=0; i<n*m; ++i)
		c2[i] = 0.;
}

void tensor2::clear() {
	delete[] c2;
	delete[] c;
	c2=NULL;
	c =NULL;
	n=0;
	m=0;
}

tensor2 & tensor2::operator=(const tensor2 & st) {
    if (this == &st) return *this;

    setsize(st.n, st.m);

    //copia los valores
    for (int i=0;i<n*m;++i)
        c2[i] = st.c2[i];

    return *this;
}

tensor2 & tensor2::operator=(const symtensor2 & st) {

    setsize(st.n, st.n);

    //copia los valores
    for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            c[i][j] = st(i,j);

    return *this;
}


const tensor2 tensor2::operator*(const tensor2 & rhs) const {
    tensor2 result;
    result.setsize(n, rhs.m);

    #pragma omp parallel for
    for (int i=0; i<n; ++i) {
        for (int j=0; j<rhs.m; ++j) {
            double sum = 0;
            for (int k=0; k<m; ++k) sum += (*this)(i,k) * rhs(k,j);
            result(i,j) = sum;
        }
    }

    return result;
}


static inline UI32 Tpos (UI32 M1, UI32 M2, UI32 v1, UI32 v2) {
    return v2*M1 + v1;
}

tensor2 & tensor2::operator=(const Sparse & rhs) {
    if (n!=rhs.ntotal) {
        if (n>0)
            clear();
        setsize(rhs.ntotal);
    }

    if (rhs.values==NULL) {
        return *this;
    }


    //foreach atom
    for (UI32 ati=0; ati<rhs.natoms; ++ati) {
        for (UI32 atj=0; atj<rhs.natoms; ++atj) {

            UI16 id1   = rhs.ids[ati];
            UI16 id2   = rhs.ids[atj];

            double * offset1 = rhs.values + rhs.a2pos[ati][atj];

            UI32 p1 = rhs.a1pos[ati];
            UI32 q1 = rhs.a1pos[atj];


            for (int b1=0; b1<rhs.nfs[id1]; ++b1) {
                for (int b2=0; b2<rhs.nfs[id2]; ++b2) {
                    double * offset2 = offset1 + rhs.f2pos[id1][id2][b1][b2];

                    UI32 p2 = p1 + rhs.f1pos[id1][b1];
                    UI32 q2 = q1 + rhs.f1pos[id2][b2];

                    UI32 mj1 = rhs.js[id1][b1];
                    UI32 mj2 = rhs.js[id2][b2];
                    UI32 mm1 = rhs.ms[id1][b1];
                    UI32 mm2 = rhs.ms[id2][b2];

                    for (UI8 j1=0; j1<mj1; ++j1) {
                        for (UI8 j2=0; j2<mj2; ++j2) {
                            double * offset3 = offset2  +  (j1*mj2+j2)*(mm1*mm2);

                            UI32 p3 = p2 + j1*mm1;
                            UI32 q3 = q2 + j2*mm2;


                            for (UI8 m1=0; m1<mm1; ++m1) {
                                for (UI8 m2=0; m2<mm2; ++m2) {
                                    UI32 Ap = Tpos(mm1,mm2, m1,m2);

                                    UI32 p4 = p3 + m1;
                                    UI32 q4 = q3 + m2;

                                    (*this)(p4, q4) = rhs(ati,atj, b1,b2, j1,j2)[Ap];
                                }
                            }

                        }
                    }

                }
            }


        }
    }

    return *this;
}

std::ostream & operator<<(std::ostream & os, const tensor2 & T) {
    for (int i=0; i<T.n; ++i) {
        for (int j=0; j<T.n; ++j) {
            os << T(i,j) << " ";
        }
        os << std::endl << std::endl;
    }
    return os;
}


void symtensor2::setsize(int m) {
    if (n==m) return;

    if (n!=0) {
        delete[] c2;
        delete[] c;
    }

	n = m;
	n2 = (n*n+n)/2;

	c2 = new double[n2];
	c  = new double*[n];

	//asigna punteros
	int p2=0;

	for (int i=0; i<n; ++i) {
		c[i] = c2 + p2;
		p2 += i+1;
	}
}

void symtensor2::zeroize() {
	for (int i=0; i<n2; ++i)
		c2[i] = 0.;
}

void symtensor2::clear() {
	delete[] c;
	delete[] c2;
	c2=NULL;
	c=NULL;
	n=0;
}

symtensor2 & symtensor2::operator=(const symtensor2 & st) {
    if (this == &st) return *this;

    setsize(st.n);

    //copia los valores
    for (int i=0;i<n2;++i)
        c2[i] = st.c2[i];

    return *this;
}

symtensor2 & symtensor2::operator=(const tensor2 & st) {

    setsize(st.n);

    //copia los valores
    for (int i=0;i<n;++i)
        for (int j=0; j<=i; ++j)
            c[i][j] = st.c[i][j];

    return *this;
}

symtensor2 & symtensor2::operator+=(const symtensor2 & st) {
    #pragma omp parallel for
    for (int i=0;i<n2;++i)
        c2[i] += st.c2[i];
    return *this;
}

symtensor2 & symtensor2::operator-=(const symtensor2 & st) {
    #pragma omp parallel for
    for (int i=0;i<n2;++i)
        c2[i] -= st.c2[i];
    return *this;
}

symtensor2 & symtensor2::operator=(double s) {
    zeroize();
    for (int i=0;i<n;++i)
        (*this)(i,i) = s;
    return *this;
}

symtensor2 & symtensor2::operator+=(double s) {
    for (int i=0;i<n;++i)
        (*this)(i,i) += s;
    return *this;
}

symtensor2 & symtensor2::operator*=(double s) {
    #pragma omp parallel for
    for (int i=0;i<n2;++i)
        c2[i] *= s;
    return *this;
}

const symtensor2 symtensor2::operator*(const symtensor2 & rhs) const {
    symtensor2 result;

    result.setsize(n);

    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<n; ++i) {
        for (int j=0; j<=i; ++j) {
            double sum = 0;
            for (int k=0; k<n; ++k)
                sum += (*this)(i,k)*rhs(k,j);
            result(i,j) = sum;
        }
    }

    return result;
}

std::ostream & operator<<(std::ostream & os, const symtensor2 & T) {
    for (int i=0; i<T.n; ++i) {
        for (int j=0; j<T.n; ++j) {
            os << T(i,j) << " ";
        }
        os << std::endl << std::endl;
    }
    return os;
}

