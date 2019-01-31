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


#ifndef __TENSORS__
#define __TENSORS__

#include <iostream>
#include <fstream>

struct tensor;
struct symtensor;

struct symtensor2;
struct Sparse;

//matriz de n x n  o  n x m ->empieza en 0 ; los indices van de 0 a n no inclusive
struct tensor2 {
    int n; //primer indice del tensor
	int m; //segundo indice del tensor

	double ** c;
	double * c2;

	tensor2() {n=0; m=0; c=NULL; c2=NULL;}
	tensor2(int p) {n=0; m=0; setsize(p);}
	tensor2(int p, int q) {n=0; m=0; setsize(p,q);}

	void setsize(int p);
	void setsize(int p, int q);
	void zeroize();
	void clear();

	void Transpose() {
	    if (n==m)
	    for (int p=0; p<n; ++p) {
	        for (int q=0; q<p; ++q) {
	            std::swap(c[p][q], c[q][p]);
	        }
	    }
	}

	void Symmetrize() {
	    if (n==m)
	    for (int p=0; p<n; ++p) {
	        for (int q=0; q<=p; ++q) {
	            c[p][q] += c[q][p];
	            c[q][p] =  c[p][q];
	        }
	    }
	}

	inline double & operator()(int i, int j) {
		return c[i][j];
	}

	inline double operator()(int i, int j) const {
		return c[i][j];
	}

    tensor2 & operator+=(const tensor2 & s) {
        for (int i=0;i<n*m;++i)
            c2[i] += s.c2[i];
        return *this;
    }

    tensor2 & operator-=(const tensor2 & s) {
        for (int i=0;i<n*m;++i)
            c2[i] -= s.c2[i];
        return *this;
    }

    tensor2 & operator*=(double s) {
        for (int i=0;i<n*m;++i)
            c2[i] *= s;
        return *this;
    }

    const tensor2 operator*(const tensor2 & rhs) const;
	tensor2 & operator=(const tensor2 & st);
	tensor2 & operator=(const tensor  & st);
	tensor2 & operator=(const symtensor  & st);

	tensor2 & operator=(const symtensor2 & st);
	tensor2 & operator=(const Sparse & rhs);
};

std::ostream & operator<<(std::ostream & os, const tensor2 & T);

//tensor triangular de rango 2; cuidado: los indices van de 0 a n no inclusive
struct symtensor2 {
	int n;
	int n2;
	double **c;
	double *c2;

	symtensor2() {n=0; c=NULL; c2=NULL;}
	symtensor2(int m) {n=0; setsize(m);}

	void setsize(int m);
	void zeroize();
	void clear();

	inline double & operator()(int i, int j) {
		if(i>=j) return c[i][j];
		else return c[j][i];
	}

	inline double operator()(int i, int j) const {
		if(i>=j) return c[i][j];
		else return c[j][i];
	}

    symtensor2 & operator= (double rhs);
    symtensor2 & operator+=(double rhs);
    symtensor2 & operator*=(double rhs);
    symtensor2 & operator+=(const symtensor2 & st);
    symtensor2 & operator-=(const symtensor2 & st);

    symtensor2 & operator= (const tensor2 & st);
    symtensor2 & operator= (const symtensor2 & st);
    symtensor2 & operator= (const symtensor  & st);
    symtensor2 & operator+=(const symtensor  & st);
    symtensor2 & operator+=(const tensor  & st);

    const symtensor2 operator*(const symtensor2 & rhs) const;
};

std::ostream & operator<<(std::ostream & os, const symtensor2 & T);

#endif
