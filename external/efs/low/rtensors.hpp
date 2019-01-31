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

#ifndef __RTENSORS__
#define __RTENSORS__

#include <queue>
#include <cstdlib>
#include <algorithm>

template<class T> class r3tensor {
    private:
        T * v3;
        T ** v2;
        T *** v;
        int Nx;
        int Ny;
        int Nz;

    public:
        r3tensor() {
            Nx = Ny = Nz = 0;
            v3 = NULL;
            v2 = NULL;
            v  = NULL;
        }

        int N1() const {
            return Nx;
        }

        int N2() const {
            return Ny;
        }

        int N3() const {
            return Nz;
        }



        void set(int nx, int ny, int nz) {
            if (v3!=NULL) delete[] v3;
            if (v2!=NULL) delete[] v2;
            if (v !=NULL) delete[] v;

            v3 = new T[nx*ny*nz];
            v2 = new T*[nx*ny];
            v  = new T**[nx];

            Nx = nx;
            Ny = ny;
            Nz = nz;

            for (int i=0; i<Nx; ++i) {
                v[i] = v2 + i*Ny;
            }

            for (int i=0; i<Nx; ++i) {
                for (int j=0; j<Ny; ++j) {
                    v[i][j] = v3 + (i*Ny+j)*Nz;
                }
            }

        }

        r3tensor(int nx, int ny, int nz) {
            v3 = NULL;
            v2 = NULL;
            v  = NULL;
            set(nx,ny,nz);
        }

        inline const T & operator()(int i, int j, int k) const {
            return v[i][j][k];
        }

        inline T & operator()(int i, int j, int k) {
            return v[i][j][k];
        }

        ~r3tensor() {
            delete[] v3;
            delete[] v2;
            delete[] v;
        }
};

template<class T> class r2tensor {
  protected:
    T * v2;
    T ** v;
    int n;
    int m;
    int * mm;

  public:
    r2tensor() {
        n = m = 0;
        v2 = NULL;
        v  = NULL;
        mm = NULL;
    }

    int N() const {
        return n;
    }

    int M() const {
        return m;
    }

    int M(int i) const {
        return mm[i];
    }

    void set(unsigned short int * lens, int nn) {
        delete[] v;
        delete[] v2;
        delete[] mm;

        unsigned int tlen=0;
        for (int i=0; i<nn; ++i) tlen += lens[i];

        v2 = new T[tlen];
        v  = new T*[nn];
        mm = new int[nn];

        v[0] = v2;
        for (int i=0; i<nn-1; ++i) {
            v[i+1] = v[i] + lens[i];
        }

        n  = nn;
        m=0;
        for (int i=0; i<nn; ++i) {
            mm[i] = lens[i];
            m = std::max(m, mm[i]);
        }
    }

    void set(int nnn, int mmm) {
        delete[] v;
        delete[] v2;
        delete[] mm;
        mm = NULL;

        v2 = new T[nnn*mmm];
        v  = new T*[nnn];
        n  = nnn;
        m  = mmm;

        for (int i=0; i<nnn; ++i) {
            v[i] = v2 + i*mmm;
        }
    }

    void set(int nnn, int mmm, T**w, T*w2) {
        v2 = w2;
        v  = w;
        n  = nnn;
        m  = mmm;

        for (int i=0; i<nnn; ++i) {
            v[i] = v2 + i*mmm;
        }
    }

    r2tensor(int nnn, int mmm) {
        v  = NULL;
        v2 = NULL;
        mm = NULL;
        set(nnn,mmm);
    }

    inline const T & operator()(int i, int j) const {
        return v[i][j];
    }

    inline T & operator()(int i, int j) {
        return v[i][j];
    }

    inline const T * operator[](int i) const {
        return v[i];
    }

    inline T * operator[](int i) {
        return v[i];
    }

    void Nullify() {
        v2 = NULL;
        v  = NULL;
    }

    ~r2tensor() {
        delete[] v2;
        delete[] v;
        delete[] mm;
    }

};

template<class T> class r1tensor {
  private:
    T * v;
    int n;
  public:
    r1tensor() {
        v  = NULL;
        n  = 0;
    }

    void set(int nn) {
        delete[] v;
        v  = new T[nn];
        n  = nn;
    }

    r1tensor(int nn) {
        v = NULL;
        set(nn);
    }

    int N() const {
        return n;
    }

    inline const T & operator[](int i) const {
        return v[i];
    }

    inline T & operator[](int i) {
        return v[i];
    }

    inline const T & operator()(int i) const {
        return v[i];
    }

    inline T & operator()(int i) {
        return v[i];
    }

    ~r1tensor() {
        delete[] v;
    }


    const r2tensor<T> operator*(const r1tensor<T> & rhs) const {
        r2tensor<T> ret;
        ret.set(n, rhs.n);

        for (int i=0; i<n; ++i)
            for (int j=0; j<n; ++j)
                ret(i,j) = v[i] * rhs.v[j];

        return ret;
    }
};

template<class T> class r2tensor<r2tensor<T> > {
  protected:

    r2tensor<T> * v2;
    r2tensor<T> ** v;
    T * w2;
    T ** w;

    int n;
    int m;

  public:
    r2tensor() {
        n = m = 0;
        v2 = NULL;
        v  = NULL;
        w2 = NULL;
        w  = NULL;
    }

    int N() const {
        return n;
    }

    int M() const {
        return m;
    }

    void set(int nnn, int mmm) {
        delete[] v;
        delete[] v2;

        v2 = new r2tensor<T>[nnn*mmm];
        v  = new r2tensor<T>*[nnn];
        n  = nnn;
        m  = mmm;

        for (int i=0; i<nnn; ++i) {
            v[i] = v2 + i*mmm;
        }
    }

    r2tensor(int nn, int mm) {

        v2 = NULL;
        v  = NULL;
        w2 = NULL;
        w  = NULL;
        set(nn,mm);
    }

    inline const r2tensor<T> & operator()(int i, int j) const {
        return v[i][j];
    }

    inline r2tensor<T> & operator()(int i, int j) {
        return v[i][j];
    }

    inline const r2tensor<T> * operator[](int i) const {
        return v[i];
    }

    inline r2tensor<T> * operator[](int i) {
        return v[i];
    }

    void Nullify() {
        v2 = NULL;
        v  = NULL;
    }

    ~r2tensor() {

        for (int nn=0; nn<n;  ++nn) {
            for (int mm=0; mm<m;  ++mm) {
                v[nn][mm].Nullify(); //otherwise it will try to deallocate a pointer that is not owned by the object
            }
        }

        delete[] v2;
        delete[] v;
        delete[] w2;
        delete[] w;
    }

    void set(const r1tensor<int> & ns, const r1tensor<int> & ms) {
        set(ns.N(), ns.N());

        int nw2 = 0;
        int nw  = 0;

        for (int nn=0; nn<n;  ++nn) {
            for (int mm=0; mm<m;  ++mm) {
                nw2 += ns(nn) * ms(mm);
                nw  += ns(nn);
            }
        }

        w2 = new T[nw2];
        w  = new T*[nw];

        T *   ww2 = w2;
        T **  ww  = w;


        for (int nn=0; nn<n;  ++nn) {
            for (int mm=0; mm<m;  ++mm) {
                v[nn][mm].set(ns(nn), ms(mm), ww, ww2);

                ww2 += ns(nn) * ms(mm);
                ww  += ns(nn);
            }
        }


    }

};


template <class T> struct Array {
    int n;
    T * keys;

    Array() {
        n = 0;
        keys = NULL;
    }

    Array(int nelem) {
        n = nelem;
        keys = new T[n];
    }

    void set(int nelem) {
        //if (n!=0) delete[] keys;
        n = nelem;
        keys = new T[n];
    }

    int search(const T & val) const {
        for (int i=0; i<n; ++i)
            if (keys[i] == val) return i;
        return -1; // should not happen
    }

    bool isin(const T & val, int nmax) const {
        for (int i=0; i<nmax; ++i)
            if (keys[i] == val) return true;
        return false;
    }

	T & operator[] (int i) {
		return keys[i];
	}

    const T & operator[] (int i) const {
		return keys[i];
	}

    /*
	Array<T> & operator=(const std::set<T> & TheSet) {

	    typedef std::set<T> set_T;

	    set_T::iterator it;
	    int i=0;

	    n = TheSet.size();
	    keys = new T[n];

	    for (it=TheSet.begin(); it!=TheSet.end(); ++it) {
	        keys[i] = *it;
	        ++i;
	    }

	    return *this;
	}
	*/


    void clear() {
        if (n!=0) delete[] keys;
        n = 0;
    }

    ~Array() {
        if (n!=0) delete[] keys;
    }

    Array<T> & operator=(const Array<T> & rhs) {
        set(rhs.n);
        for (int i=0; i<n; ++i)
            keys[i] = rhs.keys[i];
        return *this;
    }
};

#endif
