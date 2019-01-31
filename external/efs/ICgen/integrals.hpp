#ifndef __INTEGRALS__
#define __INTEGRALS__
#include <set>
#include <map>
#include <iostream>
#include "../defs.hpp"

typedef UI32 IDINT;
typedef UI8  AM;

class integralCart;

class ikernel {
  public:
    UI8 u;
    UI8 v;
    UI8 s;
    UI8 t;
    UI8 m;

    UI8 K4L; // loop
    UI8 FE;  // type of kernel (0: regular  1: auxiliary  2: second auxiliary  -1: don't compute)
    UI8 align[1];

    ikernel();
    ikernel(int ru, int rv, int rs, int rt, int rm);
    ikernel & operator=(const integralCart & rhs);
    ikernel & operator+=(const ikernel & rhs);
    ikernel & operator-=(const ikernel & rhs);
    ikernel & operator|=(const ikernel & rhs);
    ikernel operator|(const ikernel & rhs) const;
    ikernel operator+(const ikernel & rhs) const;
    ikernel operator-(const ikernel & rhs) const;

    bool operator<(const ikernel & rhs) const;
};

class integralCart {
    public:
    UI8 ex; UI8 ey; UI8 ez;
    UI8 bx; UI8 by; UI8 bz;
    UI8 fx; UI8 fy; UI8 fz;
    UI8 dx; UI8 dy; UI8 dz;
    UI8 px; UI8 py; UI8 pz;
    UI8 qx; UI8 qy; UI8 qz;
    UI8 rx; UI8 ry; UI8 rz;
    UI8 a; UI8 b; UI8 p; UI8 c; UI8 d; UI8 q; UI8 m;
    UI8 la; UI8 lb; UI8 lc; UI8 ld;

    integralCart();
    integralCart(int nn);
    integralCart(int rla, int rlb, int rlc, int rld);    //(4)
    integralCart(int rax, int ray, int raz, int rbx, int rby, int rbz, int rcx, int rcy, int rcz, int rdx, int rdy, int rdz);    //(12)
    integralCart(int rex, int rey, int rez, int rfx, int rfy, int rfz, int rrx, int rry, int rrz, int ra, int rb, int rp, int rc, int rd, int rq, int rm); //(16)
    integralCart(int rex, int rey, int rez, int rfx, int rfy, int rfz, int rpx, int rpy, int rpz, int rqx, int rqy, int rqz, int rrx, int rry, int rrz,
                 int ra, int rb, int rp, int rc, int rd, int rq, int rm); //(22)
    integralCart(int rex, int rey, int rez, int rbx, int rby, int rbz, int rfx, int rfy, int rfz, int rdx, int rdy, int rdz, int rpx, int rpy, int rpz,
                 int rqx, int rqy, int rqz, int rrx, int rry, int rrz, int ra, int rb, int rp, int rc, int rd, int rq, int rm, int rla, int rlb, int rlc, int rld); //(32)
    integralCart & operator=(const ikernel & rhs);
    integralCart & operator|=(const integralCart & rhs);
    integralCart operator|(const integralCart & rhs) const;
    integralCart & operator+=(const integralCart & rhs);
    integralCart operator+(const integralCart & rhs) const;
    virtual bool operator<(const integralCart & rhs) const;
    bool operator==(const integralCart & rhs) const;
    bool operator!=(const integralCart & rhs) const;
} __attribute__ ((aligned(32)));

class integralSph {
  public:
    SI8 la;
    SI8 ma;
    SI8 lb;
    SI8 mb;
    SI8 lc;
    SI8 mc;
    SI8 ld;
    SI8 md;

    SI8 b;
    SI8 p;
    SI8 d;
    SI8 q;
    SI8 m;

    integralSph();
    integralSph(int rla, int rma, int rlb, int rmb, int rlc, int rmc, int rld, int rmd);
    integralSph(int rla, int rma, int rlb, int rmb, int rlc, int rmc, int rld, int rmd,
                int rb, int rp, int rd, int rq, int rm);
    integralSph(int rla, int rma, int rlb, int rmb, int rlc, int rmc, int rld, int rmd,
                int rb, int rp, int rd, int rq, int rm, bool rab);
    integralSph(int rla, int rma, int rlb, int rmb, int rlc, int rmc, int rld, int rmd,
                int rb, int rp, int rd, int rq, int rm, bool rab, bool rcdab);

    integralSph operator|(const integralSph & rhs) const;
    integralSph & operator+=(const integralSph & rhs);
    integralSph operator+(const integralSph & rhs) const;
    integralSph & operator-=(const integralSph & rhs);
    integralSph operator-(const integralSph & rhs) const;
    bool operator<(const integralSph & rhs) const;
};


std::ostream & operator<< (std::ostream & os, const ikernel & rhs);
std::ostream & operator<< (std::ostream & os, const integralCart & rhs);
std::ostream & operator<< (std::ostream & os, const integralSph & rhs);

#endif
