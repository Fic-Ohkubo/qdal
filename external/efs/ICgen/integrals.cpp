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


#include "integrals.hpp"
using namespace std;

ikernel::ikernel() {
    u = v = s = t = m = K4L = FE = align[0] = 0;
}

ikernel::ikernel(int ru, int rv, int rs, int rt, int rm) {
    u = ru;
    v = rv;
    s = rs;
    t = rt;
    m = rm;

    K4L = 0;
    FE  = 0;

    align[0]  = 0;
}

ikernel & ikernel::operator=(const integralCart & rhs) {
    u = rhs.b;
    v = rhs.p;
    s = rhs.d;
    t = rhs.q;
    m = rhs.m;

    K4L = 0;
    FE  = 0;

    return *this;
}

ikernel & ikernel::operator|=(const ikernel & rhs) {
    u *= rhs.u;
    v *= rhs.v;
    s *= rhs.s;
    t *= rhs.t;
    m *= rhs.m;

    return *this;
}

ikernel & ikernel::operator+=(const ikernel & rhs) {
    u += rhs.u;
    v += rhs.v;
    s += rhs.s;
    t += rhs.t;
    m += rhs.m;

    return *this;
}

ikernel & ikernel::operator-=(const ikernel & rhs) {
    u -= rhs.u;
    v -= rhs.v;
    s -= rhs.s;
    t -= rhs.t;
    m -= rhs.m;

    return *this;
}

ikernel ikernel::operator|(const ikernel & rhs) const {
    ikernel ret;
    ret = *this;
    ret |= rhs;
    return ret;
}

ikernel ikernel::operator+(const ikernel & rhs) const {
    ikernel ret;
    ret = *this;
    ret += rhs;
    return ret;
}

ikernel ikernel::operator-(const ikernel & rhs) const {
    ikernel ret;
    ret = *this;
    ret -= rhs;
    return ret;
}

bool ikernel::operator<(const ikernel & rhs) const {
    if (K4L!=rhs.K4L) return (K4L<rhs.K4L);
    if (FE !=rhs.FE)  return (FE >rhs.FE); // first SAEs, then AEs, then kernels

    if (m!=rhs.m) return (m<rhs.m);
    if (v!=rhs.v) return (v<rhs.v);
    if (u!=rhs.u) return (u<rhs.u);
    if (t!=rhs.t) return (t<rhs.t);
    if (s!=rhs.s) return (s<rhs.s);

    return false;
}





integralCart::integralCart() {
    ex = ey = ez = bx = by = bz = fx = fy = fz = dx = dy = dz = px = py = pz = qx = qy = qz = rx = ry = rz = a = b = p = c = d = q = m = la = lb = lc = ld = 0;
}

integralCart::integralCart(int nn) {
    ex = ey = ez = bx = by = bz = fx = fy = fz = dx = dy = dz = px = py = pz = qx = qy = qz = rx = ry = rz = a = b = p = c = d = q = m = la = lb = lc = ld = nn;
}


integralCart::integralCart(int rla, int rlb, int rlc, int rld) {
    la = rla;
    lb = rlb;
    lc = rlc;
    ld = rld;
    ex = ey = ez = bx = by = bz = fx = fy = fz = dx = dy = dz = px = py = pz = qx = qy = qz = rx = ry = rz = a = b = p = c = d = q = m = 0;
}

//(12)
integralCart::integralCart(int rax, int ray, int raz, int rbx, int rby, int rbz, int rcx, int rcy, int rcz, int rdx, int rdy, int rdz) {
    ex = rax;
    ey = ray;
    ez = raz;
    bx = rbx;
    by = rby;
    bz = rbz;
    fx = rcx;
    fy = rcy;
    fz = rcz;
    dx = rdx;
    dy = rdy;
    dz = rdz;
    px = py = pz = qx = qy = qz = rx = ry = rz = b = p = d = q = m = 0;
    la = lb = lc = ld = 0; //never used
    a = c = 0; //never used
}
//(16)
integralCart::integralCart(int rex, int rey, int rez, int rfx, int rfy, int rfz, int rrx, int rry, int rrz, int ra, int rb, int rp, int rc, int rd, int rq, int rm) {
    ex = rex;
    ey = rey;
    ez = rez;
    fx = rfx;
    fy = rfy;
    fz = rfz;
    rx = rrx;
    ry = rry;
    rz = rrz;
    a = ra;
    b = rb;
    p = rp;
    c = rc;
    d = rd;
    q = rq;
    m = rm;
    la = lb = lc = ld = 0; //never used
    bx = 0;
    by = 0;
    bz = 0;
    dx = 0;
    dy = 0;
    dz = 0;
    px = py = pz = qx = qy = qz = 0;
}

//(22)
integralCart::integralCart(int rex, int rey, int rez,
             int rfx, int rfy, int rfz,
             int rpx, int rpy, int rpz,
             int rqx, int rqy, int rqz,
             int rrx, int rry, int rrz,
             int ra, int rb, int rp, int rc, int rd, int rq, int rm) {
    ex = rex;
    ey = rey;
    ez = rez;
    fx = rfx;
    fy = rfy;
    fz = rfz;
    rx = rrx;
    ry = rry;
    rz = rrz;
    a = ra;
    b = rb;
    p = rp;
    c = rc;
    d = rd;
    q = rq;
    m = rm;
    la = lb = lc = ld = 0; //never used
    bx = 0;
    by = 0;
    bz = 0;
    dx = 0;
    dy = 0;
    dz = 0;
    px = rpx;
    py = rpy;
    pz = rpz;
    qx = rqx;
    qy = rqy;
    qz = rqz;
}


//(32)
integralCart::integralCart(int rex, int rey, int rez,
             int rbx, int rby, int rbz,
             int rfx, int rfy, int rfz,
             int rdx, int rdy, int rdz,
             int rpx, int rpy, int rpz,
             int rqx, int rqy, int rqz,
             int rrx, int rry, int rrz,
             int ra, int rb, int rp, int rc, int rd, int rq, int rm,
             int rla, int rlb, int rlc, int rld
             ) {
    ex = rex;
    ey = rey;
    ez = rez;
    bx = rbx;
    by = rby;
    bz = rbz;
    fx = rfx;
    fy = rfy;
    fz = rfz;
    dx = rdx;
    dy = rdy;
    dz = rdz;

    px = rpx;
    py = rpy;
    pz = rpz;
    qx = rqx;
    qy = rqy;
    qz = rqz;
    rx = rrx;
    ry = rry;
    rz = rrz;

    a = ra;
    b = rb;
    p = rp;
    c = rc;
    d = rd;
    q = rq;
    m = rm;

    la = rla;
    lb = rlb;
    lc = rlc;
    ld = rld;
}




integralCart & integralCart::operator=(const ikernel & rhs) {
    ex = ey = ez = bx = by = bz = fx = fy = fz = dx = dy = dz = 0;
    px = py = pz = qx = qy = qz = rx = ry = rz = 0;
    a = 0;
    b = rhs.u;
    p = rhs.v;
    c = 0;
    d = rhs.s;
    q = rhs.t;
    m = rhs.m;
    la = lb = lc = ld = 0;
    return *this;
}

integralCart & integralCart::operator|=(const integralCart & rhs) {
    ex *= rhs.ex;
    ey *= rhs.ey;
    ez *= rhs.ez;
    bx *= rhs.bx;
    by *= rhs.by;
    bz *= rhs.bz;
    fx *= rhs.fx;
    fy *= rhs.fy;
    fz *= rhs.fz;
    dx *= rhs.dx;
    dy *= rhs.dy;
    dz *= rhs.dz;

    px *= rhs.px;
    py *= rhs.py;
    pz *= rhs.pz;
    qx *= rhs.qx;
    qy *= rhs.qy;
    qz *= rhs.qz;

    rx *= rhs.rx;
    ry *= rhs.ry;
    rz *= rhs.rz;
    //ret.a  *= rhs.a;
    b  *= rhs.b;
    p  *= rhs.p;
    //ret.c  *= rhs.c;
    d  *= rhs.d;
    q  *= rhs.q;
    m  *= rhs.m;
    //ret.la *= rhs.la;
    //ret.lb *= rhs.lb;
    //ret.lc *= rhs.lc;
    //ret.ld *= rhs.ld;
    return *this;
}

integralCart integralCart::operator|(const integralCart & rhs) const {
    integralCart r = *this;
    r |= rhs;
    return r;
}

//in order to apply RRs in a more sensible way
integralCart & integralCart::operator+=(const integralCart & rhs) {
    ex += rhs.ex;
    ey += rhs.ey;
    ez += rhs.ez;
    bx += rhs.bx;
    by += rhs.by;
    bz += rhs.bz;
    fx += rhs.fx;
    fy += rhs.fy;
    fz += rhs.fz;
    dx += rhs.dx;
    dy += rhs.dy;
    dz += rhs.dz;
    px += rhs.px;
    py += rhs.py;
    pz += rhs.pz;
    qx += rhs.qx;
    qy += rhs.qy;
    qz += rhs.qz;
    rx += rhs.rx;
    ry += rhs.ry;
    rz += rhs.rz;
    a  += 0;
    b  += rhs.b;
    p  += rhs.p;
    c  += 0;
    d  += rhs.d;
    q  += rhs.q;
    m  += rhs.m;
    la += 0;
    lb += 0;
    lc += 0;
    ld += 0;
    return *this;
}

integralCart integralCart::operator+(const integralCart & rhs) const{
    integralCart r = *this;
    r += rhs;
    return r;
}

bool integralCart::operator<(const integralCart & rhs) const {
    /*
    nopn-portable (endianness)
    if (*(unsignd long long int*)(&ex) != *(unsigned long long int*)(&rhs.ex)) return *(unsigned long long int*)(&ex) < *(unsigned long long int*)(&rhs.ex);
    if (*(unsigned long long int*)(&fz) != *(unsigned long long int*)(&rhs.fz)) return *(unsigned long long int*)(&fz) < *(unsigned long long int*)(&rhs.fz);
    if (*(unsigned long long int*)(&qy) != *(unsigned long long int*)(&rhs.qy)) return *(unsigned long long int*)(&qy) < *(unsigned long long int*)(&rhs.qy);
    if (*(unsigned long long int*)(&c)  != *(unsigned long long int*)(&rhs.c))  return *(unsigned long long int*)(&c)  < *(unsigned long long int*)(&rhs.c);
    return false;
    */

    if (ex!=rhs.ex) return (ex<rhs.ex);
    if (ey!=rhs.ey) return (ey<rhs.ey);
    if (ez!=rhs.ez) return (ez<rhs.ez);

    if (bx!=rhs.bx) return (bx<rhs.bx);
    if (by!=rhs.by) return (by<rhs.by);
    if (bz!=rhs.bz) return (bz<rhs.bz);

    if (fx!=rhs.fx) return (fx<rhs.fx);
    if (fy!=rhs.fy) return (fy<rhs.fy);
    if (fz!=rhs.fz) return (fz<rhs.fz);

    if (dx!=rhs.dx) return (dx<rhs.dx);
    if (dy!=rhs.dy) return (dy<rhs.dy);
    if (dz!=rhs.dz) return (dz<rhs.dz);

    if (px!=rhs.px) return (px<rhs.px);
    if (py!=rhs.py) return (py<rhs.py);
    if (pz!=rhs.pz) return (pz<rhs.pz);
    if (qx!=rhs.qx) return (qx<rhs.qx);
    if (qy!=rhs.qy) return (qy<rhs.qy);
    if (qz!=rhs.qz) return (qz<rhs.qz);
    if (rx!=rhs.rx) return (rx<rhs.rx);
    if (ry!=rhs.ry) return (ry<rhs.ry);
    if (rz!=rhs.rz) return (rz<rhs.rz);

    if (a!=rhs.a) return (a<rhs.a);
    if (b!=rhs.b) return (b<rhs.b);
    if (p!=rhs.p) return (p<rhs.p);
    if (c!=rhs.c) return (c<rhs.c);
    if (d!=rhs.d) return (d<rhs.d);
    if (q!=rhs.q) return (q<rhs.q);
    if (m!=rhs.m) return (m<rhs.m);

    if (la!=rhs.la) return (la<rhs.la);
    if (lb!=rhs.lb) return (lb<rhs.lb);
    if (lc!=rhs.lc) return (lc<rhs.lc);
    if (ld!=rhs.ld) return (ld<rhs.ld);

    return false;
}

bool integralCart::operator==(const integralCart & rhs) const {
    return (ex==rhs.ex && ey==rhs.ey && ez==rhs.ez && bx==rhs.bx && by==rhs.by && bz==rhs.bz &&
             fx==rhs.fx && fy==rhs.fy && fz==rhs.fz && dx==rhs.dx && dy==rhs.dy && dz==rhs.dz &&
             px==rhs.px && py==rhs.py && pz==rhs.pz && qx==rhs.qx && qy==rhs.qy && qz==rhs.qz &&
             rx==rhs.rx && ry==rhs.ry && rz==rhs.rz &&
             a==rhs.a && b==rhs.b && p==rhs.p && c==rhs.c && d==rhs.d && q==rhs.q && m==rhs.m &&
             la==rhs.la && lb==rhs.lb &&  lc==rhs.lc && ld==rhs.ld);

    //return ((*(long long int*)(&ex)) == (*(long long int*)(&rhs.ex))) && ((*(long long int*)(&fz)) == (*(long long int*)(&rhs.fz))) && ((*(long long int*)(&qy)) == (*(long long int*)(&rhs.qy))) && ((*(long long int*)(&c)) == (*(long long int*)(&rhs.c)));
}

bool integralCart::operator!=(const integralCart & rhs) const {
    return !(*this == rhs);
}





integralSph::integralSph() {
    la = ma = lb = mb = lc = mc = ld = md = b = p = d = q = m = 0;
}

integralSph::integralSph(int rla, int rma, int rlb, int rmb, int rlc, int rmc, int rld, int rmd) {
    la = rla;
    ma = rma;
    lb = rlb;
    mb = rmb;
    lc = rlc;
    mc = rmc;
    ld = rld;
    md = rmd;
    b = p = d = q = m = 0;
}

integralSph::integralSph(int rla, int rma, int rlb, int rmb, int rlc, int rmc, int rld, int rmd,
            int rb, int rp, int rd, int rq, int rm) {
    la = rla;
    ma = rma;
    lb = rlb;
    mb = rmb;
    lc = rlc;
    mc = rmc;
    ld = rld;
    md = rmd;
    b = rb;
    p = rp;
    d = rd;
    q = rq;
    m = rm;
}

integralSph::integralSph(int rla, int rma, int rlb, int rmb, int rlc, int rmc, int rld, int rmd,
            int rb, int rp, int rd, int rq, int rm, bool rab) {
    la = rla;
    ma = rma;
    lb = rlb;
    mb = rmb;
    lc = rlc;
    mc = rmc;
    ld = rld;
    md = rmd;
    b = rb;
    p = rp;
    d = rd;
    q = rq;
    m = rm;

    if (rab) {
        swap (la, lb);
        swap (ma, mb);
    }
}

integralSph::integralSph(int rla, int rma, int rlb, int rmb, int rlc, int rmc, int rld, int rmd,
            int rb, int rp, int rd, int rq, int rm, bool rab, bool rcdab) {
    la = rla;
    ma = rma;
    lb = rlb;
    mb = rmb;
    lc = rlc;
    mc = rmc;
    ld = rld;
    md = rmd;
    b = rb;
    p = rp;
    d = rd;
    q = rq;
    m = rm;

    if (rab) {
        swap (la, lb);
        swap (ma, mb);
    }

    if (rcdab) {
        swap(la,lc);
        swap(ma,mc);

        swap(lb,ld);
        swap(mb,md);

        swap(b, d);
        swap(p, q);
    }
}

integralSph integralSph::operator|(const integralSph & rhs) const {
    integralSph ret;
    ret = *this;
    ret.la *= rhs.la;
    ret.ma *= rhs.ma;
    ret.lb *= rhs.lb;
    ret.mb *= rhs.mb;
    ret.lc *= rhs.lc;
    ret.mc *= rhs.mc;
    ret.ld *= rhs.ld;
    ret.md *= rhs.md;

    ret.b *= rhs.b;
    ret.p *= rhs.p;
    ret.d *= rhs.d;
    ret.q *= rhs.q;
    ret.m *= rhs.m;
    return ret;
}

//in order to apply RRs in a more sensible way
integralSph & integralSph::operator+=(const integralSph & rhs) {
    la += rhs.la;
    ma += rhs.ma;
    lb += rhs.lb;
    mb += rhs.mb;
    lc += rhs.lc;
    mc += rhs.mc;
    ld += rhs.ld;
    md += rhs.md;

    b += rhs.b;
    p += rhs.p;
    d += rhs.d;
    q += rhs.q;
    m += rhs.m;
    return *this;
}

integralSph integralSph::operator+(const integralSph & rhs) const{
    integralSph r = *this;
    r += rhs;
    return r;
}

integralSph & integralSph::operator-=(const integralSph & rhs) {
    la -= rhs.la;
    ma -= rhs.ma;
    lb -= rhs.lb;
    mb -= rhs.mb;
    lc -= rhs.lc;
    mc -= rhs.mc;
    ld -= rhs.ld;
    md -= rhs.md;

    b -= rhs.b;
    p -= rhs.p;
    d -= rhs.d;
    q -= rhs.q;
    m -= rhs.m;
    return *this;
}

integralSph integralSph::operator-(const integralSph & rhs) const {
    integralSph r = *this;
    r -= rhs;
    return r;
}

bool integralSph::operator<(const integralSph & rhs) const {
    if (la!=rhs.la) return (la<rhs.la);
    if (ma!=rhs.ma) return (ma<rhs.ma);

    if (lb!=rhs.lb) return (lb<rhs.lb);
    if (mb!=rhs.mb) return (mb<rhs.mb);

    if (lc!=rhs.lc) return (lc<rhs.lc);
    if (mc!=rhs.mc) return (mc<rhs.mc);

    if (ld!=rhs.ld) return (ld<rhs.ld);
    if (md!=rhs.md) return (md<rhs.md);

    if (b!=rhs.b) return (b<rhs.b);
    if (p!=rhs.p) return (p<rhs.p);

    if (d!=rhs.d) return (d<rhs.d);
    if (q!=rhs.q) return (q<rhs.q);

    if (m!=rhs.m) return (m<rhs.m);

    return false;
}




ostream & operator<< (ostream & os, const ikernel & rhs) {
    os << int(rhs.u) << " " << int(rhs.v) << " " << int(rhs.s) << " " << int(rhs.t) << " " << int(rhs.m) << "  " << int(rhs.K4L) << " " << int(rhs.FE);
    return os;
}

ostream & operator<< (ostream & os, const integralCart & rhs) {
    os
     << int(rhs.ex) << int(rhs.ey) << int(rhs.ez) << " "
     << int(rhs.bx) << int(rhs.by) << int(rhs.bz) << " "
     << int(rhs.fx) << int(rhs.fy) << int(rhs.fz) << " "
     << int(rhs.dx) << int(rhs.dy) << int(rhs.dz) << "  "
     << int(rhs.px) << int(rhs.py) << int(rhs.pz) << " "
     << int(rhs.qx) << int(rhs.qy) << int(rhs.qz) << " "
     << int(rhs.rx) << " " << int(rhs.ry) << " " << int(rhs.rz) << "  "
     << int(rhs.a) << " " << int(rhs.b) << " " << int(rhs.p) << " "
     << int(rhs.c) << " " << int(rhs.d) << " " << int(rhs.q) << " "
     << int(rhs.m)  << "   "
     << int(rhs.la) << " " << int(rhs.lb) << " " << int(rhs.lc) << " " << int(rhs.ld);

    return os;
}

ostream & operator<< (ostream & os, const integralSph & rhs) {
    os
     << int(rhs.la) << " " << int(rhs.lb) << " " << int(rhs.lc) << " " << int(rhs.ld) << "  "
     << int(rhs.ma) << " " << int(rhs.mb) << " " << int(rhs.mc) << " " << int(rhs.md) << "  "
     << int(rhs.b) << " " << int(rhs.p) << " "
     << int(rhs.d) << " " << int(rhs.q) << " "
     << int(rhs.m);
    return os;
}


