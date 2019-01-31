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

#ifndef __AFFINE__
#define __AFFINE__

#include <cmath>
#include <iostream>


class vector3;
class point;


class vector3 {
    friend class point;


  public:
    double x;
    double y;
    double z;

    vector3() {}

    vector3(double xx, double yy, double zz) {
        x = xx;
        y = yy;
        z = zz;
    }

    inline void operator()(double xx, double yy, double zz) {
        x = xx;
        y = yy;
        z = zz;
    }

    inline vector3 & operator=(const vector3 & rhs) {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    inline vector3 & operator+=(const vector3 & rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    };

    inline vector3 & operator-=(const vector3 & rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    };

    inline vector3 & operator*=(double rhs) {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    };

    inline vector3 & operator/=(double rhs) {
        double irhs = 1./rhs;
        x *= irhs;
        y *= irhs;
        z *= irhs;
        return *this;
    };

    inline const vector3 operator+(const vector3 & rhs) const {
        vector3 v = *this;
        v += rhs;
        return v;
    };

    inline const vector3 operator-(const vector3 & rhs) const {
        vector3 v = *this;
        v -= rhs;
        return v;
    };

    inline const vector3 operator*(double rhs) const {
        vector3 v = *this;
        v *= rhs;
        return v;
    };

    inline const vector3 operator/(double rhs) const {
        vector3 v = *this;
        v /= rhs;
        return v;
    };

    //inner product
    inline double operator*(const vector3 & rhs) const {
        double r;
        r  = x*rhs.x;
        r += y*rhs.y;
        r += z*rhs.z;

        return r;
    };

    //cross product (operator* cannot be overloaded for different return types)
    inline const vector3 operator|(const vector3 & rhs) const {
        vector3 v;
        v.x = y*rhs.z - z*rhs.y;
        v.y = z*rhs.x - x*rhs.z;
        v.z = x*rhs.y - y*rhs.x;

        return v;
    };

};

class point {
  public:
    double x;
    double y;
    double z;

  public:
    point() {
        x = y = z = 0;
    }

    point(double ix, double iy, double iz) {
        x = ix;
        y = iy;
        z = iz;
    }

    inline point & operator=(const point & rhs) {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    inline point & operator+=(const vector3 & rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    inline point & operator-=(const vector3 & rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }


    inline const point operator+(const vector3 & rhs) const {
        point v = *this;
        v += rhs;
        return v;
    }

    inline const point operator-(const vector3 & rhs) const {
        point v = *this;
        v -= rhs;
        return v;
    };

    inline const vector3 operator-(const point & rhs) const {
        vector3 r;
        r.x = x - rhs.x;
        r.y = y - rhs.y;
        r.z = z - rhs.z;

        return r;
    }
};

inline vector3 Xprod (const vector3 & v1, const vector3 & v2) {
	vector3 vt;
	vt.x = v1.y*v2.z - v1.z*v2.y;
	vt.y = v1.z*v2.x - v1.x*v2.z;
	vt.z = v1.x*v2.y - v1.y*v2.x;
	return vt;
}

inline const vector3 operator* (double d, const vector3 & v) {
    vector3 r;
    r.x = d*v.x;
    r.y = d*v.y;
    r.z = d*v.z;

    return r;
}

inline double norm(const vector3 & v) {
    return sqrt(v*v);
}

inline void normalize(vector3 & v) {
    v /= norm(v);
}

inline vector3 FrameRot (const vector3 & x, const vector3 & y, const vector3 & z, const vector3 & v) {
    vector3 r;
    r.x = v*x;
    r.y = v*y;
    r.z = v*z;
    return r;
}

inline vector3 FrameRotT(const vector3 & x, const vector3 & y, const vector3 & z, const vector3 & v) {
    vector3 r;
    r.x = v.x*x.x + v.y*y.x + v.z*z.x;
    r.y = v.x*x.y + v.y*y.y + v.z*z.y;
    r.z = v.x*x.z + v.y*y.z + v.z*z.z;
    return r;
}


#endif
