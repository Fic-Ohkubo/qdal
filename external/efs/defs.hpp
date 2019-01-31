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


#ifndef __DEFS__ //se asegura de que no se dupliquen declaraciones
#define __DEFS__

// program limits
// **************

const unsigned long long int   MAX_TILES_PER_BLOCK = 4096; //1024;
const unsigned char             LMAX                =    5;  //up to h functions are permitted at the moment (still only up to g qic functions are present)
const unsigned char             maxK                =   32;  //maximum contraction degree allowed (for now)
const unsigned char             maxJ                =   16;  //maximum number of contraction functions allowed


// machine-specific constants
// **************************
#define CACHE_LINE_SIZE 64
#define AVXD_PER_CACHE_LINE    2
#define MM128_PER_CACHE_LINE  4
#define DOUBLES_PER_CACHE_LINE 8
#define FLOATS_PER_CACHE_LINE 16

typedef unsigned long long int  UI64;  //64 bits
typedef unsigned long int       UI32;  //32 bits
typedef unsigned short int      UI16;  //16 bits
typedef unsigned char            UI8;   // 8 bits
typedef           char            SI8;   // 8 bits, signed

typedef UI8   BYTE;   // minimum addressable memory unit
typedef UI64  PSIZE;  // integer with the same size of a pointer


// actual constants
// ****************

const UI64 Mword = 1024*1024;
const UI64 Gword = 1024*Mword;

const double PI   = 3.1415926535897932384626433832795; //M_PI is not compiler-independant

const double BOHR = 0.529177249240; // one bohr en amstrongs (10^-10 m)
const double AMS  = 1/BOHR;         // one amstrong in bohrs


// stuf which is computed from the previous definitions
// ****************************************************
const UI8    MMAX = (LMAX+1)*(LMAX+2)/2; //(h+1)*;
const UI8    LSP  = LMAX + 1; //code for SP GTOs

//in order to declare array sizes
const UI8    LM1  =  1*LMAX+1;
const UI8    LM2  =  2*LMAX+1;
const UI8    LM4  =  4*LMAX+1;
const UI8    LM6  =  6*LMAX+1; //needed
const UI8    LM8  =  8*LMAX+1; //needed
const UI8    LM12 = 12*LMAX+1; //needed


const UI16 maxK2 = maxK*maxK;
const UI32 maxK4 = maxK2*maxK2;
const UI16 maxJ2 = maxJ*maxJ;
const UI32 maxJ4 = maxJ2*maxJ2;

//enumerate atom 4-tuple geometries
enum GEOM {ABCD, AACD, ABAD, AACC, ABAB, AAAD, AAAA, NONE};

// do not attempt QR factorization
//#define __QR_GC__ // perform QR factorization of general contraction functions

#endif

