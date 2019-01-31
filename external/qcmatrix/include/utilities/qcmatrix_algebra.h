/* QcMatrix: an abstract matrix library
   Copyright 2012-2015 Bin Gao

   QcMatrix is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   QcMatrix is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with QcMatrix. If not, see <http://www.gnu.org/licenses/>.

   This file defines some basic algebraic functions.

   2012-04-04, Bin Gao:
   * first version
*/

#if !defined(QCMATRIX_ALGEBRA_H)
#define QCMATRIX_ALGEBRA_H

/* uses the function rand() */
#include <stdlib.h>
/* uses QReal */
#include "types/qcmatrix_basic_types.h"

/* absolute value of a number */
#define QAbs(a) (((a)>=0) ? (a) : -(a))
/* maximum of two numbers */
#define QMax(a,b) (((a)>b) ? (a) : (b))
/* minimum of two numbers */
#define QMin(a,b) (((a)<b) ? (a) : (b))
/* random integer on [a,b] */
#define QRandInt(a,b) ((rand()%(b-a+1))+a)
/* random real number on [a,b] */
#define QRandReal(a,b) ((rand()/(QReal)RAND_MAX)*(b-a)+a)

#endif
