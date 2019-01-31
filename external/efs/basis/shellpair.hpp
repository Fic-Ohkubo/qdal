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


#ifndef __SHELL_PAIR__
#define __SHELL_PAIR__

#include "../defs.hpp"

class ShellPairPrototype;
class point;

struct ShellPair {
    const ShellPairPrototype * GP;                                  // 8 bytes
    double                   * WW;                                  // 8 bytes

	UI32 ata;       // number of atom A                              // 4 bytes
	UI32 atb;       // number of atom B                              // 4 bytes

	float logCS;   // Cauchy-Schwarz parameter logarithm            // 4 bytes

	UI16 nK2;       // number of significant gaussians of the pair   // 2 bytes
	UI16 flags;     // samefunctionm sameatom, invert prototype      // 2 bytes


	void Form (const ShellPairPrototype & prototype, const point & c1, const point & c2, double * ww, double lmo);
	void Form (const ShellPairPrototype & prototype, const point & c1, double * ww);

    double  getW(int a, int b) const;
    int     getKb() const;
    int     getKa(int Kb) const;
    int     getLa() const;
    int     getLb() const;
    int     getFa() const;
    int     getFb() const;
    bool    sameF() const;
    bool    sameAtom() const;
    bool    Inverted() const;

} __attribute__((aligned(32)));

#endif
