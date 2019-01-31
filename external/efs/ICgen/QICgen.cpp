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




    Interpreted code and Inner Contraction routines generator for the Echidna Fock Solver

    Jaime Axel Rosal Sandberg, August 2013

    Input required:
    L      : maximum angular momentum
    QIC_DIR: directory where interpreted code files are to be stored
    EFS_DIR: directory where the specialized routines for the Echidna Fock Solver will be written

*/


#include <iostream>
#include <queue>
#include <string>
#include "../2eints/IICinit.hpp"
#include "../2eints/quimera.hpp"
#include "../math/angular.hpp"


#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

using namespace std;
using namespace LibAngular;


std::string QIC_DIR;
std::string K2C_DIR;


int main(int argc, char* argv[]) {

    if (argc<2) {
        cout << "Angular momentum not specified" << endl;
        return 1;
	}

	char ll = argv[1][0];

	if ( ll<'0' || ll>'4') {
        cout << "Error: angular momentum out of bounds (0-4)" << endl;
        return 1;
	}

	int L = int(ll-'0');


    //chek for QICDIR; issue a message if not set
	{
        const char * val =  argv[2]; //::getenv("QIC_DIR");
        QIC_DIR = val;
        if (val==0) cout << "Directory for QIC files not specified" << endl;
        else        cout << "QIC directory:   " << string(val) << endl;
    }

    //chek for EFSDIR; issue a message if not set
	{
        const char * val = argv[3]; //::getenv("EFS_DIR");
        K2C_DIR = val;

        if (val==0) cout << "Directory for K2C files not specified" << endl;
        else        cout << "K2C directory:   " << string(val) << endl;
    }

    InitCartList();
    InitSHList();



    std::priority_queue <ERItype> SetLater;

    for (int la=0; la<=L; ++la) {
        for (int lb=0; lb<=la; ++lb) {
            for (int lc=0; lc<=L; ++lc) {
                for (int ld=0; ld<=lc; ++ld) {

                    ERItype type;

                    // 4-center and degenerate 4-center integrals
                    if (la>lc || la==lc && lb>=ld) {
                        type(ABCD, la, lb, lc, ld, true, false, false);
                        SetLater.push(type);

                        type(ABAD, la, lb, lc, ld, true, false, false);
                        SetLater.push(type);

                        type(ABAB, la, lb, lc, ld, true, false, false);
                        SetLater.push(type);
                    }

                    // 3-center integrals
                    type(AACD, la, lb, lc, ld, true, false, false);
                    SetLater.push(type);

                    // 2-center integrals
                    if (la>lc || la==lc && lb>=ld) {
                        type(AACC, la, lb, lc, ld, true, false, false);
                        SetLater.push(type);
                    }

                    // 1-center integrals
                    if ((la+lb+lc+ld)%2)                continue;
                    if ((lb==0) && (ld==0) && (la!=lc)) continue;

                    type(AAAA, la, lb, lc, ld, true, false, false);
                    SetLater.push(type);
                }
            }
        }
    }


    cout << "Writing interpreted code routines" << endl;

    int n = SetLater.size();

    #pragma omp parallel for schedule(dynamic)
    for (UI32 i=0; i<n; ++i) {

        ERItype type;

        #pragma omp critical
        {
            type = SetLater.top();
            SetLater.pop();
            cout << ".";
            cout.flush();
        }

        ERIroutine IC;

        IC.Set(type.La, type.Lb, type.Lc, type.Ld, type.geometry, type.isCDR);
        IC.Write();
    }
    cout << endl;




    cout << "Writing headers and loading routines" << endl;

    {
        string k2cfile= K2C_DIR + "/K2C.cpp";
        ofstream file;

        file.open(k2cfile.c_str());
        file.precision(16);


        //file << "#ifndef __K2C__" << endl;
        //file << "#define __K2C__" << endl;
        //file << endl;
        file << "class cacheline64;" << endl;
        file << "class ERIgeometries64;" << endl;
        file << "class PrimitiveSet;" << endl;
        file << "class ShellPairPrototype;" << endl;
        file << "class p_ERIbuffer;" << endl;
        file << endl;

        //
        for (int la=0; la<=L; ++la) {
            for (int lb=0; lb<=la; ++lb) {
                for (int lcd=0; lcd<=la+lb; ++lcd) {
                    file << "void K2C_ABCD_" << L2S(la) << L2S(lb) << L2S(lcd) <<"(cacheline64 * (&F0), const ERIgeometries64 & vars8, const PrimitiveSet & PSab, const ShellPairPrototype & ABp, double ikcd, double rcd, p_ERIbuffer & buffer, cacheline64 & AB2, cacheline64 & X2);" << endl;
                }
            }
        }

        file << endl;
        file << "#include \"2eints/IICinit.hpp\"" << endl;
        file << endl;

        file << "void ERIroutine::LoadK2C() {" << endl;

        file << "  int lcd = lc+ld;"  << endl;
        file << "  if (geometry==ABCD) {" << endl;


        for (int la=0; la<=L; ++la) {
            for (int lb=0; lb<=la; ++lb) {
                for (int lcd=0; lcd<=la+lb; ++lcd) {
                    file << "    if (la=="<<la<<" && lb=="<<lb<<" && lcd=="<<lcd<<") InnerContractionRoutine = K2C_ABCD_"<<L2S(la)<<L2S(lb)<<L2S(lcd)<<";" << endl;
                }
            }
        }

        file << "  }" << endl;
        file << "}" << endl;


        //file << "#endif" << endl;
        file << endl;

    }


    return 0;
}
