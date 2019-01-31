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


#ifndef __PLOTTERS__
#define __PLOTTERS__

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stack>
#include <cmath>
#include "../low/chrono.hpp"
#include "../low/MPIwrap.hpp"

//skip printing messages in DALTON
#ifdef PRG_DALTON
#define _NO_PLOT true
#else
#define _NO_PLOT false
#endif


class MessagePlotter {
    template<class T>
    friend MessagePlotter & operator<< (MessagePlotter & MP, const T & out);
  protected:
    std::string module;
    MessagePlotter * Link;
    std::stack<Chronometer*> Chronos;
    Chronometer * LastChrono;
    std::ostream * os;
    std::ofstream file;
    bool openfile;

    int Depth;
    int MaxDepth;
    int Precision;
    int Width;

    int count;
    int LastPush;

    bool status;
    bool newline;
    bool freshlevel;

  public:

    MessagePlotter();
    //MessagePlotter(const std::string & filename);
    ~MessagePlotter();
    void SetFile(const std::string & filename);
    void Init();

    void SetMaxDepth(int iDepth);
    void Push();
    void Pop(bool same=false);
    int TotDepth() const;
    void precision(int p);
    void width(int p);
    void Out(const char * c) const;
    void SetModule(const std::string & rhs);
    void SetParent(MessagePlotter * parent);
    void Disable();
    void Enable();
    bool CheckStatus() const;

    void NewLine();

    template<class T> void Print(const T & out) {
        if (_NO_PLOT) return;

        NewLine();
        if (newline) for (int i=0; i<TotDepth(); ++i) {(*os) << "  "; ++count;} newline = false;

        int oPrecision = (*os).precision();
        int oWidth     = (*os).width();

        (*os).precision(Precision);
        (*os).width(Width);

        (*os) << out;

        (*os).precision(oPrecision);
        (*os).width(oWidth);
    }

    void Print(const std::string & out);
    void Print(const char * out);
    void Print(double out);

    MessagePlotter & operator<<(std::ostream & (*f)(std::ostream & ) );
};

template<class T> MessagePlotter & operator<< (MessagePlotter & MP, const T & out) {
    if (_NO_PLOT) return MP;
    // only print messages to screen and/or file if it is the master MPI process
    if (ThisNode.IsMaster()) {
        if (MP.CheckStatus())
            MP.Print(out);

        if (MP.Link!=NULL && MP.Link->os != MP.os) *(MP.Link) << out;
    }

    return MP;
}

#endif
