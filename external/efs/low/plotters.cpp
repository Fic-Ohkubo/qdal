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


#include <iostream>
#include "plotters.hpp"
using namespace std;


MessagePlotter::MessagePlotter() {
    if (_NO_PLOT) return;
    os = &(std::cout);
    openfile = false;

    os->setf(std::ios::fixed);
    Precision = os->precision();

    Init();
}

MessagePlotter::~MessagePlotter() {
    if (_NO_PLOT) return;
    if (openfile) file.close();
}

void MessagePlotter::SetFile(const std::string & filename) {
    if (_NO_PLOT) return;
    file.open(filename.c_str());
    openfile = true;

    os = &(file);

    os->setf(std::ios::fixed);
    Precision = os->precision();
}

void MessagePlotter::Init() {
    if (_NO_PLOT) return;
    MaxDepth = 100;
    Depth = 0;
    status = true;
    Link = NULL;
    Width     = 6; //std::cout.width();
    LastChrono = new Chronometer;
    LastChrono->Start();
    newline = true;
    count = 0;
    LastPush = 0;
    freshlevel = true;
}

void MessagePlotter::SetMaxDepth(int iDepth) {
    if (_NO_PLOT) return;
    MaxDepth = iDepth;
}

void MessagePlotter::Push() {
    if (_NO_PLOT) return;
    if (!ThisNode.IsMaster()) return;

    os->flush();

    if (CheckStatus()) {
        LastPush = TotDepth();
        freshlevel = true;
        if (Link!=NULL) Link->LastPush = LastPush;

        Chronos.push(LastChrono);
        LastChrono = new Chronometer;
        LastChrono->Start();
    }

    ++Depth;

}


void MessagePlotter::Pop(bool same) {
    if (_NO_PLOT) return;
    if (!ThisNode.IsMaster()) return;

    --Depth;

    // nothing was printed between last push and this pop
    if (CheckStatus()) {
        same = (LastPush==TotDepth() && freshlevel);

        LastChrono->Stop();

        if (!same) {
            for (int i=0; i<TotDepth(); ++i) *os << "  ";
            *os << " Time: " << *LastChrono << endl;
        }
        else {
            for (int i=count; i<80; ++i)    *os << " ";
            *os << *LastChrono << endl;
        }
        newline = true;

        delete LastChrono;
        LastChrono = Chronos.top();
        Chronos.pop();
        freshlevel = false;
    }
}

//maybe it has to break
void MessagePlotter::NewLine() {
    if (_NO_PLOT) return;
    if (!ThisNode.IsMaster()) return;

    if (Link!=NULL && Depth==0 && Link->freshlevel) {
        Link->freshlevel ;
        freshlevel = true;
        LastPush = Link->LastPush;
    }

    // if last action was a push
    if (freshlevel) {
        *os << endl;
        newline = true;
    }

    freshlevel=false;


    if (newline) {
        count = TotDepth();
        for (int i=0; i<TotDepth(); ++i) (*os) << "  ";
        newline = false;
    }
}

int MessagePlotter::TotDepth() const {
    if (_NO_PLOT) return 0;
    int pDepth = Depth;
    if (Link!=NULL) pDepth += Link->TotDepth();

    return pDepth;
}

void MessagePlotter::precision(int p) {
    Precision = p;
}

void MessagePlotter::width(int p) {
    Width = p;
}

void MessagePlotter::Out(const char * c) const {
    if (_NO_PLOT) return;
    if (!ThisNode.IsMaster()) return;

    if (CheckStatus()) {
        *os << *c;
        (*os).flush();
    }
}

void MessagePlotter::SetModule(const std::string & rhs) {
    if (_NO_PLOT) return;
    module = rhs;
}

void MessagePlotter::SetParent(MessagePlotter * parent) {
    if (_NO_PLOT) return;
    Link = parent;
}

void MessagePlotter::Disable() {
    status = false;
}

void MessagePlotter::Enable() {
    status = true;
}

bool MessagePlotter::CheckStatus() const {
    if (_NO_PLOT) return true;
    bool ret = status;
    if (Link!=NULL) ret &= Link->CheckStatus();
    ret &= (TotDepth()<MaxDepth);
    return ret;
}




void MessagePlotter::Print(const std::string & out) {
    if (_NO_PLOT) return;
    if (!ThisNode.IsMaster()) return;

    NewLine();
    (*os) << out;

    count += out.length();
}

void MessagePlotter::Print(const char * out) {
    if (_NO_PLOT) return;
    if (!ThisNode.IsMaster()) return;

    std::string str = out;

    Print(str);
}

#include "float.h"

bool IsNumber(double x) {
    return (x == x);
}

bool IsFinite(double x) {
    return (x<=DBL_MAX && x>=-DBL_MAX);
}

void MessagePlotter::Print(double out) {
    if (_NO_PLOT) return;
    if (!ThisNode.IsMaster()) return;

    NewLine();

    int oPrecision = (*os).precision();
    int oWidth     = (*os).width();

    (*os).precision(Precision);

    //align, etc
    if (IsNumber(out) && IsFinite(out)) {
        double cc = out;
        if (cc>=0) (*os) << " "; //for the minus sign
        cc = fabs(cc);
        int w = Width;
        for (;cc>=1.; cc*=0.1) --w;
        if (fabs(out) < 1.0) --w;
        for (;w>0; --w) (*os) << " ";
    }

    (*os) << out;

    (*os).precision(oPrecision);
}

//solo para endl, que es template <typename T> T& endl (T&) y provoca llamadas ambiguas
MessagePlotter & MessagePlotter::operator<<(std::ostream & (*f)(std::ostream & ) ) {
    if (_NO_PLOT) return *this;
    if (!ThisNode.IsMaster()) return *this;

    if (CheckStatus()) {
        f( (*os) );
        newline = true;
        count = 0;
    }

    if (Link!=NULL && Link->os != os) *(Link) << f;

    return *this;
}

