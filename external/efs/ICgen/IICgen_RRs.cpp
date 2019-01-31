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


#include <fstream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <cstdlib>
#include <string>
#include <queue>

#include "IICgen.hpp"

#include "../math/angular.hpp"
#include "../low/chrono.hpp"
#include "../low/cache.hpp"

using namespace std;
using namespace LibAngular;

enum CENTER {AAA, BBB, CCC, DDD};

static inline UI32 min(UI32 a, UI32 b) {
    return (a<b)?a:b;
}

static inline UI32 max(UI32 a, UI32 b) {
    return (a>b)?a:b;
}



//for CDR/AERR routines
void EvaluationScheme::SieveK0K1(const std::set<ikernel> & setK0, std::set<ikernel> & setK1, bool ABt, bool CDt, bool ACt) {

    //                   u v  s t  m
    const ikernel FR1   (0,0, 0,0, 1);
    const ikernel FR2   (1,0, 0,0, 1);
    const ikernel FR3   (0,0, 1,0, 1);
    const ikernel FR4   (2,0, 0,0, 1);
    const ikernel FR5   (1,0, 1,0, 1);
    const ikernel FR6   (0,0, 2,0, 1);

    const ikernel fmask (0,1, 0,1, 0);


    const ikernel ER1   (0,1, 0,0, 1);
    const ikernel ER2   (0,0, 0,1, 1);

    const ikernel emask (1,0, 1,0, 0);


    set<ikernel> sete;

    //CDR
    //===

    UI8 totcnt = 0;

    if (ACt)        ++totcnt;
    if (ABt && ACt) ++totcnt;
    if (CDt && ACt) ++totcnt;
    if (ABt)        ++totcnt;
    if (ABt && CDt) ++totcnt;
    if (CDt)        ++totcnt;


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > fgroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setK0.begin(); it!=setK0.end(); ++it) {
        ikernel index = *it | fmask;
        fgroupset[index].insert(*it);
    }



    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=fgroupset.begin(); git!=fgroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            //count number of matching parent kernels in the set
            UI8 cnt = 0;

            if (ACt)        if (git->second.count(R+FR1) == 1) ++cnt;
            if (ABt && ACt) if (git->second.count(R+FR2) == 1) ++cnt;
            if (CDt && ACt) if (git->second.count(R+FR3) == 1) ++cnt;
            if (ABt)        if (git->second.count(R+FR4) == 1) ++cnt;
            if (ABt && CDt) if (git->second.count(R+FR5) == 1) ++cnt;
            if (CDt)        if (git->second.count(R+FR6) == 1) ++cnt;


            //iff the kernel can be computed from the ones present in the group, register the contracted exponential and
            if (cnt==totcnt) {
                if (ACt)        imap[R].VD1 = R+FR1;
                if (ABt && ACt) imap[R].VD2 = R+FR2;
                if (CDt && ACt) imap[R].VD3 = R+FR3;
                if (ABt)        imap[R].VD4 = R+FR4;
                if (ABt && CDt) imap[R].VD5 = R+FR5;
                if (CDt)        imap[R].VD6 = R+FR6;

                ikernel E = R;
                E.FE = 1;
                imap[R].VDE = E;
                imap[R].ETYPE = CDR0;
                imap[R].aux   = R.m;

                sete.insert(E);
            }
            //otherwise register the integral as necessary
            else {
                imap[R].VD1 = R;
                setK1.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }


    //EARR
    //====


    map<ikernel, set<ikernel> > egroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=sete.begin(); it!=sete.end(); ++it) {
        ikernel index = *it | emask;
        egroupset[index].insert(*it);
    }


    //loop for each group
    for (git=egroupset.begin(); git!=egroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            UI8 cnt = 0;

            if (git->second.count(R+ER1)==1) ++cnt;
            if (git->second.count(R+ER2)==1) ++cnt;

            if (cnt==2) {
                imap[R].VD1 = R+ER1;
                imap[R].VD2 = R+ER2;

                imap[R].ETYPE = AERR0;
            }
            else {
                imap[R].VD1 = R;
                setK1.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }



}

void EvaluationScheme::SieveK1K2(const std::set<ikernel> & setK1, std::set<ikernel> & setK2, bool ABt, bool CDt, bool ACt) {

    //                    u v  0 q  m
    const ikernel FR1   (0,0, 0,0, 1);
    const ikernel FR2   (1,0, 0,0, 1);
    const ikernel FR3   (0,0, 0,1, 1);
    const ikernel FR4   (2,0, 0,0, 1);
    const ikernel FR5   (1,0, 0,1, 1);
    const ikernel FR6   (0,0, 0,2, 1);

    const ikernel fmask (0,1, 1,0, 0);


    const ikernel ER1   (0,1, 0,0, 1);
    const ikernel ER2   (0,0, 0,1, 1);

    const ikernel emask (1,0, 1,0, 0);


    set<ikernel> sete;

    //CDR
    //===


    UI8 totcnt = 0;

    if (ACt)        ++totcnt;
    if (ABt && ACt) ++totcnt;
    if (CDt && ACt) ++totcnt;
    if (ABt)        ++totcnt;
    if (ABt && CDt) ++totcnt;
    if (CDt)        ++totcnt;


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > fgroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setK1.begin(); it!=setK1.end(); ++it) {
        ikernel index = *it | fmask;
        fgroupset[index].insert(*it);
    }



    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=fgroupset.begin(); git!=fgroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            //count number of matching parent kernels in the set
            UI8 cnt = 0;

            if (ACt)        if (git->second.count(R+FR1) == 1) ++cnt;
            if (ABt && ACt) if (git->second.count(R+FR2) == 1) ++cnt;
            if (CDt && ACt) if (git->second.count(R+FR3) == 1) ++cnt;
            if (ABt)        if (git->second.count(R+FR4) == 1) ++cnt;
            if (ABt && CDt) if (git->second.count(R+FR5) == 1) ++cnt;
            if (CDt)        if (git->second.count(R+FR6) == 1) ++cnt;


            //iff the kernel can be computed from the ones present in the group, register the contracted exponential and
            if (cnt==totcnt) {
                if (ACt)        imap[R].VD1 = R+FR1;
                if (ABt && ACt) imap[R].VD2 = R+FR2;
                if (CDt && ACt) imap[R].VD3 = R+FR3;
                if (ABt)        imap[R].VD4 = R+FR4;
                if (ABt && CDt) imap[R].VD5 = R+FR5;
                if (CDt)        imap[R].VD6 = R+FR6;

                ikernel E = R;
                E.FE = 1;
                imap[R].VDE = E;
                imap[R].ETYPE = CDR1;
                imap[R].aux   = R.m;

                sete.insert(E);
            }
            //otherwise register the integral as necessary
            else {
                imap[R].VD1 = R;
                setK2.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }


    //EARR
    //====

    map<ikernel, set<ikernel> > egroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=sete.begin(); it!=sete.end(); ++it) {
        ikernel index = *it | emask;
        egroupset[index].insert(*it);
    }


    //loop for each group
    for (git=egroupset.begin(); git!=egroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            UI8 cnt = 0;

            if (git->second.count(R+ER1)==1) ++cnt;
            if (git->second.count(R+ER2)==1) ++cnt;

            if (cnt==2) {
                imap[R].VD1 = R+ER1;
                imap[R].VD2 = R+ER2;

                imap[R].ETYPE = AERR1;
            }
            else {
                imap[R].VD1 = R;
                setK2.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }
}

void EvaluationScheme::SieveK2K3(const std::set<ikernel> & setK2, std::set<ikernel> & setK3, bool ABt, bool CDt, bool ACt) {

    //                    u v  0 0  m
    const ikernel FR1   (0,0, 0,0, 1);
    const ikernel FR2   (1,0, 0,0, 1);
    const ikernel FR3   (2,0, 0,0, 1);

    const ikernel fmask (0,1, 1,1, 0);


    const ikernel ER1   (0,1, 0,0, 1);
    const ikernel ER2   (0,0, 0,0, 1);

    const ikernel emask (1,0, 1,1, 0);


    set<ikernel> sete;

    //CDR
    //===


    UI8 totcnt = 0;

    bool AQt = ACt || CDt;

    if (AQt)        ++totcnt;
    if (AQt && ABt) ++totcnt;
    if (ABt)        ++totcnt;


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > fgroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setK2.begin(); it!=setK2.end(); ++it) {
        ikernel index = *it | fmask;
        fgroupset[index].insert(*it);
    }



    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=fgroupset.begin(); git!=fgroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            //count number of matching parent kernels in the set
            UI8 cnt = 0;

            if (AQt)        if (git->second.count(R+FR1) == 1) ++cnt;
            if (ABt && AQt) if (git->second.count(R+FR2) == 1) ++cnt;
            if (ABt)        if (git->second.count(R+FR3) == 1) ++cnt;

            //iff the kernel can be computed from the ones present in the group, register the contracted exponential and
            if (cnt==totcnt) {
                if (AQt)        imap[R].VD1 = R+FR1;
                if (ABt && AQt) imap[R].VD2 = R+FR2;
                if (ABt)        imap[R].VD3 = R+FR3;

                ikernel E = R;
                E.FE = 1;
                imap[R].VDE = E;
                imap[R].ETYPE = CDR2;
                imap[R].aux   = R.m;

                sete.insert(E);
            }
            //otherwise register the integral as necessary
            else {
                imap[R].VD1 = R;
                setK3.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }


    //EARR
    //====

    map<ikernel, set<ikernel> > egroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=sete.begin(); it!=sete.end(); ++it) {
        ikernel index = *it | emask;
        egroupset[index].insert(*it);
    }


    //loop for each group
    for (git=egroupset.begin(); git!=egroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            UI8 cnt = 0;

            if (git->second.count(R+ER1)==1) ++cnt;
            if (git->second.count(R+ER2)==1) ++cnt;

            if (cnt==2) {
                imap[R].VD1 = R+ER1;
                imap[R].VD2 = R+ER2;

                imap[R].ETYPE = AERR2;
            }
            else {
                imap[R].VD1 = R;
                setK3.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }

}

void EvaluationScheme::SieveK3K4(const std::set<ikernel> & setK3, std::set<ikernel> & setK4, bool ABt, bool CDt, bool ACt) {

    //                   0 p  0 0  m
    const ikernel FR1   (0,0, 0,0, 1);
    const ikernel FR2   (0,1, 0,0, 1);
    const ikernel FR3   (0,2, 0,0, 1);

    const ikernel fmask (1,0, 1,1, 0);


    const ikernel ER1   (0,1, 0,0, 1);
    const ikernel ER2   (0,0, 0,0, 1);

    const ikernel emask (1,0, 1,1, 0);


    set<ikernel> sete;


    //CDR
    //===

    UI8 totcnt = 0;

    bool AQt = ACt || CDt;

    if (AQt)        ++totcnt;
    if (AQt && ABt) ++totcnt;
    if (ABt)        ++totcnt;


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > fgroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setK3.begin(); it!=setK3.end(); ++it) {
        ikernel index = *it | fmask;
        fgroupset[index].insert(*it);
    }



    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=fgroupset.begin(); git!=fgroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            //count number of matching parent kernels in the set
            UI8 cnt = 0;

            if (AQt)        if (git->second.count(R+FR1) == 1) ++cnt;
            if (ABt && AQt) if (git->second.count(R+FR2) == 1) ++cnt;
            if (ABt)        if (git->second.count(R+FR3) == 1) ++cnt;

            //iff the kernel can be computed from the ones present in the group, register the contracted exponential and
            if (cnt==totcnt) {
                if (AQt)        imap[R].VD1 = R+FR1;
                if (ABt && AQt) imap[R].VD2 = R+FR2;
                if (ABt)        imap[R].VD3 = R+FR3;

                ikernel E = R;
                E.FE = 1;
                imap[R].VDE = E;
                imap[R].ETYPE = CDR3;
                imap[R].aux   = R.m;

                sete.insert(E);
            }
            //otherwise register the integral as necessary
            else {
                imap[R].VD1 = R;
                setK4.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }


    //EARR
    //====

    map<ikernel, set<ikernel> > egroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=sete.begin(); it!=sete.end(); ++it) {
        ikernel index = *it | emask;
        egroupset[index].insert(*it);
    }


    //loop for each group
    for (git=egroupset.begin(); git!=egroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            UI8 cnt = 0;

            if (git->second.count(R+ER1)==1) ++cnt;
            if (git->second.count(R+ER2)==1) ++cnt;

            if (cnt==2) {
                imap[R].VD1 = R+ER1;
                imap[R].VD2 = R+ER2;

                imap[R].ETYPE = AERR3;
            }
            else {
                imap[R].VD1 = R;
                setK4.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }

    }

}





//for CDR/AERR routines
void EvaluationScheme::UseCDR0(const std::set<ikernel> & setF0, std::set<ikernel> & setE0, std::set<ikernel> & setF1, bool ABt, bool CDt, bool ACt) {

    //                   u v  s t  m
    const ikernel FR1   (0,0, 0,0, 1);
    const ikernel FR2   (1,0, 0,0, 1);
    const ikernel FR3   (0,0, 1,0, 1);
    const ikernel FR4   (2,0, 0,0, 1);
    const ikernel FR5   (1,0, 1,0, 1);
    const ikernel FR6   (0,0, 2,0, 1);

    const ikernel fmask (0,1, 0,1, 0);


    UI8 totcnt = 0;

    if (ACt)        ++totcnt;
    if (ABt && ACt) ++totcnt;
    if (CDt && ACt) ++totcnt;
    if (ABt)        ++totcnt;
    if (ABt && CDt) ++totcnt;
    if (CDt)        ++totcnt;


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > fgroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setF0.begin(); it!=setF0.end(); ++it) {
        ikernel index = *it | fmask;
        fgroupset[index].insert(*it);
    }



    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=fgroupset.begin(); git!=fgroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            //count number of matching parent kernels in the set
            UI8 cnt = 0;

            if (ACt)        if (git->second.count(R+FR1) == 1) ++cnt;
            if (ABt && ACt) if (git->second.count(R+FR2) == 1) ++cnt;
            if (CDt && ACt) if (git->second.count(R+FR3) == 1) ++cnt;
            if (ABt)        if (git->second.count(R+FR4) == 1) ++cnt;
            if (ABt && CDt) if (git->second.count(R+FR5) == 1) ++cnt;
            if (CDt)        if (git->second.count(R+FR6) == 1) ++cnt;


            //iff the kernel can be computed from the ones present in the group, register the contracted exponential and
            if (cnt==totcnt) {
                if (ACt)        imap[R].VD1 = R+FR1;
                if (ABt && ACt) imap[R].VD2 = R+FR2;
                if (CDt && ACt) imap[R].VD3 = R+FR3;
                if (ABt)        imap[R].VD4 = R+FR4;
                if (ABt && CDt) imap[R].VD5 = R+FR5;
                if (CDt)        imap[R].VD6 = R+FR6;

                ikernel E = R;
                E.FE = 1;
                imap[R].VDE = E;
                imap[R].ETYPE = CDR0;
                imap[R].aux   = R.m;

                setE0.insert(E);
            }
            //otherwise register the integral as necessary
            else {
                imap[R].VD1 = R;
                setF1.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }
}

void EvaluationScheme::UseAERR0(const std::set<ikernel> & setE0, std::set<ikernel> & setE1) {

    const ikernel ER1   (0,1, 0,0, 1);
    const ikernel ER2   (0,0, 0,1, 1);

    const ikernel emask (1,0, 1,0, 0);


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > egroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setE0.begin(); it!=setE0.end(); ++it) {
        ikernel index = *it | emask;
        egroupset[index].insert(*it);
    }

    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=egroupset.begin(); git!=egroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            UI8 cnt = 0;

            if (git->second.count(R+ER1)==1) ++cnt;
            if (git->second.count(R+ER2)==1) ++cnt;

            if (cnt==2) {
                imap[R].VD1 = R+ER1;
                imap[R].VD2 = R+ER2;

                imap[R].ETYPE = AERR0;
            }
            else {
                imap[R].VD1 = R;
                setE1.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }
}


void EvaluationScheme::UseCDR1(const std::set<ikernel> & setF1, std::set<ikernel> & setE1, std::set<ikernel> & setF2, bool ABt, bool CDt, bool ACt) {

    //                    u v  0 q  m
    const ikernel FR1   (0,0, 0,0, 1);
    const ikernel FR2   (1,0, 0,0, 1);
    const ikernel FR3   (0,0, 0,1, 1);
    const ikernel FR4   (2,0, 0,0, 1);
    const ikernel FR5   (1,0, 0,1, 1);
    const ikernel FR6   (0,0, 0,2, 1);

    const ikernel fmask (0,1, 1,0, 0);

    //CDR
    //===

    UI8 totcnt = 0;

    if (ACt)        ++totcnt;
    if (ABt && ACt) ++totcnt;
    if (CDt && ACt) ++totcnt;
    if (ABt)        ++totcnt;
    if (ABt && CDt) ++totcnt;
    if (CDt)        ++totcnt;


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > fgroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setF1.begin(); it!=setF1.end(); ++it) {
        ikernel index = *it | fmask;
        fgroupset[index].insert(*it);
    }


    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=fgroupset.begin(); git!=fgroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            //count number of matching parent kernels in the set
            UI8 cnt = 0;

            if (ACt)        if (git->second.count(R+FR1) == 1) ++cnt;
            if (ABt && ACt) if (git->second.count(R+FR2) == 1) ++cnt;
            if (CDt && ACt) if (git->second.count(R+FR3) == 1) ++cnt;
            if (ABt)        if (git->second.count(R+FR4) == 1) ++cnt;
            if (ABt && CDt) if (git->second.count(R+FR5) == 1) ++cnt;
            if (CDt)        if (git->second.count(R+FR6) == 1) ++cnt;


            //iff the kernel can be computed from the ones present in the group, register the contracted exponential and
            if (cnt==totcnt) {
                if (ACt)        imap[R].VD1 = R+FR1;
                if (ABt && ACt) imap[R].VD2 = R+FR2;
                if (CDt && ACt) imap[R].VD3 = R+FR3;
                if (ABt)        imap[R].VD4 = R+FR4;
                if (ABt && CDt) imap[R].VD5 = R+FR5;
                if (CDt)        imap[R].VD6 = R+FR6;

                ikernel E = R;
                E.FE = 1;
                imap[R].VDE = E;
                imap[R].ETYPE = CDR1;
                imap[R].aux   = R.m;

                setE1.insert(E);
            }
            //otherwise register the integral as necessary
            else {
                imap[R].VD1 = R;
                setF2.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }
}

void EvaluationScheme::UseAERR1(const std::set<ikernel> & setE1, std::set<ikernel> & setE2) {

    const ikernel ER1   (0,1, 0,0, 1);
    const ikernel ER2   (0,0, 0,1, 1);

    const ikernel emask (1,0, 1,0, 0);

    //EARR
    //====

    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > egroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setE1.begin(); it!=setE1.end(); ++it) {
        ikernel index = *it | emask;
        egroupset[index].insert(*it);
    }

    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=egroupset.begin(); git!=egroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            UI8 cnt = 0;

            if (git->second.count(R+ER1)==1) ++cnt;
            if (git->second.count(R+ER2)==1) ++cnt;

            if (cnt==2) {
                imap[R].VD1 = R+ER1;
                imap[R].VD2 = R+ER2;

                imap[R].ETYPE = AERR1;
            }
            else {
                imap[R].VD1 = R;
                setE2.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }
}


void EvaluationScheme::UseCDR2(const std::set<ikernel> & setF2, std::set<ikernel> & setE2, std::set<ikernel> & setF3, bool ABt, bool CDt, bool ACt) {

    //                    u v  0 0  m
    const ikernel FR1   (0,0, 0,0, 1);
    const ikernel FR2   (1,0, 0,0, 1);
    const ikernel FR3   (2,0, 0,0, 1);

    const ikernel fmask (0,1, 1,1, 0);

    //CDR
    //===


    UI8 totcnt = 0;

    bool AQt = ACt || CDt;

    if (AQt)        ++totcnt;
    if (AQt && ABt) ++totcnt;
    if (ABt)        ++totcnt;


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > fgroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setF2.begin(); it!=setF2.end(); ++it) {
        ikernel index = *it | fmask;
        fgroupset[index].insert(*it);
    }



    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=fgroupset.begin(); git!=fgroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            //count number of matching parent kernels in the set
            UI8 cnt = 0;

            if (AQt)        if (git->second.count(R+FR1) == 1) ++cnt;
            if (ABt && AQt) if (git->second.count(R+FR2) == 1) ++cnt;
            if (ABt)        if (git->second.count(R+FR3) == 1) ++cnt;

            //iff the kernel can be computed from the ones present in the group, register the contracted exponential and
            if (cnt==totcnt) {
                if (AQt)        imap[R].VD1 = R+FR1;
                if (ABt && AQt) imap[R].VD2 = R+FR2;
                if (ABt)        imap[R].VD3 = R+FR3;

                ikernel E = R;
                E.FE = 1;
                imap[R].VDE = E;
                imap[R].ETYPE = CDR2;
                imap[R].aux   = R.m;

                setE2.insert(E);
            }
            //otherwise register the integral as necessary
            else {
                imap[R].VD1 = R;
                setF3.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }

}

void EvaluationScheme::UseAERR2(const std::set<ikernel> & setE2, std::set<ikernel> & setE3) {

    const ikernel ER1   (0,1, 0,0, 1);
    const ikernel ER2   (0,0, 0,0, 1);

    const ikernel emask (1,0, 1,1, 0);

    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > egroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setE2.begin(); it!=setE2.end(); ++it) {
        ikernel index = *it | emask;
        egroupset[index].insert(*it);
    }

    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=egroupset.begin(); git!=egroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            UI8 cnt = 0;

            if (git->second.count(R+ER1)==1) ++cnt;
            if (git->second.count(R+ER2)==1) ++cnt;

            if (cnt==2) {
                imap[R].VD1 = R+ER1;
                imap[R].VD2 = R+ER2;

                imap[R].ETYPE = AERR2;
            }
            else {
                imap[R].VD1 = R;
                setE3.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }
}


// given a set of K3 integrals, generates a set with the K3 contracted exponentials needed and a set with the integrals that have to be accumulated in K4
void EvaluationScheme::UseCDR3 (const std::set<ikernel> & setF3, std::set<ikernel> & setE3, std::set<ikernel> & setF4, bool ABt, bool CDt, bool ACt) {

    //                   0 p  0 0  m
    const ikernel FR1   (0,0, 0,0, 1);
    const ikernel FR2   (0,1, 0,0, 1);
    const ikernel FR3   (0,2, 0,0, 1);

    const ikernel fmask (1,0, 1,1, 0);


    //CDR
    //===

    UI8 totcnt = 0;

    bool AQt = ACt || CDt;

    if (AQt)        ++totcnt;
    if (AQt && ABt) ++totcnt;
    if (ABt)        ++totcnt;


    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > fgroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setF3.begin(); it!=setF3.end(); ++it) {
        ikernel index = *it | fmask;
        fgroupset[index].insert(*it);
    }



    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=fgroupset.begin(); git!=fgroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            //count number of matching parent kernels in the set
            UI8 cnt = 0;

            if (AQt)        if (git->second.count(R+FR1) == 1) ++cnt;
            if (ABt && AQt) if (git->second.count(R+FR2) == 1) ++cnt;
            if (ABt)        if (git->second.count(R+FR3) == 1) ++cnt;

            //iff the kernel can be computed from the ones present in the group, register the contracted exponential and
            if (cnt==totcnt) {
                if (AQt)        imap[R].VD1 = R+FR1;
                if (ABt && AQt) imap[R].VD2 = R+FR2;
                if (ABt)        imap[R].VD3 = R+FR3;

                ikernel E = R;
                E.FE = 1;
                imap[R].VDE = E;
                imap[R].ETYPE = CDR3;
                imap[R].aux   = R.m;

                setE3.insert(E);
            }
            //otherwise register the integral as necessary
            else {
                imap[R].VD1 = R;
                setF4.insert(R);
            }
        }

        //append the recurrence relations
        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }
}
// given a set of K3 contracted exponentials, finds the recurrence relations to compute those which are linear dependant and generates a set with the irreducible ones
void EvaluationScheme::UseAERR3(const std::set<ikernel> & setE3, std::set<ikernel> & setE4) {

    const ikernel ER1   (0,1, 0,0, 1);
    const ikernel ER2   (0,0, 0,0, 1);

    const ikernel emask (1,0, 1,1, 0);

    //EARR
    //====

    set<ikernel>::const_iterator it;

    map<ikernel, set<ikernel> > egroupset;

    //classify kernels in "orthogonal" groupsets
    for (it=setE3.begin(); it!=setE3.end(); ++it) {
        ikernel index = *it | emask;
        egroupset[index].insert(*it);
    }


    map<ikernel, set<ikernel> >::iterator git;

    //loop for each group
    for (git=egroupset.begin(); git!=egroupset.end(); ++git) {
        set<ikernel>::const_iterator it;

        map<ikernel,VarDepK4> imap;

        for (it=git->second.begin(); it!=git->second.end(); ++it) {
            ikernel R = *it;

            UI8 cnt = 0;

            if (git->second.count(R+ER1)==1) ++cnt;
            if (git->second.count(R+ER2)==1) ++cnt;

            //can be computed from other e's
            if (cnt==2) {
                imap[R].VD1 = R+ER1;
                imap[R].VD2 = R+ER2;

                imap[R].ETYPE = AERR3;
            }
            //must be accumulated
            else {
                imap[R].VD1 = R;
                setE4.insert(R);
            }
        }

        map<ikernel,VarDepK4>::const_iterator itm;

        for (itm=imap.begin(); itm!=imap.end();++itm) {
            K4Var KV;

            KV.K = itm->first;
            KV.V = itm->second;

            AppendK4(KV);
        }
    }

}


void EvaluationScheme::UseCDR4 (const std::set<ikernel> & setF4, std::set<ikernel> & setE4, int totL) {

    const ikernel FR1   (0,0, 0,0, 1);

    // find lowest gamma and save RRs to obtain them all
    set<ikernel>::const_iterator it=setF4.begin();

    int lowL = it->m;
    for (it=setF4.begin(); it!=setF4.end(); ++it) lowL = min(lowL, it->m);


    map<ikernel,VarDepK4> imap;

    for (int m=lowL; m<totL; ++m) {
        ikernel R(0,0,0,0, m); R.K4L=4;

        ikernel E = R; E.FE = 1;
        imap[R].VD1   = R+FR1;
        imap[R].VDE   = E;
        imap[R].ETYPE = CDR4;
        imap[R].aux   = m;

        setE4.insert(E); // add the auxiliary exponential to the list
    }
    ikernel R(0,0,0,0, totL); R.K4L=4;
    imap[R].ETYPE = CDR4;


    //append the recurrence relations
    map<ikernel,VarDepK4>::const_iterator itm;

    for (itm=imap.begin(); itm!=imap.end();++itm) {
        K4Var KV;

        KV.K = itm->first;
        KV.V = itm->second;

        AppendK4(KV);
    }
}

void EvaluationScheme::UseAERR4(const std::set<ikernel> & setE4, int totL) {

    const ikernel ER1   (0,0, 0,0, 1);

    // find the lowest m exponential and save RRs to obtain them
    set<ikernel>::const_iterator it=setE4.begin();

    int lowL = it->m;
    for (it=setE4.begin(); it!=setE4.end(); ++it) lowL = min(lowL, it->m);


    map<ikernel,VarDepK4> imap;

    for (int m=lowL; m<totL; ++m) {
        ikernel E(0,0,0,0, m); E.FE=1; E.K4L=4;

        imap[E].VD1 = E + ER1;
        imap[E].ETYPE = AERR4;
    }

    ikernel E(0,0,0,0, totL); E.FE=1; E.K4L=4;
    imap[E].ETYPE = AERR4;


    map<ikernel,VarDepK4>::const_iterator itm;

    for (itm=imap.begin(); itm!=imap.end();++itm) {
        K4Var KV;

        KV.K = itm->first;
        KV.V = itm->second;

        AppendK4(KV);
    }

}






//link K4 semicontracted variables to the ones in the immediate previous loop

void EvaluationScheme::addS(const set<ikernel> & set0, set<ikernel> & setf) {

    set<ikernel>::const_reverse_iterator rit;

    //copy elements to evaluation list in REVERSE ORDER
    for (rit=set0.rbegin(); rit!=set0.rend(); ++rit) {
        ikernel K1 = (*rit);

        K1.t += K1.s;
        K1.s   = 0;
        K1.K4L = 1;
        setf.insert(K1);

        //registers the current set of variables with a link to the required kernel
        K4Var KV;
        KV.K = *rit;
        if (K1.FE==0) KV.V.ETYPE = K1F;
        if (K1.FE==1) KV.V.ETYPE = K1E;
        if (K1.FE==2) KV.V.ETYPE = K1D;
        KV.V.VD1   = K1;
        KV.V.aux   = rit->s;

        AppendK4(KV);
    }
}

void EvaluationScheme::addT(const set<ikernel> & set0, set<ikernel> & setf) {

    set<ikernel>::const_reverse_iterator rit;

    //copy elements to evaluation list in REVERSE ORDER
    for (rit=set0.rbegin(); rit!=set0.rend(); ++rit) {
        ikernel K2 = (*rit);
        K2.t   = 0;
        K2.K4L = 2;
        setf.insert(K2);


        K4Var KV;
        KV.K = *rit;
        if (K2.FE==0) KV.V.ETYPE = K2F;
        if (K2.FE==1) KV.V.ETYPE = K2E;
        if (K2.FE==2) KV.V.ETYPE = K2D;
        KV.V.VD1   = K2;
        KV.V.aux   = rit->t;

        AppendK4(KV);
    }
}

void EvaluationScheme::addU(const set<ikernel> & set0, set<ikernel> & setf) {

    set<ikernel>::const_reverse_iterator rit;

    //copy elements to evaluation list in REVERSE ORDER
    for (rit=set0.rbegin(); rit!=set0.rend(); ++rit) {
        ikernel K3 = (*rit);
        K3.v += K3.u;
        K3.u   =  0;
        K3.K4L =  3;
        setf.insert(K3);


        K4Var KV;
        KV.K = *rit;
        if (K3.FE==0) KV.V.ETYPE = K3F;
        if (K3.FE==1) KV.V.ETYPE = K3E;
        if (K3.FE==2) KV.V.ETYPE = K3D;
        KV.V.VD1   = K3;
        KV.V.aux   = rit->u;

        AppendK4(KV);
    }
}

void EvaluationScheme::addV(const set<ikernel> & set0, set<ikernel> & setf) {

    set<ikernel>::const_reverse_iterator rit;

    //copy elements to evaluation list in REVERSE ORDER
    for (rit=set0.rbegin(); rit!=set0.rend(); ++rit) {
        ikernel K4 = (*rit);
        K4.v   =  0;
        K4.K4L =  4;
        setf.insert(K4);

        K4Var KV;
        KV.K = *rit;
        if (K4.FE==0) KV.V.ETYPE = K4F;
        if (K4.FE==1) KV.V.ETYPE = K4E;
        if (K4.FE==2) KV.V.ETYPE = K4D;
        KV.V.VD1   = K4;
        KV.V.aux   = rit->v;

        AppendK4(KV);
    }
}


void EvaluationScheme::addBoys(int totL) {

    {
        ikernel F0m(0,0,0,0,totL); F0m.K4L = 4;

        K4Var KV;
        KV.K = F0m;
        KV.V.ETYPE = BOYS;

        AppendK4(KV);
    }

    {
        ikernel F0m(0,0,0,0,totL); F0m.K4L = 4; F0m.FE  = 1;

        K4Var KV;
        KV.K = F0m;
        KV.V.ETYPE = BOYS;

        AppendK4(KV);
    }
}


void EvaluationScheme::add0(int maxm) {
    //registers ALL variables; otherwise F0m generators' list is inconsistent
    for (int m=maxm; m>=0; --m) {
        ikernel F0m(0,0,0,0,m);
        F0m.K4L = 4;

        K4Var KV;
        KV.K = F0m;
        KV.V.ETYPE = BOYS;

        AppendK4(KV);
    }
}

void EvaluationScheme::add00(int maxm) {

    //registers ALL variables; otherwise F0m generators' list is inconsistent
    for (int m=maxm; m>=0; --m) {
        ikernel F0m(0,0,0,0,m);
        F0m.K4L = 4;

        K4Var KV;
        KV.K = F0m;
        KV.V.ETYPE = BOYS;

        AppendK4(KV);
    }

    for (int m=maxm-1; m>=0; --m) {
        ikernel F0m(0,0,0,0,m);
        F0m.K4L = 4;

        F0m.FE  = 1;

        K4Var KV;
        KV.K = F0m;
        KV.V.ETYPE = BOYS;

        AppendK4(KV);
    }
}

void EvaluationScheme::add000(int maxm) {

    //registers ALL variables; otherwise F0m generators' list is inconsistent
    for (int m=maxm; m>=0; --m) {
        ikernel F0m(0,0,0,0,m);
        F0m.K4L = 4;

        K4Var KV;
        KV.K = F0m;
        KV.V.ETYPE = BOYS;

        AppendK4(KV);
    }

    for (int m=maxm-1; m>=0; --m) {
        ikernel F0m(0,0,0,0,m);
        F0m.K4L = 4;

        F0m.FE  = 1;

        K4Var KV;
        KV.K = F0m;
        KV.V.ETYPE = BOYS;

        AppendK4(KV);
    }

    for (int m=maxm-1; m>=0; --m) {
        ikernel F0m(0,0,0,0,m);
        F0m.K4L = 4;

        F0m.FE  = 2;

        K4Var KV;
        KV.K = F0m;
        KV.V.ETYPE = BOYS;

        AppendK4(KV);
    }
}



void EvaluationScheme::MakeKernels(const set<integralCart> & set0, set<ikernel> & kset) {
    set<integralCart>::const_iterator it;

    for (it=set0.begin(); it!=set0.end(); ++it) {
        int u = it->b;
        int v = it->p;
        int s = it->d;
        int t = it->q;
        int m = it->m;

        ikernel k; //(u, v, s, t, m);
        k = *it;

        kset.insert(k);
    }


    //appends the kernels to the MIRROR evaluation scheme in the order the K4 contraction evaluated them
    set<ikernel>::const_reverse_iterator kit;

    for (kit=kset.rbegin(); kit!=kset.rend(); ++kit) {
        IntVar IV;
        IV.I = *kit;
        IV.V.ETYPE = KERNELS;
        IV.V.VD1 = *kit;

        Append(IV);
    }

}

void EvaluationScheme::AppendSet(const set<integralCart> & set0) {
    set<integralCart>::const_iterator it;

    //copy elements to evaluation list
    for (it=set0.begin(); it!=set0.end(); ++it) {
        IntVar IV;
        IV.I = *it;
        Append(IV);
    }
}



//1CRR

integralCart Rxcc(integralCart R, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABx, bool CDx, bool ACx) {
    //                        ex ey ez fx fy fz   x  y  z  a b p  c d q  m
    const integralCart MMDX1( 0, 0, 0, 0, 0, 0, -1, 0, 0, 0,1,0, 0,0,0, 1);
    const integralCart MMDX2( 0, 0, 0, 0, 0, 0, -1, 0, 0, 0,0,0, 0,1,0, 1);
    const integralCart MMDX3( 0, 0, 0, 0, 0, 0, -1, 0, 0, 0,0,0, 0,0,0, 1);
    const integralCart MMDX4( 0, 0, 0, 0, 0, 0, -2, 0, 0, 0,0,0, 0,0,0, 1);

    int rx = R.rx;

    if (iset[rx].count(R) == 1) return R;

    //x>0
    if (rx>0) {
        if (ABx)
        iset[rx][R].VD1 = Rxcc(R+MMDX1, iset, fset, ABx, CDx, ACx);
        if (CDx)
        iset[rx][R].VD2 = Rxcc(R+MMDX2, iset, fset, ABx, CDx, ACx);
        if (ACx)
        iset[rx][R].VD3 = Rxcc(R+MMDX3, iset, fset, ABx, CDx, ACx);
        if (rx>1)
        iset[rx][R].VD4 = Rxcc(R+MMDX4, iset, fset, ABx, CDx, ACx);

        iset[rx][R].ETYPE = MMDX;
        iset[rx][R].aux   = rx-1;
    }
    else {
        fset.insert(R);
        iset[0][R].aux   = 0;
    }

    return R;
}

integralCart Rycc(integralCart R, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABy, bool CDy, bool ACy) {
    //                        ex ey ez fx fy fz   x  y  z  a b p  c d q  m
    const integralCart MMDY1( 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,1,0, 0,0,0, 1);
    const integralCart MMDY2( 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,0,0, 0,1,0, 1);
    const integralCart MMDY3( 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,0,0, 0,0,0, 1);
    const integralCart MMDY4( 0, 0, 0, 0, 0, 0, 0, -2, 0, 0,0,0, 0,0,0, 1);

    int ry = R.ry;
    if (iset[ry].count(R) == 1) return R;

    if (ry>0) {
        if (ABy)
        iset[ry][R].VD1 = Rycc(R+MMDY1, iset, fset, ABy, CDy, ACy);
        if (CDy)
        iset[ry][R].VD2 = Rycc(R+MMDY2, iset, fset, ABy, CDy, ACy);
        if (ACy)
        iset[ry][R].VD3 = Rycc(R+MMDY3, iset, fset, ABy, CDy, ACy);
        if (ry>1)
        iset[ry][R].VD4 = Rycc(R+MMDY4, iset, fset, ABy, CDy, ACy);

        iset[ry][R].ETYPE = MMDY;
        iset[ry][R].aux   = ry-1;
    }
    else {
        fset.insert(R);
        iset[0][R].aux   = 0;
    }

    return R;
}

integralCart Rzcc(integralCart R, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABz, bool CDz, bool ACz) {
    //                        ex ey ez fx fy fz   x  y  z  a b p  c d q  m
    const integralCart MMDZ1( 0, 0, 0, 0, 0, 0,  0, 0,-1, 0,1,0, 0,0,0, 1);
    const integralCart MMDZ2( 0, 0, 0, 0, 0, 0,  0, 0,-1, 0,0,0, 0,1,0, 1);
    const integralCart MMDZ3( 0, 0, 0, 0, 0, 0,  0, 0,-1, 0,0,0, 0,0,0, 1);
    const integralCart MMDZ4( 0, 0, 0, 0, 0, 0,  0, 0,-2, 0,0,0, 0,0,0, 1);

    int rz = R.rz;
    if (iset[rz].count(R) == 1) return R;

    if (rz>0) {
        if (ABz)
        iset[rz][R].VD1 = Rzcc(R+MMDZ1, iset, fset, ABz, CDz, ACz);
        if (CDz)
        iset[rz][R].VD2 = Rzcc(R+MMDZ2, iset, fset, ABz, CDz, ACz);
        if (ACz)
        iset[rz][R].VD3 = Rzcc(R+MMDZ3, iset, fset, ABz, CDz, ACz);
        if (rz>1)
        iset[rz][R].VD4 = Rzcc(R+MMDZ4, iset, fset, ABz, CDz, ACz);

        iset[rz][R].ETYPE = MMDZ;
        iset[rz][R].aux   = rz-1;
    }
    else {
        fset.insert(R);
        iset[0][R].aux   = 0;
    }

    return R;
}


void EvaluationScheme::addRxcc(const set<integralCart> & set0, set<integralCart> & setf, bool ABx, bool CDx, bool ACx) {
    integralCart mask( 1,1,1, 1,1,1, 0,1,1, 1,1,1, 1,1,1, 0);
    if (ABx) mask.b = 0;
    if (CDx) mask.d = 0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            Rxcc(*it3, imap, setf, ABx, CDx, ACx);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addRycc(const set<integralCart> & set0, set<integralCart> & setf, bool ABy, bool CDy, bool ACy) {
    integralCart mask( 1,1,1, 1,1,1, 1,0,1, 1,1,1, 1,1,1, 0);
    if (ABy) mask.b = 0;
    if (CDy) mask.d = 0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            Rycc(*it3, imap, setf, ABy, CDy, ACy);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addRzcc(const set<integralCart> & set0, set<integralCart> & setf, bool ABz, bool CDz, bool ACz) {
    integralCart mask( 1,1,1, 1,1,1, 1,1,0, 1,1,1, 1,1,1, 0);
    if (ABz) mask.b = 0;
    if (CDz) mask.d = 0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            Rzcc(*it3, imap, setf, ABz, CDz, ACz);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}


//CTEs

integralCart CTEbx(integralCart CTE,  map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABx) {
    //                          ex ey ez fx fy fz                  rx ry rz  a b p  c d q  m
    const integralCart CTERR1( -1, 0, 0, 0, 0, 0, -1,0,0, 0,0,0,  -1, 0, 0, 0,0,0, 0,0,0, 0);
    const integralCart CTERR2( -1, 0, 0, 0, 0, 0,  0,0,0, 0,0,0,   0, 0, 0, 0,1,0, 0,0,0, 0);
    const integralCart CTERR3( -1, 0, 0, 0, 0, 0,  1,0,0, 0,0,0,   1, 0, 0, 0,0,1, 0,0,0, 0);

    int ex = CTE.ex;

    if (ex==0) CTE.px = 0; //the extra var doesn't matter

    if (iset[ex].count(CTE) == 1) return CTE;

    //ex>0
    if (ex>0) {
        int px = CTE.px;
        if (px>0)
        iset[ex][CTE].VD1 = CTEbx(CTE+CTERR1, iset, fset, ABx);
        if (ABx)
        iset[ex][CTE].VD2 = CTEbx(CTE+CTERR2, iset, fset, ABx);
        iset[ex][CTE].VD3 = CTEbx(CTE+CTERR3, iset, fset, ABx);

        iset[ex][CTE].ETYPE = CTEBX;
        iset[ex][CTE].aux = px;
    }
    else {
        fset.insert(CTE);
        iset[0][CTE].aux = 0;
    }

    return CTE;
}

integralCart CTEby(integralCart CTE,  map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABy) {
    //                          ex ey ez fx fy fz                  rx ry rz  a b p  c d q  m
    const integralCart CTERR1(  0,-1, 0, 0, 0, 0,  0,-1,0, 0,0,0,   0,-1, 0, 0,0,0, 0,0,0, 0);
    const integralCart CTERR2(  0,-1, 0, 0, 0, 0,  0, 0,0, 0,0,0,   0, 0, 0, 0,1,0, 0,0,0, 0);
    const integralCart CTERR3(  0,-1, 0, 0, 0, 0,  0, 1,0, 0,0,0,   0, 1, 0, 0,0,1, 0,0,0, 0);

    int ey = CTE.ey;

    if (ey==0) CTE.py = 0; //the extra var doesn't matter

    if (iset[ey].count(CTE) == 1) return CTE;

    //ey>0
    if (ey>0) {
        int py = CTE.py;
        if (py>0)
        iset[ey][CTE].VD1 = CTEby(CTE+CTERR1, iset, fset, ABy);
        if (ABy)
        iset[ey][CTE].VD2 = CTEby(CTE+CTERR2, iset, fset, ABy);
        iset[ey][CTE].VD3 = CTEby(CTE+CTERR3, iset, fset, ABy);

        iset[ey][CTE].ETYPE = CTEBY;
        iset[ey][CTE].aux = py;
    }
    else {
        fset.insert(CTE);
        iset[0][CTE].aux = 0;
    }

    return CTE;
}

integralCart CTEbz(integralCart CTE,  map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABz) {
    //                          ex ey ez fx fy fz                  rx ry rz  a b p  c d q  m
    const integralCart CTERR1( 0, 0,-1, 0, 0, 0,  0,0,-1, 0,0,0,   0, 0,-1, 0,0,0, 0,0,0, 0);
    const integralCart CTERR2( 0, 0,-1, 0, 0, 0,  0,0, 0, 0,0,0,   0, 0, 0, 0,1,0, 0,0,0, 0);
    const integralCart CTERR3( 0, 0,-1, 0, 0, 0,  0,0, 1, 0,0,0,   0, 0, 1, 0,0,1, 0,0,0, 0);

    int ez = CTE.ez;

    if (ez==0) CTE.pz = 0; //the extra var doesn't matter

    if (iset[ez].count(CTE) == 1) return CTE;

    //ez>0
    if (ez>0) {
        int pz = CTE.pz;
        if (pz>0)
        iset[ez][CTE].VD1 = CTEbz(CTE+CTERR1, iset, fset, ABz);
        if (ABz)
        iset[ez][CTE].VD2 = CTEbz(CTE+CTERR2, iset, fset, ABz);
        iset[ez][CTE].VD3 = CTEbz(CTE+CTERR3, iset, fset, ABz);

        iset[ez][CTE].ETYPE = CTEBZ;
        iset[ez][CTE].aux = pz;
    }
    else {
        fset.insert(CTE);
        iset[0][CTE].aux = 0;
    }

    return CTE;
}


integralCart CTEkx(integralCart CTE,  map<integralCart,VarDep> * iset, set<integralCart> & fset, bool CDx) {
    //                          ex ey ez fx fy fz                  rx ry rz  a b p  c d q  m
    const integralCart CTERR1(  0, 0, 0, -1, 0, 0, 0,0,0, -1,0,0, -1, 0, 0, 0,0,0, 0,0,0, 0);
    const integralCart CTERR2(  0, 0, 0, -1, 0, 0, 0,0,0,  0,0,0,  0, 0, 0, 0,0,0, 0,1,0, 0);
    const integralCart CTERR3(  0, 0, 0, -1, 0, 0, 0,0,0,  1,0,0,  1, 0, 0, 0,0,0, 0,0,1, 0);

    int fx = CTE.fx;

    if (fx==0) CTE.qx = 0; //the extra var doesn't matter

    if (iset[fx].count(CTE) == 1) return CTE;

    //fx>0
    if (fx>0) {
        int qx = CTE.qx;
        if (qx>0)
        iset[fx][CTE].VD1 = CTEkx(CTE+CTERR1, iset, fset, CDx);
        if (CDx)
        iset[fx][CTE].VD2 = CTEkx(CTE+CTERR2, iset, fset, CDx);
        iset[fx][CTE].VD3 = CTEkx(CTE+CTERR3, iset, fset, CDx);

        iset[fx][CTE].ETYPE = CTEKX;
        iset[fx][CTE].aux = qx;
    }
    else {
        fset.insert(CTE);
        iset[0][CTE].aux = 0;
    }

    return CTE;
}

integralCart CTEky(integralCart CTE,  map<integralCart,VarDep> * iset, set<integralCart> & fset, bool CDy) {
    //                          ex ey ez fx fy fz                  rx ry rz  a b p  c d q  m
    const integralCart CTERR1( 0, 0, 0, 0,-1, 0, 0,0,0,  0,-1,0,   0,-1, 0, 0,0,0, 0,0,0, 0);
    const integralCart CTERR2( 0, 0, 0, 0,-1, 0, 0,0,0,  0, 0,0,   0, 0, 0, 0,0,0, 0,1,0, 0);
    const integralCart CTERR3( 0, 0, 0, 0,-1, 0, 0,0,0,  0, 1,0,   0, 1, 0, 0,0,0, 0,0,1, 0);

    int fy = CTE.fy;

    if (fy==0) CTE.qy = 0; //the extra var doesn't matter

    if (iset[fy].count(CTE) == 1) return CTE;


    //fy>0
    if (fy>0) {
        int qy = CTE.qy;
        if (qy>0)
        iset[fy][CTE].VD1 = CTEky(CTE+CTERR1, iset, fset, CDy);
        if (CDy)
        iset[fy][CTE].VD2 = CTEky(CTE+CTERR2, iset, fset, CDy);
        iset[fy][CTE].VD3 = CTEky(CTE+CTERR3, iset, fset, CDy);

        iset[fy][CTE].ETYPE = CTEKY;
        iset[fy][CTE].aux = qy;
    }
    else {
        fset.insert(CTE);
        iset[0][CTE].aux = 0;
    }

    return CTE;
}

integralCart CTEkz(integralCart CTE,  map<integralCart,VarDep> * iset, set<integralCart> & fset, bool CDz) {
    //                         ex ey ez fx fy fz                 rx ry rz  a b p  c d q  m
    const integralCart CTERR1( 0, 0, 0, 0, 0,-1, 0,0,0, 0,0,-1,  0, 0,-1, 0,0,0, 0,0,0, 0);
    const integralCart CTERR2( 0, 0, 0, 0, 0,-1, 0,0,0, 0,0, 0,  0, 0, 0, 0,0,0, 0,1,0, 0);
    const integralCart CTERR3( 0, 0, 0, 0, 0,-1, 0,0,0, 0,0, 1,  0, 0, 1, 0,0,0, 0,0,1, 0);

    int fz = CTE.fz;

    if (fz==0) CTE.qz = 0; //the extra var doesn't matter

    if (iset[fz].count(CTE) == 1) return CTE;

    //fz>0
    if (fz>0) {
        int qz = CTE.qz;
        if (qz>0)
        iset[fz][CTE].VD1 = CTEkz(CTE+CTERR1, iset, fset, CDz);
        if (CDz)
        iset[fz][CTE].VD2 = CTEkz(CTE+CTERR2, iset, fset, CDz);
        iset[fz][CTE].VD3 = CTEkz(CTE+CTERR3, iset, fset, CDz);

        iset[fz][CTE].ETYPE = CTEKZ;
        iset[fz][CTE].aux = qz;
    }
    else {
        fset.insert(CTE);
        iset[0][CTE].aux = 0;
    }

    return CTE;
}


void EvaluationScheme::addCTEbx(const set<integralCart> & set0, set<integralCart> & setf, bool ABx) {
    //                       exeyez fxfyfz                rxryrz  a b p  c d q  m
    integralCart mask( 0,1,1, 1,1,1, 0,1,1, 1,1,1,  0,1,1, 1,1,0, 1,1,1, 1);
    if (ABx) mask.b=0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            CTEbx(*it3, imap, setf, ABx);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addCTEby(const set<integralCart> & set0, set<integralCart> & setf, bool ABy) {
    //                       exeyez fxfyfz                rxryrz  a b p  c d q  m
    integralCart mask( 1,0,1, 1,1,1, 1,0,1, 1,1,1,  1,0,1, 1,1,0, 1,1,1, 1);
    if (ABy) mask.b=0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            CTEby(*it3, imap, setf, ABy);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addCTEbz(const set<integralCart> & set0, set<integralCart> & setf, bool ABz) {
    //                       exeyez fxfyfz                rxryrz  a b p  c d q  m
    integralCart mask( 1,1,0, 1,1,1, 1,1,0, 1,1,1,  1,1,0, 1,1,0, 1,1,1, 1);
    if (ABz) mask.b=0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            CTEbz(*it3, imap, setf, ABz);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}


void EvaluationScheme::addCTEkx(const set<integralCart> & set0, set<integralCart> & setf, bool CDx) {
    //                       exeyez fxfyfz                rxryrz  a b p  c d q  m
    integralCart mask( 1,1,1, 0,1,1, 1,1,1, 0,1,1,  0,1,1, 1,1,1, 1,1,0, 1);
    if (CDx) mask.d=0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            CTEkx(*it3, imap, setf, CDx);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addCTEky(const set<integralCart> & set0, set<integralCart> & setf, bool CDy) {
    //                       exeyez fxfyfz                rxryrz  a b p  c d q  m
    integralCart mask( 1,1,1, 1,0,1, 1,1,1, 1,0,1,  1,0,1, 1,1,1, 1,1,0, 1);
    if (CDy) mask.d=0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            CTEky(*it3, imap, setf, CDy);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addCTEkz(const set<integralCart> & set0, set<integralCart> & setf, bool CDz) {
    //                exeyez fxfyfz                rxryrz  a b p  c d q  m
    integralCart mask( 1,1,1, 1,1,0, 1,1,1, 1,1,0,  1,1,0, 1,1,1, 1,1,0, 1);
    if (CDz) mask.d=0;

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            CTEkz(*it3, imap, setf, CDz);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}


//HRRs

integralCart HRRbx(integralCart HRR, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABx) {
    const integralCart HRR1(1,0,0,-1,0,0,0,0,0,0,0,0);
    const integralCart HRR2(0,0,0,-1,0,0,0,0,0,0,0,0);

    int bx = HRR.bx;

    //already initialized, return id of element
    if (iset[bx].count(HRR) == 1) return HRR;


    //ax - bx
    if (bx>0) {
        //if (!ABx) return HRRbx(HRR + HRR1,  iset, fset, ABx);

        iset[bx][HRR].VD1 = HRRbx(HRR+HRR1,  iset, fset, ABx);
        if (ABx)
        iset[bx][HRR].VD2 = HRRbx(HRR+HRR2,  iset, fset, ABx);

        iset[bx][HRR].ETYPE = HRRBX;
    }
    else {
        fset.insert(HRR);
        iset[0][HRR].aux = 0;
    }

    return HRR;
}

integralCart HRRby(integralCart HRR, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABy) {
    const integralCart HRR1(0,1,0,0,-1,0,0,0,0,0,0,0);
    const integralCart HRR2(0,0,0,0,-1,0,0,0,0,0,0,0);

    int by = HRR.by;

    //already initialized, return id of element
    if (iset[by].count(HRR) == 1) return HRR;


    //cx - dx
    if (by>0) {
        //if (!ABy) return HRRby(HRR + HRR1,  iset, fset, ABy);

        iset[by][HRR].VD1 = HRRby(HRR+HRR1,  iset, fset, ABy);
        if (ABy)
        iset[by][HRR].VD2 = HRRby(HRR+HRR2,  iset, fset, ABy);

        iset[by][HRR].ETYPE = HRRBY;
    }
    else {
        fset.insert(HRR);
        iset[0][HRR].aux = 0;
    }

    return HRR;
}

integralCart HRRbz(integralCart HRR, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool ABz) {
    const integralCart HRR1(0,0,1,0,0,-1,0,0,0,0,0,0);
    const integralCart HRR2(0,0,0,0,0,-1,0,0,0,0,0,0);

    int bz = HRR.bz;

    //already initialized, return id of element
    if (iset[bz].count(HRR) == 1) return HRR;

    //az - bz
    if (bz>0) {
        //if (!ABz) return HRRbz(HRR + HRR1,  iset, fset, ABz);

        iset[bz][HRR].VD1 = HRRbz(HRR+HRR1,  iset, fset, ABz);
        if (ABz)
        iset[bz][HRR].VD2 = HRRbz(HRR+HRR2,  iset, fset, ABz);

        iset[bz][HRR].ETYPE = HRRBZ;
    }
    else {
        fset.insert(HRR);
        iset[0][HRR].aux = 0;
    }

    return HRR;
}



integralCart HRRkx(integralCart HRR, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool CDx) {
    const integralCart HRR1(0,0,0,0,0,0,1,0,0,-1,0,0);
    const integralCart HRR2(0,0,0,0,0,0,0,0,0,-1,0,0);

    int dx = HRR.dx;

    //already initialized, return id of element
    if (iset[dx].count(HRR) == 1) return HRR;

    //cx - dx
    if (dx>0) {
        //if (!CDx) return HRRkx(HRR + HRR1,  iset, fset, CDx);

        iset[dx][HRR].VD1 = HRRkx(HRR + HRR1,  iset, fset, CDx);
        if (CDx)
        iset[dx][HRR].VD2 = HRRkx(HRR + HRR2,  iset, fset, CDx);

        iset[dx][HRR].ETYPE = HRRKX;
    }
    else {
        fset.insert(HRR);
        iset[0][HRR].aux = 0;
    }

    return HRR;
}

integralCart HRRky(integralCart HRR, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool CDy) {
    const integralCart HRR1(0,0,0,0,0,0,0,1,0,0,-1,0);
    const integralCart HRR2(0,0,0,0,0,0,0,0,0,0,-1,0);

    int dy = HRR.dy;

    //already initialized, return id of element
    if (iset[dy].count(HRR) == 1) return HRR;

    if (dy>0) {
        //if (!CDy) return HRRky(HRR + HRR1,  iset, fset, CDy);

        iset[dy][HRR].VD1 = HRRky(HRR + HRR1, iset, fset, CDy);
        if (CDy)
        iset[dy][HRR].VD2 = HRRky(HRR + HRR2, iset, fset, CDy);

        iset[dy][HRR].ETYPE = HRRKY;
    }
    else {
        fset.insert(HRR);
        iset[0][HRR].aux = 0;
    }

    return HRR;
}

integralCart HRRkz(integralCart HRR, map<integralCart,VarDep> * iset, set<integralCart> & fset, bool CDz) {
    const integralCart HRR1(0,0,0,0,0,0,0,0,1,0,0,-1);
    const integralCart HRR2(0,0,0,0,0,0,0,0,0,0,0,-1);

    int dz = HRR.dz;

    //already initialized, return id of element
    if (iset[dz].count(HRR) == 1) return HRR;


    //cz - dz
    if (dz>0) {
        //if (!CDz) return HRRkz(HRR + HRR1,  iset, fset, CDz);

        iset[dz][HRR].VD1 = HRRkz(HRR + HRR1, iset, fset, CDz);
        if (CDz)
        iset[dz][HRR].VD2 = HRRkz(HRR + HRR2, iset, fset, CDz);

        iset[dz][HRR].ETYPE = HRRKZ;
    }
    else {
        fset.insert(HRR);
        iset[0][HRR].aux = 0;
    }

    return HRR;
}



void EvaluationScheme::addHRRbx(const set<integralCart> & set0, set<integralCart> & setf, bool ABx) {
    const integralCart mask(0,1,1, 0,1,1, 1,1,1, 1,1,1);

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            HRRbx(*it3, imap, setf, ABx);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addHRRby(const set<integralCart> & set0, set<integralCart> & setf, bool ABy) {
    const integralCart mask(1,0,1, 1,0,1, 1,1,1, 1,1,1);

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            HRRby(*it3, imap, setf, ABy);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addHRRbz(const set<integralCart> & set0, set<integralCart> & setf, bool ABz) {
    const integralCart mask(1,1,0, 1,1,0, 1,1,1, 1,1,1);

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            HRRbz(*it3, imap, setf, ABz);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}


void EvaluationScheme::addHRRkx(const set<integralCart> & set0, set<integralCart> & setf, bool CDx) {
    const integralCart mask(1,1,1, 1,1,1, 0,1,1, 0,1,1);

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            HRRkx(*it3, imap, setf, CDx);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addHRRky(const set<integralCart> & set0, set<integralCart> & setf, bool CDy) {
    const integralCart mask(1,1,1, 1,1,1, 1,0,1, 1,0,1);

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            HRRky(*it3, imap, setf, CDy);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}

void EvaluationScheme::addHRRkz(const set<integralCart> & set0, set<integralCart> & setf, bool CDz) {
    const integralCart mask(1,1,1, 1,1,1, 1,1,0, 1,1,0);

    set<integralCart>::const_iterator it;

    map<integralCart, set<integralCart> > groupset;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        integralCart index = *it | mask;
        groupset[index].insert(*it);
    }

    map<integralCart, set<integralCart> >::iterator it2;

    //loop for each group
    for (it2=groupset.begin(); it2!=groupset.end(); ++it2) {
        set<integralCart>::iterator it3;

        map<integralCart,VarDep> imap[1+4*LMAX+4];

        //loop for all elements in the group
        for (it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
            HRRkz(*it3, imap, setf, CDz);
        }

        AddGroup(imap, it2->second);
    }
    AddGroups();
}




void intSphLink(integralCart Sph, VarDep & VD, set<integralCart> & fset, CENTER center, int l)  {

    int m;
    if(center==AAA) m = Sph.la;
    if(center==BBB) m = Sph.lb;
    if(center==CCC) m = Sph.lc;
    if(center==DDD) m = Sph.ld;

    VD.aux = m;

    for (int p=0; p<SHList[l][m].nps; ++p) {
        int x = SHList[l][m].T[p].nx;
        int y = SHList[l][m].T[p].ny;
        int z = SHList[l][m].T[p].nz;

        integralCart Sph2 = Sph;
        if(center==AAA) {
            Sph2.la = 0; //-1;
            Sph2.ex += x;
            Sph2.ey += y;
            Sph2.ez += z;
        }
        if(center==BBB) {
            Sph2.lb = 0; //-1;
            Sph2.bx += x;
            Sph2.by += y;
            Sph2.bz += z;
        }
        if(center==CCC) {
            Sph2.lc = 0; //-1;
            Sph2.fx += x;
            Sph2.fy += y;
            Sph2.fz += z;
        }
        if(center==DDD) {
            Sph2.ld = 0; //-1;
            Sph2.dx += x;
            Sph2.dy += y;
            Sph2.dz += z;
        }

        if      (p==0) VD.VD1 = Sph2;
        else if (p==1) VD.VD2 = Sph2;
        else if (p==2) VD.VD3 = Sph2;
        else if (p==3) VD.VD4 = Sph2;
        else if (p==4) VD.VD5 = Sph2;
        else if (p==5) VD.VD6 = Sph2;

        fset.insert(Sph2);
    }

    if(center==AAA) VD.ETYPE = SPHA;
    if(center==BBB) VD.ETYPE = SPHB;
    if(center==CCC) VD.ETYPE = SPHC;
    if(center==DDD) VD.ETYPE = SPHD;
}


void EvaluationScheme::addSphA (const set<integralCart> & set0, set<integralCart> & fset) {
    set<integralCart>::const_iterator it;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        VarDep VD;
        intSphLink(*it, VD ,fset, AAA, La);
        IntVar IV;
        IV.I = *it;
        IV.V = VD;
        Append(IV);
    }
}

void EvaluationScheme::addSphB (const set<integralCart> & set0, set<integralCart> & fset) {
    set<integralCart>::const_iterator it;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        VarDep VD;
        intSphLink(*it, VD ,fset, BBB, Lb);
        IntVar IV;
        IV.I = *it;
        IV.V = VD;
        Append(IV);
    }
}

void EvaluationScheme::addSphC (const set<integralCart> & set0, set<integralCart> & fset) {
    set<integralCart>::const_iterator it;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        VarDep VD;
        intSphLink(*it, VD ,fset, CCC, Lc);
        IntVar IV;
        IV.I = *it;
        IV.V = VD;
        Append(IV);
    }
}

void EvaluationScheme::addSphD (const set<integralCart> & set0, set<integralCart> & fset) {
    set<integralCart>::const_iterator it;

    //classify integrals in "orthogonal" groupset
    for (it=set0.begin(); it!=set0.end(); ++it) {
        VarDep VD;
        intSphLink(*it, VD ,fset, DDD, Ld);
        IntVar IV;
        IV.I = *it;
        IV.V = VD;
        Append(IV);
    }
}

static inline int Tpos (int M1, int M2, int M3, int M4, int v1, int v2, int v3, int v4) {
    return v4 * M3*M2*M1 + v3*M2*M1 + v2*M1 + v1;
}


void EvaluationScheme::addInitial (set<integralCart> & fset) {


    //REORDER
    for (int d=0; d<nmS[Ld]; ++d) {
        for (int c=0; c<nmS[Lc]; ++c) {
            for (int b=0; b<nmS[Lb]; ++b) {
                for (int a=0; a<nmS[La]; ++a) {

                    int pos = Tpos (nmS[La], nmS[Lb], nmS[Lc], nmS[Ld], a, b, c, d);

                    integralCart ERI0(a,b,c,d); ERI0.m = -1; //mark it
                    integralCart ERIf(a,b,c,d);

                    VarDep V;
                    V.VD1 = ERIf;
                    V.aux = pos;   //this will be the position in the final tensor
                    V.ETYPE = REORDER;

                    IntVar IV;

                    IV.I = ERI0;
                    IV.V = V;

                    EVList.push_back(IV);

                    fset.insert(ERIf);
                }
            }
        }
    }
}




//FOR THE SKS ALGORITHM
//*********************


/*

int SRRa(integralSph SRR, map< integralSph, RRid > * mapRR) {
    const integralSph SRR1(0,0,1, 0,0,0, -1,0,0,0,  0,0,0,0);
    const integralSph SRR2(1,0,0, 0,0,0, -1,0,0,0, -1,0,0,0);
    const integralSph SRR3(0,1,0, 0,0,0, -1,0,0,0, -1,0,0,0);

    //SRR.ma = abs(SRR.ma);

    int la = SRR.la;
    int ma = SRR.ma;

    if (la<ma) return 0;

    //already initialized, return id of element
    if (mapRR[la].count(SRR) == 1) {
        mapRR[la][SRR].ref = true;
        ++mapRR[la][SRR].count;
        return mapRR[la][SRR].id;
    }

    RRid ids;

    //cx - dx
    if (la>0) {
        ids.pid1 = SRRa(SRR + SRR1,  mapRR);
        ids.pid2 = SRRa(SRR + SRR2,  mapRR);
        ids.pid3 = SRRa(SRR + SRR3,  mapRR);
    }

    //set id of element
    ids.id = mapRR[la].size() +1;
    mapRR[la][SRR] = ids;
    mapRR[la][SRR].ref = true;
    return ids.id;
}

int SRRb(integralSph SRR, map< integralSph, RRid > * mapRR, bool ABx, bool ABy, bool ABz) {
    const integralSph SRR1(0,0,1, 0,0,0, 0,-1,0,0, 0, 0,0,0);
    const integralSph SRR2(1,0,0, 0,0,0, 0,-1,0,0, 0,-1,0,0);
    const integralSph SRR3(0,1,0, 0,0,0, 0,-1,0,0, 0,-1,0,0);

    const integralSph SRR4(0,0,0, 0,0,0, 0,-1,0,0, 0, 0,0,0);
    const integralSph SRR5(0,0,0, 0,0,0, 0,-1,0,0, 0,-1,0,0);
    const integralSph SRR6(0,0,0, 0,0,0, 0,-1,0,0, 0,-1,0,0);

    //SRR.mb = abs(SRR.mb);

    int lb = SRR.lb;
    int mb = SRR.mb;

    if (lb<mb) return 0;

    //cout << SRR << endl;

    //already initialized, return id of element
    if (mapRR[lb].count(SRR) == 1) {
        mapRR[lb][SRR].ref = true;
        ++mapRR[lb][SRR].count;
        return mapRR[lb][SRR].id;
    }

    RRid ids;

    //cx - dx
    if (lb>0) {
        ids.pid1 = SRRb(SRR + SRR1,  mapRR, ABx, ABy, ABz);
        ids.pid2 = SRRb(SRR + SRR2,  mapRR, ABx, ABy, ABz);
        ids.pid3 = SRRb(SRR + SRR3,  mapRR, ABx, ABy, ABz);
        if (ABx) ids.pid4 = SRRb(SRR + SRR4,  mapRR, ABx, ABy, ABz);
        if (ABy) ids.pid5 = SRRb(SRR + SRR5,  mapRR, ABx, ABy, ABz);
        if (ABz) ids.pid6 = SRRb(SRR + SRR6,  mapRR, ABx, ABy, ABz);
    }

    //set id of element
    ids.id = mapRR[lb].size() +1;
    mapRR[lb][SRR] = ids;
    mapRR[lb][SRR].ref = true;
    return ids.id;
}

int SRRc(integralSph SRR, map< integralSph, RRid > * mapRR) {
    const integralSph SRR1(0,0,0, 0,0,1, 0,0,-1,0, 0,0, 0,0);
    const integralSph SRR2(0,0,0, 1,0,0, 0,0,-1,0, 0,0,-1,0);
    const integralSph SRR3(0,0,0, 0,1,0, 0,0,-1,0, 0,0,-1,0);

    //SRR.mc = abs(SRR.mc);

    int lc = SRR.lc;
    int mc = SRR.mc;

    if (lc<mc) return 0;

    //cout << SRR << endl;

    //already initialized, return id of element
    if (mapRR[lc].count(SRR) == 1) {
        mapRR[lc][SRR].ref = true;
        ++mapRR[lc][SRR].count;
        return mapRR[lc][SRR].id;
    }

    RRid ids;

    //cx - dx
    if (lc>0) {
        ids.pid1 = SRRc(SRR + SRR1,  mapRR);
        ids.pid2 = SRRc(SRR + SRR2,  mapRR);
        ids.pid3 = SRRc(SRR + SRR3,  mapRR);
    }

    //set id of element
    ids.id = mapRR[lc].size() +1;
    mapRR[lc][SRR] = ids;
    mapRR[lc][SRR].ref = true;
    return ids.id;
}

int SRRd(integralSph SRR, map< integralSph, RRid > * mapRR, bool CDx, bool CDy, bool CDz) {
    const integralSph SRR1(0,0,0, 0,0,1, 0,0,0,-1, 0,0,0, 0);
    const integralSph SRR2(0,0,0, 1,0,0, 0,0,0,-1, 0,0,0,-1);
    const integralSph SRR3(0,0,0, 0,1,0, 0,0,0,-1, 0,0,0,-1);

    const integralSph SRR4(0,0,0, 0,0,0, 0,0,0,-1, 0,0,0, 0);
    const integralSph SRR5(0,0,0, 0,0,0, 0,0,0,-1, 0,0,0,-1);
    const integralSph SRR6(0,0,0, 0,0,0, 0,0,0,-1, 0,0,0,-1);

    //SRR.md = abs(SRR.md);

    int ld = SRR.ld;
    int md = SRR.md;

    if (ld<md) return 0;

    //cout << SRR << endl;

    //already initialized, return id of element
    if (mapRR[ld].count(SRR) == 1) {
        mapRR[ld][SRR].ref = true;
        ++mapRR[ld][SRR].count;
        return mapRR[ld][SRR].id;
    }

    RRid ids;

    //cx - dx
    if (ld>0) {
        ids.pid1 = SRRd(SRR + SRR1,  mapRR, CDx, CDy, CDz);
        ids.pid2 = SRRd(SRR + SRR2,  mapRR, CDx, CDy, CDz);
        ids.pid3 = SRRd(SRR + SRR3,  mapRR, CDx, CDy, CDz);
        if (CDx) ids.pid4 = SRRd(SRR + SRR4,  mapRR, CDx, CDy, CDz);
        if (CDy) ids.pid5 = SRRd(SRR + SRR5,  mapRR, CDx, CDy, CDz);
        if (CDz) ids.pid6 = SRRd(SRR + SRR6,  mapRR, CDx, CDy, CDz);
    }

    //set id of element
    ids.id = mapRR[ld].size() +1;
    mapRR[ld][SRR] = ids;
    mapRR[ld][SRR].ref = true;
    return ids.id;
}


int CTEkx(integralSph CTE, map<integralSph, RRid> * mapRR, bool CDx) {
    //                         ex ey ez  fx fy fz                 rx ry rz  a b p  c d q  m
    const integralSph CTERR1(  0, 0, 0, -1, 0, 0, 0,0,0, -1,0,0, -1, 0, 0, 0,0,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR2(  0, 0, 0, -1, 0, 0, 0,0,0,  0,0,0,  0, 0, 0, 0,0,0, 0,1,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR3(  0, 0, 0, -1, 0, 0, 0,0,0,  1,0,0,  1, 0, 0, 0,0,0, 0,0,1, 0, 0,0,0,0, 0,0,0,0);

    int fx = CTE.fx;

    if (fx==0) CTE.qx = 0; //the extra var doesn't matter

    //already initialized, return id of element
    if (mapRR[fx].count(CTE) == 1) {
        mapRR[fx][CTE].ref = true;
        ++mapRR[fx][CTE].count;
        return mapRR[fx][CTE].id;
    }

    RRid ids;

    //apply RR further still
    if (fx>0) {
        int qx = CTE.qx;
        if (qx>0)
        ids.pid1 = CTEkx(CTE+CTERR1, mapRR, CDx);
        if (CDx) ids.pid2 = CTEkx(CTE+CTERR2, mapRR, CDx);
        ids.pid3 = CTEkx(CTE+CTERR3, mapRR, CDx);
    }

    //set id of element
    ids.id = mapRR[fx].size() +1;
    mapRR[fx][CTE] = ids;
    mapRR[fx][CTE].ref = true;
    return ids.id;
}

int CTEky(integralSph CTE, map<integralSph, RRid> * mapRR, bool CDy) {
    //                        ex ey ez fx fy fz                   rx ry rz  a b p  c d q  m
    const integralSph CTERR1( 0, 0, 0, 0,-1, 0, 0,0,0,  0,-1,0,   0,-1, 0, 0,0,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR2( 0, 0, 0, 0,-1, 0, 0,0,0,  0, 0,0,   0, 0, 0, 0,0,0, 0,1,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR3( 0, 0, 0, 0,-1, 0, 0,0,0,  0, 1,0,   0, 1, 0, 0,0,0, 0,0,1, 0, 0,0,0,0, 0,0,0,0);

    int fy = CTE.fy;

    if (fy==0) CTE.qy = 0; //the extra var doesn't matter

    if (mapRR[fy].count(CTE) == 1) {
        mapRR[fy][CTE].ref = true;
        ++mapRR[fy][CTE].count;
        return mapRR[fy][CTE].id;
    }

    RRid ids;

    //fy>0
    if (fy>0) {
        int qy = CTE.qy;
        if (qy>0)
        ids.pid1 = CTEky(CTE+CTERR1, mapRR, CDy);
        if (CDy) ids.pid2 = CTEky(CTE+CTERR2, mapRR, CDy);
        ids.pid3 = CTEky(CTE+CTERR3, mapRR, CDy);
    }

    //set id of element
    ids.id = mapRR[fy].size() +1;
    mapRR[fy][CTE] = ids;
    mapRR[fy][CTE].ref = true;
    return ids.id;
}

int CTEkz(integralSph CTE, map<integralSph, RRid> * mapRR, bool CDz) {
    //                        ex ey ez fx fy fz                 rx ry rz  a b p  c d q  m
    const integralSph CTERR1( 0, 0, 0, 0, 0,-1, 0,0,0, 0,0,-1,  0, 0,-1, 0,0,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR2( 0, 0, 0, 0, 0,-1, 0,0,0, 0,0, 0,  0, 0, 0, 0,0,0, 0,1,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR3( 0, 0, 0, 0, 0,-1, 0,0,0, 0,0, 1,  0, 0, 1, 0,0,0, 0,0,1, 0, 0,0,0,0, 0,0,0,0);

    int fz = CTE.fz;

    if (fz==0) CTE.qz = 0; //the extra var doesn't matter

    if (mapRR[fz].count(CTE) == 1) {
        mapRR[fz][CTE].ref = true;
        ++mapRR[fz][CTE].count;
        return mapRR[fz][CTE].id;
    }

    RRid ids;

    //fz>0
    if (fz>0) {
        int qz = CTE.qz;
        if (qz>0)
        ids.pid1 = CTEkz(CTE+CTERR1, mapRR, CDz);
        if (CDz) ids.pid2 = CTEkz(CTE+CTERR2, mapRR, CDz);
        ids.pid3 = CTEkz(CTE+CTERR3, mapRR, CDz);
    }

    //set id of element
    ids.id = mapRR[fz].size() +1;
    mapRR[fz][CTE] = ids;
    mapRR[fz][CTE].ref = true;
    return ids.id;
}

int CTEbx(integralSph CTE, map<integralSph, RRid> * mapRR, bool ABx) {
    //                         ex ey ez fx fy fz                  rx ry rz  a b p  c d q  m
    const integralSph CTERR1( -1, 0, 0, 0, 0, 0, -1,0,0, 0,0,0,  -1, 0, 0, 0,0,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR2( -1, 0, 0, 0, 0, 0,  0,0,0, 0,0,0,   0, 0, 0, 0,1,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR3( -1, 0, 0, 0, 0, 0,  1,0,0, 0,0,0,   1, 0, 0, 0,0,1, 0,0,0, 0, 0,0,0,0, 0,0,0,0);

    int ex = CTE.ex;

    if (ex==0) CTE.px = 0; //the extra var doesn't matter

    if (mapRR[ex].count(CTE) == 1) {
        mapRR[ex][CTE].ref = true;
        ++mapRR[ex][CTE].count;
        return mapRR[ex][CTE].id;
    }

    RRid ids;

    //ex>0
    if (ex>0) {
        int px = CTE.px;
        if (px>0)
        ids.pid1 = CTEbx(CTE+CTERR1, mapRR, ABx);
        if (ABx) ids.pid2 = CTEbx(CTE+CTERR2, mapRR, ABx);
        ids.pid3 = CTEbx(CTE+CTERR3, mapRR, ABx);
    }

    //set id of element
    ids.id = mapRR[ex].size() +1;
    mapRR[ex][CTE] = ids;
    mapRR[ex][CTE].ref = true;
    return ids.id;
}

int CTEby(integralSph CTE, map<integralSph, RRid> * mapRR, bool ABy) {
    //                         ex ey ez fx fy fz                   rx ry rz  a b p  c d q  m
    const integralSph CTERR1(  0,-1, 0, 0, 0, 0,  0,-1,0, 0,0,0,   0,-1, 0, 0,0,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR2(  0,-1, 0, 0, 0, 0,  0, 0,0, 0,0,0,   0, 0, 0, 0,1,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR3(  0,-1, 0, 0, 0, 0,  0, 1,0, 0,0,0,   0, 1, 0, 0,0,1, 0,0,0, 0, 0,0,0,0, 0,0,0,0);

    int ey = CTE.ey;

    if (ey==0) CTE.py = 0; //the extra var doesn't matter

    if (mapRR[ey].count(CTE) == 1) {
        mapRR[ey][CTE].ref = true;
        ++mapRR[ey][CTE].count;
        return mapRR[ey][CTE].id;
    }

    RRid ids;


    //ey>0
    if (ey>0) {
        int py = CTE.py;
        if (py>0)
        ids.pid1 = CTEby(CTE+CTERR1, mapRR, ABy);
        if (ABy) ids.pid2 = CTEby(CTE+CTERR2, mapRR, ABy);
        ids.pid3 = CTEby(CTE+CTERR3, mapRR, ABy);
    }

    //set id of element
    ids.id = mapRR[ey].size() +1;
    mapRR[ey][CTE] = ids;
    mapRR[ey][CTE].ref = true;
    return ids.id;
}

int CTEbz(integralSph CTE, map<integralSph, RRid> * mapRR, bool ABz) {
    //                        ex ey ez fx fy fz                   rx ry rz  a b p  c d q  m
    const integralSph CTERR1( 0, 0,-1, 0, 0, 0,  0,0,-1, 0,0,0,   0, 0,-1, 0,0,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR2( 0, 0,-1, 0, 0, 0,  0,0, 0, 0,0,0,   0, 0, 0, 0,1,0, 0,0,0, 0, 0,0,0,0, 0,0,0,0);
    const integralSph CTERR3( 0, 0,-1, 0, 0, 0,  0,0, 1, 0,0,0,   0, 0, 1, 0,0,1, 0,0,0, 0, 0,0,0,0, 0,0,0,0);

    int ez = CTE.ez;

    if (ez==0) CTE.pz = 0; //the extra var doesn't matter

    if (mapRR[ez].count(CTE) == 1)  {
        mapRR[ez][CTE].ref = true;
        ++mapRR[ez][CTE].count;
        return mapRR[ez][CTE].id;
    }

    RRid ids;

    //ez>0
    if (ez>0) {
        int pz = CTE.pz;
        if (pz>0)
        ids.pid1 = CTEbz(CTE+CTERR1, mapRR, ABz);
        if (ABz) ids.pid2 = CTEbz(CTE+CTERR2, mapRR, ABz);
        ids.pid3 = CTEbz(CTE+CTERR3, mapRR, ABz);
    }

    //set id of element
    ids.id = mapRR[ez].size() +1;
    mapRR[ez][CTE] = ids;
    mapRR[ez][CTE].ref = true;
    return ids.id;
}


//contracted2 Rxcc
int Rxcc(integralSph R, map<integralSph, RRid> * mapRR, bool ABx, int CDx, bool ACx) {
    //                       ex ey ez fx fy fz                x  y  z  a b p  c d q  m
    const integralSph MMDX1( 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0, -1, 0, 0, 0,1,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDX2( 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0, -1, 0, 0, 0,0,0, 0,1,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDX3( 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0, -1, 0, 0, 0,0,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDX4( 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0, -2, 0, 0, 0,0,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);

    int rx = R.rx;

    if (mapRR[rx].count(R) == 1) {
        mapRR[rx][R].ref = true;
        ++mapRR[rx][R].count;
        return mapRR[rx][R].id;
    }

    RRid ids;

    //x>0
    if (rx>0) {
        if (ABx)  ids.pid1 = Rxcc(R+MMDX1, mapRR, ABx, CDx, ACx);
        if (CDx)  ids.pid2 = Rxcc(R+MMDX2, mapRR, ABx, CDx, ACx);
        if (ACx)  ids.pid3 = Rxcc(R+MMDX3, mapRR, ABx, CDx, ACx);

        if (rx>1) ids.pid4 = Rxcc(R+MMDX4, mapRR, ABx, CDx, ACx);
    }

    //set id of element
    ids.id = mapRR[rx].size() +1;
    mapRR[rx][R] = ids;
    mapRR[rx][R].ref = true;
    return ids.id;
}

//contracted2 Rycc
int Rycc(integralSph R, map<integralSph, RRid> * mapRR, bool ABy, bool CDy, bool ACy) {
    //                       ex ey ez fx fy fz               x  y  z  a b p  c d q  m
    const integralSph MMDY1( 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0, 0,-1, 0, 0,1,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDY2( 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0, 0,-1, 0, 0,0,0, 0,1,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDY3( 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0, 0,-1, 0, 0,0,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDY4( 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0, 0,-2, 0, 0,0,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);

    int ry = R.ry;
    if (mapRR[ry].count(R) == 1) {
        mapRR[ry][R].ref = true;
        ++mapRR[ry][R].count;
        return mapRR[ry][R].id;
    }

    RRid ids;

    //y>0
    if (ry>0) {
        if (ABy)  ids.pid1 = Rycc(R+MMDY1, mapRR, ABy, CDy, ACy);
        if (CDy)  ids.pid2 = Rycc(R+MMDY2, mapRR, ABy, CDy, ACy);
        if (ACy)  ids.pid3 = Rycc(R+MMDY3, mapRR, ABy, CDy, ACy);

        if (ry>1) ids.pid4 = Rycc(R+MMDY4, mapRR, ABy, CDy, ACy);
    }


    ids.id = mapRR[ry].size() +1;
    mapRR[ry][R] = ids;
    mapRR[ry][R].ref = true;
    return ids.id;
}

//contracted2 Rzcc
int Rzcc(integralSph R, map<integralSph, RRid> * mapRR, bool ABz, bool CDz, bool ACz) {
    //                       ex ey ez fx fy fz                x  y  z  a b p  c d q  m
    const integralSph MMDZ1( 0, 0, 0, 0, 0, 0,  0,0,0,0,0,0, 0, 0,-1, 0,1,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDZ2( 0, 0, 0, 0, 0, 0,  0,0,0,0,0,0, 0, 0,-1, 0,0,0, 0,1,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDZ3( 0, 0, 0, 0, 0, 0,  0,0,0,0,0,0, 0, 0,-1, 0,0,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);
    const integralSph MMDZ4( 0, 0, 0, 0, 0, 0,  0,0,0,0,0,0, 0, 0,-2, 0,0,0, 0,0,0, 1, 0,0,0,0, 0,0,0,0);

    int rz = R.rz;
    if (mapRR[rz].count(R) == 1) {
        mapRR[rz][R].ref = true;
        ++mapRR[rz][R].count;
        return mapRR[rz][R].id;
    }

    RRid ids;

    //y>0
    if (rz>0) {
        if (ABz)  ids.pid1 = Rzcc(R+MMDZ1, mapRR, ABz, CDz, ACz);
        if (CDz)  ids.pid2 = Rzcc(R+MMDZ2, mapRR, ABz, CDz, ACz);
        if (ACz)  ids.pid3 = Rzcc(R+MMDZ3, mapRR, ABz, CDz, ACz);

        if (rz>1) ids.pid4 = Rzcc(R+MMDZ4, mapRR, ABz, CDz, ACz);
    }


    ids.id = mapRR[rz].size() +1;
    mapRR[rz][R] = ids;
    mapRR[rz][R].ref = true;
    return ids.id;
}



// RR wrap functions


void SRRa(map< integralSph, RRid > & map0, map< integralSph, RRid > * mapSRR) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = SRRa(it->first, mapSRR);
        it->second.pid2 = it->first.la;
        mapSRR[it->first.la][it->first].mark = true;
        mapSRR[it->first.la][it->first].ref  = false;
        mapSRR[it->first.la][it->first].pidf = it->second.id;
        //cout << endl;
    }
}

void SRRb(map< integralSph, RRid > & map0, map< integralSph, RRid > * mapSRR, bool ABx, bool ABy, bool ABz) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = SRRb(it->first, mapSRR, ABx, ABy, ABz);
        it->second.pid2 = it->first.lb;
        mapSRR[it->first.lb][it->first].mark = true;
        mapSRR[it->first.lb][it->first].ref  = false;
        mapSRR[it->first.lb][it->first].pidf = it->second.id;
    }
}

void SRRc(map< integralSph, RRid > & map0, map< integralSph, RRid > * mapSRR) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = SRRc(it->first, mapSRR);
        it->second.pid2 = it->first.lc;
        mapSRR[it->first.lc][it->first].mark = true;
        mapSRR[it->first.lc][it->first].ref  = false;
        mapSRR[it->first.lc][it->first].pidf = it->second.id;
    }
}

void SRRd(map< integralSph, RRid > & map0, map< integralSph, RRid > * mapSRR, bool CDx, bool CDy, bool CDz) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = SRRd(it->first, mapSRR, CDx, CDy, CDz);
        it->second.pid2 = it->first.ld;
        mapSRR[it->first.ld][it->first].mark = true;
        mapSRR[it->first.ld][it->first].ref  = false;
        mapSRR[it->first.ld][it->first].pidf = it->second.id;
    }
}


void CTEkx(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapCTE, bool CDx) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = CTEkx(it->first, mapCTE, CDx);
        it->second.pid2 = it->first.fx;
        mapCTE[it->first.fx][it->first].mark = true;
        mapCTE[it->first.fx][it->first].ref  = false;
        mapCTE[it->first.fx][it->first].pidf = it->second.id;
    }
}

void CTEky(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapCTE, bool CDy) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = CTEky(it->first, mapCTE, CDy);
        it->second.pid2 = it->first.fy;
        mapCTE[it->first.fy][it->first].mark = true;
        mapCTE[it->first.fy][it->first].ref  = false;
        mapCTE[it->first.fy][it->first].pidf = it->second.id;
    }
}

void CTEkz(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapCTE, bool CDz) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = CTEkz(it->first, mapCTE, CDz);
        it->second.pid2 = it->first.fz;
        mapCTE[it->first.fz][it->first].mark = true;
        mapCTE[it->first.fz][it->first].ref  = false;
        mapCTE[it->first.fz][it->first].pidf = it->second.id;
    }
}

void CTEbx(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapCTE, bool ABx) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = CTEbx(it->first, mapCTE, ABx);
        it->second.pid2 = it->first.ex;
        mapCTE[it->first.ex][it->first].mark = true;
        mapCTE[it->first.ex][it->first].ref  = false;
        mapCTE[it->first.ex][it->first].pidf = it->second.id;
    }
}

void CTEby(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapCTE, bool ABy) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = CTEby(it->first, mapCTE, ABy);
        it->second.pid2 = it->first.ey;
        mapCTE[it->first.ey][it->first].mark = true;
        mapCTE[it->first.ey][it->first].ref  = false;
        mapCTE[it->first.ey][it->first].pidf = it->second.id;
    }
}

void CTEbz(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapCTE, bool ABz) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = CTEbz(it->first, mapCTE, ABz);
        it->second.pid2 = it->first.ez;
        mapCTE[it->first.ez][it->first].mark = true;
        mapCTE[it->first.ez][it->first].ref  = false;
        mapCTE[it->first.ez][it->first].pidf = it->second.id;
    }
}


void Rxcc(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapR, bool ABx, bool CDx, bool ACx) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = Rxcc(it->first, mapR, ABx, CDx, ACx);
        it->second.pid2 = it->first.rx;
        mapR[it->first.rx][it->first].mark = true;
        mapR[it->first.rx][it->first].ref  = false;
        mapR[it->first.rx][it->first].pidf = it->second.id;
    }
}

void Rycc(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapR, bool ABy, bool CDy, bool ACy) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = Rycc(it->first, mapR, ABy, CDy, ACy);
        it->second.pid2 = it->first.ry;
        mapR[it->first.ry][it->first].mark = true;
        mapR[it->first.ry][it->first].ref  = false;
        mapR[it->first.ry][it->first].pidf = it->second.id;
    }
}

void Rzcc(map<integralSph, RRid> & map0, map<integralSph, RRid> * mapR, bool ABz, bool CDz, bool ACz) {
    map<integralSph, RRid>::iterator it;

    for (it=map0.begin(); it!=map0.end(); ++it) {
        it->second.pid1 = Rzcc(it->first, mapR, ABz, CDz, ACz);
        it->second.pid2 = it->first.rz;
        mapR[it->first.rz][it->first].mark = true;
        mapR[it->first.rz][it->first].ref  = false;
        mapR[it->first.rz][it->first].pidf = it->second.id;
    }
}

*/
