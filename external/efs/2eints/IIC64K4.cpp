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



    This is an implementation of the K4+MIRROR algorithm to compute 2-electron integrals.
    If you use it in your research, please cite the original  papers:

    "An algorithm for the efficient evaluation of two-electron repulsion integrals over contracted Gaussian-type basis functions",
    Jaime Axel Rosal Sandberg, Zilvinas Rinkevicius, J. Chem. Phys. 137, 234105 (2012); http://dx.doi.org/10.1063/1.4769730

    "New recurrence relations for analytic evaluation of two-electron repulsion integrals over highly contracted gaussian-type orbitals",
    Jaime Axel Rosal Sandberg, Zilvinas Rinkevicius, In preparation
*/

#include <string.h>
#include "../defs.hpp"

#include "../2eints/IICinit.hpp"
#include "../2eints/quimera.hpp"
#include "../libquimera/libquimera.hpp"
#include "../basis/SPprototype.hpp"
#include "../integrals/rotations.hpp"
#include "../integrals/atomprod.hpp"
#include "../math/gamma.hpp"

using namespace LibIGamma;


static const UI8 PF = 8; //prefetch distance
static const int DPC = DOUBLES_PER_CACHE_LINE;


static const double im2[32] = {
    1.,     1./3.,  1./5.,  1./7.,  1./9.,  1./11., 1./13., 1./15.,
    1./17., 1./19., 1./21., 1./23., 1./25., 1./27., 1./29., 1./31.,
    1./33., 1./35., 1./37., 1./39., 1./41., 1./43., 1./45., 1./47.,
    1./49., 1./51., 1./53., 1./55., 1./57., 1./59., 1./61., 1./63.,
};

// Gamma[n+1/2]/Sqrt[PI]
static const double iii[] =
{0.5, 0.75, 1.875, 6.5625, 29.53125, 162.421875, 1055.7421875, 7918.06640625, 67303.564453125,
639383.8623046875, 6.7135305541992188e6, 7.7205601373291016e7,
 9.6507001716613770e8, 1.3028445231742859e10,
 1.8891245586027145e11, 2.9281430658342075e12,
 4.8314360586264424e13, 8.4550131025962743e14,
 1.5641774239803107e16, 3.0501459767616059e17,
 6.2527992523612922e18, 1.3443518392576778e20,
 3.0247916383297751e21, 7.1082603500749715e22,
 1.7415237857683680e24, 4.4408856537093384e25,
 1.1768346982329747e27, 3.2362954201406804e28,
 9.2234419474009391e29, 2.7209153744832770e31,
 8.2987918921739949e32, 2.6141194460348084e34};


template <int labcd> inline void GammaD(double & e, double & m, double KpqR2) {

    double ** gamma = IncompleteGammas[labcd].gamma_table;

    if (KpqR2+vg_step>vg_max) {
        double ir2 = 1/(KpqR2);
        m = PI3 * sqrt(ir2);

        if (labcd>0) {
            double ir4  = ir2  * ir2;
            double ir8  = ir4  * ir4;
            double ir16 = ir8  * ir8;
            double ir32 = ir16 * ir16;
            double ir64 = ir32 * ir32;

            if (labcd&1)  m *= ir2;
            if (labcd&2)  m *= ir4;
            if (labcd&4)  m *= ir8;
            if (labcd&8)  m *= ir16;
            if (labcd&16) m *= ir32;
            if (labcd&32) m *= ir64;

            m *= iii[labcd-1];
        }

        e = 0;
    }
    else {
        double p = ivg_step * KpqR2;
        int pos = int(p+0.5);
        double x0 = vg_step*double(pos);

        double Ax1 = x0 - KpqR2;
        double Ax2 = 0.5 * Ax1*Ax1;
        double Ax3 = 0.33333333333333333333 * Ax1*Ax2;

        //m = PI54 * (gamma[pos][labcd+1] + gamma[pos][labcd+2] * Ax1 + gamma[pos][labcd+3] * Ax2 + gamma[pos][labcd+4] * Ax3);
        m = PI54 * (gamma[pos][1] + gamma[pos][2] * Ax1 + gamma[pos][3] * Ax2 + gamma[pos][4] * Ax3);
        e = PI54 * gamma[pos][0] * (1+Ax1+Ax2+Ax3);
    }
}


template <int Labcd> void Gamma(cacheline64 * e, cacheline64 * m, double iKpq, const cacheline64 & R2, const cacheline64 & W) {

    double Kpq = 0.5/iKpq;

    cacheline64 KpqR2;
    KpqR2 = R2 * Kpq;

    //only non-vector part of the algorithm
    for (int k=0; k<DPC; ++k) {
        double ed, md;
        GammaD<Labcd>(ed, md, KpqR2(k));
        (*e)(k) = ed;
        (*m)(k) = md;
    }

    //this could be precomputed and stored
    double Kpqn = sqrt(Kpq); {
        Kpq += Kpq;

        double Kpq2  = Kpq   * Kpq;
        double Kpq4  = Kpq2  * Kpq2;
        double Kpq8  = Kpq4  * Kpq4;
        double Kpq16 = Kpq8  * Kpq8;
        double Kpq32 = Kpq16 * Kpq16;

        if (Labcd&1)  Kpqn *= Kpq;
        if (Labcd&2)  Kpqn *= Kpq2;
        if (Labcd&4)  Kpqn *= Kpq4;
        if (Labcd&8)  Kpqn *= Kpq8;
        if (Labcd&16) Kpqn *= Kpq16;
        if (Labcd&32) Kpqn *= Kpq32;
    }

    cacheline64 WW;
    WW = W * Kpqn;

    if (Labcd>0) *e *= WW;
    *m *= WW;

    /*
    //AER downward recursion
    for (int i=Labcd-1; i>=0; --i) e[i] = e[i+1] * iKpq;

    //downward recursion
    for (int i=Labcd-1; i>=0; --i) m[i] = (R2 * m[i+1] + e[i]) * (1./double(2*i+1));
    */
}

//gamma for 1-center integrals
void gamma(int Labcd, double * e, double * m, double Kpq) {

    double Kpqn = PI54 * sqrt(Kpq);
    Kpq += Kpq;

    for (int l=0; l<Labcd; ++l) {
        e[l] =                       Kpqn; //this is the exponential
        m[l] = (1./double(2*l+1)) * Kpqn; //this is gamma
        Kpqn *= Kpq;
    }

    //last gamma
    m[Labcd] = (1./double(2*Labcd+1)) * Kpqn;
}


template <bool SAB, bool SCD> void ERIroutine::CalcGammas (const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, p_ERIbuffer & buffer, bool OnlyJ) const {

    void (*cGamma)(cacheline64 * e, cacheline64 * m, double Kpq, const cacheline64 & R2, const cacheline64 & W) ;

    if (Lt== 0) cGamma = Gamma< 0>;
    if (Lt== 1) cGamma = Gamma< 1>;
    if (Lt== 2) cGamma = Gamma< 2>;
    if (Lt== 3) cGamma = Gamma< 3>;
    if (Lt== 4) cGamma = Gamma< 4>;
    if (Lt== 5) cGamma = Gamma< 5>;
    if (Lt== 6) cGamma = Gamma< 6>;
    if (Lt== 7) cGamma = Gamma< 7>;
    if (Lt== 8) cGamma = Gamma< 8>;
    if (Lt== 9) cGamma = Gamma< 9>;
    if (Lt==10) cGamma = Gamma<10>;
    if (Lt==11) cGamma = Gamma<11>;
    if (Lt==12) cGamma = Gamma<12>;
    if (Lt==13) cGamma = Gamma<13>;
    if (Lt==14) cGamma = Gamma<14>;
    if (Lt==15) cGamma = Gamma<15>;
    if (Lt==16) cGamma = Gamma<16>;
    if (Lt==17) cGamma = Gamma<17>;
    if (Lt==18) cGamma = Gamma<18>;
    if (Lt==19) cGamma = Gamma<19>;
    if (Lt==20) cGamma = Gamma<20>;


    const int JA = ABp.Ja;
    const int JB = ABp.Jb;
    const int JC = CDp.Ja;
    const int JD = CDp.Jb;

    const PrimitiveSet & PSab = ABp.Psets[ET.nKab];
    const PrimitiveSet & PScd = CDp.Psets[ET.nKcd];

    cacheline64 * KAB = (cacheline64*)buffer.KAB;
    cacheline64 * KCD = (cacheline64*)buffer.KCD;

    if (SAB) {
        int ab=0;
        int addr0 = ET.ap12[0]*ET.nK2ab;
        int addr1 = ET.ap12[1]*ET.nK2ab;
        int addr2 = ET.ap12[2]*ET.nK2ab;
        int addr3 = ET.ap12[3]*ET.nK2ab;
        int addr4 = ET.ap12[4]*ET.nK2ab;
        int addr5 = ET.ap12[5]*ET.nK2ab;
        int addr6 = ET.ap12[6]*ET.nK2ab;
        int addr7 = ET.ap12[7]*ET.nK2ab;

        for (int b=0; b<PSab.nKb; ++b) {
            for (int a=0; a<PSab.nKa[b]; ++a) {
                int ba = b*ABp.Ka + a;
                KAB[ab].set( ET.wSP12[addr0+ba], ET.wSP12[addr1+ba], ET.wSP12[addr2+ba], ET.wSP12[addr3+ba],
                             ET.wSP12[addr4+ba], ET.wSP12[addr5+ba], ET.wSP12[addr6+ba], ET.wSP12[addr7+ba] );

                ++ab;
            }
        }

    }


    if (SCD) {
        int cd = 0;
        int addr0 = ET.ap34[0]*ET.nK2cd;
        int addr1 = ET.ap34[1]*ET.nK2cd;
        int addr2 = ET.ap34[2]*ET.nK2cd;
        int addr3 = ET.ap34[3]*ET.nK2cd;
        int addr4 = ET.ap34[4]*ET.nK2cd;
        int addr5 = ET.ap34[5]*ET.nK2cd;
        int addr6 = ET.ap34[6]*ET.nK2cd;
        int addr7 = ET.ap34[7]*ET.nK2cd;

        for (int d=0; d<PScd.nKb; ++d) {
            for (int c=0; c<PScd.nKa[d]; ++c) {
                int dc = d*CDp.Ka + c;

                KCD[cd].set( ET.wSP34[addr0+dc], ET.wSP34[addr1+dc], ET.wSP34[addr2+dc], ET.wSP34[addr3+dc],
                             ET.wSP34[addr4+dc], ET.wSP34[addr5+dc], ET.wSP34[addr6+dc], ET.wSP34[addr7+dc] );
                ++cd;
            }
        }
    }




    const cacheline64 & ABz = vars8.ABz;
    const cacheline64 & CDy = vars8.CDy;
    const cacheline64 & CDz = vars8.CDz;

    const cacheline64 & ACx = vars8.ACx;
    const cacheline64 & ACy = vars8.ACy;
    const cacheline64 & ACz = vars8.ACz;


    cacheline64 AC2;


    if (SAB && SCD) {
        AC2  = ACx*ACx;
    }
    else if (!SAB && SCD) {
        AC2  = ACx*ACx + ACy*ACy;
    }
    else if (!SAB && !SCD) {
        AC2  = ACz*ACz;
    }


    cacheline64 * F0 = (cacheline64*)buffer.F0;


    int dc = 0;

    for (int d=0; d<PScd.nKb; ++d) {
        for (int c=0; c<PScd.nKa[d]; ++c) {

            double rcd  = CDp.BCP.r[d][c];

            cacheline64 AQy;
            cacheline64 AQz;
            cacheline64 AQ2;

            if (SAB && SCD) {
                AQy = ACy + CDy * rcd;
                AQz = ACz + CDz * rcd;
                AQ2 = AC2 + AQy*AQy;
            }
            else if (!SAB && SCD) {
                AQz = ACz + CDz * rcd;
                AQ2 = AC2 + AQz*AQz;
            }


            int ba = 0;
            for (int b=0; b<PSab.nKb; ++b) {
                for (int a=0; a<PSab.nKa[b]; ++a) {
                    double rab  = ABp.BCP.r[b][a];

                    cacheline64 ww;

                    ww.set(ABp.BCP.K[b][a]*CDp.BCP.K[d][c]);

                    if (SAB && SCD) {
                        cacheline64 PQz;
                        cacheline64 PQ2;

                        PQz = AQz - ABz*rab;
                        PQ2 = AQ2 + PQz*PQz; //gPQ2 = gAQ2 - gAQAB*rab + gAB2*(rab*rab);
                        ww *=  KAB[ba] * KCD[dc];

                        double iKpq  = (ABp.BCP.k[b][a] + CDp.BCP.k[d][c]);
                        cGamma (F0, F0+1, iKpq, PQ2, ww);
                    }
                    else if (!SAB && SCD) {
                        ww *=  KCD[dc];

                        double iKpq  = (ABp.BCP.k[b][a] + CDp.BCP.k[d][c]);
                        cGamma (F0, F0+1, iKpq, AQ2, ww);
                    }
                    else if (!SAB && !SCD) {

                        double iKpq  = (ABp.BCP.k[b][a] + CDp.BCP.k[d][c]);
                        cGamma (F0, F0+1, iKpq, AC2, ww);
                    }


                    F0 += 2;

                    ++ba;
                }
            }

            ++dc;
        }
    }

}


//for AAAA
void ERIroutine::ContractCDR_GC(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, double * __restrict__ uv_m_st, p_ERIbuffer & buffer, bool OnlyJ) const {

    const int JA = ABp.Ja;
    const int JB = ABp.Jb;
    const int JC = CDp.Ja;
    const int JD = CDp.Jb;


    for (int i=0; i<JA*JB*JC*JD*memK1J4e; ++i) buffer.k1j4e[i] = 0;
    for (int i=0; i<JA*JB*JC*JD*memK1J4f; ++i) buffer.k1j4f[i] = 0;

    for (int d=0; d<CDp.Kb; ++d) {

        for (int i=0; i<JA*JB*JC*memK2J3e; ++i) buffer.k2j3e[i] = 0;
        for (int i=0; i<JA*JB*JC*memK2J3f; ++i) buffer.k2j3f[i] = 0;

        for (int c=0; c<CDp.Ka; ++c) {

            for (int i=0; i<JA*JB*memK3J2e; ++i) buffer.k3j2e[i] = 0;
            for (int i=0; i<JA*JB*memK3J2f; ++i) buffer.k3j2f[i] = 0;

            for (int b=0; b<ABp.Kb; ++b) {

                for (int i=0; i<JA*memK4J1e; ++i) buffer.k4j1e[i] = 0;
                for (int i=0; i<JA*memK4J1f; ++i) buffer.k4j1f[i] = 0;

                for (int a=0; a<ABp.Ka; ++a) {

                    double iKpq  = (ABp.BCP.k[b][a] + CDp.BCP.k[d][c]);
                    double Kpq = 0.5/iKpq;

                    double e[Lt];
                    double m[Lt+1];

                    double em[2*Lt+1];

                    gamma(Lt, em, em+Lt, Kpq);

                    /*
                    if (UseCASE && OnlyJ) {
                        double Kpqw2 = (Kpq*case_w2)/(Kpq+case_w2);

                        double e[Lt];
                        double m[Lt+1];

                        gamma(Lt, buffer.f0, Kpqw2);

                        for (int l=0; l<2*Lt+1; ++l) {
                            buffer.j0[l] -= buffer.f0[l];
                        }
                    }
                    */

                    double KK = ABp.BCP.K[b][a]*CDp.BCP.K[d][c];

                    for (int i=0; i<2*Lt+1; ++i) em[i] *= KK;


                    //K4J1
                    double * E = buffer.k4j1e;

                    for (int ja=0; ja<JA; ++ja) {
                        double Vv[32];
                        Vv[0] = ABp.BCP.Na[a][ja];
                        for (int i=1; i<=maxV; ++i)
                            Vv[i] = Vv[i-1] * ABp.BCP.k[b][a];

                        //K4 CS
                        for (int i=0; i<ninstrK4[K4E]; ++i) {
                            UI32 op1  = nseqK4[K4E][i].op1;
                            UI32 dest = nseqK4[K4E][i].dest;

                            UI16 aux  = nseqK4[K4E][i].aux;
                            //double V = ABp.BCP.V[b][a][aux] * ABp.BCP.Na[a][ja];

                            E[dest] += em[op1] * Vv[aux];
                        }
                        E += memK4J1e;
                    }

                    double * F = buffer.k4j1f;

                    for (int ja=0; ja<JA; ++ja) {
                        double Vv[32];
                        Vv[0] = ABp.BCP.Na[a][ja];
                        for (int i=1; i<=maxV; ++i)
                            Vv[i] = Vv[i-1] * ABp.BCP.k[b][a];

                        //K4 CS
                        for (int i=0; i<ninstrK4[K4F]; ++i) {
                            UI32 op1  = nseqK4[K4F][i].op1;
                            UI32 dest = nseqK4[K4F][i].dest;

                            UI16 aux  = nseqK4[K4F][i].aux;
                            //double V = ABp.BCP.V[b][a][aux] * ABp.BCP.Na[a][ja];

                            F[dest] += em[op1] * Vv[aux];
                        }
                        F += memK4J1f;
                    }

                }

                double * e = buffer.k4j1e;
                double * f = buffer.k4j1f;

                double * E = buffer.k3j2e;
                double * F = buffer.k3j2f;

                //K3J1
                for (int ja=0; ja<JA; ++ja) {

                    //AERR K3
                    for (int i=0; i<ninstrK4[AERR3]; ++i) {
                        UI32 dest = nseqK4[AERR3][i].dest;
                        UI32 op1  = nseqK4[AERR3][i].op1;
                        UI32 op2  = nseqK4[AERR3][i].op2;

                        e[dest] = e[op1] + e[op2] * CDp.BCP.k[d][c];
                    }

                    //CDR K3
                    for (int i=0; i<ninstrK4[CDR3]; ++i) {
                        UI32 dest = nseqK4[CDR3][i].dest;
                        UI32 ope  = nseqK4[CDR3][i].ope;
                        UI16 aux  = nseqK4[CDR3][i].aux;
                        double i2m1 = 1./double(2*aux + 1);

                        f[dest] = (e[ope]) * i2m1;
                    }


                    //K3J2
                    for (int jb=0; jb<JB; ++jb) {
                        double Uu[32];
                        Uu[0] = ABp.BCP.Nb[b][jb];
                        for (int i=1; i<=maxU; ++i)
                            Uu[i] = Uu[i-1] * ABp.BCP.b1[b];

                        for (int i=0; i<ninstrK4[K3E]; ++i) {
                            UI32 op1  = nseqK4[K3E][i].op1;
                            UI32 dest = nseqK4[K3E][i].dest;
                            UI16 aux  = nseqK4[K3E][i].aux;
                            //double U = ABp.BCP.U[b][aux] * ABp.BCP.Nb[b][jb];
                            E[dest] += e[op1] * Uu[aux];
                        }
                        E += memK3J2e;
                    }

                    //K3J2
                    for (int jb=0; jb<JB; ++jb) {
                        double Uu[32];
                        Uu[0] = ABp.BCP.Nb[b][jb];
                        for (int i=1; i<=maxU; ++i)
                            Uu[i] = Uu[i-1] * ABp.BCP.b1[b];

                        for (int i=0; i<ninstrK4[K3F]; ++i) {
                            UI32 op1  = nseqK4[K3F][i].op1;
                            UI32 dest = nseqK4[K3F][i].dest;
                            UI16 aux  = nseqK4[K3F][i].aux;
                            //double U = ABp.BCP.U[b][aux] * ABp.BCP.Nb[b][jb];
                            F[dest] += f[op1] * Uu[aux];
                        }
                        F += memK3J2f;
                    }

                    e += memK4J1e;
                    f += memK4J1f;
                }
            }

            //K2J2
            double * e = buffer.k3j2e;
            double * f = buffer.k3j2f;

            double * E = buffer.k2j3e;
            double * F = buffer.k2j3f;

            for (int jab=0; jab<JA*JB; ++jab) {

                //AERR K2
                for (int i=0; i<ninstrK4[AERR2]; ++i) {
                    UI32 dest = nseqK4[AERR2][i].dest;
                    UI32 op1  = nseqK4[AERR2][i].op1;
                    UI32 op2  = nseqK4[AERR2][i].op2;

                    e[dest] = e[op1] + e[op2] * CDp.BCP.k[d][c];
                }

                //CDR K3
                for (int i=0; i<ninstrK4[CDR2]; ++i) {
                    UI32 dest = nseqK4[CDR2][i].dest;
                    UI32 ope  = nseqK4[CDR2][i].ope;
                    UI16 aux  = nseqK4[CDR2][i].aux;
                    double i2m1 = 1./double(2*aux + 1);

                    f[dest] = (e[ope]) * i2m1;
                }


                //K2J3
                for (int jc=0; jc<JC; ++jc) {
                    double Tt[32];
                    Tt[0] = CDp.BCP.Na[c][jc];
                    for (int i=1; i<=maxT; ++i)
                        Tt[i] = Tt[i-1] * CDp.BCP.k[d][c];

                    for (int i=0; i<ninstrK4[K2E]; ++i) {
                        UI32 op1  = nseqK4[K2E][i].op1;
                        UI32 dest = nseqK4[K2E][i].dest;
                        UI16 aux  = nseqK4[K2E][i].aux;
                        //double T = CDp.BCP.V[d][c][aux] * CDp.BCP.Na[c][jc];
                        E[dest] += e[op1] * Tt[aux];
                    }
                    E += memK2J3e;
                }

                //K2J3
                for (int jc=0; jc<JC; ++jc) {
                    double Tt[32];
                    Tt[0] = CDp.BCP.Na[c][jc];
                    for (int i=1; i<=maxT; ++i)
                        Tt[i] = Tt[i-1] * CDp.BCP.k[d][c];

                    for (int i=0; i<ninstrK4[K2F]; ++i) {
                        UI32 op1  = nseqK4[K2F][i].op1;
                        UI32 dest = nseqK4[K2F][i].dest;
                        UI16 aux  = nseqK4[K2F][i].aux;
                        //double T = CDp.BCP.V[d][c][aux] * CDp.BCP.Na[c][jc];
                        F[dest] += f[op1] * Tt[aux];
                    }
                    F += memK2J3f;
                }

                e += memK3J2e;
                f += memK3J2f;
            }

        }

        //K1J3
        double * e = buffer.k2j3e;
        double * f = buffer.k2j3f;

        double * E = buffer.k1j4e;
        double * F = buffer.k1j4f;

        for (int jabc=0; jabc<JA*JB*JC; ++jabc) {

            //AERR K1
            for (int i=0; i<ninstrK4[AERR1]; ++i) {
                UI32 dest = nseqK4[AERR1][i].dest;
                UI32 op1  = nseqK4[AERR1][i].op1;
                UI32 op2  = nseqK4[AERR1][i].op2;

                e[dest] = e[op1] + e[op2];
            }

            //CDR K1
            for (int i=0; i<ninstrK4[CDR1]; ++i) {
                UI32 dest = nseqK4[CDR1][i].dest;
                UI32 ope  = nseqK4[CDR1][i].ope;
                UI16 aux  = nseqK4[CDR1][i].aux;
                double i2m1 = 1./double(2*aux + 1);

                f[dest] = (e[ope]) * i2m1;
            }


            //K1J4
            for (int jd=0; jd<JD; ++jd) {
                double Ss[32];
                Ss[0] = CDp.BCP.Nb[d][jd];
                for (int i=1; i<=maxS; ++i)
                    Ss[i] = Ss[i-1] * CDp.BCP.b1[d];

                for (int i=0; i<ninstrK4[K1E]; ++i) {
                    UI32 op1  = nseqK4[K1E][i].op1;
                    UI32 dest = nseqK4[K1E][i].dest;
                    UI16 aux  = nseqK4[K1E][i].aux;
                    //double S = CDp.BCP.U[d][aux] * CDp.BCP.Nb[d][jd];
                    E[dest] += e[op1] * Ss[aux];
                }
                E += memK1J4e;
            }

            //K1J4
            for (int jd=0; jd<JD; ++jd) {
                double Ss[32];
                Ss[0] = CDp.BCP.Nb[d][jd];
                for (int i=1; i<=maxS; ++i)
                    Ss[i] = Ss[i-1] * CDp.BCP.b1[d];

                for (int i=0; i<ninstrK4[K1F]; ++i) {
                    UI32 op1  = nseqK4[K1F][i].op1;
                    UI32 dest = nseqK4[K1F][i].dest;
                    UI16 aux  = nseqK4[K1F][i].aux;
                    //double S = CDp.BCP.U[d][aux] * CDp.BCP.Nb[d][jd];
                    F[dest] += f[op1] * Ss[aux];
                }
                F += memK1J4f;
            }

            e += memK2J3e;
            f += memK2J3f;
        }

    }

    //K0J4
    double * e = buffer.k1j4e;
    double * f = buffer.k1j4f;

    for (int jabcd=0; jabcd<JA*JB*JC*JD; ++jabcd) {

        //AERR K0
        for (int i=0; i<ninstrK4[AERR0]; ++i) {
            UI32 dest = nseqK4[AERR0][i].dest;
            UI32 op1  = nseqK4[AERR0][i].op1;
            UI32 op2  = nseqK4[AERR0][i].op2;

            e[dest] = e[op1] + e[op2];
        }

        //CDR K1
        for (int i=0; i<ninstrK4[CDR0]; ++i) {
            UI32 dest = nseqK4[CDR0][i].dest;
            UI32 ope  = nseqK4[CDR0][i].ope;
            UI16 aux  = nseqK4[CDR0][i].aux;
            double i2m1 = 1./double(2*aux + 1);

            f[dest] = (e[ope]) * i2m1;
        }

        e += memK1J4e;
        f += memK1J4f;
    }


    //copy to the kernel buffer

    double * ff = buffer.k1j4f;
    double * dd = uv_m_st;

    for (int jabcd=0; jabcd<JA*JB*JC*JD; ++jabcd) {
        register Op_MIRROR * s2 = eseq;

        for (int i=0; s2<nseq[KERNELS]; ++i,++s2) {
            UI32 op1  = s2->op1;
            UI32 dest = s2->dest;

            if (s2->aux==0) dd[dest] = ff[op1];
            else            dd[dest] = 0;
        }

        ff += memK1J4f;
        dd += nKernels;
    }

}


inline void LoopsContraction2(cacheline64 * __restrict__ D, const cacheline64 * __restrict__ S, const double * __restrict__ Vs,
                                const VILfma * __restrict__ IIC, int ninstr) {

    for (int i=0; i<ninstr; ++i) {
        cacheline64       * dest = D + IIC[i].dest;
        const cacheline64 * op1  = S + IIC[i].op1;
        UI16 aux                 =     IIC[i].aux;

        //__builtin_prefetch(D + s2[PF].dest, 1, 1);

        store(dest, load(dest) + load(op1) * Vs[aux]);
    }

}



template <bool SAB, bool SCD> void ERIroutine::InnerContractCDR_GC(cacheline64 * (&F0), const ERIgeometries64 & vars8, const PrimitiveSet & PSab,
                                      const ShellPairPrototype & ABp, double ikcd, double rcd, p_ERIbuffer & buffer,
                                      cacheline64 & AB2, cacheline64 & X2) const {


    const cacheline64 & ABz = vars8.ABz;
    const cacheline64 & CDy = vars8.CDy;
    const cacheline64 & CDz = vars8.CDz;

    const cacheline64 & ACx = vars8.ACx;
    const cacheline64 & ACy = vars8.ACy;
    const cacheline64 & ACz = vars8.ACz;

    const int JA = ABp.Ja;
    const int JB = ABp.Jb;

    const int nJ1 = JA;
    const int nJ2 = JA*JB;

    const int mK3f = nJ2*memK3J2f;
    const int mK4f = nJ1*memK4J1f;

    const int mK3e = nJ2*memK3J2e;
    const int mK4e = nJ1*memK4J1e;

    {
        memset(buffer.K3J2e, 0, mK3e*sizeof(cacheline64));
        memset(buffer.K3J2f, 0, mK3f*sizeof(cacheline64));

        cacheline64 AQy;
        cacheline64 AQz;

        cacheline64 X2Y2;
        cacheline64 AQ2;

        cacheline64 AQAB;


        if (SAB && SCD) {
            AQy  = ACy + CDy * rcd;
            AQz  = ACz + CDz * rcd;

            X2Y2 = X2   + AQy*AQy;
            AQ2  = X2Y2 + AQz*AQz;

            AQAB = (AQz*ABz); AQAB += AQAB;
        }
        else if (!SAB && SCD) {
            AQz  = ACz + CDz * rcd;
            AQ2  = X2 + AQz*AQz;
        }

        for (int b=0; b<PSab.nKb; ++b) {
            memset(buffer.K4J1e, 0, mK4e*sizeof(cacheline64));
            memset(buffer.K4J1f, 0, mK4f*sizeof(cacheline64));

            for (int a=0; a<PSab.nKa[b]; ++a) {

                double ikab = ABp.BCP.k[b][a];

                double iKpq  = (ikab + ikcd);

                cacheline64 PQz;
                cacheline64 R2;

                if (SAB && SCD) {
                    double rab  = ABp.BCP.r[b][a];
                    PQz = AQz - ABz*rab;
                    R2 = X2Y2 + PQz*PQz;
                }

                cacheline64 * F0e = (cacheline64*)buffer.F0e;
                cacheline64 * F0f = (cacheline64*)buffer.F0f;

                //copy the gamma of higher L and the contracted exponential
                F0e[0] = F0[0];
                F0f[0] = F0[1];

                //AERR K4 (uncontracted)
                {

                    for (int i=1; i<K4VILcode.nbi[AERR4]; ++i) {
                        cacheline64       * dest = F0e + K4VILcode.pAERR4[i].dest;
                        const cacheline64 * op1  = F0e + K4VILcode.pAERR4[i].op1;

                        *dest = *op1 * iKpq;
                    }

                }

                //CDR K4 (uncontracted downward recursion)
                {

                    for (int i=1; i<K4VILcode.nbi[CDR4]; ++i) {
                        cacheline64       * dest = F0f + K4VILcode.pCDR4[i].dest;
                        const cacheline64 * op1  = F0f + K4VILcode.pCDR4[i].op1;
                        const cacheline64 * ope  = F0e + K4VILcode.pCDR4[i].ope;
                        UI16 aux                 = K4VILcode.pCDR4[i].aux;

                        if      ( SAB &&  SCD) *dest = ( R2 * *op1 + *ope) * im2[aux];
                        else if (!SAB &&  SCD) *dest = (AQ2 * *op1 + *ope) * im2[aux];
                        else if (!SAB && !SCD) *dest = ( X2 * *op1 + *ope) * im2[aux];
                    }

                }

                //K4 E + F
                {
                    const cacheline64 * J0e = (cacheline64*)F0e;
                    cacheline64 * J1e = (cacheline64*)buffer.K4J1e;

                    const cacheline64 * J0f = (cacheline64*)F0f;
                    cacheline64 * J1f = (cacheline64*)buffer.K4J1f;

                    double Vv[maxJ][32] __attribute__((aligned(CACHE_LINE_SIZE)));

                    for (int j=0; j<JA; ++j) {
                        Vv[j][0] = ABp.BCP.Na[a][j];

                        for (int v=1; v<=maxV; ++v)
                            Vv[j][v] = Vv[j][v-1] * ABp.BCP.k[b][a];
                    }

                    #ifdef __QR_GC__
                    for (int j=0; j<min(JA, a+1); ++j) {
                        LoopsContraction2(J1e+j*memK4J1e, J0e, Vv[j], K4VILcode.pK4E, K4VILcode.nbi[K4E]);
                        LoopsContraction2(J1f+j*memK4J1f, J0f, Vv[j], K4VILcode.pK4F, K4VILcode.nbi[K4F]);
                    }
                    #else
                    for (int j=0; j<JA; ++j) {
                        LoopsContraction2(J1e+j*memK4J1e, J0e, Vv[j], K4VILcode.pK4E, K4VILcode.nbi[K4E]);
                        LoopsContraction2(J1f+j*memK4J1f, J0f, Vv[j], K4VILcode.pK4F, K4VILcode.nbi[K4F]);
                    }
                    #endif
                }

                F0 += 2;
            }

            //K3J1 //AERR K3
            {
                cacheline64 * Jf = (cacheline64*)buffer.K4J1f;
                cacheline64 * Je = (cacheline64*)buffer.K4J1e;

                //__m128d ikcd_v = _mm_load1_pd(&ikcd);
                //u_atom ikcd_v(ikcd);

                for (int nj1=0; nj1<nJ1; ++nj1) {
                    for (int i=0; i<K4VILcode.nbi[AERR3]; ++i) {
                        cacheline64       * dest = Je + K4VILcode.pAERR3[i].dest;
                        const cacheline64 * op1  = Je + K4VILcode.pAERR3[i].op1;
                        const cacheline64 * op2  = Je + K4VILcode.pAERR3[i].op2;

                        store(dest, load(op1) + load(op2) * ikcd);
                    }

                    Je += memK4J1e;
                }
            }

            //CDR3
            {

                cacheline64 AQABb;
                cacheline64 AB2b2;

                if (SAB && SCD) {
                    AQABb = AQAB * ABp.BCP.b1[b];
                    AB2b2 = AB2  * ABp.BCP.b2[b];
                }

                cacheline64 * Jf = (cacheline64*)buffer.K4J1f;
                cacheline64 * Je = (cacheline64*)buffer.K4J1e;

                for (int nj1=0; nj1<nJ1; ++nj1) {
                    for (int i=0; i<K4VILcode.nbi[CDR3]; ++i) {
                        cacheline64       * dest = Jf + K4VILcode.pCDR3[i].dest;
                        const cacheline64 * op1  = Jf + K4VILcode.pCDR3[i].op1;
                        const cacheline64 * op2  = Jf + K4VILcode.pCDR3[i].op2;
                        const cacheline64 * op3  = Jf + K4VILcode.pCDR3[i].op3;
                        const cacheline64 * ope  = Je + K4VILcode.pCDR3[i].ope;
                        UI16 aux  = K4VILcode.pCDR3[i].aux;

                        if (SAB && SCD) {
                            store(dest, (AQ2 * load(op1) - AQABb * load(op2) + AB2b2 * load(op3) + load(ope)) * im2[aux]);
                        }
                        else if (!SAB && SCD) {
                            store(dest, (AQ2 * load(op1)                                         + load(ope)) * im2[aux]);
                        }
                        else if (!SAB && !SCD) {
                            store(dest, (X2 * load(op1)                                         + load(ope)) * im2[aux]);
                        }
                    }


                    Jf += memK4J1f;
                    Je += memK4J1e;
                }
            }

            //K3 E + F
            {
                const cacheline64 * J1e = (cacheline64*)buffer.K4J1e;
                cacheline64 * J2e = (cacheline64*)buffer.K3J2e;

                const cacheline64 * J1f = (cacheline64*)buffer.K4J1f;
                cacheline64 * J2f = (cacheline64*)buffer.K3J2f;

                double Uu[maxJ][32] __attribute__((aligned(CACHE_LINE_SIZE)));

                for (int j=0; j<JB; ++j) {
                    Uu[j][0] = ABp.BCP.Nb[b][j];

                    for (int u=1; u<=maxU; ++u)
                        Uu[j][u] = Uu[j][u-1] * ABp.BCP.b1[b];
                }


                for (int nj1=0; nj1<nJ1; ++nj1) {
                    #ifdef __QR_GC__
                    for (int j=0; j<min(JB, b+1); ++j) {
                        LoopsContraction2(J2e+j*memK3J2e, J1e, Uu[j], K4VILcode.pK3E, K4VILcode.nbi[K3E]);
                        LoopsContraction2(J2f+j*memK3J2f, J1f, Uu[j], K4VILcode.pK3F, K4VILcode.nbi[K3F]);
                    }
                    #else
                    for (int j=0; j<JB; ++j) {
                        LoopsContraction2(J2e+j*memK3J2e, J1e, Uu[j], K4VILcode.pK3E, K4VILcode.nbi[K3E]);
                        LoopsContraction2(J2f+j*memK3J2f, J1f, Uu[j], K4VILcode.pK3F, K4VILcode.nbi[K3F]);
                    }
                    #endif

                    J2e += JB*memK3J2e;
                    J1e += memK4J1e;

                    J2f += JB*memK3J2f;
                    J1f += memK4J1f;
                }
            }

        }

        //AERR K2
        {
            cacheline64 * J2 = (cacheline64*)buffer.K3J2e;

            //__m128d ikcd_v = _mm_load1_pd(&ikcd);
            //u_atom ikcd_v(ikcd);

            //K2J2
            for (int nj2=0; nj2<nJ2; ++nj2) {
                for (int i=0; i<K4VILcode.nbi[AERR2]; ++i) {
                    cacheline64       * dest = J2 + K4VILcode.pAERR2[i].dest;
                    const cacheline64 * op1  = J2 + K4VILcode.pAERR2[i].op1;
                    const cacheline64 * op2  = J2 + K4VILcode.pAERR2[i].op2;

                    store(dest, load(op1) + load(op2) * ikcd);
                }

                J2 += memK3J2e;
            }
        }

        //CDR K2
        {
            cacheline64 * J2 = (cacheline64*)buffer.K3J2f;
            cacheline64 * Je = (cacheline64*)buffer.K3J2e;

            for (int nj2=0; nj2<nJ2; ++nj2) {
                for (int i=0; i<K4VILcode.nbi[CDR2]; ++i) {
                    cacheline64       * dest = J2 + K4VILcode.pCDR2[i].dest;
                    const cacheline64 * op1  = J2 + K4VILcode.pCDR2[i].op1;
                    const cacheline64 * op2  = J2 + K4VILcode.pCDR2[i].op2;
                    const cacheline64 * op3  = J2 + K4VILcode.pCDR2[i].op3;
                    const cacheline64 * ope  = Je + K4VILcode.pCDR2[i].ope;
                    UI16 aux  = K4VILcode.pCDR2[i].aux;

                    if (SAB && SCD) {
                        store(dest, (AQ2 * load(op1) - AQAB * load(op2) + AB2 * load(op3) + load(ope)) * im2[aux]);
                    }
                    else if (!SAB && SCD) {
                        store(dest, (AQ2 * load(op1)                                      + load(ope)) * im2[aux]);
                    }
                    else if (!SAB && !SCD) {
                        store(dest, (X2 * load(op1)                                      + load(ope)) * im2[aux]);
                    }
                }

                J2 += memK3J2f;
                Je += memK3J2e;
            }
        }

    }

}


template <bool SAB, bool SCD> void ERIroutine::ContractCDR_GC(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * __restrict__ uv_m_st8, p_ERIbuffer & buffer) const {

    const PrimitiveSet & PSab = ABp.Psets[ET.nKab];
    const PrimitiveSet & PScd = CDp.Psets[ET.nKcd];

    const cacheline64 & ABz = vars8.ABz;
    const cacheline64 & CDy = vars8.CDy;
    const cacheline64 & CDz = vars8.CDz;

    const cacheline64 & ACx = vars8.ACx;
    const cacheline64 & ACy = vars8.ACy;
    const cacheline64 & ACz = vars8.ACz;

    const int JA = ABp.Ja;
    const int JB = ABp.Jb;
    const int JC = CDp.Ja;
    const int JD = CDp.Jb;

    const int nJ1 = JA;
    const int nJ2 = JA*JB;
    const int nJ3 = JA*JB*JC;
    const int nJ4 = JA*JB*JC*JD;


    const int mK1f = nJ4*memK1J4f;
    const int mK2f = nJ3*memK2J3f;
    const int mK3f = nJ2*memK3J2f;
    const int mK4f = nJ1*memK4J1f;

    const int mK1e = nJ4*memK1J4e;
    const int mK2e = nJ3*memK2J3e;
    const int mK3e = nJ2*memK3J2e;
    const int mK4e = nJ1*memK4J1e;


    cacheline64 CD2;
    cacheline64 ACCD;
    cacheline64 ABCD;

    cacheline64 AC2;
    cacheline64 ACAB;
    cacheline64 AB2;

    cacheline64 X2;

    cacheline64 X2Y2;


    if (SAB && SCD) {
        CD2  =        CDy*CDy + CDz*CDz;
        ACCD =        ACy*CDy + ACz*CDz; ACCD += ACCD;
        ABCD =                  ABz*CDz; ABCD += ABCD;

        X2   = ACx*ACx;

        AC2  =      X2+ACy*ACy+ACz*ACz;
        ACAB =                 ABz*ACz; ACAB += ACAB;
        AB2  =                 ABz*ABz;
    }
    else if (!SAB && SCD) {
        CD2  =                  CDz*CDz;
        ACCD =                  ACz*CDz; ACCD += ACCD;
        X2   = ACx*ACx;
        X2Y2 = X2 + ACy*ACy;
        AC2  =      X2Y2+ACz*ACz;
    }
    else if (!SAB && !SCD) {
        AC2  =      ACx*ACx+ACy*ACy+ACz*ACz;
    }



    cacheline64 * F0 = (cacheline64*)buffer.F0;


    memset(buffer.K1J4e, 0, mK1e*sizeof(cacheline64));
    memset(buffer.K1J4f, 0, mK1f*sizeof(cacheline64));

    for (int d=0; d<PScd.nKb; ++d) {
        memset(buffer.K2J3e, 0, mK2e*sizeof(cacheline64));
        memset(buffer.K2J3f, 0, mK2f*sizeof(cacheline64));

        for (int c=0; c<PScd.nKa[d]; ++c) {

            double ikcd  = CDp.BCP.k[d][c];
            double rcd   = CDp.BCP.r[d][c];

            if      (InnerContractionRoutine!=NULL &&  SAB &&  SCD)
                InnerContractionRoutine      ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2);
            else if (InnerContractionRoutine!=NULL && !SAB &&  SCD)
                InnerContractionRoutine      ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2Y2);
            else if (InnerContractionRoutine!=NULL && !SAB &&  !SCD)
                InnerContractionRoutine      ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, AC2);

            else if (SAB && SCD)
                InnerContractCDR_GC<SAB,SCD> ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2);
            else if (!SAB && SCD)
                InnerContractCDR_GC<SAB,SCD> ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, X2Y2);
            else if (!SAB && !SCD)
                InnerContractCDR_GC<SAB,SCD> ( F0, vars8, PSab, ABp, ikcd, rcd, buffer, AB2, AC2);

            {
                const cacheline64 * J2e = (cacheline64*)buffer.K3J2e;
                cacheline64 * J3e = (cacheline64*)buffer.K2J3e;

                const cacheline64 * J2f = (cacheline64*)buffer.K3J2f;
                cacheline64 * J3f = (cacheline64*)buffer.K2J3f;

                double Tt[maxJ][32] __attribute__((aligned(CACHE_LINE_SIZE)));

                for (int j=0; j<JC; ++j) {
                    Tt[j][0] = CDp.BCP.Na[c][j];

                    for (int t=1; t<=maxT; ++t)
                        Tt[j][t] = Tt[j][t-1] * CDp.BCP.k[d][c];
                }

                //K2J3e+f contraction
                for (int nj2=0; nj2<nJ2; ++nj2) {
                  #ifdef __QR_GC__
                    for (int j=0; j<min(JC, c+1); ++j) {
                        LoopsContraction2(J3e+j*memK2J3e, J2e, Tt[j], K4VILcode.pK2E, K4VILcode.nbi[K2E]);
                        LoopsContraction2(J3f+j*memK2J3f, J2f, Tt[j], K4VILcode.pK2F, K4VILcode.nbi[K2F]);
                    }
                  #else
                    for (int j=0; j<JC; ++j) {
                        LoopsContraction2(J3e+j*memK2J3e, J2e, Tt[j], K4VILcode.pK2E, K4VILcode.nbi[K2E]);
                        LoopsContraction2(J3f+j*memK2J3f, J2f, Tt[j], K4VILcode.pK2F, K4VILcode.nbi[K2F]);
                    }
                  #endif

                    J3e += JC*memK2J3e;
                    J2e += memK3J2e;

                    J3f += JC*memK2J3f;
                    J2f += memK3J2f;
                }

            }


        }

        //AERR K1
        {
            cacheline64 * J3 = (cacheline64*)buffer.K2J3e;

            //K1J3
            for (int nj3=0; nj3<nJ3; ++nj3) {
                for (int i=0; i<K4VILcode.nbi[AERR1]; ++i) {
                    cacheline64       * dest = J3 + K4VILcode.pAERR1[i].dest;
                    const cacheline64 * op1  = J3 + K4VILcode.pAERR1[i].op1;
                    const cacheline64 * op2  = J3 + K4VILcode.pAERR1[i].op2;

                    store(dest, load(op1) + load(op2));
                }


                J3 += memK2J3e;
            }
        }

        cacheline64 CD2d2;
        cacheline64 ACCDd;
        cacheline64 ABCDd;

        if (SAB && SCD) {
            CD2d2 = CD2  * CDp.BCP.b2[d];
            ACCDd = ACCD * CDp.BCP.b1[d];
            ABCDd = ABCD * CDp.BCP.b1[d];
        }
        else if (!SAB && SCD) {
            CD2d2 = CD2  * CDp.BCP.b2[d];
            ACCDd = ACCD * CDp.BCP.b1[d];
        }

        //CDR K1
        {
            cacheline64 * J3 = (cacheline64*)buffer.K2J3f;
            cacheline64 * Je = (cacheline64*)buffer.K2J3e;

            for (int nj3=0; nj3<nJ3; ++nj3) {
                for (int i=0; i<K4VILcode.nbi[CDR1]; ++i) {
                    cacheline64       * dest = J3 + K4VILcode.pCDR1[i].dest;
                    const cacheline64 * op1  = J3 + K4VILcode.pCDR1[i].op1;
                    const cacheline64 * op2  = J3 + K4VILcode.pCDR1[i].op2;
                    const cacheline64 * op3  = J3 + K4VILcode.pCDR1[i].op3;
                    const cacheline64 * op4  = J3 + K4VILcode.pCDR1[i].op4;
                    const cacheline64 * op5  = J3 + K4VILcode.pCDR1[i].op5;
                    const cacheline64 * op6  = J3 + K4VILcode.pCDR1[i].op6;
                    const cacheline64 * ope  = Je + K4VILcode.pCDR1[i].ope;
                    UI16 aux  = K4VILcode.pCDR1[i].aux;

                    if (SAB && SCD) {
                        store(dest, (AC2 * load(op1) - ACAB * load(op2) + ACCDd * load(op3) +
                         AB2 * load(op4) - ABCDd * load(op5) + CD2d2 * load(op6) +
                        load(ope)) * im2[aux]);
                    }
                    else if (!SAB && SCD) {
                        store(dest, (AC2 * load(op1)                    + ACCDd * load(op3) +
                                                               CD2d2 * load(op6) +
                        load(ope)) * im2[aux]);
                    }
                    else if (!SAB && !SCD) {
                        store(dest, (AC2 * load(op1) +
                        load(ope)) * im2[aux]);
                    }
                }

                J3 += memK2J3f;
                Je += memK2J3e;
            }

        }

        //K1J4 E + F
        {
            const cacheline64 * J3e = (cacheline64*)buffer.K2J3e;
            cacheline64 * J4e = (cacheline64*)buffer.K1J4e;

            const cacheline64 * J3f = (cacheline64*)buffer.K2J3f;
            cacheline64 * J4f = (cacheline64*)buffer.K1J4f;

            double Ss[maxJ][32] __attribute__((aligned(CACHE_LINE_SIZE)));

            for (int j=0; j<JD; ++j) {
                Ss[j][0] = CDp.BCP.Nb[d][j];

                for (int s=1; s<=maxS; ++s)
                    Ss[j][s] = Ss[j][s-1] * CDp.BCP.b1[d];
            }


            for (int nj3=0; nj3<nJ3; ++nj3) {

              #ifdef __QR_GC__
                for (int j=0; j<min(JD, d+1); ++j) {
                    LoopsContraction2(J4e+j*memK1J4e, J3e, Ss[j], K4VILcode.pK1E, K4VILcode.nbi[K1E]);
                    LoopsContraction2(J4f+j*memK1J4f, J3f, Ss[j], K4VILcode.pK1F, K4VILcode.nbi[K1F]);
                }
              #else
                for (int j=0; j<JD; ++j) {
                    LoopsContraction2(J4e+j*memK1J4e, J3e, Ss[j], K4VILcode.pK1E, K4VILcode.nbi[K1E]);
                    LoopsContraction2(J4f+j*memK1J4f, J3f, Ss[j], K4VILcode.pK1F, K4VILcode.nbi[K1F]);
                }
              #endif

                J4e += JD*memK1J4e;
                J3e += memK2J3e;

                J4f += JD*memK1J4f;
                J3f += memK2J3f;
            }
        }
    }


    //AERR K0
    {
        cacheline64 * J4 = (cacheline64*)buffer.K1J4e;

        //K0J4
        for (int nj4=0; nj4<nJ4; ++nj4) {
            for (int i=0; i<K4VILcode.nbi[AERR0]; ++i) {
                cacheline64       * dest = J4 + K4VILcode.pAERR0[i].dest;
                const cacheline64 * op1  = J4 + K4VILcode.pAERR0[i].op1;
                const cacheline64 * op2  = J4 + K4VILcode.pAERR0[i].op2;

                store(dest, load(op1) + load(op2));
            }

            J4 += memK1J4e;
        }
    }

    //CDR K0
    {
        cacheline64 * J4 = (cacheline64*)buffer.K1J4f;
        cacheline64 * Je = (cacheline64*)buffer.K1J4e;

        for (int nj4=0; nj4<nJ4; ++nj4) {
            for (int i=0; i<K4VILcode.nbi[CDR0]; ++i) {
                cacheline64       * dest = J4 + K4VILcode.pCDR0[i].dest;
                const cacheline64 * op1  = J4 + K4VILcode.pCDR0[i].op1;
                const cacheline64 * op2  = J4 + K4VILcode.pCDR0[i].op2;
                const cacheline64 * op3  = J4 + K4VILcode.pCDR0[i].op3;
                const cacheline64 * op4  = J4 + K4VILcode.pCDR0[i].op4;
                const cacheline64 * op5  = J4 + K4VILcode.pCDR0[i].op5;
                const cacheline64 * op6  = J4 + K4VILcode.pCDR0[i].op6;
                const cacheline64 * ope  = Je + K4VILcode.pCDR0[i].ope;
                UI16 aux  = K4VILcode.pCDR0[i].aux;

                if (SAB && SCD) {
                    store(dest, (AC2 * load(op1) - ACAB * load(op2) + ACCD * load(op3) +
                     AB2 * load(op4) - ABCD * load(op5) + CD2 * load(op6) +
                    load(ope)) * im2[aux]);
                }
                else if (!SAB && SCD) {
                    store(dest, (AC2 * load(op1)                    + ACCD * load(op3) +
                                                          CD2 * load(op6) +
                    load(ope)) * im2[aux]);
                }
                else if (!SAB && !SCD) {
                    store(dest, (AC2 * load(op1) +
                    load(ope)) * im2[aux]);
                }
            }


            J4 += memK1J4f;
            Je += memK1J4e;
        }
    }

    //copy to the kernel buffer
    {
        const cacheline64 * J4 = (cacheline64*)buffer.K1J4f;
        cacheline64 * I = uv_m_st8;

        for (int nj4=0; nj4<nJ4; ++nj4) {
            register Op_MIRROR * s2 = eseq;

            for (int i=0; s2<nseq[KERNELS]; ++i,++s2) {
                cacheline64       * dest = I + s2->dest;
                const cacheline64 * op1 = J4 + s2->op1;

                __builtin_prefetch(J4 + s2[PF].op1);
                __builtin_prefetch(I  + s2[PF].dest, 1, 1);

                if (s2->aux == 0) store(dest, load(op1));
                else              memset(dest, 0, sizeof(cacheline64)); // zero the kernel
            }

            J4 += memK1J4f;
            I += nKernels;
        }
    }


}


void p_Qalgorithm::K4bench(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * uv_m_st8, p_ERIbuffer & buffer, bool OnlyJ) const {


    K4benchmark * bench = p_Q.Benchmarks[eritype];

    bench->AddCall(ET,ABp,CDp);

    //has to initialize buffer first
    //IC->CalcGammas(vars8, ET, ABp, CDp, buffer, OnlyJ);

    if      (eritype.geometry==ABCD) IC->CalcGammas<true,  true>  (vars8, ET, ABp, CDp, buffer, OnlyJ);
    else if (eritype.geometry==AACD) IC->CalcGammas<false, true>  (vars8, ET, ABp, CDp, buffer, OnlyJ);
    else if (eritype.geometry==AACC) IC->CalcGammas<false, false> (vars8, ET, ABp, CDp, buffer, OnlyJ);
    else                             IC->CalcGammas<true,  true>  (vars8, ET, ABp, CDp, buffer, OnlyJ);


    Chronometer chronoC;
    chronoC.Start();

    if (!IsDynamic) (*SC)(buffer,vars8, ET, ABp, CDp, uv_m_st8);
    else if (eritype.geometry==ABCD) IC->ContractCDR_GC<true,  true>  (vars8, ET, ABp, CDp, uv_m_st8, buffer);
    else if (eritype.geometry==AACD) IC->ContractCDR_GC<false, true>  (vars8, ET, ABp, CDp, uv_m_st8, buffer);
    else if (eritype.geometry==AACC) IC->ContractCDR_GC<false, false> (vars8, ET, ABp, CDp, uv_m_st8, buffer);
    else                             IC->ContractCDR_GC<true,  true>  (vars8, ET, ABp, CDp, uv_m_st8, buffer);

    chronoC.Stop();

    double At = chronoC.GetTotalTime();

    #pragma omp atomic
    bench->deltaK4 += At;
}

void p_Qalgorithm::K4(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * uv_m_st8, p_ERIbuffer & buffer, bool OnlyJ) const {

    if      (eritype.geometry==ABCD) IC->CalcGammas<true,  true>  (vars8, ET, ABp, CDp, buffer, OnlyJ);
    else if (eritype.geometry==AACD) IC->CalcGammas<false, true>  (vars8, ET, ABp, CDp, buffer, OnlyJ);
    else if (eritype.geometry==AACC) IC->CalcGammas<false, false> (vars8, ET, ABp, CDp, buffer, OnlyJ);
    else                             IC->CalcGammas<true,  true>  (vars8, ET, ABp, CDp, buffer, OnlyJ);

    if      (eritype.geometry==ABCD) IC->ContractCDR_GC<true,  true>  (vars8, ET, ABp, CDp, uv_m_st8, buffer);
    else if (eritype.geometry==AACD) IC->ContractCDR_GC<false, true>  (vars8, ET, ABp, CDp, uv_m_st8, buffer);
    else if (eritype.geometry==AACC) IC->ContractCDR_GC<false, false> (vars8, ET, ABp, CDp, uv_m_st8, buffer);
    else                             IC->ContractCDR_GC<true,  true>  (vars8, ET, ABp, CDp, uv_m_st8, buffer);

    //IC->CalcGammas     (vars8, ET, ABp, CDp, buffer, OnlyJ);
    //IC->ContractCDR_GC (vars8, ET, ABp, CDp, uv_m_st8, buffer);

    //if (IsDynamic) IC->ContractCDR_GC(vars8,ET,ABp,CDp,uv_m_st8,buffer);
    //else          (*SC)(buffer,vars8, ET, ABp, CDp, uv_m_st8);
}


