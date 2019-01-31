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


#include <set>
#include "libquimera.hpp"
#include "../2eints/quimera.hpp"
#include "../math/angular.hpp"
#include "../math/gamma.hpp"


LibQuimera::ERIbuffer::ERIbuffer() {
    pERIbuffer = new p_ERIbuffer;
}

LibQuimera::ERIbuffer::~ERIbuffer() {
    delete pERIbuffer;
}

size_t LibQuimera::ERIbuffer::SetBuffers(cacheline * mem, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp,  const Qalgorithm & routine) {
    return pERIbuffer->SetBuffers( mem, ABp, CDp,  *(routine.pQalgorithm));
}


LibQuimera::Qalgorithm::Qalgorithm() {
    pQalgorithm = new p_Qalgorithm();
}

LibQuimera::Qalgorithm::~Qalgorithm() {
    delete pQalgorithm;
}

UI32 LibQuimera::Qalgorithm::GetNKernels() const {
    return pQalgorithm->nKernels;
}

size_t LibQuimera::Qalgorithm::MemSize(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp) const {
    return pQalgorithm->MemSize(ABp, CDp);
}

//execute
void LibQuimera::Qalgorithm::K4bench(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * uv_m_st8, ERIbuffer & buffer, bool UseCase) const {
    pQalgorithm->K4bench(vars8, ET, ABp, CDp, uv_m_st8, *(buffer.pERIbuffer), UseCase);
}

void LibQuimera::Qalgorithm::K4(const ERIgeometries64 & vars8, const ERITile64 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline64 * uv_m_st8, ERIbuffer & buffer, bool UseCase) const {
    pQalgorithm->K4(vars8, ET, ABp, CDp, uv_m_st8, *(buffer.pERIbuffer), UseCase);
}

void LibQuimera::Qalgorithm::K4(const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, double * uv_m_st, ERIbuffer & buffer, bool UseCase) const {
    pQalgorithm->K4(ABp, CDp, uv_m_st, *(buffer.pERIbuffer), UseCase);
}


void LibQuimera::Qalgorithm::MIRRORbench(const cacheline64 *   mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, ERIbuffer & buffer) const {
    pQalgorithm->MIRRORbench(mem, vars8, ERI8, *(buffer.pERIbuffer));
}

void LibQuimera::Qalgorithm::MIRROR(const cacheline64 *   mem, const ERIgeometries64 & vars8, cacheline64 * ERI8, ERIbuffer & buffer) const {
    pQalgorithm->MIRROR(mem, vars8, ERI8, *(buffer.pERIbuffer));
}

void LibQuimera::Qalgorithm::MIRROR(const double *  uv_m_st, double * W, ERIbuffer & buffer) const {
    pQalgorithm->MIRROR(uv_m_st, W, *(buffer.pERIbuffer));
}



void LibQuimera::Qalgorithm::K4(const ERIgeometries32 & vars16, const ERITile32 & ET, const ShellPairPrototype & ABp, const ShellPairPrototype & CDp, cacheline32 * uv_m_st16, ERIbuffer & buffer, bool UseCase) const {
    pQalgorithm->K4(vars16, ET, ABp, CDp, uv_m_st16, *(buffer.pERIbuffer), UseCase);
}

void LibQuimera::Qalgorithm::MIRROR(const cacheline32 *   mem, const ERIgeometries32 & vars16, cacheline32 * ERI16, ERIbuffer & buffer) const {
    pQalgorithm->MIRROR(mem, vars16, ERI16, *(buffer.pERIbuffer));
}



MessagePlotter LibQuimera::Quimera::QMessenger;
MessagePlotter LibQuimera::Quimera::QBenchmarker;

LibQuimera::Quimera::Quimera() {
    QMessenger.SetModule  ("Quimera module");
    QBenchmarker.SetModule ("Quimera module benchmark");

    pQuimera = &p_Q;
}

LibQuimera::Quimera::~Quimera() {
    pQuimera = NULL;
    //delete pQuimera;
}

const Qalgorithm * LibQuimera::Quimera::SelectAlgorithm (GEOM geometry, UI8  La, UI8  Lb, UI8  Lc, UI8  Ld) const {
    return pQuimera->SelectAlgorithm(geometry, La, Lb, Lc, Ld);
}

void LibQuimera::Quimera::ListNeeded     (GEOM geometry, UI8 La, UI8  Lb, UI8  Lc, UI8  Ld) {
    pQuimera->ListNeeded(geometry, La,  Lb, Lc, Ld);
}

void LibQuimera::Quimera::GenerateNeeded (bool overwrite) {
    pQuimera->GenerateNeeded (overwrite);
}

void LibQuimera::Quimera::Statistics() const {
    pQuimera->Statistics();
}



namespace LibQuimera __attribute__ ((visibility ("default"))) {
    Quimera Q;

    double  case_w2 = 0.25;
    bool    UseCASE = false;
    bool    UseCDR  = true;
    bool    UseGC   = true;
}

bool LibQuimera::Initialized = false;
std::string LibQuimera::IC_DIR = "";

// initializes some computing constants, including the gamma function and spherical harmonics
int LibQuimera::InitLibQuimera() {

    if (!Initialized) {
        LibAngular::InitCartList();
        LibAngular::InitSHList();
        LibIGamma::InitIncompleteGammas();
        Quimera::QBenchmarker.SetFile("quimera.benchmark.out");
        Initialized = true;
    }

    return 0;
}



