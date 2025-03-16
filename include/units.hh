// This file is part of CCQE, a quasi-elastic reaction modeling module.
// Portions of this file are derived from the NuWro project.
//
// Copyright (C) 2025 Anna Ershova
// Copyright (C) NuWro Developers
//
// CCQE is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// CCQE is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CCQE. If not, see <https://www.gnu.org/licenses/>.

#ifndef _units_
#define _units_

static const double c = 1;

static const double MeV = 1;
static const double GeV = 1000 * MeV;

const double G = 1.16639 * 1E-5 / (GeV * GeV);

const double sin2thetaW = 0.2312215;  // weak mixing angle

static const double Mv = 0.84 * GeV;
static const double p_mass = 938.272;
static const double n_mass = 939.565;
static const double M = (p_mass + n_mass)/2;
static const double ga = 1.267;
static const double gas = 0;
static const double Ma = 1.03 * GeV;
static const double Mpi = 134.9768 * MeV;
#endif
