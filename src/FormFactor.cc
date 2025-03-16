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

#include "FormFactor.hh"
#include <iostream>

FormFactor::FormFactor(particle part_, bool ifCC_):part(part_),
                                                   ifCC(ifCC_)
{
  if(ifCC)
  {
    ptrGa = &FormFactor::GaCC;
    ptrFv12 = &FormFactor::Fv12CC;
  }
  else
  {
    ptrGa = &FormFactor::GaNC;
    ptrFv12 = &FormFactor::Fv12NC;
  }
}

// nuwro heritage
void FormFactor::bbba05_FF(const double Q2)
{
  double q2 = -Q2;
  double tau = -q2 / 4.0 / pow2(M);
  emff.GEp = (1.0 - tau * 0.0578) /
           (1.0 + tau * (11.1 + tau * (13.6 + tau * 33.0)));
  emff.GEn =
      tau * (1.25 + tau * 1.30) /
      (1.0 +
       tau * (-9.86 +
              tau * (305.0 + tau * (-758.0 + tau * 802.0))));
  emff.GMp = 2.792847351 * (1.0 + tau * 0.15) /
           (1.0 + tau * (11.1 + tau * (19.6 + tau * 7.54)));
  emff.GMn = -1.91304273 * (1.0 + tau * 1.81) /
           (1.0 + tau * (14.1 + tau * (20.7 + tau * 68.7)));
}

double FormFactor::Fv12pn(double Q2, particle p, int num)
{
  bbba05_FF(Q2);
  double q2 = -Q2;
  double z = Q2/pow2(M);

  double Ge = 0, Gm = 0;

  if(ifCC)
  {
    Ge = emff.GEp - emff.GEn;
    Gm = emff.GMp - emff.GMn;
  }
  else
  {
    if(p == particle::proton)
    {
      Ge = 0.5 * (emff.GEp - emff.GEn) - 2 * sin2thetaW * emff.GEp;
      Gm = 0.5 * (emff.GMp - emff.GMn) - 2 * sin2thetaW * emff.GMp;
    }
    else if(p == particle::neutron)
    {
      Ge = 0.5 * (emff.GEn - emff.GEp) - 2 * sin2thetaW * emff.GEn;
      Gm = 0.5 * (emff.GMn - emff.GMp) - 2 * sin2thetaW * emff.GMn;
    }
  }

  if(num == 1)
  {
    return (Ge + z / 4 * Gm) / (1. + z / 4);
  }
  else if (num == 2)
  {
    return (Gm - Ge) / (1. + z / 4);
  }
}

double FormFactor::Fv12CC(double Q2, particle p, int num)
{
  return Fv12pn(Q2, particle::proton, num);
}

double FormFactor::Fv12NC(double Q2, particle p, int num)
{
  return Fv12pn(Q2, p, num);
}

double FormFactor::GaCC(double Q2)
{
  return ga / pow2( 1 + Q2 / pow2(Ma) );
}

double FormFactor::GaNC(double Q2)
{
  double g = 0;
  if(part == particle::proton) g = ga;
  else if(part == particle::neutron) g = -ga;
  return 0.5 / pow2(1 + Q2 / pow2(Ma)) * (g - gas);
}

double FormFactor::Fs12(double Q2, double num)
{
  double F10 = 0;
  if(num == 1) F10 = 0.53;
  else if(num == 2) F10 = -0.4;

  return F10 * 1 / (1 + Q2 / (4 * M)) * 1 / pow2(1 + Q2 / (4 * Mv));
}

double FormFactor::Fp(double Q2)
{
  double q2 = -Q2;
	return 2 * M * M * GaCC(Q2) / (Mpi * Mpi + Q2);
}