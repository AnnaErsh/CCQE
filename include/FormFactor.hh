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

#ifndef INCLUDE_FORMFACTOR_HH_
#define INCLUDE_FORMFACTOR_HH_ 1

#include <functional>

#include "units.hh"

struct EMFFparam
{
  double GEp, GEn;
  double GMp, GMn;
  EMFFparam() : GEp(0), GEn(0), GMp(0), GMn(0) {}
};

enum particle {proton, neutron};

double inline pow2(double x) { return x * x; }

class FormFactor
{
  public:
    FormFactor(particle, bool);
    ~FormFactor();
    /**
		* pseudoscalar axial form factor
		*/
    double Fp       (double Q2)
    void   bbba05_FF(double Q2)

    /**
     * Axial form factor, ptrGa will be GaCC or GaNC depending on the channel
    */
    std::function<double(FormFactor*,double)> ptrGa;
    double Ga  (double x) {ptrGa(this,x);}
    /**
     * Vector form factor, ptrGa will be GaCC or GaNC depending on the channel
    */
    std::function<double(FormFactor*,double,particle,int)> ptrFv12;
    double Fv12(double x, particle p, int y) {ptrFv12(this,x,p,y);}

  private:
    double Fs12  (double Q2, double num);
    /**
     * Proton and neutron electromagnetic form factors
    */
    double Fv12pn(double Q2, particle p, int num);
    /**
     * Vector form factor for the CC channel
    */
    double Fv12CC(double Q2, particle p, int num);
    /**
     * Vector form factor for the NC channel
    */
    double Fv12NC(double Q2, particle p, int num);
    /**
     * Axial form factor for the CC channel
    */
    double GaCC  (double Q2);
    /**
     * Axial form factor for the NC channel
    */
    double GaNC  (double Q2);
    EMFFparam emff;
    particle part;
    bool ifCC;
};

#endif // INCLUDE_FORMFACTOR_HH_