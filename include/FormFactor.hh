#ifndef G4INCLFormFactor_HH_
#define G4INCLFormFactor_HH_ 1

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
    double Fp       (double Q2);
    void   bbba05_FF(double Q2);

    /**
     * Axial form factor, ptrGa will be GaCC or GaNC depending on the channel
    */
    std::function<double(FormFactor*,double)> ptrGa;
    double Ga  (double x) {ptrGa(this,x);};
    /**
     * Vector form factor, ptrGa will be GaCC or GaNC depending on the channel
    */
    std::function<double(FormFactor*,double,particle,int)> ptrFv12;
    double Fv12(double x, particle p, int y) {ptrFv12(this,x,p,y);};

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

#endif