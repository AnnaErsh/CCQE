#ifndef G4INCLFormFactor_HH_
#define G4INCLFormFactor_HH_ 1

#include <functional>

#include "jednostki.h"

struct EMFFparam
{
  double GEp, GEn;
  double GMp, GMn;
  EMFFparam() : GEp(0), GEn(0), GMp(0), GMn(0) {}
};

enum particle {proton, neutron};

double inline pow2(double x) { return x * x; }

static const double Mv = 0.84 * GeV;
static const double p_mass = 938.272;
static const double n_mass = 939.565;
static const double M = (p_mass + n_mass)/2;
static const double ga = 1.267;
static const double gas = 0;
static const double Ma = 1.03 * GeV;
static const double Mpi = 134.9768 * MeV;

class FormFactor
{
  public:
    FormFactor(particle, bool);
    ~FormFactor();
    double Fp       (double Q2);
    void   bbba05_FF(double Q2);

    std::function<double(FormFactor*,double)>       ptrGa;
    double Ga  (double x)        {ptrGa(this,x);};
    std::function<double(FormFactor*,double,particle,int)> ptrFv12;
    double Fv12(double x, particle p, int y) {ptrFv12(this,x,p,y);};

  private:
    double Fs12  (double Q2, double num);
    double Fv12pn(double Q2, particle p, int num);
    double Fv12CC(double Q2, particle p, int num);
    double Fv12NC(double Q2, particle p, int num);
    double GaCC  (double Q2);
    double GaNC  (double Q2);
    EMFFparam emff;
    particle part;
    bool ifCC;
};

#endif