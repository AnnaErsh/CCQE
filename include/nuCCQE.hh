#ifndef G4INCLnuQELChannel_HH_
#define G4INCLnuQELChannel_HH_ 1

#include "jednostki.h"
#include "FormFactor.hh"
#include <fstream>

namespace G4INCL {
	class nuQELChannel {

	public:
		nuQELChannel(bool, bool, bool);
		virtual ~nuQELChannel();
		void fillFinalState();
		double QELAmplitude(double q2, double Enu, bool ifnu);
	private:
		double Fp(double Q2);
		double z(double q2);
		double A(double q2);
		double B(double q2);
		double C(double q2);
		void ABCtest(double q2, double Enu);

		std::ofstream out;
		std::ofstream out1;

		// to change later:
		// double p_mass = 938.272;
		// double n_mass = 939.565;

		double m = 105.658; // muon mass
		const double cosThetac = 0.974213;

		FormFactor *ff;
		particle part;

		double f = 0;
		double ml = 0; // mass of the outgoing lepton

		// double (FormFactor::*Ga)(double);
    // double (FormFactor::*FV12)(double, int);

		// int type;
	};
}
#endif