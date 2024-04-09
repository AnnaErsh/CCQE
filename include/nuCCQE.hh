#ifndef G4INCLnuQELChannel_HH_
#define G4INCLnuQELChannel_HH_ 1

#include "units.hh"
#include "FormFactor.hh"
#include <fstream>

namespace G4INCL {
	class nuQELChannel {

	public:
		/**
		 * Constructor of the main Llewellyn-Smith calculation
		 * @param ifCC if the channel is charged current or neutral current
		 * @param ifproton if the outgoing hadron is a proton or neutron
		 * @param ifmuon if the outgoing lepton is a muon
		*/
		nuQELChannel(bool, bool, bool);
		virtual ~nuQELChannel();
		/**
		 * Should provide the final state particles in INCL
		 * @note not yet implemented
		*/
		void fillFinalState();
		/**
		 * Main function that calculates the amplitude of the QEL channel 
		 * @param q2 square of four-momentum transfer
		 * @param Enu neutrino energy
		 * @param ifnu neutrino or antineutrino reaction
		 * @return d|\sigma| / d|q^2|
		*/
		double QELAmplitude(double q2, double Enu, bool ifnu);
	private:
		/**
		 * Calculates z, z = q^2/M^2
		*/
		double z(double q2);
		/**
		 * A function of the Llewellyn-Smith formula
		*/
		double A(double q2);
		/**
		 * B function of the Llewellyn-Smith formula
		*/
		double B(double q2);
		/**
		 * C function of the Llewellyn-Smith formula
		*/
		double C(double q2);
		/**
		 * Debug function to validate amplitude calculation
		 * @param q2 square of four-momentum transfer
		 * @param Enu neutrino energy
		*/
		void ABCtest(double q2, double Enu);

		double m = 105.658; // muon mass
		const double cosThetac = 0.974213;

		FormFactor *ff;
		particle part;

		double f = 0;
		double ml = 0; // mass of the outgoing lepton
	};
}
#endif