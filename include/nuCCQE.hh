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

#ifndef INCLUDE_NUCCQE_HH_
#define INCLUDE_NUCCQE_HH_ 1

#include <fstream>
#include "units.hh"
#include "FormFactor.hh"

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
} // namespace G4INCL
#endif // INCLUDE_NUCCQE_HH_