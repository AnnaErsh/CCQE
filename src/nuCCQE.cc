#include "nuCCQE.hh"
#include <iostream>

# define M_PI 3.14159265358979323846

namespace G4INCL {

	nuQELChannel::nuQELChannel(bool ifCC, bool ifproton = true, bool ifmuon = true)
	{
    out.open("ABC300anu_ncn.dat", std::fstream::out);
    out1.open("INCLABC300NCn_antinu.dat", std::fstream::out);
    if(ifproton) part = particle::proton;
    else part = particle::neutron;
    
    ff = new FormFactor(part, ifCC);

    if(ifmuon) ml = 105.65 * MeV;
    else ml = 0.511 * MeV;

    if(ifCC)
    {
      f = cosThetac;
      m = 105.658 * MeV;
    }
    else 
    {
      f = 1;
      m = 0;
      ml = 0;
    }
	}

	nuQELChannel::~nuQELChannel()
  {
  }

  void nuQELChannel::fillFinalState() {};

  double nuQELChannel::z(double q2)
  {
  	return q2 / pow2(M);
  }

  double nuQELChannel::A(double q2)
  {
    double Q2 = -q2;

    return 1./4. * (pow2(m/M) - z(q2)) * 
      (
        (4 - z(q2)) * pow2((ff->Ga)(Q2)) -
        (4 + z(q2)) * pow2((ff->Fv12)(Q2, part, 1)) -
        z(q2) * pow2((ff->Fv12)(Q2, part, 2)) * (1 + 1./4. * z(q2)) - 
        4 * (ff->Fv12)(Q2, part, 1) * (ff->Fv12)(Q2, part, 2) * z(q2) - 
        (pow2(m/M)) * 
        (
          pow2((ff->Fv12)(Q2, part, 1) + (ff->Fv12)(Q2, part, 2)) +
          pow2((ff->Ga)(Q2) + 2 * ff->Fp(Q2)) +
          (z(q2) - 4) * pow2(ff->Fp(Q2))
        )
      );
  }

	double nuQELChannel::B(double q2)
	{
    double Q2 = -q2;
    std::cout<<"B: "<<q2<<" "<<z(q2)<<" "<<(ff->Ga)(Q2)<<" "<<(ff->Fv12)(Q2,part, 1)<<" "<<(ff->Fv12)(Q2,part, 2)<<" "<<-z(q2) * (ff->Ga)(Q2) * ((ff->Fv12)(Q2,part, 1) + (ff->Fv12)(Q2,part, 2))<<std::endl;
    return -z(q2) * (ff->Ga)(Q2) * ((ff->Fv12)(Q2, part, 1) + (ff->Fv12)(Q2, part, 2));
	}

	double nuQELChannel::C(double q2)
	{
    double Q2 = -q2;
    // std::cout<<q2<<" "<<1./4. * (pow2((ff->Ga)(Q2)) + pow2((ff->Fv12)(Q2, 1)) - z(q2) * pow2((ff->Fv12)(Q2, 2)/2))<<" "<<pow2((ff->Ga)(Q2))<<" "<<pow2((ff->Fv12)(Q2, 1))<<" "<<z(q2)<<" "<<pow2((ff->Fv12)(Q2, 2))<<std::endl;
    return 1./4. * (pow2((ff->Ga)(Q2)) + pow2((ff->Fv12)(Q2, part, 1)) - z(q2) * pow2((ff->Fv12)(Q2, part, 2)/2));
	}

  void nuQELChannel::ABCtest(double q2, double Enu)
  {
    ff->bbba05_FF(-q2);

    // ff->bbba05_FF(244324);
    // A(-364984);
    double s_u = 4*Enu*M + q2 - ml;
    double suM2 = s_u / M / M;
    double k = -1;

    std::cout<<q2<<" amp = "<<pow2(M * G * f / Enu) / 8 / M_PI * (A(q2) + suM2 * (B(q2) + suM2 * C(q2)))<<std::endl;
  }

  double nuQELChannel::QELAmplitude(double q2, double Enu, bool ifnu)
  {
    // std::cout<<"QELAmplitude q2 = "<<q2<<std::endl;
    ff->bbba05_FF(-q2);

    // ff->bbba05_FF(244324);
    B(-52968.5);
  	double s_u = 4*Enu*M + q2 - ml;
    double k = 0;
    if(ifnu) k = 1;
    else k = -1;
    double suM2 = s_u / M / M;

    out<<q2<<" "<<A(q2)<<" "<<B(q2)<<" "<<C(q2)<<std::endl;
    // ABCtest(-130143, 1100);
    // std::cout<<"!"<<-364984<<" "<<A(-364984)<<" "<<B(-364984)<<" "<<C(-364984)<<std::endl;

    // out1<<(A(q2) + suM2 * (B(q2) + suM2 * C(q2)))<<" "<<-q2<<std::endl;
    // std::cout<<(A(q2) + suM2 * (B(q2) + suM2 * C(q2)))<<" "<<-q2<<std::endl;
    out1<<pow2(M * G * f / Enu) / 8 / M_PI * 
           (A(q2) + suM2 * (k * B(q2) + suM2 * C(q2)))<<" "<<-q2<<std::endl;
  	return pow2(M * G * f / Enu) / 8 / M_PI * 
  	  		 (A(q2) + suM2 * (k * B(q2) + suM2 * C(q2)));
  }

}