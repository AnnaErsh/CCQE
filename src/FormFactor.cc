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

  // std::cout<<"bbba05_FF: "<<Q2<<" "<<tau<<" "<<emff.GMp<<" "<<emff.GMn<<std::endl;
}

double FormFactor::Fv12pn(double Q2, particle p, int num)
{
  bbba05_FF(Q2);
  double q2 = -Q2;
  double z = Q2/pow2(M);

  double Ge = 0, Gm = 0;

  if(ifCC)
  {
    std::cout<<"hey!CC"<<std::endl;
    Ge = emff.GEp - emff.GEn;
    Gm = emff.GMp - emff.GMn;
  }
  else
  {
    if(p == particle::proton)
    {
      std::cout<<"hey!!!"<<std::endl;
      Ge = 0.5 * (emff.GEp - emff.GEn) - 2 * sin2thetaW * emff.GEp;
      Gm = 0.5 * (emff.GMp - emff.GMn) - 2 * sin2thetaW * emff.GMp;
    }
    else if(p == particle::neutron)
    {
      std::cout<<"hey!"<<std::endl;
      Ge = 0.5 * (emff.GEn - emff.GEp) - 2 * sin2thetaW * emff.GEn;
      Gm = 0.5 * (emff.GMn - emff.GMp) - 2 * sin2thetaW * emff.GMn;
    }
  }

  // if(p == particle::proton)
  // {
  //   // Ge = emff.GEp;
  //   // Gm = emff.GMp;
  //   Ge = emff.GEp - emff.GEn;
  //   Gm = emff.GMp - emff.GMn;
  // }
  // else if(p == particle::neutron)
  // {
  //   Ge = emff.GEn;
  //   Gm = emff.GMn;
  //   // Ge = emff.GEp - emff.GEn;
  //   // Gm = emff.GMp - emff.GMn;
  // }

  if(Q2 == 52968.5)
  {
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"Q2 = "<<Q2<<"z/4 = "<<z/4<<std::endl;
  std::cout<<"Ge = "<<Ge<<" GEp = "<<emff.GEp<<" GEn = "<<emff.GEn<<std::endl;
  std::cout<<"Gm = "<<Gm<<" GMp = "<<emff.GMp<<" GMn = "<<emff.GMn<<std::endl;
  std::cout<<"GEp - GEn = "<<emff.GEp - emff.GEn<<" GMp - GMn = "<<emff.GMp - emff.GMn<<std::endl;  
  std::cout<<" num = "<<num<<" f1 = "<<(Ge + z/4 * Gm)/(1. + z/4)<<" f2 = "<<(Gm - Ge)/(1. + z/4)<<std::endl; 
  std::cout<<"1/(1 + z/4) = "<<1/(1 + z/4)<<std::endl;

  }


  if(num == 1)
  {
    // std::cout<<"num = "<<num<<" "<<(Ge + z/4 * Gm)/(1. + z/4)<<std::endl;
    return (Ge + z/4 * Gm)/(1. + z/4);
  }
  else if (num == 2)
  { 
      // std::cout<<q2<<" Gm = "<<Gm<<" Ge = "<<Ge<<" z = "<<z/4<<" F22 = "<<(Gm - Ge)/(1. + z/4) * (Gm - Ge)/(1. + z/4)<<std::endl;

    // std::cout<<"num = "<<num<<" "<<(Gm - Ge)/(1. + z/4)<<std::endl;
    return (Gm - Ge)/(1. + z/4);
  }
}

double FormFactor::Fv12CC(double Q2, particle p, int num)
{
  // std::cout<<"Fv12CC = "<<Fv12pn(Q2, particle::proton, 1)<<" "<<Fv12pn(Q2, particle::neutron, 1)<<std::endl;
  // if(num == 1) return Fv12pn(Q2, particle::proton, 1) - Fv12pn(Q2, particle::neutron, 2);
  // else if (num == 2) return Fv12pn(Q2, particle::proton, 1) - Fv12pn(Q2, particle::neutron, 2);
  return Fv12pn(Q2, particle::proton, num);
}

double FormFactor::Fv12NC(double Q2, particle p, int num)
{
  // double k = 0;
  // if(part == particle::proton) k = 1;
  // else if(part == particle::neutron) k = -1;

  // return 0.5 * (Fv12pn(Q2, particle::proton, num) - Fv12pn(Q2, particle::neutron, num)) -
  //          2 * sin2thetaW * Fv12pn(Q2, part, num) - 
  //        0.5 * Fs12(Q2, num);
  return Fv12pn(Q2, p, num);
}

double FormFactor::GaCC(double Q2) { return ga / pow2( 1 + Q2 / pow2(Ma) ); }

double FormFactor::GaNC(double Q2)
{
  double g = 0;
  if(part == particle::proton) g = ga;
  else if(part == particle::neutron) g = -ga;
  return 0.5/pow2(1 + Q2/pow2(Ma)) * (g - gas);
}

double FormFactor::Fs12(double Q2, double num)
{
  double F10 = 0;
  if(num == 1) F10 = 0.53;
  else if(num == 2) F10 = -0.4;

  return F10 * 1/(1+Q2/(4 * M)) * 1/pow2(1+Q2/(4 * Mv));
}

double FormFactor::Fp(double Q2)
{
  double q2 = -Q2;
    // std::cout<<"Fp: "<<2*M*M*GaCC(Q2)/(Mpi*Mpi + Q2)<<" M = "<<M*M<<" Mpi = "<<Mpi*Mpi<<" GaCC(Q2) = "<<GaCC(Q2)<<" Mpi*Mpi + Q2 = "<<Mpi*Mpi + Q2<<std::endl;

	return 2*M*M*GaCC(Q2)/(Mpi*Mpi + Q2);
  // return -1.267 / pow2(1 - q2 / pow2(Ma));
}