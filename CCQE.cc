#include <stdlib.h>
#include <time.h>
#include <vector>

#include <TGraph.h>
#include <TCanvas.h>
#include <TH1D.h>

#include "units.hh"
#include "nuCCQE.hh"
#include "vect.h"

using namespace std;

double jakob(vect nu,vect N0,vect lepton)
{
  vec vcms = (vect(nu) + vect(N0)).v ();
  if(vcms*vcms>=1)
     return 0;
  else
  {
    nu.boost (-vcms);
    lepton.boost (-vcms);	
    return vec (nu).length () * vec (lepton).length () * 4;
  }
}

int main(int argc, char const *argv[])
{
  if (argc != 5) 
  {
    std::cerr << "Usage: " << argv[0] << " ifCCQE if_proton if_muon if_neutrino" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  bool ifCC = argv[1];
  bool if_proton = argv[2];
  bool if_muon = argv[3];
  bool if_neutrino = argv[4];
	G4INCL::nuQELChannel qel(ifCC, if_proton, if_muon);

	srand(time(NULL));

	// double Enu=((double)rand()/(double)RAND_MAX) * 10000 * MeV; // random double between 0 and 10 GeV
	// double Eb = 34 * MeV; 																			// binding energy in C
	// const double s = std::pow(initial_nucleon_mass,2) + 2 * (initial_nucleon_Energy*neutrino_Energy - initial_nucleon_Momentum.dot(neutrino_Momentum)); 

  const double e_nu = 300 * MeV;
  vector<double> v_QQ;
  
  for (int i = 1; i <= 10; ++i)
  {
    v_QQ.push_back((i * 100) * (i * 100));
  }

  vector<double> v_points;

  for(int q = 0; q < v_QQ.size(); q++)
  {
    v_points.push_back(qel.QELAmplitude(-v_QQ[q], e_nu, if_neutrino));
    cout<<v_QQ[q]<<" "<<qel.QELAmplitude(-v_QQ[q], e_nu, if_neutrino)<<endl;
  }

	return 0;
}