#include <stdlib.h>
#include <time.h>
#include <vector>

#include <TGraph.h>
#include <TCanvas.h>
#include <TH1D.h>

#include "jednostki.h"
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

void GraphStyle(TGraph* gr, int marker, int color/*, int line*/)
{
  gr->SetMarkerStyle(marker);
  gr->SetDrawOption("PC");
  gr->SetLineColorAlpha(color, 0.5);
  gr->SetMarkerColor(color);
  gr->SetLineWidth(4);
  // gr->SetLineStyle(line);
  gr->SetFillStyle(0);
}

int main(int argc, char const *argv[])
{
	G4INCL::nuQELChannel qel(false, false, false);

	srand(time(NULL));

	// double Enu=((double)rand()/(double)RAND_MAX) * 10000 * MeV; // random double between 0 and 10 GeV
	// double Eb = 34 * MeV; 																			// binding energy in C
	
	// const double s = std::pow(initial_nucleon_mass,2) + 2 * (initial_nucleon_Energy*neutrino_Energy - initial_nucleon_Momentum.dot(neutrino_Momentum)); 

  const double e_nu = 300 * MeV;
  vector<double> v_QQ;
  
  for (int i = 1; i <= 10; ++i)
  {
    v_QQ.push_back((i*100)*(i*100));
  }

  const int n_trials = 1;
  vector<double> v_points;

  for(int q = 0; q < v_QQ.size(); q++)
  {
    double sum = 0;
    for (int n = 0; n < n_trials; n++)
    {
      double ampl = qel.QELAmplitude(-v_QQ[q], e_nu, true);
      sum += ampl;
      // cout<<"-v_QQ[q] = "<<-v_QQ[q]<<endl;
    }
    v_points.push_back(sum / n_trials);
  }

  TGraph* gr = new TGraph(v_QQ.size(), v_QQ.data(), v_points.data());
  TCanvas* c1 = new TCanvas();
  GraphStyle(gr, 20, kBlue);
  gr->Draw("APC");
  c1->Print(("qel_test_"+to_string(e_nu)+".pdf").c_str());

  double x[100] = {0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 
                0.43, 0.45, 0.47, 0.49, 0.51, 0.53, 0.55, 0.57, 0.59, 0.61, 0.63, 0.65, 0.67, 0.69, 0.71, 0.73, 0.75, 0.77, 0.79, 0.81, 0.83, 
                0.85, 0.87, 0.89, 0.91, 0.93, 0.95, 0.97, 0.99, 1.01, 1.03, 1.05, 1.07, 1.09, 1.11, 1.13, 1.15, 1.17, 1.19, 1.21, 1.23, 1.25, 
                1.27, 1.29, 1.31, 1.33, 1.35, 1.37, 1.39, 1.41, 1.43, 1.45, 1.47, 1.49, 1.51, 1.53, 1.55, 1.57, 1.59, 1.61, 1.63, 1.65, 1.67, 
                1.69, 1.71, 1.73, 1.75, 1.77, 1.79, 1.81, 1.83, 1.85, 1.87, 1.89, 1.91, 1.93, 1.95, 1.97, 1.99};

  double y[100] = {2.07915, 2.07744, 2.10644, 2.0738, 2.00386, 2.03852, 1.99474, 1.96889, 1.88606, 1.87903, 1.78467, 1.79862, 1.73328, 1.68337, 
                   1.65979, 1.55518, 1.56537, 1.49821, 1.42726, 1.44424, 1.37551, 1.34936, 1.30625, 0.583502, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  double x_incl[10] = {0.01, 0.04, 0.09, 0.16, 0.25, 0.36, 0.49};

  TGraph* gr_nuwro_test = new TGraph(100, x, y);
  GraphStyle(gr_nuwro_test, 20, kRed);
  gr_nuwro_test->Draw("APC");
  c1->Print("nuwro_test.pdf");

  TH1D* h_nuwro = new TH1D("h_nuwro", "", 100, x);
  for (int i = 0; i < 100; ++i)
  {
    h_nuwro->SetBinContent(i, y[i]);
  }
  h_nuwro->SetLineColor(kRed);

  for (int i = 0; i < v_QQ.size(); ++i) v_QQ[i]/=1000000;

  TH1D* h_incl = new TH1D("h_incl", "", 6, x_incl);

  for (int i = 0; i < 7; ++i)
  {
    cout<<v_QQ.size()<<" "<<v_QQ[i]<<endl;
    h_incl->SetBinContent(i, v_points[i]);
  }

  h_nuwro->Scale(1./h_nuwro->Integral());
  h_incl->Scale(1./h_incl->Integral());

  h_nuwro->Scale(1.,"WIDTH");
  h_incl->Scale(1.,"WIDTH");

  h_incl->Draw("HIST");
  h_nuwro->Draw("HIST SAME");

  c1->Print("test.pdf");
	return 0;
}