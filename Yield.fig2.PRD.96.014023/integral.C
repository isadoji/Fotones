//-----------------------------------------------------------
// Description:
//
//       Vegas integal macro. This macro is to integrate  eq 14 of PhysRevD.96.014023.pdf numerical calculation.
//
// Environment:
//      ROOT
//
// Author List:
//       Isabel Domínguez Jiménez
//       isadoji@uas.edu.mx
//-----------------------------------------------------------

#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include <iostream>
#include "TStopwatch.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include <limits>

const Int_t kMax=10000000;
Int_t nT;
Double_t eqT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

double eq =0.3; //{2/3, 1/3, 1/3}
float beta=0.0;
double eta=3;
double betaf=TMath::Power(1-beta,2)*TMath::Power(1-TMath::Power(beta,2),-1);
double omegaq;
double B = 0.0249;
double gammas=2;
int nc = 0;
double n1func( double x){
   nc++;
   return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(x,2)*betaf)*
				      TMath::Power(gammas,-1))-1,-1);
}
double n2func( double x){
  nc++;
  return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(0.1,2)+
						 TMath::Power(omegaq-x,2))*betaf*
				     TMath::Power(gammas,-1))-1,-1);
}
double argb( double x){
  nc++;
  return betaf*((TMath::Power(x,2))-(x*(omegaq))+(TMath::Power(omegaq,2)))
    *TMath::Power(2*B*eq,-1);
}
double argB1( double x){
  nc++;
  return (1+betaf*argb(x))*TMath::Power(betaf*argb(x),-1);
}

double pol( double x){
  nc++;
  return betaf*(2*TMath::Power(x,2)-(x*omegaq)+TMath::Power(omegaq,2))*TMath::Power(omegaq,-1);
}
double func( double x){
  nc++;
  return TMath::Exp(-argb(x))*n1func(x)*n2func(x)*pol(x)*(TMath::BesselI0(argb(x))-TMath::BesselI1(argb(x)));
}

void  testIntegPerf(double x1, double x2, int n=1){
  double dx = (x2-x1)/double(n);
  gsT->Branch("nT",&nT,"nT/I");
  gsT->Branch("eqT",eqT,"eqT[nT]/D");

  gsT->Branch("omegaT",omegaT,"omegaT[nT]/D");
  gsT->Branch("integralT",integralT,"integralT[nT]/D");

  double eqa[2] = {0.3,0.6};

  double omegaqa[15] = {0.01,0.02,0.03,0.04,0.05,0.2,0.3,0.4,0.5,1,1.5,2,2.5,3,3.2};
  double s1;
  for (Int_t k=0; k<2; k++) {
      double eq = eqa[k];
	    for (Int_t j=0; j<15; j++) {
          omegaq = omegaqa[j];
	s1= 0.0;
	ROOT::Math::WrappedFunction<> f1(func);
	ROOT::Math::Integrator ig(f1 );
	nc = 0;
//	for (int i = 0; i < n; ++i) {
	 //double x = x1+ dx*i;
	double x=omegaq;
	  s1 = ig.Integral(x1,x);
	  printf("%8f,%8f,%8f. %8f\n",x1,x,omegaq,s1);
	  //}
	  //	printf("%8f, %8f\n",omegaq,s1);
	  eqT[nT]=eq;
	 	omegaT[nT]=omegaq;
	  integralT[nT]=s1;
	  nT++;
  }
  }

gsT->Fill();
char name[60];
  sprintf(name,"yield.root");
  cout << " Writtiing file ->" << name << endl;
  TFile fOut(name, "recreate");
  fOut.cd();
  gsT->Write();

}
void integral(double a = 0.001, double b = 3.2)
{
  testIntegPerf(a, b);
}
