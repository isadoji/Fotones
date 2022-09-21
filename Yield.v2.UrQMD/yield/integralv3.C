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

const Int_t kMax=1000000;
Int_t nT;
Double_t eqT[kMax];
Double_t tT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");


float mq = 0.0001;
int g = 2;
double alphas = TMath::Power(TMath::Pi(),-1);
double alphaem = TMath::Power(137,-1);
double gammas=2;

float beta=0.15;
double eta=3;
double betaf=TMath::Power(1-beta,2)*TMath::Power(1-TMath::Power(beta,2),-1);
//double tau = 0.3;//0.06 fm
float tau;//0.03 fm
double omegaq;
double B;
double eq=0.3;
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
double pol( double x){
  nc++;
  return betaf*(2*TMath::Power(x,2)-(x*omegaq)+TMath::Power(omegaq,2))*
    TMath::Power(omegaq,-1);
}
double func( double x){
  nc++;
   return TMath::Exp(-argb(x))*n1func(x)*n2func(x)*pol(x)*
     (TMath::BesselI0(argb(x))-TMath::BesselI1(argb(x)));
}

void  testIntegPerf(double x1, double x2, int n = 10){
  double dx = (x2-x1)/double(n);
  
  gsT->Branch("nT",&nT,"nT/I");
  gsT->Branch("eqT",eqT,"eqT[nT]/D");
  gsT->Branch("tT",tT,"tT[nT]/D");
  gsT->Branch("omegaT",omegaT,"omegaT[nT]/D");
  gsT->Branch("integralT",integralT,"integralT[nT]/D");

  double eqa[2] = {0.3,0.6};
  double omegaqa[8] = {0.1,0.5,1,1.5,2,2.5,3,3.5};
  double s1;

  // FILE *fp = fopen("datos20-40.t0.5.txt","r");
  //Float_t t,eBxm,eBy,eBym,eBzm, eBmpi,eB,V;
  FILE *fp = fopen("0.40.Cu.txt","r");
  Float_t t,Npart,eBmpi,eB,V;
  Int_t ncols;
  Int_t nlines = 0;
  while (1) {
    //  ncols = fscanf(fp,"%f %f %f %f %f %f %f %f",&t, &eBxm, &eBy,
    //		   &eBym, &eBzm, &eBmpi, &eB, &V);
    ncols = fscanf(fp,"%f %f %f %f %f",&t,&Npart,&eBmpi,&eB,&V);
    if (ncols < 0) break;    
    B = eB;
    tau = t*TMath::Power(0.2,-1);
    
    for (Int_t k=0; k<2; k++) {
      double eq = eqa[k];
      double Vol= V*alphaem*alphas*alphas*eq*eq*tau*TMath::Power(2,-1)*
	TMath::Power(2*TMath::Pi(),-6);
      for (Int_t j=0; j<8; j++) {
      s1= 0.0;
      ROOT::Math::WrappedFunction<> f1(func);
      ROOT::Math::Integrator ig(f1 );
      nc = 0;
      for (int i = 0; i < n; ++i) {
	double x = x1 + dx*i;
	s1 += ig.Integral(x1,x);
	}

      printf("%8f %8f %8f %8f\n",eq, tau, omegaq, tau*Vol*s1 );
      eqT[nT]=eq;
      tT[nT]=tau;
      omegaT[nT]=omegaq;
      integralT[nT]=tau*Vol*s1;
      nT++;
    }
    nlines++;
    }
    }
  fclose(fp);

  gsT->Fill();
  char name[60];
  sprintf(name,"prueba.root");
  cout << " Writtiing file ->" << name << endl;
  TFile fOut(name, "recreate");
  fOut.cd();
  gsT->Write();
}
void integralv3(double a = 0.1, double b = 3)
{
  testIntegPerf(a, b);
}

