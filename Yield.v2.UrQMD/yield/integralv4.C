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

const Int_t kMax=100000;
Int_t nT;
Double_t eqT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

double eq;
float beta=0.15;
double eta=3;
double betaf=TMath::Power(1-beta,2)*TMath::Power(1-TMath::Power(beta,2),-1);
//double betaf=1;
double omegaq;
double B;
double tau;
double Vol;
double lambdas=2;
int nc = 0;

double n1func( double x){
   nc++;
   return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(x,2)*betaf)*
				      TMath::Power(lambdas,-1))-1,-1);
}
double n2func( double x){
  nc++;
  return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(0.1,2)+
						 TMath::Power(omegaq-x,2))*betaf*
				     TMath::Power(lambdas,-1))-1,-1);
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

void  testIntegPerf(double x1, double x2, int n = 10000){

  gsT->Branch("nT",&nT,"nT/I");
  gsT->Branch("eqT",eqT,"eqT[nT]/D");
  gsT->Branch("tT",tT,"tT[nT]/D");
  gsT->Branch("VT",VT,"VT[nT]/D");
  gsT->Branch("omegaT",omegaT,"omegaT[nT]/D");
  gsT->Branch("integralT",integralT,"integralT[nT]/D");

  double eqa[2] = {0.3,0.6};
  double omegaqa[13] = {0.01,0.05,0.1,0.2,0.3,0.4,0.5,1,1.5,2,2.5,3,3.2};
  double s1;

  FILE *fp = fopen("40.60.txt","r");
  Float_t t,eBxm,eBy,eBym,eBzm, eBmpi,eB,V;
  //FILE *fp = fopen("0.40.Cu.txt","r");
  //Float_t t,Npart,eBmpi,eB,V;
  Int_t ncols;
  Int_t nlines = 0;
  while (1) {
       ncols = fscanf(fp,"%f %f %f %f %f %f %f %f",&t, &eBxm, &eBy,
    		   &eBym, &eBzm, &eBmpi, &eB, &V);
       //ncols = fscanf(fp,"%f %f %f %f %f",&t,&Npart,&eBmpi,&eB,&V);
    if (ncols < 0) break;
    B = eB;
    tau=t;
    Vol=V;
    if (tau >0 && tau < 0.12){
      for (Int_t k=0; k<2;k++) {
      eq = eqa[k];
      for (Int_t j=0; j<13; j++) {
      omegaq = omegaqa[j];
      s1= 0.0;
      ROOT::Math::WrappedFunction<> f1(func);
      ROOT::Math::Integrator ig(f1 );
      nc = 0;
      //double x = omegaq;
      for (int i = 0; i < n+1; ++i) {
	double dx = (omegaq-x1)/double(n);
	double x = x1 + dx*i;
	s1 += ig.Integral(x1,x);
	//printf("%8f %8f %8f\n",x1, x, omegaq);
      }
      printf("%8f %8f %8f %8f\n",tau,omegaq,Vol,s1);
      
      eqT[nT]=eq;
      tT[nT]=tau;
      VT[nT]=Vol;
      omegaT[nT]=omegaq;
      integralT[nT]=s1;
      nT++;
      }
      }
    }
  }
    nlines++;
  fclose(fp);

  gsT->Fill();
  char name[60];
  sprintf(name,"40.60.Au.Au.beta0.15.root");
  cout << " Writtiing file ->" << name << endl;
  TFile fOut(name, "recreate");
  fOut.cd();
  gsT->Write();

}
void integralv4(double a = 0.001, double b = 3.51)
{
  testIntegPerf(a, b);
}
