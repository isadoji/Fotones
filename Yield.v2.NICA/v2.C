#include "TMath.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/Functor.h"
#include "TTree.h"
#include "TFile.h"

const Int_t kMax=10000000;
Int_t nT;
Double_t eqT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
Double_t qintT[kMax];
Double_t integralv2T[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

double eq;
double omegaq,qint;
float beta=0.15;
double eta=3;
double betaf=TMath::Power(1-beta,2)*TMath::Power(1-TMath::Power(beta,2),-1);
double B;
double tau;
double Vol;
double lambdas=2;
int nc = 0;

double n1funcv2( const double *  x){
   nc++;
   return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(x[0],2)*betaf)*
                                      TMath::Power(lambdas,-1))-1,-1);
}
double n2funcv2( const double *  x){
  nc++;
  return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(x[1]-x[0],2)*betaf)*
                                     TMath::Power(lambdas,-1))-1,-1);
}
double argbv2( const double * x){
  nc++;
  return betaf*((TMath::Power(x[0],2))-(x[0]*(x[1]))+(TMath::Power(x[1],2)))
    *TMath::Power(2*B*eq,-1);
}
double argB1v2( const double * x){
  nc++;
  return  (1+argbv2(x))*TMath::Power(argbv2(x),-1);
}
double polv2( const double *  x){
  nc++;
  return betaf*(2*TMath::Power(x[0],2)-(x[0]*x[1])+TMath::Power(x[1],2));

}
double f2v2(const double *  x){
  nc++;
  return TMath::Exp(-argbv2(x))*((2*n1funcv2(x)*n2funcv2(x))+n1funcv2(x))*polv2(x)*
    (TMath::BesselI0(argbv2(x))-argB1v2(x)*TMath::BesselI1(argbv2(x)));
}

int v2()
{
  gsT->Branch("nT",&nT,"nT/I");
  gsT->Branch("eqT",eqT,"eqT[nT]/D");
  gsT->Branch("tT",tT,"tT[nT]/D");
  gsT->Branch("VT",VT,"VT[nT]/D");
  gsT->Branch("omegaT",omegaT,"omegaT[nT]/D");
  gsT->Branch("integralT",integralT,"integralT[nT]/D");
  gsT->Branch("qintT",qintT,"qintT[nT]/D");
  gsT->Branch("integralv2T",integralv2T,"integralv2T[nT]/D");
  
  const double RESULT = 1.0;
  const double ERRORLIMIT = 1E-3;
  int status = 0;
  double eqa[2] = {0.3,0.6};

  FILE *fp = fopen("20.40.txt","r");
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
    std::cout << "B: " << B << " tau:  " << t << " Vol:  " << V << std::endl;
   
  for (Int_t k=0; k<2;k++) {
    eq = eqa[k];
     for (Int_t i=0; i<35; i++) {
      qint = 0.1*i;
      for (Int_t j=0; j<35; j++) {
        omegaq = 0.1*j;
      ROOT::Math::Functor wf(&f2v2,2);
      double a[2] = {0,0};
      double b[2] = {omegaq,qint};
      
      ROOT::Math::GSLMCIntegrator * ig2 = new ROOT::Math::GSLMCIntegrator( ROOT::Math::IntegrationMultiDim::kVEGAS, 1.E-6, 1.E-4, 100);
      ig2->SetFunction(wf);
      double val = ig2->Integral(a,b);
      //   std::cout << "integral result: " << val << " " << ig2->ChiSqr() << std::endl;
   std::cout << "integral result: " << val << " " << omegaq << std::endl;
   status += std::fabs(val-RESULT) > ERRORLIMIT;
   eqT[nT]=eq;
   tT[nT]=tau;
   VT[nT]=Vol;
   omegaT[nT]=omegaq;
   integralT[nT]=val;
   qintT[nT]=qint; 
   integralv2T[nT]=val;
  nT++;
      }
  }
    }
  }
  
  gsT->Fill();
  char name[60];
  sprintf(name,"20.40.v2.root");
  cout << " Writtiing file ->" << name << endl;
  TFile fOut(name, "recreate");
  fOut.cd();
  gsT->Write();
  return status;

}
