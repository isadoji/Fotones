#include "TMath.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/Functor.h"
#include "TTree.h"
#include "TFile.h"

const Int_t kMax=100000;
Int_t nT;
Double_t eqT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

double eq;
double beta=0.0;
double eta=3;
double betaf=TMath::Power(1-beta,2)*TMath::Power(1-TMath::Power(beta,2),-1);
//double betaf=1;
double omegaq;
double B;
double tau;
double Vol;
double lambdas=2;
double alphas = TMath::Power(TMath::Pi(),-1);
double mu;
int nc = 0;
double n1func( const double *  x){
   nc++;
   return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(mu,2)+
						  (TMath::Power(x[0],2)*betaf))*
                                      TMath::Power(lambdas,-1))-1,-1);
}
double n2func( const double *  x){
  nc++;
  return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(mu,2)+
                                                 (TMath::Power(omegaq-x[0],2)*betaf))*
                                     TMath::Power(lambdas,-1))-1,-1);
}
double argb( const double * x){
  nc++;
  return betaf*((TMath::Power(x[0],2))-(x[0]*(omegaq))+(TMath::Power(omegaq,2)))
    *TMath::Power(2*B*eq,-1);
}
double argB1( const double * x){
  nc++;
  return  (1+betaf*argb(x))*TMath::Power(betaf*argb(x),-1);
}
double pol( const double *  x){
  nc++;
  return betaf*(2*TMath::Power(x[0],2)-(x[0]*omegaq)+TMath::Power(omegaq,2))*
    TMath::Power(omegaq,-1);

}

double f2(const double *  x){
  nc++;
  return TMath::Exp(-argb(x))*((2*n1func(x)*n2func(x))+n2func(x))*
    (TMath::BesselI0(argb(x))-TMath::BesselI1(argb(x)));
  //pol(x)*
}

int yield()
{
  gsT->Branch("nT",&nT,"nT/I");
  gsT->Branch("eqT",eqT,"eqT[nT]/D");
  gsT->Branch("tT",tT,"tT[nT]/D");
  gsT->Branch("VT",VT,"VT[nT]/D");
  gsT->Branch("omegaT",omegaT,"omegaT[nT]/D");
  gsT->Branch("integralT",integralT,"integralT[nT]/D");
  
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
    //mu=4*alphas*B*TMath::Power(3*TMath::Pi(),-1);
    mu=0;
    std::cout << "B: " << B << " tau:  " << t << " Vol:  " << V << " mu:  " << mu << std::endl;
   
  for (Int_t k=0; k<2;k++) {
    eq = eqa[k];
    for (Int_t j=0; j<35; j++) {
        omegaq = 0.1*j;
	//    double omegaqa[13] = {0.01,0.05,0.1,0.2,0.3,0.4,0.5,1,1.5,2,2.5,3,3.2};
	//for (Int_t j=0; j<13; j++) {
	//omegaq = omegaqa[j];
      ROOT::Math::Functor wf(&f2,1);
      double a[2] = {0};
      double b[2] = {omegaq};
      
      ROOT::Math::GSLMCIntegrator * ig2 = new ROOT::Math::GSLMCIntegrator( ROOT::Math::IntegrationMultiDim::kVEGAS, 1.E-6, 1.E-4, 5000);
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
   nT++;

  }
    }
  }
  gsT->Fill();
  char name[60];
  sprintf(name,"20.40.beta.0.prueba.root");
  cout << " Writtiing file ->" << name << endl;
  TFile fOut(name, "recreate");
  fOut.cd();
  gsT->Write();
  return status;

}
