//-----------------------------------------------------------
// Description:
//
//       Vegas integal macro. This macro is to integrate  V2mag in 1D of eq 18 of PhysRevD.96.014023.pdf numerical calculation.
//
// Environment:
//      ROOT
//
// Author List:
//       Isabel Domínguez Jiménez
//       isadoji@uas.edu.mx
//-----------------------------------------------------------
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

double eq =0.3; //{2/3, 1/3, 1/3}
float beta=0.25;
double eta=3;
double betaf=TMath::Power(1-beta,2)*TMath::Power(1-TMath::Power(beta,2),-1);
double omegaq;
//double B = 0.074;//t=0.6 fm
//double B = 0.055;//t = 1 fm
//double B = 0.006;//t=2 fm
double B = 0.02;//Jorge

double lambdas =2;
int nc = 0;

double n1func( const double *  x){
   nc++;
   return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(x[0],2)*betaf)*
                                      TMath::Power(lambdas,-1))-1,-1);
}
double n2func( const double *  x){
  nc++;
  return eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(0.1,2)+
                                                 TMath::Power(omegaq-x[0],2))*betaf*
                                     TMath::Power(lambdas,-1))-1,-1);
}
double argb( const double * x){
  nc++;
  return betaf*((TMath::Power(x[0],2))-(x[0]*(omegaq))+(TMath::Power(omegaq,2)))
    *TMath::Power(2*B*eq,-1);
}
double argB1( const double * x){
  nc++;
  return  (1+argb(x))*TMath::Power(argb(x),-1);
}
double pol( const double *  x){
  nc++;
  return betaf*(2*TMath::Power(x[0],2)-(x[0]*omegaq)+TMath::Power(omegaq,2));

}

double f2(const double *  x){
  nc++;
  return TMath::Exp(-argb(x))*((1*n1func(x)*n2func(x)))*pol(x)*
    (TMath::BesselI0(argb(x))-argB1(x)*TMath::BesselI1(argb(x)));
}

int v2_1D()
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

  for (Int_t k=0; k<2;k++) {
    eq = eqa[k];
    for (Int_t j=0; j<350; j++) {
      omegaq = 0.01*j;
      ROOT::Math::Functor wf(&f2,1);
      double a[1] = {0};
      double b[1] = {omegaq};

      ROOT::Math::GSLMCIntegrator * ig2 = new ROOT::Math::GSLMCIntegrator( ROOT::Math::IntegrationMultiDim::kVEGAS, 1.E-6, 1.E-4, 10000);
      ig2->SetFunction(wf);
      double val = ig2->Integral(a,b);
      //   std::cout << "integral result: " << val << " " << ig2->ChiSqr() << std::endl;
   std::cout << "integral result: " << val << " " << omegaq << std::endl;
   status += std::fabs(val-RESULT) > ERRORLIMIT;
   eqT[nT]=eq;
   //   tT[nT]=tau;
   //VT[nT]=Vol;
   omegaT[nT]=omegaq;
   integralT[nT]=val;
   nT++;

  }
    }

  gsT->Fill();
  char name[60];
  sprintf(name,"prueba1D.root");
  cout << " Writtiing file ->" << name << endl;
  TFile fOut(name, "recreate");
  fOut.cd();
  gsT->Write();
  return status;

}
