#include "v2.h"
#include "Math/GSLMCIntegrator.h"
#include "TMath.h"
#include "TF2.h"
#include "TFile.h"

#include "TTree.h"
#include "TCanvas.h"

#include <Riostream.h>
#include "TLegend.h"
#include "TLegendEntry.h"

#include "Math/Functor.h"
#include <cmath>
#include "TSystem.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "Math/SpecFuncMathMore.h"

const Int_t kMax=100000;
Int_t nT;
Double_t eqT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

void mcintegralv2() {
  gsT->Branch("nT",&nT,"nT/I");
  gsT->Branch("eqT",eqT,"eqT[nT]/D");
  gsT->Branch("tT",tT,"tT[nT]/D");
  gsT->Branch("VT",VT,"VT[nT]/D");
  gsT->Branch("omegaT",omegaT,"omegaT[nT]/D");
  gsT->Branch("integralT",integralT,"integralT[nT]/D");

  double omegaq;
  double eq;
  double xx;
  v2 * f1 = new v2();

  double omegaqa[13] = {0.01,0.05,0.1,0.2,0.3,0.4,0.5,1,1.5,2,2.5,3,3.2};
    double eqa[2] = {0.3,0.6};
    for (Int_t k=0; k<2;k++) {
      eq = eqa[k];
      for (Int_t j=0; j<13; j++) {
  
  // Obtener un metodo de integracion de MC
  ROOT::Math::GSLMCIntegrator * nminteg = new ROOT::Math::
    GSLMCIntegrator( ROOT::Math::IntegrationMultiDim::kVEGAS, 1.E-6, 1.E-4, 50000);
  
  // Pasar al integrador, la funcion
  nminteg->SetFunction( *(ROOT::Math::IMultiGenFunction*)f1);
    
  // Limites de integracion 
    double xmin[2];double xmax[2];
	omegaq = omegaqa[j];
	xmin[0] = 0.0; // minimo variable 1
	xmax[0] = omegaq; // maximo variable 1
	
	xmin[1] = 0.0; // minimo variable 2
	xmax[1] = omegaq; // maximo variable 2
    
    double result = nminteg->Integral( xmin , xmax );
    
    std::cout << xmin[0] << "  " << omegaq << "  " <<  "integral result: "
	      << result << " " << nminteg->ChiSqr() << std::endl;

    eqT[nT]=eq;
    omegaT[nT]=omegaq;
    integralT[nT]=result;
    nT++;
      
       
 // Limpiar la memoria   
  delete f1;
  delete nminteg;
      }
    }
  gsT->Fill();
  char name[60];
  sprintf(name,"prueba.root");
  cout << " Writtiing file ->" << name << endl;
  TFile fOut(name, "recreate");
  fOut.cd();
  gsT->Write();


}
