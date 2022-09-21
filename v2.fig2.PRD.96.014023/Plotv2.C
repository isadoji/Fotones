//-----------------------------------------------------------
// Description:
//
//       Plot macro. This macro is to plot Vegas integal v2mag of PhysRevD.96.014023.pdf.
//
// Environment:
//      ROOT
//
// Author List:
//       Isabel Domínguez Jiménez
//       isadoji@uas.edu.mx
//-----------------------------------------------------------
TH1D* v2H = new TH1D("v2H","",70,0,3.5);
const Int_t kMax=1000000;
Int_t nT;
Double_t eqT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t omegaT[kMax];
Double_t qintT[kMax];
Double_t integralT[kMax];

TTree* gsT = new  TTree("gsT","g's tree");

double alphaem = TMath::Power(137,-1);
double alphas = TMath::Power(TMath::Pi(),-1);
float mq = 0.0001;
int g = 2;
/*double tau = 0.06;
double VolT = 450; //{450.47, 383.875, 310.03, 256.823, 165.179};
double AnchoDisco = 0.25*TMath::Power(0.197,-1);
// (* Esfera de diametro 12fm con contraccion de Lorentz*)
double Vol = VolT*AnchoDisco;
*/
void Plotv2(){
  char name[50];
  sprintf(name,"v2.root");
  cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);
  TTree *gsT = (TTree*) file->Get("gsT");
  gsT->SetBranchAddress("nT",&nT);
  gsT->SetBranchAddress("eqT",eqT);
  gsT->SetBranchAddress("tT",tT);
  gsT->SetBranchAddress("VT",VT);
  gsT->SetBranchAddress("omegaT",omegaT);
  gsT->SetBranchAddress("qintT",qintT);
  gsT->SetBranchAddress("integralT",integralT);

  Float_t n=0;
  Int_t entries = (Int_t) gsT->GetEntries();
  for (Int_t i =0; i< entries;i++){
    gsT->GetEntry(i);
    for (Int_t j =0; j< nT;j++){
      gsT->GetEntry(j);
      Double_t eq = eqT[j];
      Double_t omega = omegaT[j];
      Double_t Vol = VT[j];
      Double_t tau = tT[j];
      Double_t integral = TMath::Power(TMath::Pi(),-6)*TMath::Power(2,-1)*
      TMath::Pi()*Vol*tau*TMath::Power(200,-1)*eq*eq*alphaem*alphas*alphas*
	integralT[j];

      //      if(integral >0){
      if(eq == 0.3){
	n++;
	v2H->Fill(omega,2*integral);
      }
      if(eq == 0.6){
	v2H->Fill(omega,integral);
      //}
    }
  }
  }

  TCanvas *pt = new TCanvas("pt", "pt particles distribution",50,50,500,500);
  gStyle->SetOptStat(false);
  // gStyle->SetOptStat("me");
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
  pt->SetRightMargin(0.0465116);
  pt->SetTopMargin(0.39916);
  pt->SetFillColor(0);

  TPad *pt_1 = new TPad("pt_1", "pt_1",0.0,0.0,0.98,0.98);
  pt_1->Draw();
  pt_1->cd();
  pt_1->Range(-20.6795,-0.0133333,105.496,0.062549);
  pt_1->SetFillColor(0);
  pt_1->SetBorderSize(2);
  pt_1->SetLeftMargin(0.17);
  pt_1->SetRightMargin(0.05);
  pt_1->SetTopMargin(0.1);
  pt_1->SetBottomMargin(0.15);
  pt_1->SetFrameLineWidth(2);
  pt_1->SetFrameLineWidth(2);
  //pt_1->SetLogy();

  v2H->GetXaxis()->SetTitle("#omega_{q} [GeV]");
  v2H->GetXaxis()->CenterTitle(true);
  v2H->GetXaxis()->SetTitleSize(0.04);
  v2H->GetXaxis()->SetLabelSize(0.03);
  v2H->GetXaxis()->SetTitleOffset(1.4);
  v2H->GetYaxis()->SetTitle("v_{2}^{mag}");
  v2H->GetYaxis()->CenterTitle(true);
  v2H->GetYaxis()->SetTitleSize(0.04);
  v2H->GetYaxis()->SetLabelSize(0.03);
  v2H->GetYaxis()->SetTitleOffset(2);

  Float_t norm = TMath::Sqrt(n);
  //  v2H->Scale(norm*2);
  v2H->Draw("hist p");
  //  v2H->GetYaxis()->SetRangeUser(0,2);
  //v2H->GetXaxis()->SetRangeUser(0.0,3.5);
  v2H->SetMarkerColor(1);
  v2H->SetMarkerStyle(20);
  v2H->SetMarkerSize(0.7);


  TLegend *leg2 = new TLegend(0.5,0.7,0.9,0.89);
  // leg2->SetTextFont(62);
  //leg2->SetTextSize(0.04);
  // leg2->SetLineColor(0);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(0);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry("","Au-Au #sqrt{s_{NN}} = 200 GeV","");
  //  leg2->AddEntry("v2H ","0-20 %","p");
  leg2->Draw();
  pt->SaveAs("v2mag.pdf");
  pt->SaveAs("v2mag.eps");
  //pt->SaveAs("v2mag.root");



    }
