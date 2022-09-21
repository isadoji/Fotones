TH1D* v2H = new TH1D("v2H","",35,0,3.5);
TH1D* dataH = new TH1D("dataH","",35,0,35);
TH1D* tpc1H = new TH1D("tpc1H","",35,0,3.5);

const Int_t kMax=10000000;
Int_t nT;
Double_t eqT[kMax];
Double_t omegaT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

double alphaem = TMath::Power(137,-1);
double alphas = TMath::Power(TMath::Pi(),-1);

void Plotv2(){
  char name[50];
  sprintf(name,"prueba.root");
  cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);
  TTree *gsT = (TTree*) file->Get("gsT");
  gsT->SetBranchAddress("nT",&nT);
  gsT->SetBranchAddress("eqT",eqT);
  gsT->SetBranchAddress("tT",tT);
  gsT->SetBranchAddress("VT",VT);
  gsT->SetBranchAddress("omegaT",omegaT);
  gsT->SetBranchAddress("integralT",integralT);

  Int_t entries = (Int_t) gsT->GetEntries();
  for (Int_t i =0; i< entries;i++){
    gsT->GetEntry(i);
    for (Int_t j =0; j< nT;j++){
      gsT->GetEntry(j);
      Double_t eq = eqT[j];
      Double_t omega = omegaT[j];
      Double_t tau = tT[j];
      Double_t Vol = VT[j];
      //cout << nT << endl;
      if(tau> 0.05 && tau < 0.07){
	//Double_t integral = TMath::Power(2*TMath::Pi(),-5)*TMath::Power(2,-1)*
	//  TMath::Pi()*Vol*tau*TMath::Power(200,-1)*eq*eq*alphaem*alphas*alphas*
	 Double_t integral = integralT[j];
	if (integral > -10000){
	  cout << integral << " " << tau << endl;
	  if(eq == 0.3) v2H->Fill(omega,2*integral);
	  if(eq == 0.6) v2H->Fill(omega,integral);
	  }
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
  v2H->GetXaxis()->SetTitle("#v");
  v2H->GetXaxis()->CenterTitle(true);
  v2H->GetXaxis()->SetTitleSize(0.04);
  v2H->GetXaxis()->SetLabelSize(0.03);
  v2H->GetXaxis()->SetTitleOffset(1.4);
  v2H->GetYaxis()->SetTitle("(dN/dv");
  v2H->GetYaxis()->CenterTitle(true);
  v2H->GetYaxis()->SetTitleSize(0.04);
  v2H->GetYaxis()->SetLabelSize(0.03);
  v2H->GetYaxis()->SetTitleOffset(2);

  Float_t scale = TMath::Power(2*nT, -1);
  //v2H = (TH1D*)v2H->DrawNormalized("hist p",1/v2H->Integral());
  //v2H->Scale(scale);
  v2H->Draw("hist p");
  //  v2H->GetYaxis()->SetRangeUser(0.00001,20);
  //v2H->GetXaxis()->SetRangeUser(0.0,3.5);
  v2H->SetMarkerColor(1);
  v2H->SetMarkerStyle(20);
  v2H->SetMarkerSize(0.7);

  /*dataH->Add(tpc1H,-1);
  dataH->Draw("sames hist p ");
  dataH->SetMarkerColor(1);                                                                  
  dataH->SetMarkerStyle(24);
  dataH->SetMarkerSize(0.7);
  */
  
  TLegend *leg2 = new TLegend(0.5,0.7,0.9,0.89);
  // leg2->SetTextFont(62);
  //leg2->SetTextSize(0.04);
  // leg2->SetLineColor(0);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(0);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry("","UrQMD, Au-Au #sqrt{s_{NN}} = 200 GeV","");
  leg2->AddEntry("v2H","0-20 %","p");
  //leg2->Draw();
  //  pt->SaveAs("V2.200GeV.old.pdf");
  //pt->SaveAs("V2.200GeV.old.eps");

    }
