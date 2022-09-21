TH1D* b0H = new TH1D("b0H","",40,0,4);
TH1D* b4H = new TH1D("b4H","",40,0,4);
TH1D* b9H = new TH1D("b9H","",40,0,4);
TH1D* b14H = new TH1D("b14H","",40,0,4);


const Int_t kMax=1000000;
Int_t nT;
Double_t eqT[kMax];
Double_t bT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

void Plot(){
  char name[50];
  //  sprintf(name,"integral0.06fm.root");
    sprintf(name,"prueba.root");
  cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);

  TTree *gsT = (TTree*) file->Get("gsT");
  gsT->SetBranchAddress("nT",&nT);
  gsT->SetBranchAddress("eqT",eqT);
  gsT->SetBranchAddress("bT",bT);
  gsT->SetBranchAddress("omegaT",omegaT);
  gsT->SetBranchAddress("integralT",integralT);


  Int_t entries = (Int_t) gsT->GetEntries();
  for (Int_t i =0; i< entries;i++){
    gsT->GetEntry(i);


  for (Int_t j =0; j< 272;j++){
    gsT->GetEntry(j);
    //cout << "entries=  " << j << endl;
    Double_t eq = eqT[j];
      Double_t b = bT[j];
      Double_t omega = omegaT[j];
      Double_t integral = integralT[j];

      if(b == 0 && eq == 0.3) {
	b0H->Fill(omega,2*integral);
      }
      if(b== 0 && eq == 0.6) {
	b0H->Fill(omega,integral);
      }
      if(b == 4 && eq == 0.3) {
	b4H->Fill(omega,2*integral);
	}
      if(b== 4 && eq == 0.6) {
	  b4H->Fill(omega,integral);
      }
      if(b == 9 && eq == 0.3) {
	b9H->Fill(omega,2*integral);
	}
      if(b== 9 && eq == 0.6) {
	  b9H->Fill(omega,integral);
      }
      if(b == 14 && eq == 0.3) {
	b14H->Fill(omega,2*integral);
	}
      if(b== 14 && eq == 0.6) {
	  b14H->Fill(omega,integral);
      }
  }
  }


  TCanvas *pt = new TCanvas("pt", "pt particles distribution",50,50,500,500);
  gStyle->SetOptStat(false);
  //  gStyle->SetOptStat("me");
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
  pt_1->SetLogy();

  b0H->GetXaxis()->SetTitle("#omega_{q} [GeV]");
  b0H->GetXaxis()->CenterTitle(true);
  b0H->GetXaxis()->SetTitleSize(0.04);
  b0H->GetXaxis()->SetLabelSize(0.03);

  b0H->GetXaxis()->SetTitleOffset(1.4);
  b0H->GetYaxis()->SetTitle("(1/2#pi#omega_{q}) dN/d#omega_{q} [GeV^{-2}]");
  b0H->GetYaxis()->CenterTitle(true);
  b0H->GetYaxis()->SetTitleSize(0.04);
  b0H->GetYaxis()->SetLabelSize(0.03);
  b0H->GetYaxis()->SetTitleOffset(2);


  b0H->Draw("hist p");
  b0H->SetMarkerColor(1);
  b0H->SetMarkerStyle(3);
  b4H->Draw("same hist p");
  b4H->SetMarkerColor(1);
  b4H->SetMarkerStyle(21);
  b9H->Draw("same hist p");
  b9H->SetMarkerColor(1);
  b9H->SetMarkerStyle(4);
  b14H->Draw("same hist p");
  b14H->SetMarkerColor(1);
  b14H->SetMarkerStyle(23);

TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);
 leg->AddEntry("","#tau = 0.06 fm","");
 leg->AddEntry("b0H","b = 0","p");
 leg->AddEntry("b4H","b = 4","p");
 leg->AddEntry("b9H","b = 9","p");
 leg->AddEntry("b14H","b = 14","p");
 leg->Draw();
}
