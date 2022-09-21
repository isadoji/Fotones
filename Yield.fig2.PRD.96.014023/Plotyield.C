//-----------------------------------------------------------
// Description:
//
//       Plot macro. This macro is to plot Vegas integal of Fig2.a of PhysRevD.96.014023.pdf.
//
// Environment:
//      ROOT
//
// Author List:
//       Isabel Domínguez Jiménez
//       isadoji@uas.edu.mx
//-----------------------------------------------------------
TH1D* b0H = new TH1D("b0H","",35,0,3.5);
TH1D* sumH = new TH1D("sumH","",35,0,3.5);
TH1D* dataH = new TH1D("dataH","",35,0,3.5);
TH1D* tpc1H = new TH1D("tpc1H","",35,0,3.5);
TH1D* tpc2H = new TH1D("tpc2H","",35,0,3.5);

double AnchoDisco = 0.25*TMath::Power(0.197,-1);
double VolT = 450; //{450.47, 383.875, 310.03, 256.823, 165.179};
double alphaem = TMath::Power(137,-1);
double alphas = TMath::Power(TMath::Pi(),-2);
double tau=0.6/200; //fm->Gev^-1

float mq = 0.0001;
int g = 2;
////REVISAR CONSTANTES///////
double Vol= AnchoDisco*VolT*alphaem*alphas*tau
  *TMath::Pi()*TMath::Power(2,-1)*TMath::Power(2*TMath::Pi(),-6);



const Int_t kMax=10000000;
Int_t nT;
Double_t eqT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

void Plotyield(){
  char name[50];
  sprintf(name,"yield.root");
  cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);

  TTree *gsT = (TTree*) file->Get("gsT");
  gsT->SetBranchAddress("nT",&nT);
  gsT->SetBranchAddress("eqT",eqT);
  gsT->SetBranchAddress("omegaT",omegaT);
  gsT->SetBranchAddress("integralT",integralT);

  FILE *fp = fopen("dataexp/phenix.txt","r");
  FILE *fp1 = fopen("dataexp/tpnc1.txt","r");
  Float_t x,y,ex,ey;
  Float_t x0,y0,ex0,ey0;
  Float_t x1,y1;
  Float_t x2,y2;

  Int_t ncols, ncols0, ncols1, ncols2;
  Int_t nlines = 0;

 while (1) {
    ncols = fscanf(fp,"%f %f %f %f",&x, &y, &ex, &ey);
    if (ncols <= 0) break;
    dataH->Fill(x,y);
    tpc1H->Fill(x1,y1);
  }

  Int_t entries = (Int_t) gsT->GetEntries();
  for (Int_t i =0; i< entries;i++){
    gsT->GetEntry(i);
    cout << nT << endl;
    for (Int_t j =0; j< nT;j++){
    gsT->GetEntry(j);
    //    cout << "entries=  " << j << endl;
    Double_t eq = eqT[j];
    Double_t omega = omegaT[j];
    Double_t integral = eq*eq*Vol*integralT[j];
    //cout << tau << endl;
    if(eq == 0.3) {
      b0H->Fill(omega,2*integral);
    }
    if(eq == 0.6) {
      b0H->Fill(omega,integral);
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


  //  b0H = (TH1D*)b0H->DrawNormalized("",10*b0H->Integral());
  b0H->Scale(nT*10/b0H->Integral());
  b0H->GetYaxis()->SetRangeUser(0.00001,10000);
  //b0H->GetXaxis()->SetRangeUser(0.0,3.2);
  b0H->Draw("histo p");
  b0H->SetMarkerColor(1);
  b0H->SetMarkerStyle(3);
  dataH->Add(tpc1H,-1);
  dataH->Draw("sames hist p ");
  dataH->SetMarkerColor(1);
  dataH->SetMarkerStyle(24);
  dataH->SetMarkerSize(0.7);


  TLegend *leg2 = new TLegend(0.6,0.7,0.89,0.88);
  leg2->SetTextFont(62);
  //leg2->SetTextSize(0.04);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(0);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);

  leg2->AddEntry("dataH","PHENIX Au+Au 20-40%  - Data","p");
  //leg2->AddEntry("b0H","UrQMD, Au-Au #sqrt{s_{NN}}=200 GeV","p");
  //leg2->AddEntry("","#tau = 0.06 fm","");
  leg2->Draw();
  pt->SaveAs("Yield.Au-Au.20-40.pdf");
}
