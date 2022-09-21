TH1D* yield03H = new TH1D("yield03H","",35,0,3.5);
TH1D* dataH = new TH1D("dataH","",35,0,3.5);
TH1D* tpc1H = new TH1D("tpc1H","",35,0,3.5);

const Int_t kMax=100000;
Int_t nT;
Double_t eqT[kMax];
Double_t omegaT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

double alphaem = TMath::Power(137,-1);
double alphas = TMath::Power(TMath::Pi(),-1);
float mq = 0.0001;
int g = 2;

void Plotyield(){
  char name[50];
  sprintf(name,"20.40.beta.0.root");
  cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);
  TTree *gsT = (TTree*) file->Get("gsT");
  gsT->SetBranchAddress("nT",&nT);
  gsT->SetBranchAddress("eqT",eqT);
  gsT->SetBranchAddress("tT",tT);
  gsT->SetBranchAddress("VT",VT);
  gsT->SetBranchAddress("omegaT",omegaT);
  gsT->SetBranchAddress("integralT",integralT);
  //   gsT->Scan("integralT");

  FILE *fp = fopen("phenix20.40.txt","r");
  Int_t ncols;
  Float_t x,y,ex,ey;
  while (1) {
    ncols = fscanf(fp,"%f %f %f %f",&x, &y, &ex, &ey);
    if (ncols <= 0) break;
    dataH->Fill(x,y);
  }
  
  FILE *fp1 = fopen("tpnc20.40.txt","r");
  Float_t x1,y1;
  Int_t ncols1;
  //  Int_t nlines = 0;
  while (1) {
    ncols1 = fscanf(fp1,"%f %f",&x1, &y1);
    if (ncols1 <= 0) break;
    tpc1H->Fill(x1,y1);
    //cout << x1 << " " << endl;
  }
  

  Int_t n=0;
  Int_t entries = (Int_t) gsT->GetEntries();
 for (Int_t i =0; i< entries;i++){
   gsT->GetEntry(i);
   for (Int_t j =0; j< nT;j++){
     gsT->GetEntry(j);
     Double_t eq = eqT[j];
     Double_t omega = omegaT[j];
     Double_t tau = tT[j];
     Double_t Vol = VT[j];
     Double_t integral = TMath::Power(2*TMath::Pi(),-6)*TMath::Power(2,-1)*
     TMath::Pi()*Vol*tau*TMath::Power(200,-1)*eq*eq*alphaem*alphas*alphas*
     integralT[j];
     // Double_t integral = integralT[j];
     
     if (integral > 0){
     n++; 
      cout << integral << " " << eq  << " " << omega  << " " << tau << " " << Vol<< endl;
      if(eq == 0.3) yield03H->Fill(omega,2*integral);
       if(eq == 0.6) yield03H->Fill(omega,integral);
      }}}

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
  pt_1->SetLogy();

  yield03H->GetXaxis()->SetTitle("#omega_{q} [GeV]");
  yield03H->GetXaxis()->CenterTitle(true);
  yield03H->GetXaxis()->SetTitleSize(0.04);
  yield03H->GetXaxis()->SetLabelSize(0.03);
  yield03H->GetXaxis()->SetTitleOffset(1.4);
  yield03H->GetYaxis()->SetTitle("(1/2#pi#omega_{q}) dN/d#omega_{q} [GeV^{-2}]");
  yield03H->GetYaxis()->CenterTitle(true);
  yield03H->GetYaxis()->SetTitleSize(0.04);
  yield03H->GetYaxis()->SetLabelSize(0.03);
  yield03H->GetYaxis()->SetTitleOffset(1.5);


  //   yield03H = (TH1D*)yield03H->DrawNormalized("hist p",yield03H->Integral());
  yield03H->Scale(nT*n/10);//Deltat = 0.01 (500 pasos)*entradas/2 histogranas (1/3,1/6)
  yield03H->Draw("hist p");
  yield03H->GetYaxis()->SetRangeUser(0.00001,5000);
  //yield03H->GetXaxis()->SetRangeUser(0.0,3.5);
  yield03H->SetMarkerColor(1);
  yield03H->SetMarkerStyle(21);
  yield03H->SetMarkerSize(0.7);

    
  dataH->Add(tpc1H,-1);
  dataH->Draw("sames hist p ");
  dataH->SetMarkerColor(1);                                                                  
  dataH->SetMarkerStyle(25);
  dataH->SetMarkerSize(0.7);
  
  /*
  tpc1H->Draw("sames hist p ");
  tpc1H->SetMarkerColor(1);                                                                  
  tpc1H->SetMarkerStyle(24);
  tpc1H->SetMarkerSize(0.7);
  
  
  TLegend *leg2 = new TLegend(0.5,0.7,0.9,0.89);
  // leg2->SetTextFont(62);
  //leg2->SetTextSize(0.04);
  // leg2->SetLineColor(0);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(0);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry("dataH","PHENIX Au+Au 20-40%  - Direct","p");
  leg2->AddEntry("yield03H","UrQMD, Au-Au #sqrt{s_{NN}} = 200 GeV","p");
  leg2->AddEntry("","#Delta#tau = 0.0-0.5 fm, #beta = 0, 20-40 %","");
  leg2->Draw();
  // pt->SaveAs("fig5_2.pdf");
  //pt->SaveAs("fig5_2.eps");

  */  
    }
