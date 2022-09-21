TH1D* b0H = new TH1D("b0H","",35,0,3.5);
TH1D* b1H = new TH1D("b1H","",35,0,3.5);
TH1D* b2H = new TH1D("b2H","",35,0,3.5);
TH1D* b3H = new TH1D("b3H","",35,0,3.5);

TH1D* sumH = new TH1D("sumH","",35,0,3.5);
TH1D* dataH = new TH1D("dataH","",35,0,3.5);
TH1D* tpc1H = new TH1D("tpc1H","",35,0,3.5);
TH1D* tpc2H = new TH1D("tpc2H","",35,0,3.5);

double alphaem = TMath::Power(137,-1);
double alphas = TMath::Power(TMath::Pi(),-1);

const Int_t kMax=100000;
Int_t nT;
Double_t eqT[kMax];
Double_t omegaT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

Int_t n1T;
Double_t eq1T[kMax];
Double_t omega1T[kMax];
Double_t t1T[kMax];
Double_t V1T[kMax];
Double_t integral1T[kMax];
TTree* gs1T = new  TTree("gs1T","g's tree");

Int_t n2T;
Double_t eq2T[kMax];
Double_t omega2T[kMax];
Double_t t2T[kMax];
Double_t V2T[kMax];
Double_t integral2T[kMax];
TTree* gs2T = new  TTree("gs2T","g's tree");

Int_t n3T;
Double_t eq3T[kMax];
Double_t omega3T[kMax];
Double_t t3T[kMax];
Double_t V3T[kMax];
Double_t integral3T[kMax];
TTree* gs3T = new  TTree("gs3T","g's tree");
Int_t n,n1,n2,n3;
void Plot4(){
  char name[50];
  sprintf(name,"0.20.beta.0.root");
  cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);
  TTree *gsT = (TTree*) file->Get("gsT");
  gsT->SetBranchAddress("nT",&nT);
  gsT->SetBranchAddress("eqT",eqT);
  gsT->SetBranchAddress("tT",tT);
  gsT->SetBranchAddress("VT",VT);
  gsT->SetBranchAddress("omegaT",omegaT);
  gsT->SetBranchAddress("integralT",integralT);
  
  char name1[50];
  sprintf(name1,"20.40.beta.0.root");
  cout << " Openning file " << name1 << endl;
  TFile *file1 = new TFile(name1);
  TTree *gs1T = (TTree*) file1->Get("gsT");
  gs1T->SetBranchAddress("nT",&n1T);
  gs1T->SetBranchAddress("eqT",eq1T);
  gs1T->SetBranchAddress("tT",t1T);
  gs1T->SetBranchAddress("VT",V1T);
  gs1T->SetBranchAddress("omegaT",omega1T);
  gs1T->SetBranchAddress("integralT",integral1T);
  
  char name2[50];
  sprintf(name2,"40.60.beta.0.root");
  cout << " Openning file " << name2 << endl;
  TFile *file2 = new TFile(name2);
  TTree *gs2T = (TTree*) file2->Get("gsT");
  gs2T->SetBranchAddress("nT",&n2T);
  gs2T->SetBranchAddress("eqT",eq2T);
  gs2T->SetBranchAddress("tT",t2T);
  gs2T->SetBranchAddress("VT",V2T);
  gs2T->SetBranchAddress("omegaT",omega2T);
  gs2T->SetBranchAddress("integralT",integral2T);

  char name3[50];
  sprintf(name3,"0.40.Cu.beta.0.root");
  cout << " Openning file " << name3 << endl;
  TFile *file3 = new TFile(name3);
  TTree *gs3T = (TTree*) file3->Get("gsT");
  gs3T->SetBranchAddress("nT",&n3T);
  gs3T->SetBranchAddress("eqT",eq3T);
  gs3T->SetBranchAddress("tT",t3T);
  gs3T->SetBranchAddress("VT",V3T);
  gs3T->SetBranchAddress("omegaT",omega3T);
  gs3T->SetBranchAddress("integralT",integral3T);
  
  
  FILE *fp = fopen("phenix20.40.txt","r");
  FILE *fp1 = fopen("tpnc20.40.txt","r");
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
    for (Int_t j =0; j< nT;j++){
    gsT->GetEntry(j);
    Double_t eq = eqT[j];
    Double_t omega = omegaT[j];
    Double_t tau = tT[j];
    Double_t Vol = VT[j];
    Double_t integral = TMath::Power(TMath::Pi(),-6)*TMath::Power(2,-1)*
      TMath::Pi()*Vol*tau*TMath::Power(200,-1)*eq*eq*alphaem*alphas*alphas*
      integralT[j];
    if(integral > -1000){
      n++;
      if(eq == 0.3) b0H->Fill(omega,2*integral);
      if(eq == 0.6) b0H->Fill(omega,integral);
    }}}
  
  Int_t entries1 = (Int_t) gs1T->GetEntries();
  for (Int_t i =0; i<entries1;i++){
    gs1T->GetEntry(i);
    for (Int_t j =0; j< n1T;j++){
      gs1T->GetEntry(j);
      Double_t eq1 = eq1T[j];
      Double_t omega1 = omega1T[j];
      Double_t tau1 = t1T[j];
      Double_t Vol1 = V1T[j];
      Double_t integral1 = TMath::Power(TMath::Pi(),-6)*TMath::Power(2,-1)*
	TMath::Pi()*Vol1*tau1*TMath::Power(200,-1)*eq1*eq1*alphaem*alphas*alphas*
	integral1T[j];
      if(integral1 > -1000){
	n1++;
	if(eq1 == 0.3) b1H->Fill(omega1,2*integral1);
	if(eq1 == 0.6) b1H->Fill(omega1,integral1);
      }}}
  
  Int_t entries2 = (Int_t) gs2T->GetEntries();
  for (Int_t i =0; i< entries2;i++){
    gs2T->GetEntry(i);
    for (Int_t j =0; j< n2T;j++){
      gs2T->GetEntry(j);
      Double_t eq2 = eq2T[j];
      Double_t omega2 = omega2T[j];
      Double_t tau2 = t2T[j];
      Double_t Vol2 = V2T[j];
      Double_t integral2 = TMath::Power(TMath::Pi(),-6)*TMath::Power(2,-1)*
	TMath::Pi()*Vol2*tau2*TMath::Power(200,-1)*eq2*eq2*alphaem*alphas*alphas*
	integral2T[j];
      if(integral2 > -1000){
	n2++;
	if(eq2 == 0.3) b2H->Fill(omega2,2*integral2);
	if(eq2 == 0.6) b2H->Fill(omega2,integral2);
      }}}

  Int_t entries3 = (Int_t) gs3T->GetEntries();
  for (Int_t i =0; i< entries3;i++){
    gs3T->GetEntry(i);
    for (Int_t j =0; j< n3T;j++){
      gs3T->GetEntry(j);
      Double_t eq3 = eq3T[j];
      Double_t omega3 = omega3T[j];
      Double_t tau3 = t3T[j];
      Double_t Vol3 = V3T[j];
      Double_t integral3 = TMath::Power(TMath::Pi(),-6)*TMath::Power(2,-1)*
	TMath::Pi()*Vol3*tau3*TMath::Power(200,-1)*eq3*eq3*alphaem*alphas*alphas*
	integral3T[j];
      if(integral3 > 0){
	n3++;
	if(eq3 == 0.3) b3H->Fill(omega3,2*integral3);
	if(eq3 == 0.6) b3H->Fill(omega3,integral3);
      }}}

  cout << n << " " <<  n1 << " " << n2 << " " << n3 << endl;    
  
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
  b0H->GetYaxis()->SetTitle("(1/2#pi#omega_{q}) dN_{mag}/d#omega_{q} [GeV^{-2}]");
  b0H->GetYaxis()->CenterTitle(true);
  b0H->GetYaxis()->SetTitleSize(0.04);
  b0H->GetYaxis()->SetLabelSize(0.03);
  b0H->GetYaxis()->SetTitleOffset(1.5);

  b0H->Scale(10*n/2);
  b0H->Draw("hist p");
  //  b0H = (TH1D*)b0H->DrawNormalized("hist p",b0H->Integral());
  b0H->GetYaxis()->SetRangeUser(0.00001,5000);
  //b0H->GetXaxis()->SetRangeUser(0.0,3.5);
  b0H->SetMarkerColor(1);
  b0H->SetMarkerStyle(21);
  b0H->SetMarkerSize(0.7);
  
  b1H->Scale(10*n1/2);
  b1H->Draw("sames hist p");
  //b1H = (TH1D*)b1H->DrawNormalized("sames hist p",b1H->Integral());
  //b1H->Draw("sames histo p");
  b1H->SetMarkerColor(1);
  b1H->SetMarkerStyle(20);
  b1H->SetMarkerSize(0.7);

  b2H->Scale(10*n2/2);
  b2H->Draw("sames hist p");
  //b2H = (TH1D*)b2H->DrawNormalized("sames hist p",b2H->Integral());
  b2H->SetMarkerColor(1);
  b2H->SetMarkerStyle(22);
  b2H->SetMarkerSize(0.7);

  b3H->Scale(10*n3/2);
  b3H->Draw("sames hist p");
  //  b3H = (TH1D*)b3H->DrawNormalized("sames hist p",b3H->Integral());
  b3H->SetMarkerColor(1);
  b3H->SetMarkerStyle(3);
  b3H->SetMarkerSize(0.7);

  /*dataH->Add(tpc1H,-1);
    dataH->Draw("sames hist p ");
  dataH->SetMarkerColor(1);
  dataH->SetMarkerStyle(24);
  dataH->SetMarkerSize(0.7);
  */

  TLegend *leg2 = new TLegend(0.5,0.6,0.9,0.89);
  // leg2->SetTextFont(62);
  //leg2->SetTextSize(0.04);
  // leg2->SetLineColor(0);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(0);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry("","#Delta#tau = 0.0-0.5 fm, #beta = 0","");
  leg2->AddEntry("","UrQMD, Au-Au #sqrt{s_{NN}} = 200 GeV","");
  leg2->AddEntry("b0H","0-20 %","p");
  leg2->AddEntry("b1H","20-40 %","p");
  leg2->AddEntry("b2H","40-60 %","p");
  leg2->AddEntry("","UrQMD, Cu-Cu #sqrt{s_{NN}} = 200 GeV","");
  leg2->AddEntry("b3H","0-40 %","p");
  leg2->AddEntry("","","");

  leg2->Draw();
  pt->SaveAs("fig3_1.pdf");
  pt->SaveAs("fig3_1.eps");

    
    }
