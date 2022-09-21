TH1D* b1H = new TH1D("b1H","",40,0,4);
TH1D* b2H = new TH1D("b2H","",40,0,4);
TH1D* b3H = new TH1D("b3H","",40,0,4);
TH1D* b4H = new TH1D("b4H","",40,0,4);
TH1D* b5H = new TH1D("b5H","",40,0,4);
TH1D* b6H = new TH1D("b6H","",40,0,4);
TH1D* b7H = new TH1D("b7H","",40,0,4);
TH1D* b8H = new TH1D("b8H","",40,0,4);
TH1D* b9H = new TH1D("b9H","",40,0,4);

TH1D* sumH = new TH1D("sumH","",40,0,4);
TH1D* dataH = new TH1D("dataH","",40,0,4);
TH1D* dataCuH = new TH1D("dataCuH","",40,0,4);
TH1D* tpc1H = new TH1D("tpc1H","",40,0,4);
TH1D* tpc2H = new TH1D("tpc2H","",40,0,4);
const Int_t kMax=100000;
Int_t nT;
Double_t eqT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

double alphaem = TMath::Power(137,-1);
double alphas = TMath::Power(TMath::Pi(),-1);
float mq = 0.0001;
int g = 2;

void Plotv3(){
  char name[50];
  sprintf(name,"20.40.Au.Au.beta0.15.J.root");
  //sprintf(name,"integral20-40.t0.5.root");
  //cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);

  TTree *gsT = (TTree*) file->Get("gsT");
  gsT->SetBranchAddress("nT",&nT);
  gsT->SetBranchAddress("eqT",eqT);
  gsT->SetBranchAddress("VT",VT);
  gsT->SetBranchAddress("tT",tT);
  gsT->SetBranchAddress("omegaT",omegaT);
  gsT->SetBranchAddress("integralT",integralT);


  FILE *fp = fopen("phenix.txt","r");
  FILE *fp0 = fopen("phenixCu.txt","r");
  FILE *fp1 = fopen("tpnc1.txt","r");
  FILE *fp2 = fopen("tpnc2.txt","r");
  Float_t x,y,ex,ey;
  Float_t x0,y0,ex0,ey0;
  Float_t x1,y1;
  Float_t x2,y2;

  Int_t ncols, ncols0, ncols1, ncols2;
  Int_t nlines = 0;


  while (1) {
    ncols = fscanf(fp,"%f %f %f %f",&x, &y, &ex, &ey);
    ncols1 = fscanf(fp1,"%f %f",&x1, &y1);
    if (ncols < 0) break;
    if (ncols1 < 0) break;
    dataH->Fill(x,y);
    tpc1H->Fill(x1,y1);
  }
  while (1) {
    ncols0 = fscanf(fp0,"%f %f %f %f",&x0, &y0, &ex0, &ey0);
    ncols2 = fscanf(fp2,"%f %f",&x2, &y2);
    if (ncols0 < 0) break;
    if (ncols2 < 0) break;
    dataCuH->Fill(x0,y0);
    tpc2H->Fill(x2,y2);

  }

    Int_t entries = (Int_t) gsT->GetEntries();
  //cout << "entries=  " << entries << endl;

  for (Int_t i =0; i< entries;i++){
    gsT->GetEntry(i);
    for (Int_t j =0; j< nT;j++){
      gsT->GetEntry(j);
      Double_t eq = eqT[j];
      Float_t t = tT[j];
      Double_t omega = omegaT[j];
      //Double_t integral = integralT[j];
      Double_t Vol = VT[j];
      
      cout << eq <<  " " << Vol << endl;
      
     Double_t integral = TMath::Pi()*TMath::Power(TMath::Pi(),-6)*TMath::Power(2,-1)*
       Vol*t*TMath::Power(200,-1)*eq*eq*alphaem*alphas*alphas*integralT[j];


    if(t>0.0 &&  t<0.05 && eq == 0.3) {
      b1H->Fill(omega,2*integral);
    }
    if(t>0.0 &&  t<0.05 && eq == 0.6) {
      b1H->Fill(omega,integral);
    }
    if(t>0.1 &&  t<0.15 && eq == 0.3) {
      b2H->Fill(omega,2*integral);
    }
    if(t>0.1 &&  t<0.15 && eq == 0.6) {
      b2H->Fill(omega,integral);
    }
    if(t>0.15 &&  t<0.2 && eq == 0.3) {
      b3H->Fill(omega,2*integral);
    }
    if(t>0.15 &&  t<0.2 && eq == 0.6) {
      b3H->Fill(omega,integral);
    }
    if(t>0.2 &&  t<0.25 && eq == 0.3) {
      b4H->Fill(omega,2*integral);
    }
    if(t>0.2 &&  t<0.25 && eq == 0.6) {
      b4H->Fill(omega,integral);
    }
    if(t>0.25 &&  t<0.3 && eq == 0.3) {
      b5H->Fill(omega,2*integral);
    }
    if(t>0.25 &&  t<0.3 && eq == 0.6) {
      b5H->Fill(omega,integral);
    }
    if(t>0.3 &&  t<0.35 && eq == 0.3) {
      b6H->Fill(omega,2*integral);
    }
    if(t>0.3 &&  t<0.35 && eq == 0.6) {
      b6H->Fill(omega,integral);
    }
    if(t>0.35 &&  t<0.4 && eq == 0.3) {
      b7H->Fill(omega,2*integral);
    }
    if(t>0.35 &&  t<0.4 && eq == 0.6) {
      b7H->Fill(omega,integral);
    }
    if(t>0.4 &&  t<0.45 && eq == 0.3) {
      b8H->Fill(omega,2*integral);
    }
    if(t>0.4 &&  t<0.45 && eq == 0.6) {
      b8H->Fill(omega,integral);
    }
    if(t>0.45 &&  t<0.5 && eq == 0.3) {
      b9H->Fill(omega,2*integral);
    }
    if(t>0.45 &&  t<0.5 && eq == 0.6) {
      b9H->Fill(omega,integral);
    }
                                                
  }
  }
  /*TCanvas *pt = new TCanvas("pt", "pt particles distribution",50,50,500,500);
  gStyle->SetOptStat(false);
  //gStyle->SetOptStat("me");
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

  b1H->GetXaxis()->SetTitle("#omega_{q} [GeV]");
  b1H->GetXaxis()->CenterTitle(true);
  b1H->GetXaxis()->SetTitleSize(0.04);
  b1H->GetXaxis()->SetLabelSize(0.03);
  b1H->GetXaxis()->SetTitleOffset(1.4);
  b1H->GetYaxis()->SetTitle("(1/2#pi#omega_{q}) dN/d#omega_{q} [GeV^{-2}]");
  b1H->GetYaxis()->CenterTitle(true);
  b1H->GetYaxis()->SetTitleSize(0.04);
  b1H->GetYaxis()->SetLabelSize(0.03);
  b1H->GetYaxis()->SetTitleOffset(2);
  b1H->GetYaxis()->SetRangeUser(0.1,1000);

  b1H->Draw("hist p");
  b1H->SetMarkerColor(1);
  b1H->SetMarkerStyle(20);
  b1H->SetMarkerSize(0.5);
  b2H->Draw("same hist p");
  b2H->SetMarkerColor(2);
  b2H->SetMarkerStyle(21);
  b3H->Draw("same hist p");
  b3H->SetMarkerColor(3);
  b3H->SetMarkerStyle(21);
  b4H->Draw("same hist p");
  b4H->SetMarkerColor(4);
  b4H->SetMarkerStyle(21);
  b5H->Draw("same hist p");
  b5H->SetMarkerColor(1);
  b5H->SetMarkerStyle(22);
  
TLegend *leg = new TLegend(0.7,0.5,0.89,0.89);
 leg->AddEntry("","b = 20 - 40","");
 leg->AddEntry("b1H","#tau = 0.01","p");
 leg->AddEntry("b2H","#tau = 0.02","p");
 leg->AddEntry("b3H","#tau = 0.03","p");
 leg->AddEntry("b4H","#tau = 0.04","p");
 leg->AddEntry("b5H","#tau = 0.05","p");
 leg->Draw();
*/
 
 sumH->Add(b1H);
 sumH->Add(b2H);
 sumH->Add(b3H);
 sumH->Add(b4H);
 sumH->Add(b5H);
 sumH->Add(b6H);
 sumH->Add(b7H);
 sumH->Add(b8H);
 sumH->Add(b9H);
 //sumH->Scale(10);//add 10 histos
 
 
 TCanvas *pt2 = new TCanvas("pt", "pt particles distribution",100,100,900,900);
 gStyle->SetOptStat(false);
 //gStyle->SetOptStat("me");
 gStyle->SetOptTitle(0);
 gStyle->SetPalette(1);
 pt2->SetRightMargin(0.0465116);
 pt2->SetTopMargin(0.39916);
 pt2->SetFillColor(0);
 TPad *pt_2 = new TPad("pt_2", "pt_2",0.0,0.0,0.98,0.98);
 pt_2->Draw();
 pt_2->cd();
 pt_2->Range(-20.6795,-0.0133333,105.496,0.062549);
 pt_2->SetFillColor(0);
 pt_2->SetBorderSize(2);
 pt_2->SetLeftMargin(0.17);
 pt_2->SetRightMargin(0.05);
 pt_2->SetTopMargin(0.1);
 pt_2->SetBottomMargin(0.15);
 pt_2->SetFrameLineWidth(2);
 pt_2->SetFrameLineWidth(2);
 pt_2->SetLogy();
 sumH->GetXaxis()->SetTitle("#omega_{q} [GeV]");
 sumH->GetXaxis()->CenterTitle(true);
 sumH->GetXaxis()->SetTitleSize(0.04);
 sumH->GetXaxis()->SetLabelSize(0.03);
 sumH->GetXaxis()->SetTitleOffset(1.4);
 sumH->GetYaxis()->SetTitle("(1/2#pi#omega_{q}) dN/d#omega_{q} [GeV^{-2}]");
 sumH->GetYaxis()->CenterTitle(true);
 sumH->GetYaxis()->SetTitleSize(0.04);
 sumH->GetYaxis()->SetLabelSize(0.03);
 sumH->GetYaxis()->SetTitleOffset(2);
 sumH->GetYaxis()->SetRangeUser(0.000001,20);
 sumH->GetXaxis()->SetRangeUser(0.0,3.5);
 sumH = (TH1D*)sumH->DrawNormalized("hist p",10*sumH->Integral());
 //sumH->Draw("hist p");
 sumH->SetMarkerColor(1);
 sumH->SetMarkerStyle(20);
 sumH->SetMarkerSize(0.7);
 dataH->Add(tpc1H,-1);
 dataH->Draw("sames hist p ");
 dataH->SetMarkerColor(1);
 dataH->SetMarkerStyle(24);
 dataH->SetMarkerSize(0.7);
 /* dataCuH->Add(tpc1H,-1);
    dataCuH->Draw("sames hist p ");
    dataCuH->SetMarkerColor(1);
    dataCuH->SetMarkerStyle(24);
    dataCuH->SetMarkerSize(0.7);
 */
 
  TLegend *leg2 = new TLegend(0.6,0.7,0.89,0.88);
  // leg2->SetTextFont(62);
  //leg2->SetTextSize(0.04);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(0);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry("dataH","PHENIX Au+Au 20-40%  - Direct","p");
  leg2->AddEntry("sumH","UrQMD, Au-Au #sqrt{s_{NN}}=200 GeV","p");
  leg2->AddEntry("","#Delta#tau = 0.0-0.5 fm, #beta = 0.15","");
  /*leg2->AddEntry("dataH","PHENIX Cu+Cu 0-40%  - Background","p");
  leg2->AddEntry("sumH","UrQMD, Cu-Cu #sqrt{s_{NN}}=200 GeV","p");
  leg2->AddEntry("","#Delta#tau = 0.02-0.5 fm","");
  */
  leg2->Draw();
  
  pt2->SaveAs("Yield.Au.Au.20-40.beta0.15.eps");
  pt2->SaveAs("Yield.Au.Au.20-40.beta0.15.pdf");

}
