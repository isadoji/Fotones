TH1F* yieldH = new TH1F("yieldH","",35,0,3.5);
TH1F* v2H = new TH1F("v2H","",35,0,3.5);
TH1F* yieldv2H = new TH1F("yieldv2H","",35,0,3.5);
TH1D* dataH = new TH1D("dataH","",35,0,3.5);
TH1D* tpc1H = new TH1D("tpc1H","",35,0,3.5);
TH1F* v2totH = new TH1F("v2totH","",35,0,3.5);

const Int_t kMax=1000000;
Int_t nT;
Double_t eqT[kMax];
Double_t tT[kMax];
Double_t VT[kMax];
Double_t omegaT[kMax];
Double_t integralT[kMax];
Double_t qintT[kMax];
Double_t integralv2T[kMax];
TTree* gsT = new  TTree("gsT","g's tree");

double alphaem = TMath::Power(137,-1);
double alphas = TMath::Power(TMath::Pi(),-1);


void Plotyieldv2(){
  char name[50];
  sprintf(name,"yieldv2.root");
  cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);
  TTree *gsT = (TTree*) file->Get("gsT");
  gsT->SetBranchAddress("nT",&nT);
  gsT->SetBranchAddress("eqT",eqT);
  gsT->SetBranchAddress("tT",tT);
  gsT->SetBranchAddress("VT",VT);
  gsT->SetBranchAddress("omegaT",omegaT);
  gsT->SetBranchAddress("integralT",integralT);
  gsT->SetBranchAddress("qintT",qintT);
  gsT->SetBranchAddress("integralv2T",integralv2T);

  FILE *fp = fopen("dataexp/yield20.40.txt","r");
  Int_t ncols;
  Float_t x,y;
  FILE *fp1 = fopen("dataexp/v220.40.txt","r");
  Float_t x1,y1,y3,ex1,ey1;

  while (1) {
    ncols = fscanf(fp,"%f %f",&x, &y);
    ncols = fscanf(fp1,"%f %f %f %f %f",&x1, &y1,&ex1, &y3,&ey1);
    if (ncols <= 0) break;
    dataH->Fill(x,y);
    tpc1H->Fill(x1,y*y1);
  }
  Int_t n=0;
  Int_t entries = (Int_t) gsT->GetEntries();
  for (Int_t i =0; i< entries;i++){
    gsT->GetEntry(i);
    for (Int_t j =0; j< nT;j++){
      gsT->GetEntry(j);
      Double_t eq = eqT[j];
      Double_t tau = tT[j];
      Double_t Vol = VT[j];
      Double_t omega = omegaT[j];
      Double_t qint = qintT[j];
      //if(tau >0 && tau < 0.6) {
	Double_t cte =  TMath::Power(2*TMath::Pi(),-5)*TMath::Power(2,-1)*
	  TMath::Pi()*Vol*tau*TMath::Power(200,-1)*eq*eq*alphaem*alphas*alphas;
	Double_t integral = cte*integralT[j]*omega;
	Double_t integralv2 = cte*integralv2T[j];
	Double_t yieldv2 = integral*integralv2;
	if (integral > 0  &&integralv2 > 0){
	  n++;
	  if(eq == 0.3){
	    yieldH->Fill(omega,2*integral);
	    v2H->Fill(omega,2*integralv2);
	  }
	  if(eq == 0.6){
	    yieldH->Fill(omega,integral);
	    v2H->Fill(omega,integralv2);
	    //  }
	}
      }
    }
  }
  TCanvas *pt = new TCanvas("pt", "pt particles distribution",10,10,700,700);
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
  //  pt_1->SetGridy();
  //pt_1->SetGridx();

  v2H->GetXaxis()->SetTitle("#omega_{q}");
  v2H->GetXaxis()->CenterTitle(true);
  v2H->GetXaxis()->SetTitleSize(0.04);
  v2H->GetXaxis()->SetLabelSize(0.03);
  v2H->GetXaxis()->SetTitleOffset(1.4);
  v2H->GetYaxis()->SetTitle("");
  v2H->GetYaxis()->CenterTitle(true);
  v2H->GetYaxis()->SetTitleSize(0.04);
  v2H->GetYaxis()->SetLabelSize(0.03);
  v2H->GetYaxis()->SetTitleOffset(2);

  v2H->Scale(n/10);
  v2H->Draw("hist p");
  v2H->GetYaxis()->SetRangeUser(1e-10,1e5);
  //v2H->GetXaxis()->SetRangeUser(0.0,3.5);
  v2H->SetMarkerColor(1);
  v2H->SetMarkerStyle(20);
  v2H->SetMarkerSize(0.7);

  yieldH->Scale(n/10);
  yieldH->Draw("sames hist p");
  yieldH->SetMarkerColor(2);
  yieldH->SetMarkerStyle(20);
  yieldH->SetMarkerSize(0.7);

  /*
  for (int i = 0; i<=v2H->GetXaxis()->GetNbins(); i++){
    Double_t temp =  v2H->GetBinContent(i)*yieldH->GetBinContent(i);
    yieldv2H->Fill(0.1*yieldH->GetBin(i)-0.05, temp);
    if (tpc1H->GetBinContent(i) != 0){
      v2totH->Fill(0.1*yieldH->GetBin(i)-0.05, temp+tpc1H->GetBinContent(i));
      cout << i << " " << temp << " " << tpc1H->GetBinContent(i) << endl;
      //temp+tpc1H->GetBinContent(i) << endl ;
    }
  }

  Double_t normv2 = 1/(yieldH->Integral());
  Double_t normdata = 1/(dataH->Integral());
  Double_t normtot = 1/(yieldH->Integral()+dataH->Integral());


  //yieldv2H->Add(tpc1H,1);

  yieldv2H->Scale(normv2);
  yieldv2H->Draw("sames hist p");
  yieldv2H->SetMarkerColor(4);
  yieldv2H->SetMarkerStyle(20);
  yieldv2H->SetMarkerSize(0.7);

  //dataH->Add(tpc1H,-1);
  dataH->Draw("sames hist p ");
  dataH->SetMarkerColor(1);
  dataH->SetMarkerStyle(25);
  dataH->SetMarkerSize(0.7);


   tpc1H->Scale(normdata);
  tpc1H->Draw("sames hist p ");
  tpc1H->SetMarkerColor(1);
  tpc1H->SetMarkerStyle(25);
  tpc1H->SetMarkerSize(0.7);

  //v2totH->Add(tpc1H,1);
  v2totH->Scale(normtot);
  v2totH->Draw("sames hist p");
  v2totH->SetMarkerColor(6);
  v2totH->SetMarkerStyle(20);
  v2totH->SetMarkerSize(0.7);
  */
  TLegend *leg2 = new TLegend(0.4,0.7,0.8,0.89);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.025);
  // leg2->SetLineColor(0);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(0);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry("","UrQMD, Au-Au #sqrt{s_{NN}} = 200 GeV, 20-40 %","");
  leg2->AddEntry("v2H","v_{2}","p");
  leg2->AddEntry("yieldH","Yield","p");
  leg2->AddEntry("yieldv2H","Yield*v_{2}/Yield_{tot}","p");

  leg2->Draw();
  //  pt->SaveAs("V2.200GeV.old.pdf");
  //pt->SaveAs("V2.200GeV.old.eps");


    }
