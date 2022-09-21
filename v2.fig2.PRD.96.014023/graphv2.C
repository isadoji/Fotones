//-----------------------------------------------------------
// Description:
//
//       Plot macro. This macro is to plot the original Fig2.b of PhysRevD.96.014023.pdf numerical calculation.
//
// Environment:
//      ROOT
//
// Author List:
//       Isabel Domínguez Jiménez
//       isadoji@uas.edu.mx
//-----------------------------------------------------------
TH1D* BH = new TH1D("BH","",350,0,3.5);
TH1D* B1H = new TH1D("B1H","",350,0,3.5);

void graphv2(){
  FILE *fp = fopen("v2.03.txt","r");
  FILE *fp1 = fopen("v2.06.txt","r");
  TGraph *plot = new TGraph("v2.03.txt");ls
  plot->Draw("AL*");

  Float_t omega03, integral03;
  Float_t omega06, integral06;

  Int_t ncols, ncols1;
  Int_t nlines = 0;

  while (1) {
    ncols = fscanf(fp,"%f %f ",&omega03,&integral03);
    ncols1 = fscanf(fp,"%f %f ",&omega06,&integral06);

    if (ncols < 0) break;
    if (ncols1 < 0) break;

    BH->Fill(omega03,integral03);
    B1H->Fill(omega06,integral06);


  }

  TCanvas *pt = new TCanvas("pt", "pt particles distribution",50,50,500,500);
  gStyle->SetOptStat(false);
  //gStyle->SetOptStat("me");

  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
  pt->SetRightMargin(0.0465116);
  pt->SetTopMargin(0.39916);
  pt->SetFillColor(0);

  //  pt->Divide(2,1);
  //pt->cd(1);

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

  BH->Draw("hist  p");
  BH->SetMarkerColor(1);
  BH->SetMarkerStyle(20);
  BH->SetMarkerSize(0.7);
  BH->GetXaxis()->SetTitle("#tau [fm]");
  BH->GetXaxis()->CenterTitle(true);
  BH->GetXaxis()->SetTitleSize(0.04);
  BH->GetXaxis()->SetLabelSize(0.03);
  BH->GetXaxis()->SetTitleOffset(1.0);
  BH->GetYaxis()->SetTitle(" eB/m_{#pi}^{2}");
  BH->GetYaxis()->CenterTitle(true);
  BH->GetYaxis()->SetTitleSize(0.04);
  BH->GetYaxis()->SetLabelSize(0.03);
  BH->GetYaxis()->SetTitleOffset(1.0);
  //  BH->GetYaxis()->SetRangeUser(0.005,10);
  B1H->Draw("same hist p ");
  B1H->SetMarkerColor(1);
  B1H->SetMarkerStyle(21);
  B1H->SetMarkerSize(0.7);

  TLegend *leg1 = new TLegend(0.5,0.7,0.9,0.89);
  // leg1->SetTextFont(62);
  //leg1->SetTextSize(0.04);
  leg1->SetLineColor(0);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(1001);

  leg1->AddEntry("","UrQMD, Au-Au #sqrt{s_{NN}}=200 GeV","");
  leg1->AddEntry("BH","0-20 %","p");
  leg1->AddEntry("B1H","20-40 %","p");
  leg1->AddEntry("B2H","40-60 %","p");
  leg1->AddEntry("","UrQMD, Cu-Cu #sqrt{s_{NN}}=200 GeV","");
  leg1->AddEntry("B3H","0-40 %","p");
  leg1->Draw();

  //  pt->SaveAs("B.eps");
  //pt->SaveAs("B.pdf");


}
