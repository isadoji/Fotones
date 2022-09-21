//-----------------------------------------------------------
// Description:
//
//       Plot macro. This macro is to plot experimental data for yield photons.
//
// Environment:
//      ROOT
//
// Author List:
//       Isabel Domínguez Jiménez
//       isadoji@uas.edu.mx
//-----------------------------------------------------------
TH1D* dataH = new TH1D("dataH","",40,0,4);
TH1D* promptH = new TH1D("promptH","",40,0,4);
TH1D* thermalH = new TH1D("thermalH","",40,0,4);
TH1D* nococktailH = new TH1D("nococktailH","",40,0,4);
TH1D* tpncH = new TH1D("tpncH","",40,0,4);
void graphP(){


  FILE *fp = fopen("../dataexp/phenix.txt","r");
  FILE *fpp = fopen("../dataexp/prompt.txt","r");
  FILE *fpt = fopen("../dataexp/thermal.txt","r");
  FILE *fnc = fopen("../dataexp/non_cocktail.txt","r");
  FILE *ftpnc = fopen("../dataexp/tpnc.txt","r");
  Float_t x,y,ex,ey;
  Float_t x1,y1,x2,y2,x3,y3,x4,y4;

  Int_t ncols, ncols1, ncols2, ncols3, ncols4;
  Int_t nlines = 0;


  while (1) {
    ncols = fscanf(fp,"%f %f %f %f",&x, &y, &ex, &ey);
    ncols1 = fscanf(fpp,"%f %f",&x1, &y1);
    ncols2 = fscanf(fpt,"%f %f",&x2, &y2);
    ncols3 = fscanf(fnc,"%f %f",&x3, &y3);
    ncols4 = fscanf(ftpnc,"%f %f",&x4, &y4);

    if (ncols < 0) break;
    if (ncols1 < 0) break;
    dataH->Fill(x,y);
    promptH->Fill(x1,y1);
    thermalH->Fill(x2,y2);
    nococktailH->Fill(x3,y3);
    tpncH->Fill(x4,y4);
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
  pt_1->SetLogy();

  dataH->Draw("hist p ");
  dataH->SetMarkerColor(1);                                                                    dataH->SetMarkerStyle(20);
  promptH->Draw("same hist p");
  promptH->SetMarkerColor(2);                                                                   promptH->SetMarkerStyle(20);
  thermalH->Draw("same hist p");
  thermalH->SetMarkerColor(3);                                                                   thermalH->SetMarkerStyle(20);
  nococktailH->Draw("same hist p");
  nococktailH->SetMarkerColor(4);                                                              nococktailH->SetMarkerStyle(20);
  tpncH->Draw("same hist p");
  tpncH->SetMarkerColor(6);
  tpncH->SetMarkerStyle(20);

  TLegend *leg1 = new TLegend(0.6,0.7,0.89,0.89);
  leg1->AddEntry("dataH","PHENIX (20-40 %)","p");
  leg1->AddEntry("promptH","PROMPT","p");
  leg1->AddEntry("thermalH","THERMAL","p");
  leg1->AddEntry("nococktailH","NON-COCKTAIL","p");
  leg1->AddEntry("tpncH","PROMPT+THERMAL+NON-COCKTAIL","p");
  leg1->Draw();

  /*

  pt->cd(2);

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



  //  dataH->Add(thermalH,1);
  dataH->Draw("hist p ");
  dataH->SetMarkerColor(1);
  dataH->SetMarkerStyle(20);

 TLegend *leg2 = new TLegend(0.6,0.7,0.89,0.89);
  leg2->AddEntry("dataH","PHENIX (20-40 %)","p");
    leg2->Draw();
  */
}
