int tI = 0.0;
const float tSteps = 0.1;
const int tTotal = 50;
const int protonNumber = 79;
double_t aTime[tTotal] = {0.};
double eB[tTotal] = {};

Int_t nbin = 1000;// #events
Int_t xmin = 0;// #events
Double_t xmax = 5;// #events
Float_t norm = nbin/(TMath::Abs(xmin)+TMath::Abs(xmax));
TH1F *eBH = new TH1F("eBH","", nbin,xmin,xmax);
  
  void Bfield(){
  for(int k=10; k<11; k++){
    gROOT->Reset();
    TChain mychain("T");
    mychain.Add(TString::Format("prueba%d.root",k));
    struct particula_t 
    {
      Float_t time,X,Y,Z,E,Px,Py,Pz,Pt,P,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx,frezePy,frezePz,b,nspec,R,PXR,eBx,eBy,eBz,eByWave,eByPoint,dCb;
    } PARTICLE;
    particula_t  particle;
    Int_t nevent = mychain.GetEntries();
    cout << nevent << endl;
    for(Int_t j=0; j<nevent; j++){
      mychain.GetEvent(j);
        Float_t R = particle.R;
      if((particle.lastcoll==0 && particle.charge==1 && R>0 && particle.id==1) || (particle.lastcoll!=0 && particle.charge==1 && R>0 && particle.id==1)){
        for( int iTime = 0; iTime<=tTotal; ++iTime){
          float_t iTC = iTime*tSteps;
          Float_t Delta = std::abs(iTC-particle.time);
          if(Delta < 0.01 && TMath::Finite(particle.eBy)==1 && particle.eBy!=0){
           cout << Delta<< endl;
            eB[iTime] += sqrt(particle.eBx*particle.eBx+particle.eBy*particle.eBy+particle.eBz*particle.eBz);
          }
        }
      }
     
     // eB[j] += sqrt(particle.eBx*particle.eBx+particle.eBy*particle.eBy+particle.eBz*particle.eBz);
  }
  TCanvas* c1 = new TCanvas("c1","UrQMD test example",800,800);
  gStyle->SetOptStat(false);
  c1->SetRightMargin(0.0465116);
  c1->SetTopMargin(0.1);
  c1->SetFillColor(0);
  c1->SetLogy();
  eBH->Scale();
  eBH->Draw("hist p");
  eBH->SetMarkerColor(1);
  eBH->SetMarkerStyle(20);
  eBH->SetXTitle("#Tau [fm]");
  eBH->GetXaxis()->CenterTitle(true);
  eBH->GetXaxis()->SetTitleSize(0.04);
  eBH->GetXaxis()->SetLabelSize(0.03);
  eBH->GetXaxis()->SetTitleOffset(1.2);
  eBH->SetYTitle("e|B|/m_{#pi}^{2}");
  eBH->GetYaxis()->CenterTitle(true);
  eBH->GetYaxis()->SetTitleSize(0.04);
  eBH->GetYaxis()->SetLabelSize(0.03);
  eBH->GetYaxis()->SetTitleOffset(1.2);
  TLegend *leg = new TLegend(0.48,0.7,0.89,0.89);
  leg->SetTextFont(62);
  //leg2->SetTextSize(0.04);                                         
  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry("","UrQMD, Au-Au #sqrt{s_{NN}} = 11 GeV","");
  leg->AddEntry("","Spectators","");
  leg->AddEntry(eBH,"0-20 %","p");
  leg->Draw();
  //c1->SaveAs("b20_entre_sigma.eps");
}    

  }