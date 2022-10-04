/*void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);

  */

void StylePlots(){
  gStyle->SetOptFit(11111);
  gStyle->Reset("Plain");
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10); //
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.18); //0.12
  gStyle->SetPadLeftMargin(0.18); // 0.12
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);

  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.06,"xyz"); //0.5
  //  gStyle->SetTitleOffset(1.5,"y");  //0.95
  //  gStyle->SetTitleOffset(1.3,"x"); //0-95

 gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);

  gStyle->SetLegendBorderSize(0); //
  
  gStyle->SetLegendFillColor(0);
    gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
  }

 int tI = 0.0;
const int tF = 50.0;
const float tSteps = 0.1;
const int tTotal = 50;
const int protonNumber1 = 166.0;
const int protonNumber2 = 158.0;
double_t aTime[tTotal] = {0.};


double eByWave04[tTotal] = {};
double VWave04[tTotal] = {};


double_t particlesTotal0 = 0.,particlesTotal04 = 0.;
double_t particlesTotal0_0 = 0.,particlesTotaleos0_0 = 0.;
double_t eByWave0 = 0.;

/*
double_t eByWaveN04 = 0.;
double_t eByWaveN07 = 0.;
double_t eByWaveN10 = 0.;
double_t eByWaveN15 = 0.;
*/
Double_t norm0_0= (100*137);
Double_t normeos0_0= (100*137);

Double_t norm0=0;
Double_t norm04=0;

void CampoB(){
  for(int k=10; k<11; k++){
    
    StylePlots();

    gROOT->Reset();

    TChain mychain04("T");

    mychain04.Add(TString::Format("prueba%d.root",k));

    cout<<"Opening 0.4: "<<TString::Format("prueba%d.root",k)<<endl;
    
    TChain mychainAuAu0("T");
    
    struct particula_t 
    {
      Float_t time,X,Y,Z,E,Px,Py,Pz,Pt,P,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx,frezePy,frezePz,b,nspec,R,PXR,eBx,eBy,eBz,eByWave,eByPoint,dCb;
    } PARTICLE;
    particula_t  particle;
    
    mychain04.SetBranchAddress("particle",&particle);
    mychainAuAu0.SetBranchAddress("particle",&particle);
    
    Int_t nevent04 = mychain04.GetEntries();
    Int_t nevent0 = mychainAuAu0.GetEntries();
    
    double_t nev=50;
    Double_t epsilon=0.000001;
    double_t r = 6.38; //Au (fm)
    for(Int_t j=0; j<nevent04; j++){
      mychain04.GetEvent(j);
      if(particle.time == tSteps && particle.lastcoll==0 && particle.charge==1 && particle.id == 1){
	particlesTotal04++; 
      }
      Float_t R = particle.R;
      if((particle.lastcoll==0 && particle.charge==1 && R>0 && particle.id==1) || (particle.lastcoll!=0 && particle.charge==1 && R>0 && particle.id==1)){
        for( int iTime = 1; iTime<tTotal; ++iTime){
          float_t iTC = iTime*tSteps;
          Float_t Delta = std::abs(iTC-particle.time);
          if(Delta<0.001 && TMath::Finite(particle.eBy)==1 && particle.eBy!=0){
            eByWave04[iTime] += sqrt(particle.eBx*particle.eBx+particle.eBy*particle.eBy+particle.eBz*particle.eBz);
            VWave04[iTime] += 2*iTime*TMath::Pi()*TMath::Power(r,2)*TMath::Power(particlesTotal04/(2*protonNumber2),2/3);
          cout << eByWave04[iTime] << endl;
          }
        }
      }
    }
    cout<<"0.4 done"<<endl;
    norm0=(particlesTotal0*137)/(protonNumber1);
    norm04=(particlesTotal04*137)/(protonNumber2);
    cout<<"normalizing done"<<endl;
  } // k loop
  /*
  eByWave04[0] = (eByWave0*139.57*139.57)/normeos0_0;
  eByWave07[0] = (eByWave0*139.57*139.57)/normeos0_0;
  eByWave10[0] = (eByWave0*139.57*139.57)/normeos0_0;
  eByWave15[0] = (eByWave0*139.57*139.57)/normeos0_0; 
  */

  eByWave04[0] = 0;
  cout<<"Filling other times..."<<endl;
  for( Int_t iTime = 1; iTime<tTotal; ++iTime){
    eByWave04[iTime] =(eByWave04[iTime])/norm04;
  }
  cout<<"Array done"<<endl;
  for( Int_t iTime = 0; iTime<tTotal; ++iTime){
    double_t iTC = iTime*tSteps;
    aTime[iTime] = iTC;
  }
  cout<<"Time array done"<<endl;

  
  /*for( Int_t iTime = 0; iTime<tTotal; ++iTime){
    //cout << "eos 0  "<<   eByLWeos0[iTime] <<"  "<<  eByWaveeos0[iTime] <<"  "<< eByPointeos0[iTime]<<" " <<   aTime[iTime]  <<endl;
    cout << aTime[iTime]  <<"    " <<eByWave04[iTime]<<"    " <<VWave04[iTime]<< " " << endl;
  }
  */
  
  cout<<"Drawing plot"<<endl;
  TCanvas *c =  new TCanvas("c","",650,650);
  //    c1->SetFillStyle(4000);    
  //  gStyle->SetOptStat(false);
  //c1->SetRightMargin(0.0465116);
  // c1->SetTopMargin(0.1);
  //c1->SetFillColor(0);
  TGraph *gr = new TGraph(tTotal,aTime,eByWave04); 
  gr->SetLineColor(2);
  gr->SetLineWidth(3);
  gr->SetMarkerColor(2);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(1.5);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Time [fm]");
  gr->GetYaxis()->SetTitle("-B_{y}(0,0,0)/m_{#pi}^{2}");
  gr->GetXaxis()->CenterTitle(true);
  gr->GetYaxis()->CenterTitle(true);
  //gr->GetXaxis()->SetLimits(0.,5);
  gr->SetMinimum(0.);
 // gr->SetMaximum(0.4);
  gr->Draw("AL");

  TLegend *leg1 = new TLegend(0.45,0.47,0.80,0.89);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.05);
  leg1->SetLineColor(0);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(2);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(1012);
  leg1->SetHeader("Bi+Bi, b = 7 fm, L-W","L");
  leg1->AddEntry(gr,"#sqrt{S_{NN}} = 9 GeV","L");
  //   leg1->AddEntry(parteos0wave," Wave participants (eos 0)","L"); 
  //leg1->Draw();
  
  TCanvas *c1 =  new TCanvas("c1","",650,650);
  TGraph *gr1 = new TGraph(tTotal,aTime,VWave04); 
  gr1->SetLineColor(2);
  gr1->SetLineWidth(3);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerSize(1.5);
  gr1->SetTitle("");
  gr1->GetXaxis()->SetTitle("Time [fm]");
  gr1->GetYaxis()->SetTitle("V [GeV-3]");
  gr1->GetXaxis()->CenterTitle(true);
  gr1->GetYaxis()->CenterTitle(true);
  //gr1->GetXaxis()->SetLimits(0.,5);
  gr1->SetMinimum(0.);
 // gr1->SetMaximum(0.4);
  gr1->Draw("AL");

   /* c1->SaveAs("BiBi_Mag_eBy_9.pdf");
    c1->SaveAs("BiBi_Mag_eBy_9.C");
    c1->SaveAs("BiBi_Mag_eBy_9.eps");
    c1->SaveAs("BiBi_Mag_eBy_9.root");
  */
}


//*/

//}
