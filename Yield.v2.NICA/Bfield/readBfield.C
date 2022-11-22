Int_t nbin = 200;// #events
Int_t xmin = 0;// #events
Double_t xmax = 0.5;// #events
Float_t norm = nbin/(TMath::Abs(xmin)+TMath::Abs(xmax));
TH1F *eBH = new TH1F("eBH","", nbin,xmin,xmax);
TH1F *VH = new TH1F("VH","", nbin,xmin,xmax);

const Int_t entriesT = 10000000;
Int_t nT;
Double_t eBT[entriesT];
Double_t VT[entriesT];
Double_t timeT[entriesT];
TTree* BfieldT = new  TTree("BfieldT","B field tree ");
Double_t eB1,V1,tau1,n1;
Double_t eB2,V2,tau2,n2;
Double_t eB3,V3,tau3,n3;
Double_t eB4,V4,tau4,n4;
Double_t eB5,V5,tau5,n5;
Double_t eB[50];
Double_t V[50];
Double_t tau[50];
Int_t n[50];

void readBfield(){
char name[50];
  sprintf(name,"Bfield.root");
  cout << " Openning file " << name << endl;
  TFile *file = new TFile(name);

  TTree *BfieldT = (TTree*) file->Get("BfieldT");
  BfieldT->SetBranchAddress("nT",&nT);
  BfieldT->SetBranchAddress("eBT",eBT);
  BfieldT->SetBranchAddress("VT",VT);
  BfieldT->SetBranchAddress("timeT",timeT);

Int_t entries = BfieldT->GetEntries();
  for (Int_t i =0; i< entries;i++){
    BfieldT->GetEntry(i);
    for (Int_t j =0; j< nT;j++){
    BfieldT->GetEntry(j);
    //cout << timeT[j] << endl; 
    
    for(Int_t k =0;k<=50;k++){
    Float_t step = k*0.01; 
    if(timeT[j] <= step){
    eB[k] += eBT[j];
    V[k] += VT[j];
    tau[k] = step;
    
    //cout << tau[k] << " " << eBT[k] << endl;
    }
    }
    

  }
     }
  for(Int_t l =0;l<=50;l++){
    cout << " " << tau[l] << " " << eB[l] << " " << V[l] << endl;

    eBH->Fill(tau[l],eB[l]);
    VH->Fill(tau[l],V[l]);
  }
  
  TCanvas* c1 = new TCanvas("c1","UrQMD test example",800,800);
  //gStyle->SetOptStat(false);
  c1->SetRightMargin(0.0465116);
  c1->SetTopMargin(0.1);
  c1->SetFillColor(0);
  c1->SetLogy();
  eBH->Draw("hist p");
  eBH->SetMarkerColor(1);
  eBH->SetMarkerStyle(20);
  eBH->SetXTitle("t [fm]");
  eBH->GetXaxis()->CenterTitle(true);
  eBH->GetXaxis()->SetTitleSize(0.04);
  eBH->GetXaxis()->SetLabelSize(0.03);
  eBH->GetXaxis()->SetTitleOffset(1.2);
  eBH->SetYTitle("e|B|/m_{#pi}^{2}");
  eBH->GetYaxis()->CenterTitle(true);
  eBH->GetYaxis()->SetTitleSize(0.04);
  eBH->GetYaxis()->SetLabelSize(0.03);
  eBH->GetYaxis()->SetTitleOffset(1.2);

  TCanvas* c2 = new TCanvas("c2","UrQMD test example",800,800);
  gStyle->SetOptStat(false);
  c2->SetRightMargin(0.0465116);
  c2->SetTopMargin(0.1);
  c2->SetFillColor(0);
  VH->Draw("hist p");
  VH->SetMarkerColor(1);
  VH->SetMarkerStyle(20);
  VH->SetXTitle("t [fm]");
  VH->GetXaxis()->CenterTitle(true);
  VH->GetXaxis()->SetTitleSize(0.04);
  VH->GetXaxis()->SetLabelSize(0.03);
  VH->GetXaxis()->SetTitleOffset(1.2);
  VH->SetYTitle("V ");
  VH->GetYaxis()->CenterTitle(true);
  VH->GetYaxis()->SetTitleSize(0.04);
  VH->GetYaxis()->SetLabelSize(0.03);
  VH->GetYaxis()->SetTitleOffset(1.2);
 
}