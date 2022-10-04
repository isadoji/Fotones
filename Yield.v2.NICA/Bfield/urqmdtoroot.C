#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TMath.h"
void urqmdtoroot() 
{
  //for (int i=0; i<=160; i++){
  //std::string test = to_string(i) + ".f14";
  //ifstream pFile(to_string(i));
  //ifstream pFile("100ev_200fm_AuAutestMAG.f14");
  //ifstream pFile("100EvAuAut0.f14");
  ifstream pFile("test.f14");
  
  Float_t  time,X,Y,Z,E,Px,Py,Pz,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx,frezePy,frezePz;
  Int_t nlines = 0;
  
    
  
  TFile hfile("prueba10.root","UPDATE","Demo ROOT file with histograms & trees");    
  //  TFile hfile("100EvAuAut0.root","UPDATE","Demo ROOT file with histograms & trees");
  
  
  struct particula_t 
  {
    Float_t time,X,Y,Z,E,Px,Py,Pz,Pt,P,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx,frezePy,frezePz,b,nspec,R,PXR,eBx,eBy,eBz,eByWave,eByPoint,dCb,nev;    
  } PARTICLE;
  
  particula_t  particle;
  
  TTree *tree = new TTree("T","An example of ROOT tree with a few branches");
  tree->Branch("particle",&particle,"time:X:Y:Z:E:Px:Py:Pz:Pt:P:m:id:isoespin:charge:lastcoll:numbercoll:history:frezetime:frezeX:frezeY:frezeZ:frezeE:frezePx:frezePy:frezePz:b:nspec:R:PXR:eBx:eBy:eBz:eByWave:eByPoint:dCb:nev");
  
  char header; 
  float myarray[23]; 
  Int_t endofevent=0;
  Int_t particlesperevent=0; 
  Int_t chargedparticles=0; 
  Int_t ns = 0;
  if(pFile.is_open()) 
    {
      while(!pFile.eof())
	{
	  std::string str; 
	  std::getline(pFile,str);
	  if(str[0] == 'U' || str[0] == 'p' || str[0] == 't')continue;
	  Float_t b=-999.00, bmin=-999.00, bmax=-999.00, sigma=-999.00, nev=-999.00;
	  if(str[0] == 'i'){
	    cout<<"b="<<str.substr(36,6)<<"\n";
	    cout<<"bmin="<<str.substr(43,5)<<"\n";
	    cout<<"bmax="<<str.substr(49,5)<<"\n";
	    cout<<"sigma="<<str.substr(85,9)<<"\n";
	    std::stringstream sb(str.substr(36,6));
	    std::stringstream sbmin(str.substr(43,5));
	    std::stringstream sbmax(str.substr(49,5));
	    std::stringstream ssigma(str.substr(85,9));
	    
	    sb>>b;
	    sbmin>>bmin;
	    sbmax>>bmax;
	    ssigma>>sigma;
	    
	    particle.b = b;
	  }
	  if(str[0] == 'e' && str[1] == 'v'){
	    cout<<"N event="<<str.substr(7,9)<<"\n";
	    std::stringstream snev(str.substr(7,9));
	    snev>>nev;
	    particle.nev = nev;
	  }
	  if(str[0] == 'e' && str[1] == 'q') continue;
	  if(str[0] == 'o') continue;
	  std::stringstream ss(str);
	  Int_t position = 0 ; 
	  while(ss>>myarray[position])
	    {
	      //printf("\n position %d =  %f \n ", position, myarray[position]);
	      position++;  
	    }
	  
	  //if(myarray[0]==200.00)
	  //{
	  ++particlesperevent;          
	  time = myarray[0];
	  X = myarray[1];
	  Y = myarray[2];
	  Z = myarray[3];
	  E = myarray[4];
	  Px = myarray[5];
	  Py = myarray[6];
	  Pz = myarray[7];
	  m = myarray[8];
	  id = myarray[9];
	  isoespin = myarray[10];
	  charge = myarray[11];
	  if(charge>0 || charge <0)++chargedparticles; 
	  lastcoll = myarray[12];
	  numbercoll = myarray[13];
	  history = myarray[14];
	  frezetime = myarray[15];
	  frezeX = myarray[16];
	  frezeY = myarray[17];
	  frezeZ = myarray[18];
	  frezeE = myarray[19];
	  frezePx = myarray[20];
	  frezePy = myarray[21];
	  frezePz = myarray[22];
	  //} 
	  
	  //Calculate  Eta, Pt, and P
	  Float_t  Pt= -999.0 ,P= -999.0 ,nspec=0,eBx=-999.0,eBy=-999.0,eBz=-999.0,R=-999.0,PXR=-999.0;
	  Float_t dCb = 0.0, eByWave = 0.0, eByPoint = 0.0;
	  Pt = TMath::Sqrt(Px*Px+ Py*Py);
	  P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
	  
	  
	  R= TMath::Sqrt(X*X+Y*Y+Z*Z);
	  PXR= TMath::Sqrt((Py*Z-Pz*Y)*(Py*Z-Pz*Y)+(X*Pz-Px*Z)*(X*Pz-Px*Z)+(Px*Y-Py*X)*(Px*Y-Py*X));
	  if(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3)!=0){
	    
	    eBx = (197/139.57)*(197/139.57)*(charge)*(E*E-P*P)*(Py*Z-Pz*Y)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
	    eBy = (197/139.57)*(197/139.57)*(E*E-P*P)*(Pz*X-Px*Z)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
	    eBz = (197/139.57)*(197/139.57)*(E*E-P*P)*(Px*Y-Py*X)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
	    //}
	    //if( R>0.3 && E>0.3){
	    //Wave
	    
	    dCb = 1.44*(2/(TMath::Pi()*R)*TMath::Exp(-0.25*R*R)-TMath::Erf(TMath::Sqrt(0.25)*R)/(pow(R,2)));
	    eByWave = -(197/139.57)*(197/139.57)*(dCb)*(Pz*X-Px*Z)/(E*R);
	    //Point
	    eByPoint = (197/139.57)*(197/139.57)*(charge)*(Pz*X-Px*Z)/(pow(R,3)*E);
	  }
	  else{	    	dCb = 0.;}
	  if(numbercoll == 0){
	    ++ns;
	    nspec = ns;
	  }
	  
	  particle.time = time; 
	  particle.X = X; 
	  particle.Y = Y; 
	  particle.Z = Z;
	  particle.E = E;
	  particle.Px = Px;
	  particle.Py = Py;
	  particle.Pz = Pz;
	  particle.Pt = Pt;
	  particle.P= P;
	  particle.m = m; 
	  particle.id = id;
	  particle.isoespin = isoespin ;
	  particle.charge = charge;
	  particle.lastcoll = lastcoll ;
	  particle.numbercoll = numbercoll;   
	  particle.history = history; 
	  particle.frezetime = frezetime;
	  particle.frezeX = frezeX ;
	  particle.frezeY = frezeY ;
	  particle.frezeZ = frezeZ;
	  particle.frezeE = frezeE;   
	  particle.frezePx = frezePx; 
	  particle.frezePy = frezePy;
	  particle.frezePz = frezePz;  
	  particle.nspec = nspec;
	  particle.R = R;
	  particle.PXR = PXR;
	  particle.eBx = eBx;
	  particle.eBy = eBy;
	  particle.eBz = eBz;
	  particle.eByWave = eByWave;
	  particle.eByPoint = eByPoint;
	  particle.dCb = dCb;
	  nlines++; 
	  //event.nparticles = nlines
	  
	  tree->Fill();
	}
    }
  
  else
    {
      printf("Unable to open file");
    }
  
  
  
  printf(" found %d points\n",nlines);
  
  pFile.close();
  
  hfile.Write(); 
  hfile.Close();
  
  //f->Write();
  //}
  
} 
