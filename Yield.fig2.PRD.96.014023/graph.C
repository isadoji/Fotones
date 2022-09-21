//-----------------------------------------------------------
// Description:
//
//       Plot macro. This macro is to plot the original Fig2.a of PhysRevD.96.014023.pdf numerical calculation.
//
// Environment:
//      ROOT
//
// Author List:
//       Isabel Domínguez Jiménez
//       isadoji@uas.edu.mx
//-----------------------------------------------------------
void graph() {

   const Int_t n = 8;
   Double_t x[n] = {0.1,0.5,1,1.5,2,2.5,3,3.5};
   Double_t y1[n] = {1.670692,0.069711,0.009130,0.002397,0.000832,0.000332,0.000142,0.000065};//eq=1/3

   Double_t y2[n] = {12.905331,0.835838,0.104421,0.027248,0.009439,0.003763,0.001614,0.000736};


   Double_t b0[n];
    for (Int_t i=0; i<n; i++) {
      b0[i] =(2*y1[i])+y2[i];
      //      cout << x[i] << "" << b0[i] << endl;
   }
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

   c1->SetGrid();
   c1->SetLogy();


   TGraph *gr0 = new TGraph(n,x,b0);
   gr0->SetLineColor(2);
   gr0->SetLineWidth(4);
   gr0->SetMarkerColor(4);
   gr0->SetMarkerStyle(21);
   //gr0->SetTitle("a simple graph");
   gr0->GetXaxis()->SetTitle("#omega_{q} [GeV]");
   gr0->GetYaxis()->SetTitle("#frac{1}{2#pi#omega_{q}} #frac{dN}{d#omega_{q}} [GeV^{-2}]");
   gr0->Draw("ACP");


   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
