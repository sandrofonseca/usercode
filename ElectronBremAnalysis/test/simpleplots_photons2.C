{
           
gStyle->SetCanvasColor(10);            
gStyle->SetPadColor(0);
gStyle->SetOptStat(kFALSE);
gStyle->SetOptFit(0111);
gStyle->SetStatH(0.1); 
gStyle->SetCanvasDefH(544); //Height of canvas
gStyle->SetCanvasDefW(567); //Width of canvas

//TFile ff("muon_brem_watcher_g4_1k.root");
TFile ff("muon_brem_watcher_g4_2evt.root");

TTree *MyTree=T1;


TCanvas *Can = new TCanvas("c1","c1",129,17,926,703);
Can->SetBorderSize(2);
Can->SetFrameFillColor(0);
Can->SetLogx(0);
Can->SetLogy(0);

Can->Divide(2,2);

Can->cd(1);
TH1F *hphotons = new TH1F("hphotons","hphotons",100,0.,1000.);
MyTree->Project("hphotons","NumberPhotonsSecondaries");
//MyTree->Project("hphotons","NumberPhotonsSecondaries");
hphotons->Draw("");


Can->cd(2);
TH1F *hKenergy = new TH1F("hKenergy","hKenergy",100,0.,500.);
MyTree->Project("hKenergy","KenergyPhotonsSecondaries");



hKenergy->Draw("");

Can->cd(3);
TH1F *hpdg = new TH1F("hpdg","hpdg",100,20.,23.);
MyTree->Project("hpdg","PDGIDSecondaries");
hpdg->Draw("");

Can->cd(4);
TH1F *heta = new TH1F("heta","heta",100,-4.,4.);
MyTree->Project("heta","eta");
heta->Draw("");



Can->Update();
Can->SaveAs("muonBrem.gif");

}
