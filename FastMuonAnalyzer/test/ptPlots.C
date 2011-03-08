{
gROOT->Reset();
gStyle->Reset();

//
// Parte da clonare su tutti i files :
//
const char* fileName = "FastMuonAnalyzer_RelVal_2_2_0_pre1_pt0100.root";
//const char* fileName = "FastMuonAnalyzer_pt0010_propagatorWithMaterial.root";
string stringPt = "_pt0100";

const char* HistoLegendFull = "CMSSW_220pre1 - FullSim"; 
const char* HistoLegendFast = "CMSSW_220pre1 - FastSim";

//
//

TFile f(fileName);
Color_t fullColor = kBlue;
Color_t fastColor = kRed;
Size_t  legTextSize = 0.040;
Width_t Width=2;
Style_t Style=1;
cout << "Reading from file: " << fileName << endl;
cout << HistoLegendFull << endl;
cout << HistoLegendFast << endl;
cout << "..." << endl;



//
//  STA muon pT resolution
//

{

  string stringName = "PlotPtSTA"+stringPt;
  const char* canvasName = stringName.c_str();
  c2 = new TCanvas(canvasName,canvasName,0,0,800,700);
  c2->Divide(1,2);
  c2->SetFillColor(0);

  TH1F* afullhisto = (TH1F*)  STATrackRecoPtFullHisto->Clone();  // FullSim
  
  TH1F* afasthisto = (TH1F*)  STATrackRecoPtFastHisto->Clone();  // FastSim

  TLegend  leg(0.58,0.72,0.88,0.84);
  leg->AddEntry(afullhisto,HistoLegendFull,"l");
  leg->AddEntry(afasthisto,HistoLegendFast,"l");
  leg->SetTextSize(legTextSize);
  leg->SetFillColor(0);
  
  c2->cd(1);
  // CMSSW FastSim
  afasthisto->SetStats(kFALSE);
  afasthisto->GetXaxis()->SetTitle("p_{T} res");
  afasthisto->GetYaxis()->SetTitle("");
  afasthisto->SetLineWidth(Width);
  afasthisto->SetLineColor(fastColor);
  afasthisto->SetLineStyle(Style);
  afasthisto->Draw("");
  // CMSSW FullSim:
  afullhisto->SetStats(kFALSE);
  afullhisto->SetLineWidth(Width);
  afullhisto->SetLineColor(fullColor);
  afullhisto->SetLineStyle(Style);
  afullhisto->Draw("same"); 
  leg->Draw();
  
  c2->cd(2);
  // CMSSW FastSim
  afasthisto->SetStats(kFALSE);
  afasthisto->GetXaxis()->SetTitle("p_{T} res");
  afasthisto->GetYaxis()->SetTitle("");
  afasthisto->SetLineWidth(Width);
  afasthisto->SetLineColor(fastColor);
  afasthisto->SetLineStyle(Style);
  afasthisto->Draw("");
  // CMSSW FullSim:
  afullhisto->SetStats(kFALSE);
  afullhisto->SetLineWidth(Width);
  afullhisto->SetLineColor(fullColor);
  afullhisto->SetLineStyle(Style);
  afullhisto->Draw("same"); 
  leg->Draw();


  c2->GetPad(2)->SetLogy();
  const char* jpgName = (stringName+".jpg").c_str();
  c2->Print(jpgName);

}


//
//  GeneralTrack pT resolution
//

{

  string stringName = "PlotPtTK"+stringPt;
  const char* canvasName = stringName.c_str();
  c2 = new TCanvas(canvasName,canvasName,0,0,800,700);
  c2->Divide(1,2);
  c2->SetFillColor(0);

  TH1F* afullhisto = (TH1F*)  TKRecoPtFullHisto->Clone();  // FullSim
  
  TH1F* afasthisto = (TH1F*)  TKRecoPtFastHisto->Clone();  // FastSim

  TLegend  leg(0.58,0.72,0.88,0.84);
  leg->AddEntry(afullhisto,HistoLegendFull,"l");
  leg->AddEntry(afasthisto,HistoLegendFast,"l");
  leg->SetTextSize(legTextSize);
  leg->SetFillColor(0);
  
  c2->cd(1);
  // CMSSW FastSim
  afasthisto->SetStats(kFALSE);
  afasthisto->GetXaxis()->SetTitle("p_{T} res");
  afasthisto->GetYaxis()->SetTitle("");
  afasthisto->SetLineWidth(Width);
  afasthisto->SetLineColor(fastColor);
  afasthisto->SetLineStyle(Style);
  afasthisto->Draw("");
  // CMSSW FullSim:
  afullhisto->SetStats(kFALSE);
  afullhisto->SetLineWidth(Width);
  afullhisto->SetLineColor(fullColor);
  afullhisto->SetLineStyle(Style);
  afullhisto->Draw("same"); 
  leg->Draw();
  
  c2->cd(2);
  // CMSSW FastSim
   afasthisto->SetStats(kFALSE);
  afasthisto->GetXaxis()->SetTitle("p_{T} res");
  afasthisto->GetYaxis()->SetTitle("");
  afasthisto->SetLineWidth(Width);
  afasthisto->SetLineColor(fastColor);
  afasthisto->SetLineStyle(Style);
  afasthisto->Draw("");
  // CMSSW FullSim:
  afullhisto->SetStats(kFALSE);
  afullhisto->SetLineWidth(Width);
  afullhisto->SetLineColor(fullColor);
  afullhisto->SetLineStyle(Style);
  afullhisto->Draw("same"); 
  leg->Draw();


  c2->GetPad(2)->SetLogy();
  const char* jpgName = (stringName+".jpg").c_str();
  c2->Print(jpgName);

}


gStyle->Reset();

}
