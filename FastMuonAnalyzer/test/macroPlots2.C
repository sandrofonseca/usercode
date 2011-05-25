{

gROOT->Reset();

/////////////////////////////////////////////////////////////////// 
 TFile *a = new TFile("TrackerMuonAnalyzer_FASTfromFILE.root");
 TFile *b = new TFile("TrackerMuonAnalyzer_FASTfromFILE_MuonBremTrue.root");	
 TFile *c = new TFile("TrackerMuonAnalyzer_FULL_RelValSingleMuPt1000.root");	
////////////////////////////////////////////////////////////////////////////////////

 TCanvas *c1 = new TCanvas("c1","Validation", 0, 0, 700, 700);
  // c1->Range(0,0,1,1);
// c1->SetFillColor(41);
   c1->Clear();
   c1->cd(0); 



	
	
	//  /////////////////////////////////////////////////////
	a.cd();
    TH1F *hptDiffa = (TH1F*)DifftrackerPoutPgenInvPgen->Clone();
 
// THStack *hs = new THStack("hs",GenDilMass0_Zp331M1000->GetTitle());
 hptDiffa->SetName("DifftrackerPoutPgenInvPgen");
 hptDiffa->SetLineColor(kBlue);
 hptDiffa->SetLineStyle(2);
 hptDiffa->SetLineWidth(2); 
//hptDiffa->SetFillColor(kBlue);
 //hptDiffa->SetMarkerStyle(21);
 
 c1->cd();hptDiffa->Draw();

//hs->Add(GenDilMass0_Zp331M1000);


 /////////////////////////////////////////////////////
	b.cd();
	
	TH1F *hptDiffb = (TH1F*)DifftrackerPoutPgenInvPgen->Clone();
	hptDiffb->SetName("DifftrackerPoutPgenInvPgen");
	hptDiffb->SetLineColor(kMagenta);
	hptDiffb->SetLineStyle(5);
	hptDiffb->SetLineWidth(4); 
	
	c1->cd();hptDiffb->Draw("same");
	
	
	
	/////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////
	// /////////////////////////////////////////////////////
// 	c.cd();
	
// 	TH1F *hptDiffc = (TH1F*)trackerPtOn->Clone();
// 	hptDiffc->SetName("trackPtOut");
// 	hptDiffc->SetLineColor(kRed);
// 	hptDiffc->SetLineStyle(2);
// 	hptDiffc->SetLineWidth(2); 
// 	c1->cd();hptDiffc->Draw("same");
	
	/////////////////////////////////////////////	
//  /////////////////////////////////////////////////////
    c.cd();
    TH1F *hptDiffc = (TH1F*)DifftrackerPoutPgenInvPgen->Clone();

// THStack *hs = new THStack("hs",GenDilMass0_Zp331M1000->GetTitle());
 hptDiffc->SetName("DifftrackerPoutPgenInvPgen");
 hptDiffc->SetLineColor(kRed);
 hptDiffc->SetLineStyle(5);
 hptDiffc->SetLineWidth(7);
//hptDiffa->SetFillColor(kBlue);
 //hptDiffa->SetMarkerStyle(21);

 c1->cd();hptDiffc->Draw("same");
	
	
	// legend
	c1->cd();
	TLegend *leg = new TLegend(0.6994286,0.7896296,0.8914286,0.922963,NULL,"brNDC");
	leg->AddEntry(hptDiffa,"FastSim Standard","f");
	leg->AddEntry(hptDiffb,"FastSim+MuonBrem","f");
	leg->AddEntry(hptDiffc,"FullSim","f");

	leg->Draw();
	c1->SaveAs("DifftrackerPoutPgenInvPgen_Full_Fast.png");

////////////////////////////////////////////////////////////////////////////


}
