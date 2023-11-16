#include "../Tools/BasicHeaders.h"

// #include "../Tools/Style/tdrStyle.C"
// #include "../Tools/Style/CMS_lumi.C"
#include "../Tools/FitShortcuts.h"

void drawCosThetaPhiMapWeightedDataset(Int_t minPt = 0, Int_t maxPt = 30, Int_t centMin = 0, Int_t centMax = 90,
						const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = false; // if extra text
	// extraText = "       Internal";

	// (cos theta, phi) 2D distribution maps for CS and HX frames

	Double_t massMin = 8., massMax = 14.;

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	//	Double_t cosThetaBinning[] = {-1, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1};
	//	nCosThetaBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	Int_t nPhiBins = 23;
	Float_t phiMin = -180, phiMax = 280;

	//Double_t phiBinning[] = {-TMath::Pi(), -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 180};
	//	nPhiBins = sizeof(cosThetaBinning) / sizeof(Double_t) - 1;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Read skimmed dataset (contains angular distributions in CS and HX after kinematic cuts due to acceptance)
	RooDataSet* allDatasetCS = (RooDataSet*)f->Get("datasetCS");
	RooDataSet* allDatasetHX = (RooDataSet*)f->Get("datasetHX");

	/// import the dataset to a workspace
	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDatasetCS);
	wspace->import(*allDatasetHX);

	wspace->Print();

	/// Reduce the dataset to an interested region (cent, pT cuts)
	RooDataSet* reducedDatasetCS = ReducedWeightedDatasetCS(allDatasetCS, wspace, minPt, maxPt, massMin, massMax, cosThetaMin, cosThetaMax, phiMin, phiMax);
	RooDataSet* reducedDatasetHX = ReducedWeightedDatasetHX(allDatasetHX, wspace, minPt, maxPt, massMin, massMax, cosThetaMin, cosThetaMax, phiMin, phiMax);
	
	/// read variables in the reduced dataset in the workspace
	RooRealVar* centrality = wspace -> var("centrality");
	RooRealVar* mass = wspace -> var("mass");
	RooRealVar* pt = wspace -> var("pt");
	RooRealVar* cosThetaCS = wspace -> var("cosThetaCS");
	RooRealVar* phiCS = wspace -> var("phiCS");
	RooRealVar* cosThetaHX = wspace -> var("cosThetaHX");
	RooRealVar* phiHX = wspace -> var("phiHX");

	// Long64_t nEntries = reducedDataset -> sumEntries();

	/// Set the plot style
	gStyle->SetPadLeftMargin(.15);
	gStyle->SetTitleYOffset(1.2);	
	gStyle->SetPadRightMargin(0.2);
	gStyle->SetTitleOffset(1.4, "z");
	gStyle->SetPalette(kRainBow);

	/// Defind legend
	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space in CS 
	TCanvas* cCS = new TCanvas("cCS", "cCS", 650, 600);
	TH2* hCS= dynamic_cast<TH2*>(reducedDatasetCS->createHistogram("hCS", *cosThetaCS, Binning(nCosThetaBins,cosThetaMin,cosThetaMax), YVar(*phiCS, Binning(nPhiBins,phiMin,phiMax))));
	hCS -> SetTitle(";cos #theta_{CS}; #varphi_{CS} (#circ)");
	hCS -> Draw("COLZ");

	legend -> DrawLatexNDC(.48, .88, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV", centMin, centMax, minPt, maxPt));
	legend -> DrawLatexNDC(.48, .8, Form("%d < m_{#mu#mu} < %d GeV", 8, 14));

	hCS -> GetYaxis() -> SetRangeUser(-190, 300);
	gPad -> Update();
	CMS_lumi(cCS, "2018 PbPb miniAOD, DoubleMuon PD");
	cCS -> SaveAs(Form("frame_distrib/WeightedCS_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, minPt, maxPt), "RECREATE");

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space in HX 
	TCanvas* cHX = new TCanvas("cHX", "cHX", 650, 600);
	TH2* hHX= dynamic_cast<TH2*>(reducedDatasetHX->createHistogram("hHX", *cosThetaHX, Binning(nCosThetaBins,cosThetaMin,cosThetaMax), YVar(*phiHX, Binning(nPhiBins,phiMin,phiMax))));
	hHX -> SetTitle(";cos #theta_{HX}; #varphi_{HX} (#circ)");
	hHX -> Draw("COLZ");

	legend -> DrawLatexNDC(.48, .88, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV", centMin, centMax, minPt, maxPt));
	legend -> DrawLatexNDC(.48, .8, Form("%d < m_{#mu#mu} < %d GeV", 8, 14));

	hHX -> GetYaxis() -> SetRangeUser(-190, 300);
	gPad -> Update();
	CMS_lumi(cHX, "2018 PbPb miniAOD, DoubleMuon PD");
	cHX -> SaveAs(Form("frame_distrib/WeightedHX_cent%dto%d_pt%dto%dGeV.png", centMin, centMax, minPt, maxPt), "RECREATE");

	/// save the results in a file for later usage
	const char* outputFileName = Form("frame_distrib/WeightedCosThetaPhiMap_pt%dto%dGeV.root", minPt, maxPt);
	TFile outputFile(outputFileName, "RECREATE");

	hHX->Write();
	hCS->Write();

	outputFile.Close();

	cout << endl
	     << "CosTheta-Phi maps saved in " << outputFileName << endl;

}

void scandrawCosThetaPhiMapWeightedDataset(){

	Int_t PtEdges[9] = {0, 2, 4, 6, 8, 12, 16, 20, 30};
	Int_t NumPtEle = sizeof(PtEdges)/sizeof(Int_t);
	for(Int_t idx =0; idx < NumPtEle-1; idx++){
		drawCosThetaPhiMapWeightedDataset(PtEdges[idx], PtEdges[idx+1], 0, 90, "../Files/WeightedUpsilonSkimmedDataset.root");
	}
}

