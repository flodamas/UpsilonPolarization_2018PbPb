#include "../Tools/FitShortcuts.h"


void drawEtaDecayMuonGraph(Int_t minPt = 0, Int_t maxPt = 30, Int_t centMin = 0, Int_t centMax = 90, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180,
						const char* filename = "../Files/upsilonSkimmedDataset.root"){

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = false; // if extra text
	// extraText = "       Internal";

	// (cos theta, phi) 2D distribution maps for CS and HX frames

	Double_t massMin = 9.1, massMax = 9.8;

	Int_t nEtaBins = 50;
	Float_t etaMin = -2.5, etaMax = 2.5;

	Int_t nCosThetaBins = 20;

	Int_t nPhiBins = 23;
	

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Read skimmed dataset (contains angular distributions in CS and HX after kinematic cuts due to acceptance)
	RooDataSet* allDataset = (RooDataSet*)f->Get("dataset");
	/// import the dataset to a workspace
	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	/// Reduce the dataset to an interested region (cent, pT cuts)
	RooDataSet* reducedDataset = ReducedDataset(allDataset, wspace, centMin, centMax, minPt, maxPt, massMin, massMax, kFALSE, cosThetaMin, cosThetaMax, phiMin, phiMax);
	
	/// read variables in the reduced dataset in the workspace
	RooRealVar* centrality = wspace -> var("centrality");
	RooRealVar* mass = wspace -> var("mass");
	RooRealVar* pt = wspace -> var("pt");
	RooRealVar* eta_mupl = wspace -> var("etaLabMupl");
	RooRealVar* eta_mumi = wspace -> var("etaLabMumi");
	RooRealVar* cosThetaCS = wspace -> var("cosThetaCS");
	RooRealVar* phiCS = wspace -> var("phiCS");
	RooRealVar* cosThetaHX = wspace -> var("cosThetaHX");
	RooRealVar* phiHX = wspace -> var("phiHX");

	// Long64_t nEntries = reducedDataset -> sumEntries();

	/// Set the plot style
	gStyle->SetPadLeftMargin(.15);
	gStyle->SetTitleYOffset(1.1);	
	gStyle->SetPadRightMargin(0.19);
	gStyle->SetTitleOffset(0.9, "z");
	gStyle->SetPalette(kRainBow);

	/// Defind legend
	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space in CS 
	TCanvas* cEta = new TCanvas("cEta", "cEta", 700, 600);
	TH2* hEta= dynamic_cast<TH2*>(reducedDataset->createHistogram("hEta", *eta_mupl, Binning(nEtaBins,etaMin,etaMax), YVar(*eta_mumi, Binning(nEtaBins,etaMin,etaMax))));
	hEta -> SetTitle(";#eta_{#mu^{+}}; #eta_{#mu^{-}}");
	hEta -> Draw("COLZ");

	legend -> DrawLatexNDC(.33, .87, Form("centrality %d-%d%%", centMin, centMax));
	legend -> DrawLatexNDC(.325, .80, Form("%d < p_{T}^{#mu#mu} < %d GeV", minPt, maxPt));
	legend -> DrawLatexNDC(.36, .73, Form("%.1f < m_{#mu#mu} < %.1f GeV", massMin, massMax));
	
	legend -> DrawLatexNDC(.67, .40, isCSframe ? "Collins-Soper frame" : "Helicity frame");
	legend -> DrawLatexNDC(.62, .33, Form("%.2f < cos #theta < %.2f", cosThetaMin, cosThetaMax));
	legend -> DrawLatexNDC(.65, .25, Form("%d#circ < #varphi < %d#circ", phiMin, phiMax));

	hEta -> GetZaxis() -> SetRangeUser(0, 10);
	gPad -> Update();
	CMS_lumi(cEta, "2018 PbPb miniAOD, DoubleMuon PD");

	cEta -> SaveAs(Form("EtaDecayMuonPlots/HX_cent%dto%d_pt%dto%dGeV_cos%.1fto%.1f_phi%dto%d.png", centMin, centMax, minPt, maxPt, cosThetaMin, cosThetaMax, phiMin, phiMax), "RECREATE");
}

void drawEtaDecayMuonGraph_scan(){
	
	Int_t ptCuts[9] = {0, 2, 4, 6, 8, 12, 16, 20, 30};
	Bool_t isCSframe = kFALSE;
	Float_t cosThetaCuts[11] = {-1, -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1};
	Int_t phiCuts[7] = {-180, -120, -60, 0, 60, 120, 180};
	Int_t NumPhiArr = sizeof(phiCuts)/sizeof(Int_t);
	for(int idx=0; idx< NumPhiArr-1; idx++){
		
		drawEtaDecayMuonGraph(ptCuts[0], ptCuts[1], 0, 90, isCSframe, cosThetaCuts[5], cosThetaCuts[6], phiCuts[idx], phiCuts[idx+1],
						"../Files/upsilonSkimmedDataset.root");

	}

}
