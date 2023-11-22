#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

// reduce the whole dataset (N dimensions) to 2D (cos theta, phi)
RooDataSet* CosThetaPhiWeightedDataset(RooDataSet* allDataset, RooWorkspace* wspace, Bool_t isCSframe = kTRUE, Int_t ptMin = 0, Int_t ptMax = 30, Double_t massMin = 8, Double_t massMax = 14) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (mass > %f && mass < %f) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace->var(isCSframe ? "cosThetaCS" : "cosThetaHX")), *(wspace->var(isCSframe ? "phiCS" : "phiHX"))), kinematicCut);

	reducedDataset->SetName(kinematicCut); // just to make it unique

	//	wspace->import(*reducedDataset);

	return reducedDataset;
}

void drawCosThetaPhiMapWeightedDataset(Int_t minPt = 0, Int_t maxPt = 30, const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	// (cos theta, phi) 2D distribution maps for CS and HX frames

	Double_t massMin = 9., massMax = 11.;

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
	auto* cosThetaPhiDatasetCS = CosThetaPhiWeightedDataset(allDatasetCS, wspace, true, minPt, maxPt, massMin, massMax);
	auto* cosThetaPhiDatasetHX = CosThetaPhiWeightedDataset(allDatasetHX, wspace, false, minPt, maxPt, massMin, massMax);

	/// read variables in the reduced dataset in the workspace
	RooRealVar* cosThetaCS = wspace->var("cosThetaCS");
	RooRealVar* phiCS = wspace->var("phiCS");
	RooRealVar* cosThetaHX = wspace->var("cosThetaHX");
	RooRealVar* phiHX = wspace->var("phiHX");

	/// Set the plot style
	gStyle->SetPadLeftMargin(.15);
	gStyle->SetTitleYOffset(1.2);
	gStyle->SetPadRightMargin(0.2);
	gStyle->SetTitleOffset(1.4, "z");
	gStyle->SetPalette(kRainBow);
	gStyle->SetNumberContours(256);

	/// Defind legend
	TLatex* legend = new TLatex();
	legend->SetTextAlign(22);
	legend->SetTextSize(0.05);

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space in CS
	TCanvas* cCS = new TCanvas("cCS", "cCS", 700, 600);
	TH2* hCS = dynamic_cast<TH2*>(cosThetaPhiDatasetCS->createHistogram("hCS", *cosThetaCS, Binning(nCosThetaBins, cosThetaMin, cosThetaMax), YVar(*phiCS, Binning(nPhiBins, phiMin, phiMax))));

	hCS->SetTitle(";cos #theta_{CS}; #varphi_{CS} (#circ)");
	hCS->Draw("COLZ");

	legend->DrawLatexNDC(.48, .88, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, minPt, maxPt));
	legend->DrawLatexNDC(.48, .8, Form("%.0f < m_{#mu#mu} < %.0f GeV/c^{2}", massMin, massMax));

	hCS->GetYaxis()->SetRangeUser(-190, 300);
	gPad->Update();

	CMS_lumi(cCS, gCMSLumiText);
	cCS->SaveAs(Form("frame_distrib/WeightedCS_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, minPt, maxPt), "RECREATE");

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space in HX
	TCanvas* cHX = new TCanvas("cHX", "cHX", 700, 600);
	TH2* hHX = dynamic_cast<TH2*>(cosThetaPhiDatasetHX->createHistogram("hHX", *cosThetaHX, Binning(nCosThetaBins, cosThetaMin, cosThetaMax), YVar(*phiHX, Binning(nPhiBins, phiMin, phiMax))));
	hHX->SetTitle(";cos #theta_{HX}; #varphi_{HX} (#circ)");
	hHX->Draw("COLZ");

	legend->DrawLatexNDC(.48, .88, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, minPt, maxPt));
	legend->DrawLatexNDC(.48, .8, Form("%.0f < m_{#mu#mu} < %.0f GeV/c^{2}", massMin, massMax));

	hHX->GetYaxis()->SetRangeUser(-190, 300);
	gPad->Update();

	CMS_lumi(cHX, gCMSLumiText);
	cHX->SaveAs(Form("frame_distrib/WeightedHX_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, minPt, maxPt), "RECREATE");

	/// save the results in a file for later usage
	const char* outputFileName = Form("frame_distrib/WeightedCosThetaPhiMap_pt%dto%dGeV.root", minPt, maxPt);
	TFile outputFile(outputFileName, "RECREATE");

	hHX->Write();
	hCS->Write();

	outputFile.Close();

	cout << endl
	     << "CosTheta-Phi maps saved in " << outputFileName << endl;
}

void scandrawCosThetaPhiMapWeightedDataset() {
	for (Int_t idx = 0; idx < NPtBins; idx++) {
		drawCosThetaPhiMapWeightedDataset(gPtBinning[idx], gPtBinning[idx + 1], "../Files/WeightedUpsilonSkimmedDataset.root");
	}
}
