#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

TH2* DrawCosThetaPhiDistribution(RooDataSet* dataset, RooWorkspace& wspace, const char* frameAcronym = "CS", Int_t ptMin = 0, Int_t ptMax = 30) {
	Double_t massMin = 9., massMax = 11.;

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	Int_t nPhiBins = 25;
	Float_t phiMin = -200, phiMax = 300;

	const char* histoName = Form("%s_cent%dto%d_pt%dto%dGeV", frameAcronym, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);

	wspace.import(*dataset);
	RooRealVar cosThetaVar = *wspace.var(Form("cosTheta%s", frameAcronym));
	RooRealVar phiVar = *wspace.var(Form("phi%s", frameAcronym));

	// reduce the whole dataset (N dimensions) to 2D (cos theta, phi)
	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (mass > %f && mass < %f) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet* reducedDataset = (RooDataSet*)dataset->reduce(RooArgSet(cosThetaVar, phiVar), kinematicCut);

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space

	TCanvas* canvas = new TCanvas(Form("canvas%s", frameAcronym), "canvas", 700, 600);

	TH2* histo = dynamic_cast<TH2*>(reducedDataset->createHistogram(histoName, cosThetaVar, RooFit::Binning(nCosThetaBins, cosThetaMin, cosThetaMax), RooFit::YVar(phiVar, RooFit::Binning(nPhiBins, phiMin, phiMax))));

	histo->SetTitle(Form(";cos #theta_{%s}; #varphi_{%s} (#circ)", frameAcronym, frameAcronym));
	histo->Draw("COLZ");

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.48, .86, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend.DrawLatexNDC(.48, .78, Form("%.0f < m_{#mu#mu} < %.0f GeV/c^{2}", massMin, massMax));

	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->SetRangeUser(-200, 300);
	histo->GetYaxis()->CenterTitle();
	histo->GetZaxis()->SetMaxDigits(3);

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	canvas->SaveAs(Form("frame_distrib/%s.png", histoName), "RECREATE");

	return histo;
}

void drawCosThetaPhiMap(Int_t ptMin = 0, Int_t ptMax = 30, const char* filename = "../Files/UpsilonSkimmedDataset.root") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = false; // if extra text
	extraText = "       Internal";

	gStyle->SetPadLeftMargin(.15);
	gStyle->SetTitleYOffset(1.2);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetTitleOffset(1.2, "z");
	gStyle->SetPalette(kRainBow);
	gStyle->SetNumberContours(256);

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space
	RooWorkspace wspace("workspace");

	TH2* histoCS = DrawCosThetaPhiDistribution((RooDataSet*)f->Get("datasetCS"), wspace, "CS", ptMin, ptMax);

	TH2* histoHX = DrawCosThetaPhiDistribution((RooDataSet*)f->Get("datasetHX"), wspace, "HX", ptMin, ptMax);

	/// save the results in a file for later usage
	gSystem->mkdir("frame_distrib", kTRUE);
	TFile outputFile("frame_distrib/CosThetaPhiMaps.root", "UPDATE");

	histoHX->Write();
	histoCS->Write();

	outputFile.Close();

	cout << endl
	     << "CosTheta-Phi maps saved in " << outputFile.GetName() << endl;
}

void scanDrawCosThetaPhiMap(const char* filename = "../Files/UpsilonSkimmedDataset.root") {
	for (Int_t idx = 0; idx < NPtBins; idx++) {
		drawCosThetaPhiMap(gPtBinning[idx], gPtBinning[idx + 1], filename);
	}
}
