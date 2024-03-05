#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"

void drawMassCosThetaDistribution(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = 0, Int_t phiMax = 180) {
	Double_t massMin = 9, massMax = 10;
	Int_t nInvMassBins = (Int_t)10 * (massMax - massMin);

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	gStyle->SetPadLeftMargin(.14);
	gStyle->SetTitleYOffset(1.1);
	gStyle->SetPadRightMargin(0.17);
	gStyle->SetTitleOffset(1., "z");
	gStyle->SetPalette(kRainBow);
	gStyle->SetNumberContours(256);

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/UpsilonSkimmedDataset.root";

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	const char* datasetName = Form("dataset%s", refFrameName);
	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	// import the dataset to a workspace
	RooWorkspace wspace(Form("workspace_%s", refFrameName));
	wspace.import(*allDataset);

	RooRealVar invMass = *wspace.var("mass");

	RooRealVar cosTheta = *wspace.var(Form("cosTheta%s", refFrameName));

	auto* data = InvMassCosThetaPhiDataset(allDataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	std::unique_ptr<RooAbsData> reducedDataset{data->reduce({cosTheta, invMass}, Form("mass > %f && mass < %f", massMin, massMax))};

	double correlationFactor = reducedDataset->correlation(cosTheta, invMass);

	cout << endl
	     << "Correlation between the invariant mass and cos theta = " << correlationFactor << endl;

	const char* histoName = Form("InvMassCosTheta%s_cent%dto%d_pt%dto%dGeV", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space

	TCanvas* canvas = new TCanvas(Form("canvas%s", refFrameName), "canvas", 700, 600);

	TH2* histo = dynamic_cast<TH2*>(reducedDataset->createHistogram(histoName, cosTheta, RooFit::Binning(nCosThetaBins, cosThetaMin, cosThetaMax), RooFit::YVar(invMass, RooFit::Binning(nInvMassBins, massMin, massMax))));

	histo->SetTitle(Form(";cos #theta_{%s}; m_{#mu^{#plus}#mu^{#minus}} (GeV/c^{2})", refFrameName));
	histo->Draw("COLZ");

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.SetTextColor(kWhite);
	legend.DrawLatexNDC(.48, .86, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend.DrawLatexNDC(.48, .78, Form("correlation = %.2f%%", 100 * correlationFactor));

	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->CenterTitle();
	//histo->GetYaxis()->SetRangeUser(massMin, massMax + 2);

	histo->GetZaxis()->SetMaxDigits(3);

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	gSystem->mkdir("CorrelationStudies", kTRUE);
	canvas->SaveAs(Form("CorrelationStudies/%s.png", histoName), "RECREATE");
}
