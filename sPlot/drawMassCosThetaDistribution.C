#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"

void drawMassCosThetaDistribution(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180) {
	Double_t massMin = 7, massMax = 13;
	Int_t nInvMassBins = (Int_t)10 * (massMax - massMin);

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	gStyle->SetPadLeftMargin(.14);
	gStyle->SetTitleYOffset(1.1);
	gStyle->SetPadRightMargin(0.17);
	gStyle->SetTitleOffset(1., "z");
	SetColorPalette("KlimtGold");

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/UpsilonSkimmedDataset.root";

	RooWorkspace wspace = SetUpWorkspace(filename, refFrameName);

	RooRealVar invMass = *wspace.var("mass");

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));

	auto data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	std::unique_ptr<RooAbsData> reducedDataset{data.reduce({cosTheta, invMass}, Form("mass > %f && mass < %f", massMin, massMax))};

	double correlationFactor = reducedDataset->correlation(cosTheta, invMass);

	cout << endl
	     << "Correlation between the invariant mass and cos theta = " << correlationFactor << endl;

	const char* histoName = Form("InvMassCosTheta%s_cent%dto%d_pt%dto%dGeV", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);

	/// Draw and save the number of events(signal+background) plots in the 2D (costheta, phi) space

	TCanvas* canvas = new TCanvas("canvas", "canvas", 700, 600);

	TH2* histo = dynamic_cast<TH2*>(reducedDataset->createHistogram(histoName, cosTheta, RooFit::Binning(nCosThetaBins, cosThetaMin, cosThetaMax), RooFit::YVar(invMass, RooFit::Binning(nInvMassBins, massMin, massMax))));
	histo->SetTitle(" ");
	histo->Draw("COLZ");

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.SetTextColor(kWhite);
	legend.DrawLatexNDC(.48, .86, Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
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
