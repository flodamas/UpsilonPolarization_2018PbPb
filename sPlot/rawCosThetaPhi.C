#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "SPlotHelpers.h"

#include "../Tools/Style/Legends.h"
#include "../Tools/Style/Figures.h"

void drawAndSaveDistribution(TH2* histo, const char* name, Int_t ptMin, Int_t ptMax, const char* legend) {
	TCanvas* canvas = new TCanvas("canvas", "canvas", 700, 600);

	histo->Draw("COLZ");

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.045);
	text.SetTextColor(kWhite);

	text.DrawLatexNDC(.48, .8, Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	text.DrawLatexNDC(.48, .72, legend);

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	gSystem->mkdir("2D", kTRUE);
	canvas->SaveAs(Form("2D/%s.png", name), "RECREATE");

	delete canvas;
}

void rawCosThetaPhi(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Float_t phiMin = -180, Float_t phiMax = 180, const char* filename = "../Files/UpsilonSkimmedDataset_TriggerAcc.root") {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 10, nPhiBins = 10;

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooWorkspace wspace = SetUpWorkspace(filename);

	auto allData = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax);

	// read variables in the reduced dataset in the workspace

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));

	RooRealVar phi = *wspace.var(PhiVarName(refFrameName));

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	const char* signalShapeName = "SymDSCB";

	// background
	//int order = 2;
	//const char* bkgShapeName = Form("ChebychevOrder%d", order);
	const char* bkgShapeName = "ExpTimesErr";

	/// SPlot time!

	auto* sData = SWeightedDataset(wspace, ptMin, ptMax, signalShapeName, bkgShapeName);

	/// Draw the (cos theta, phi) distributions with and without sWeights
	gStyle->SetPadLeftMargin(.14);
	gStyle->SetTitleYOffset(1.);
	gStyle->SetPadRightMargin(0.17);
	gStyle->SetTitleOffset(1., "z");
	SetColorPalette("TamDragon");

	// select the mass window of the Y peaks for visualisation
	Float_t lowMassCut = 8.5, highMassCut = 10.5;

	const char* massCut = Form("mass > %f && mass < %f", lowMassCut, highMassCut);

	const char* histoName = Form("rawCosThetaPhi%s_cent%dto%d_pt%dto%dGeV", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);

	auto data = *(RooDataSet*)sData->reduce(massCut);

	// all data in mass window

	const char* nameAll = Form("%s_AllData", histoName);

	auto histoAll = TH2fromRooDataSet(data, nameAll, cosTheta, nCosThetaBins, cosThetaMin, cosThetaMax, phi, nPhiBins, phiMin, phiMax);

	drawAndSaveDistribution(histoAll, nameAll, ptMin, ptMax, Form("candidates in %.1f < %s < %.1f %s", lowMassCut, gMassVarTitle, highMassCut, gMassUnit));

	//allData.plotOn(frame, DrawOption("P0Z"), Name("allData"), Cut(massCut));

	// background in mass window
	RooDataSet data_weightBkg = GetSpeciesSWeightedDataset(&data, "Bkg");

	const char* nameBkg = Form("%s_Bkg", histoName);

	auto histoBkg = TH2fromRooDataSet(data_weightBkg, nameBkg, cosTheta, nCosThetaBins, cosThetaMin, cosThetaMax, phi, nPhiBins, phiMin, phiMax);

	drawAndSaveDistribution(histoBkg, nameBkg, ptMin, ptMax, Form("background in %.1f < %s < %.1f %s", lowMassCut, gMassVarTitle, highMassCut, gMassUnit));

	// Y(1S)
	RooDataSet data_weight1S = GetSpeciesSWeightedDataset(&data, "1S");

	const char* name1S = Form("%s_1S", histoName);

	auto histo1S = TH2fromRooDataSet(data_weight1S, name1S, cosTheta, nCosThetaBins, cosThetaMin, cosThetaMax, phi, nPhiBins, phiMin, phiMax);

	drawAndSaveDistribution(histo1S, name1S, ptMin, ptMax, "#varUpsilon(1S) sWeights");
}
