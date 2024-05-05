#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "SPlotHelpers.h"

#include "../Tools/Style/Legends.h"

void rawCosTheta(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t nCosThetaBins = 10, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Int_t phiMin = -180, Int_t phiMax = 180, const char* filename = "../Files/UpsilonSkimmedDataset.root") {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooWorkspace wspace = SetUpWorkspace(filename);

	auto data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax);

	// read variables in the reduced dataset in the workspace

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	const char* signalShapeName = "SymDSCB";

	// background
	//int order = 2;
	//const char* bkgShapeName = Form("ChebychevOrder%d", order);
	const char* bkgShapeName = "ExpTimesErr";

	/// SPlot time!

	auto* sData = SWeightedDataset(wspace, ptMin, ptMax, signalShapeName, bkgShapeName);

	/// Draw the cos theta distributions with and without sWeights

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Bins(nCosThetaBins), Range(cosThetaMin, cosThetaMax));

	// select the mass window of the Y peaks for visualisation
	Float_t lowMassCut = 8.5, highMassCut = 10.5;

	const char* massCut = Form("mass > %f && mass < %f", lowMassCut, highMassCut);
	data.plotOn(frame, DrawOption("P0Z"), Name("data"), Cut(massCut));

	// create sWeighted data sets
	RooDataSet data_weightBkg = GetSpeciesSWeightedDataset(sData, "Bkg");
	RooDataSet data_weight1S = GetSpeciesSWeightedDataset(sData, "1S");
	RooDataSet data_weight2S = GetSpeciesSWeightedDataset(sData, "2S");

	data_weightBkg.plotOn(frame, DrawOption("P0Z"), MarkerColor(gColorBkg), Name("dataBkg"), Cut(massCut));

	data_weight1S.plotOn(frame, DrawOption("P0Z"), MarkerColor(gColor1S), Name("data1S"));

	data_weight2S.plotOn(frame, DrawOption("P0Z"), MarkerColor(gColor2S), Name("data2S"));

	frame->GetYaxis()->SetMaxDigits(3);
	gStyle->SetExponentOffset(-0.08, 0.005, "Y");

	frame->SetMaximum(1.6 * frame->GetMaximum());

	frame->Draw();

	gPad->RedrawAxis();

	// cosmetics

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.05);
	text.DrawLatexNDC(.55, .85, Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));

	TLegend legend(.2, .8, .45, .55);
	legend.SetTextSize(.045);

	legend.AddEntry(frame->findObject("data"), Form("candidates with %.1f < %s < %.1f %s", lowMassCut, gMassVarTitle, highMassCut, gMassUnit), "lp");
	legend.AddEntry(frame->findObject("dataBkg"), "with background sWeights", "lp");

	legend.AddEntry(frame->findObject("data1S"), "with #varUpsilon(1S) sWeights", "lp");
	legend.AddEntry(frame->findObject("data2S"), "with #varUpsilon(2S) sWeights", "lp");

	legend.DrawClone();

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	gSystem->mkdir("1D", kTRUE);
	canvas->SaveAs(Form("1D/rawCosTheta%s_cent%dto%d_pt%dto%dGeV.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
}
