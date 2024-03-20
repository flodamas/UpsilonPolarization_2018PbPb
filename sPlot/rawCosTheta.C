#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../Tools/Datasets/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"

void rawCosTheta(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t nCosThetaBins = 10, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Int_t phiMin = -180, Int_t phiMax = 180, const char* filename = "../Files/UpsilonSkimmedDataset.root") {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooWorkspace wspace = SetUpWorkspace(filename, refFrameName);

	auto* data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	// read variables in the reduced dataset in the workspace

	RooRealVar* cosTheta = wspace.var(Form("cosTheta%s", refFrameName));

	Long64_t nEntries = data->sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	const char* signalShapeName = "SymDSCB";

	// background
	//int order = 2;
	//const char* bkgShapeName = Form("ChebychevOrder%d", order);
	const char* bkgShapeName = "ExpTimesErr";

	auto* invMassModel = MassFitModel(wspace, signalShapeName, bkgShapeName, ptMin, ptMax, nEntries);

	auto* fitResult = invMassModel->fitTo(*data, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), Minos(kTRUE));

	fitResult->Print("v");

	/// Draw the invariant mass distribution, to check the fit
	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	gSystem->mkdir("InvMassFits", kTRUE);
	massCanvas->SaveAs(Form("InvMassFits/rawInvMassFit_%s_cent%dto%d_pt%dto%dGeV.png", bkgShapeName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	/// SPlot time!

	SPlot sData = CreateSPlot(wspace, data, invMassModel);

	/// Draw the cos theta distributions with and without sWeights

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta->frame(Title(" "), Bins(nCosThetaBins), Range(cosThetaMin, cosThetaMax));

	// select the mass window of the Y peaks for visualisation
	Float_t lowMassCut = 8.5, highMassCut = 10.5;
	data->plotOn(frame, DrawOption("P0Z"), Name("data"), Cut(Form("mass > %f && mass < %f", lowMassCut, highMassCut)));

	// create sWeighted data sets
	RooDataSet data_weightBkg = GetSWeightedDataset(data, "Bkg");
	RooDataSet data_weight1S = GetSWeightedDataset(data, "1S");
	RooDataSet data_weight2S = GetSWeightedDataset(data, "2S");

	data_weightBkg.plotOn(frame, DrawOption("P0Z"), MarkerColor(kGray + 2), Name("dataBkg"), Cut(Form("mass > %f && mass < %f", lowMassCut, highMassCut)));

	data_weight1S.plotOn(frame, DrawOption("P0Z"), MarkerColor(kRed), Name("data1S"));

	data_weight2S.plotOn(frame, DrawOption("P0Z"), MarkerColor(kGreen + 2), Name("data2S"));

	frame->GetYaxis()->SetMaxDigits(3);

	frame->SetMaximum(1.6 * frame->GetMaximum());

	frame->Draw();

	gPad->RedrawAxis();

	/// Polarization fit
	/*
	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -2, 2);

	auto cosThetaPDF_1S = CosThetaPolarizationPDF("cosThetaPDF_1S", " ", *cosTheta, lambdaTheta);

	auto* polarizationFitResult = cosThetaPDF_1S.fitTo(data_weight1S, Save(), Extended(kTRUE), PrintLevel(+1), NumCPU(NCPUs), Range(cosThetaMin, cosThetaMax), AsymptoticError(DoAsymptoticError));

	polarizationFitResult->Print("v");

	cosThetaPDF_1S.plotOn(frame, LineColor(kRed + 1), Name("polaResult"));

	frame->Draw();
*/
	// cosmetics

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.05);
	text.DrawLatexNDC(.55, .85, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	TLegend legend(.2, .8, .45, .55);
	legend.SetTextSize(.045);

	legend.AddEntry(frame->findObject("data"), Form("candidates with %.1f < %s < %.1f %s", lowMassCut, wspace.var("mass")->GetTitle(), highMassCut, wspace.var("mass")->getUnit()), "lp");
	legend.AddEntry(frame->findObject("dataBkg"), "with background sWeights", "lp");

	legend.AddEntry(frame->findObject("data1S"), "with #varUpsilon(1S) sWeights", "lp");
	legend.AddEntry(frame->findObject("data2S"), "with #varUpsilon(2S) sWeights", "lp");

	legend.DrawClone();

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	gSystem->mkdir("1D", kTRUE);
	canvas->SaveAs(Form("1D/rawCosTheta%s_cent%dto%d_pt%dto%dGeV.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
}
