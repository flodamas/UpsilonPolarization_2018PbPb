#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "SPlotHelpers.h"

#include "../Tools/Style/Legends.h"

// compare the sPlot and raw yield extraction methods
void compareRawCosThetaDistrib(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180, const char* filename = "../Files/UpsilonSkimmedDataset.root") { //possible refFrame names: CS or HX

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 7;
	Float_t cosThetaMin = -0.7, cosThetaMax = 0.7;

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooWorkspace wspace = SetUpWorkspace(filename);

	auto data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax);

	// read variables in the reduced dataset in the workspace
	RooRealVar invMass = *wspace.var("mass");

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

	// create weighted data sets
	RooDataSet data_weight1S = GetSpeciesSWeightedDataset(sData, "1S");

	/// Standard extraction of the raw yield per bin of cos theta
	RooRealVar* yield1S = wspace.var("yield1S");
	RooRealVar* yield2S = wspace.var("yield2S");
	// loop and fit

	TH1D yieldHist("yieldHist", " ", nCosThetaBins, cosThetaMin, cosThetaMax);

	Float_t cosThetaStep = ((cosThetaMax - cosThetaMin) / nCosThetaBins);

	Float_t maxYield = 0;

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		Float_t cosThetaVal = cosThetaMin + iCosTheta * cosThetaStep;

		cout << "Invariant mass fit for cos theta = [" << cosThetaVal << ", " << cosThetaVal + cosThetaStep << "]" << endl;

		RooDataSet reducedDataset = *(RooDataSet*)data.reduce({invMass}, Form("cosTheta%s > %f && cosTheta%s < %f", refFrameName, cosThetaVal, refFrameName, cosThetaVal + cosThetaStep));

		auto* massFitResult = RawInvariantMassFit(reducedDataset, invMassModel);

		double rawYield = yield1S->getVal();
		double errorRawYield = yield1S->getError();

		yieldHist.SetBinContent(iCosTheta + 1, rawYield);
		yieldHist.SetBinError(iCosTheta + 1, errorRawYield);

		if (yield1S->getVal() > maxYield) maxYield = rawYield;
	}

	RooDataHist yieldRooDataHist("yieldRooDataHist", " ", cosTheta, Import(yieldHist));

	/// Draw the cos theta distributions

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Range(cosThetaMin, cosThetaMax));

	yieldRooDataHist.plotOn(frame, DrawOption("P0Z"), MarkerColor(kAzure + 2), Name("rawYield"));

	data_weight1S.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(gColor1S), Name("sData"));

	frame->GetYaxis()->SetMaxDigits(3);
	frame->SetMaximum(2 * maxYield);

	frame->Draw();

	gPad->RedrawAxis();

	//frame->SetMaximum(nEntries / 15);

	TLegend legend(.22, .85, .45, .65);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("%d < p_{T}^{#mu#mu} < %d GeV/c, %d#circ < #varphi_{%s} < %d#circ", ptMin, ptMax, phiMin, refFrameName, phiMax));
	legend.AddEntry(frame->findObject("rawYield"), "#varUpsilon(1S) raw yield", "lp");
	legend.AddEntry(frame->findObject("sData"), "sWeighted data", "lp");

	legend.DrawClone();

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	gSystem->mkdir("1D", kTRUE);
	canvas->SaveAs(Form("1D/compareRawCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}
