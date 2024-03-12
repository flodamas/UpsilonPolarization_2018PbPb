#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../Tools/Datasets/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

// compare the sPlot and raw yield extraction methods
void compareRawCosThetaDistrib(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180, const char* filename = "../Files/UpsilonSkimmedDataset.root") { //possible refFrame names: CS or HX

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 10;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooWorkspace wspace = SetUpWorkspace(filename, refFrameName);

	auto* data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	// read variables in the reduced dataset in the workspace
	RooRealVar* invMass = wspace.var("mass");

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

	massCanvas->SaveAs(Form("InvMassFits/rawInvMassFit_%s_cent%dto%d_pt%dto%dGeV.png", bkgShapeName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	/// SPlot time!
	SPlot sData = CreateSPlot(wspace, data, invMassModel);

	// create weighted data sets
	RooDataSet data_weight1S = GetSWeightedDataset(data, "1S");

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

		std::unique_ptr<RooAbsData> reducedDataset{data->reduce({*invMass}, Form("cosTheta%s > %f && cosTheta%s < %f", refFrameName, cosThetaVal, refFrameName, cosThetaVal + cosThetaStep))};

		auto* massFitResult = invMassModel->fitTo(*reducedDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), Minos({*yield1S, *yield2S}));

		massFitResult->Print("v");

		double rawYield = yield1S->getVal();
		double errorRawYield = yield1S->getError();

		yieldHist.SetBinContent(iCosTheta + 1, rawYield);
		yieldHist.SetBinError(iCosTheta + 1, errorRawYield);

		if (yield1S->getVal() > maxYield) maxYield = rawYield;
	}

	RooDataHist yieldRooDataHist("yieldRooDataHist", " ", *cosTheta, Import(yieldHist));

	/// Draw the cos theta distributions

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta->frame(Title(" "), Range(cosThetaMin, cosThetaMax));

	yieldRooDataHist.plotOn(frame, DrawOption("P0Z"), MarkerColor(kAzure + 2), Name("rawYield"));

	data_weight1S.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(kRed), Name("sData"));

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
