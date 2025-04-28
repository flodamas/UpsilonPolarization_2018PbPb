#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../Tools/RooFitPDFs/InvariantMassModels.h"

#include "../MonteCarlo/AccEffHelpers.h"
#include "../Polarization/PolarFitHelpers.h"

using namespace std;

RooDataSet InvMassDataset(RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Float_t cosThetaMin = -0.1, Float_t cosThetaMax = 0.1, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax) {
	RooDataSet* allDataset = (RooDataSet*)wspace.data(RawDatasetName(""));

	if (allDataset == nullptr) {
		std::cerr << "Null RooDataSet provided to the reducer method!!" << std::endl;
		return RooDataSet();
	}

	// const char* kinematicCut = Form("centrality >= %d && centrality < %d && rapidity > %f && rapidity < %f && pt > %d && pt < %d && cosTheta%s > %f && cosTheta%s < %f && %s > %d && %s < %d", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);
	const char* kinematicCut = Form("centrality >= %d && centrality < %d && mass >= %f && mass < %f && rapidity > %f && rapidity < %f && pt > %d && pt < %d && cosTheta%s > %f && cosTheta%s < %f && fabs(%s) > %d && fabs(%s) < %d", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);

	RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass"))), kinematicCut);

	wspace.import(reducedDataset, Rename(Form("dataset_cosTheta_%.2fto%.2f_absphi_%dto%d", cosThetaMin, cosThetaMax, phiMin, phiMax)));

	return reducedDataset;
}

void nominalFit_RawDataset(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180, Bool_t isPhiFolded = kTRUE, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax, const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ExpTimesErr", const char* muonAccName = "UpsilonTriggerThresholds") {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = Form("../Files/UpsilonSkimmedDataset_%s.root", muonAccName);

	const char* refFrameName = isCSframe ? "CS" : "HX";

	RooWorkspace wspace = SetUpWorkspace(filename, "");

	RooRealVar invMass = *wspace.var("mass");
	invMass.setRange("MassFitRange", massMin, massMax);

	RooRealVar cosTheta = *wspace.var(Form("cosTheta%s", refFrameName));

	RooRealVar phi = *wspace.var(Form("phi%s", refFrameName));

	wspace.Print();

	auto allDataset = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax, "");

	RooDataSet reducedDataset = InvMassDataset(wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, refFrameName, phiMin, phiMax);

	Long64_t nEntries = reducedDataset.sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	// const char* signalShapeName = "SymDSCB";
	// // const char* signalShapeName = "Johnson";

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	// BuildInvariantMassModel(wspace, signalShapeName, bkgShapeName, fitModelName, nEntries, true);
	BuildInvariantMassModel(wspace, signalShapeName, bkgShapeName, fitModelName, nEntries, true, muonAccName); // fix sigma to MC

	RooAddPdf invMassModel = *((RooAddPdf*)wspace.pdf("invMassModel"));
	invMassModel.setNormRange("MassFitRange");

	/// fit!!!
	auto* fitResult = RawInvariantMassFit(wspace, reducedDataset);

	// save the invariant mass distribution fit for further checks
	// one pad for the invariant mass data distribution with fit components, one for the pull distribution
	auto* massCanvas = new TCanvas("massCanvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, .95);

	pad1->SetTopMargin(0.003);
	pad1->SetBottomMargin(0.03);
	pad1->SetRightMargin(0.035);

	pad1->Draw();
	pad1->cd();

	RooPlot* frame = InvariantMassRooPlot(wspace, reducedDataset);
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, 0.67, 0.70, 0.94, 0.95)); // (HX, 2to6, -0.42to-0.14, 60to120): 0.67, 0.70, 0.94, 0.95 // (HX, 12to20, 0.42to0.70, 60to120): 0.63, 0.695, 0.94, 0.95

	if (!isPhiFolded)
		frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	else
		frame->addObject(RefFrameTextPhiFolded(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax, 0.595, 0.47, 0.945, 0.66, 32)); // (HX, 2to6, -0.42to-0.14, 60to120): 0.595, 0.47, 0.945, 0.66, 32 // (HX, 12to20, 0.42to0.70, 60to120): 0.605, 0.47, 0.945, 0.66

	if (strcmp(signalShapeName, "SymDSCB") == 0)
		frame->addObject(FitResultText(*wspace.var("yield1S"), ComputeSignalSignificance(wspace, 1), *wspace.var("yield2S"), ComputeSignalSignificance(wspace, 2), 0.14, 0.07, 0.48, 0.35, 12)); // (HX, 2to6, -0.42to-0.14, 60to120): 0.14, 0.07, 0.48, 0.35, 12 // (HX, 12to20, 0.42to0.70, 60to120): 0.147, 0.70, 0.458, 0.95
	else
		frame->addObject(FitResultText(*wspace.var("yield1S"), 0, *wspace.var("yield2S"), 0));

	frame->Draw();

	gPad->RedrawAxis();

	/// add 10^3 to the top left corner of the canvas
	TLatex* latex = new TLatex();
	latex->SetTextSize(0.052); // Adjust the text size as needed
	//latex->DrawLatexNDC(0.0515, 0.955, "10^{3}#times"); // (x, y) coordinates in normalized device coordinates

	// pull distribution
	massCanvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	massCanvas->cd();
	pad1->Draw();
	pad2->Draw();

	// massCanvas->Modified();
	// massCanvas->Update();

	CMS_lumi(massCanvas, gCMSLumiText);

	// const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooArgSet* signalYields = new RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S"), *wspace.var("yield3S"));

	RooArgSet* fitParams = new RooArgSet(fitResult->floatParsFinal());

	const char* totalFitModelName = GetTotalFitModelName(bkgShapeName, signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	SaveRawSignalYields(signalYields, totalFitModelName, muonAccName);

	SaveRawFitResults(fitResult, totalFitModelName, muonAccName);

	SaveRawDataFitCanvas(massCanvas, totalFitModelName, muonAccName);

	// SaveRawSignalYields(signalYields, totalFitModelName, Form("%s_woSigmaConstraints", gMuonAccName));

	// SaveRawFitResults(fitResult, totalFitModelName, Form("%s_woSigmaConstraints", gMuonAccName));

	// SaveRawDataFitCanvas(massCanvas, totalFitModelName, Form("%s_woSigmaConstraints", gMuonAccName));

	delete massCanvas;
	delete fitResult;
}

void scanNominalFit_RawDataset(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180, const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ExpTimesErr", const char* muonAccName = "UpsilonTriggerThresholds") { // possible bkgShapeName: ExpTimesErr, ChebychevOrderN, Argus, Exponential

	std::vector<Double_t> cosThetaEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	std::vector<Double_t> phiEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// Int_t numPtEle = sizeof(ptEdges) / sizeof(Int_t);

	// for (Int_t ptIdx = 0; ptIdx < NPtBins; ptIdx++) {
	for (Int_t cosThetaIdx = 0; cosThetaIdx < nCosThetaBins; cosThetaIdx++) {
		for (Int_t idx = 0; idx < nPhiBins; idx++) {
			nominalFit_RawDataset(ptMin, ptMax, isCSframe, cosThetaEdges[cosThetaIdx], cosThetaEdges[cosThetaIdx + 1], phiEdges[idx], phiEdges[idx + 1], kTRUE, 7, 11.5, signalShapeName, bkgShapeName, muonAccName);
			// nominalFit_RawDataset(gPtBinning[ptIdx], gPtBinning[ptIdx + 1], isCSframe, cosThetaEdges[cosThetaIdx], cosThetaEdges[cosThetaIdx + 1], phiEdges[idx], phiEdges[idx + 1], kTRUE, MassBinMin, 11.5, signalShapeName, bkgShapeName);
		}
	}
	// }
}
