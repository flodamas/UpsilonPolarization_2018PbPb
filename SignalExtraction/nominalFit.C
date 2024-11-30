#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/RooFitPDFs/ErrorFuncTimesExp.h"
#include "../Tools/Style/FitDistributions.h"

void nominalFit(Int_t ptMin = 0, Int_t ptMax = gPtMax, const char* bkgShapeName = "ExpTimesErr", Float_t massMin = MassBinMin, Float_t massMax = MassBinMax, const char* refFrameName = "HX", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = Form("../Files/UpsilonSkimmedDataset%s.root", gMuonAccName);

	RooWorkspace wspace = SetUpWorkspace(filename);

	RooDataSet data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax);

	RooRealVar mass = *wspace.var("mass");
	mass.setRange("MassFitRange", massMin, massMax);

	Long64_t nEntries = data.sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	const char* signalShapeName = "SymDSCB";

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	BuildInvariantMassModel(wspace, signalShapeName, bkgShapeName, fitModelName, nEntries);

	auto* fitResult = RawInvariantMassFit(wspace, data, RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S")), massMin, massMax);

	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	//const char* fitModelName = GetSignalFitName(signalShapeName, ptMin, ptMax);

	RooArgSet* signalYields = new RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S"), *wspace.var("yield3S"));

	const char* totalFitModelName = GetTotalFitModelName(bkgShapeName, signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	SaveRawSignalYields(signalYields, totalFitModelName);

	SaveRawDataFitCanvas(massCanvas, totalFitModelName);
}
