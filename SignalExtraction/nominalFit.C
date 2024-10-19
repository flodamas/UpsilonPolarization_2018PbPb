#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/RooFitPDFs/ErrorFuncTimesExp.h"
#include "../Tools/Style/FitDistributions.h"

void nominalFit(Int_t ptMin = 0, Int_t ptMax = gPtMax, const char* bkgShapeName = "ExpTimesErr", Float_t massMin = MassBinMin, Float_t massMax = MassBinMax) {
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

	auto invMassModel = MassFitModel(wspace, signalShapeName, bkgShapeName, ptMin, ptMax, nEntries);

	auto* fitResult = RawInvariantMassFit(data, invMassModel, RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S")), massMin, massMax);

	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	const char* fitModelName = GetSignalFitName(signalShapeName, ptMin, ptMax);

	RooArgSet* signalYields = new RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S"), *wspace.var("yield3S"));

	SaveRawDataSignalYields(signalYields, bkgShapeName, fitModelName, gMuonAccName);

	SaveRawDataCanvas(massCanvas, bkgShapeName, fitModelName, gMuonAccName);
}
