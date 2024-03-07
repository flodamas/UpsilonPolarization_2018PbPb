//#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"

void weightedInvMass(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)f->Get(isCSframe ? "datasetCS" : "datasetHX");

	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooDataSet* massDataset = (isCSframe) ? ReducedWeightedMassDatasetCS(allDataset, wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax) : ReducedWeightedMassDatasetHX(allDataset, wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax);

	Long64_t nEntries = massDataset->sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	const char* signalShapeName = "SymDSCB";

	// background
	//int order = 2;
	//const char* bkgShapeName = Form("ChebychevOrder%d", order);
	const char* bkgShapeName = "ExpTimesErr";

	auto* invMassModel = MassFitModel(wspace, signalShapeName, bkgShapeName, ptMin, ptMax, nEntries);

	auto* fitResult = invMassModel->fitTo(*massDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

	fitResult->Print("v");

	// draw
	TCanvas* canvas = DrawMassFitDistributions(wspace, massDataset, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	canvas->SaveAs(Form("mass_distrib/weightedEvents_%s.png", fitModelName), "RECREATE");
}
