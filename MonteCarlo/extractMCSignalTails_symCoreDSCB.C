
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/Shortcuts.h"

// crystal ball shape with symmetric Gaussian core and asymmetric tails (just like RooDSCBShape)

void extractMCSignalTails_symCoreDSCB(Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
	const char* filename = "../Files/MCUpsilonSkimmedWeightedDataset.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}

	cout << "File " << filename << " opened" << endl;

	// we only extract the tail parameters for a specific pt bin, they do not vary significantly with cos theta or phi within this pt bin
	Bool_t isCSframe = kTRUE;
	Float_t cosThetaMin = -1, cosThetaMax = 1;
	Float_t phiMin = -180, phiMax = 180;

	Float_t massMin = 8.5, massMax = 10.5;
	Int_t nBins = 80;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)file->Get("MCdataset");

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, centMin, centMax, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar* massVar = wspace->var("mass");

	// fit
	RooRealVar mean("mean", "", 9., 10.);
	RooRealVar sigma("sigma", "", .05, .15);
	RooRealVar alphaInf("alphaInf", "", 0.1, 10);
	RooRealVar orderInf("orderInf", "", 0.1, 10);
	RooRealVar alphaSup("alphaSup", "", 0.1, 10);
	RooRealVar orderSup("orderSup", "", 0.1, 40);

	RooCrystalBall signal("CB", "", *massVar, mean, sigma, alphaInf, orderInf, alphaSup, orderSup);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with a double-sided Crystal Ball PDF made of a symmetric Gaussian core and asymmetric tail distributions..." << endl;

	bool doWeightedError = true;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true), PrintLevel(-1), Minos(!doWeightedError), NumCPU(3), Range(massMin, massMax), AsymptoticError(doWeightedError)); // quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit

	fitResult->Print("v");

	TString outputName = Form("symCoreDSCB_cent%dto%d_pt%dto%d", centMin, centMax, ptMin, ptMax);

	// save signal shape parameters in a txt file to be read for data fit

	SaveMCSignalTailParameters(RooArgSet(alphaInf, orderInf, alphaSup, orderSup), outputName.Data()); // so that we don't have to refit later

	file->Close();
}
