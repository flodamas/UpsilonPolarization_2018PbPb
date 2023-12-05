//#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"

void weightedInvMass(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* signalShapeName = "SymDSCB";

	// get the tail parameters of the signal shape first in case the MC fit is needed
	RooRealVar* alphaInf = new RooRealVar("alphaInf", "", 1);
	RooRealVar* orderInf = new RooRealVar("orderInf", "", 1);
	RooRealVar* alphaSup = new RooRealVar("alphaSup", "", 1);
	RooRealVar* orderSup = new RooRealVar("orderSup", "", 1);

	RooArgSet tailParams = GetMCSignalTailParameters(alphaInf, orderInf, alphaSup, orderSup, signalShapeName, ptMin, ptMax);

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

	RooRealVar massVar = *wspace.var("mass");

	Long64_t nEntries = massDataset->sumEntries();

	/// fitting model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	auto signalModel = NominalSignalModel(wspace, alphaInf, orderInf, alphaSup, orderSup, nEntries);

	// get variables
	RooAbsPdf* signalPDF_1S = wspace.pdf("signalPDF_1S");
	RooAbsPdf* signalPDF_2S = wspace.pdf("signalPDF_2S");
	RooAbsPdf* signalPDF_3S = wspace.pdf("signalPDF_3S");

	RooRealVar* yield1S = wspace.var("yield1S");
	RooRealVar* yield2S = wspace.var("yield2S");
	RooRealVar* yield3S = wspace.var("yield3S");

	// background: error function x exponential
	RooRealVar err_mu("err_mu", "err_mu", 0, 10);
	RooRealVar err_sigma("err_sigma", "err_sigma", 0, 10);
	RooRealVar exp_lambda("exp_lambda", "m_lambda", 0, 10);

	ErrorFuncTimesExp bkgPDF("bkgPDF", "", massVar, err_mu, err_sigma, exp_lambda);
	RooRealVar yieldBkg("yieldBkg", "N background events", 0, nEntries);

	RooAddPdf fitModel("fitModel", "", {*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, bkgPDF}, {*yield1S, *yield2S, *yield3S, yieldBkg});

	auto* fitResult = fitModel.fitTo(*massDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

	fitResult->Print("v");

	wspace.import(fitModel, RecycleConflictNodes());

	// draw
	TCanvas* canvas = DrawMassFitDistributions(wspace, massDataset, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	canvas->SaveAs(Form("mass_distrib/weightedEvents_%s.png", fitModelName), "RECREATE");
}
