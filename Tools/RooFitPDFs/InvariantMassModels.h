// when the header is included several times, to avoid the redefinition error
#ifndef Invariant_Mass_GUARD
#define Invariant_Mass_GUARD
/*

All RooFit PDFs and models used in many macros compiled in a single header!

*/
using namespace RooFit;

#include "../Parameters/PhysicsConstants.h"

#include "./ErrorFuncTimesExp.h"

#include "../FitShortcuts.h"

/// PDFs for Y signal extraction

// one symmetric double-sided Crystal Ball PDF per Y resonance with PDG mass scaling for the mean and width of the excited states
// tail parameters identical for the three resonances

RooAddPdf NominalSignalModel(RooWorkspace& wspace, Long64_t nEntries = 1e6) {
	RooRealVar mass = *wspace.var("mass");

	Long64_t initYield = nEntries / 100;

	// get the tail parameters, assuming that they have been imported to the workspace first!!
	RooRealVar alphaInf = *wspace.var("alphaInfSymDSCB");
	RooRealVar orderInf = *wspace.var("orderInfSymDSCB");
	RooRealVar alphaSup = *wspace.var("alphaSupSymDSCB");
	RooRealVar orderSup = *wspace.var("orderSupSymDSCB");

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.35, 9.55);
	RooRealVar sigma_1S("sigma_1S", "", .04, .02, .15);

	RooCrystalBall signalPDF_1S("signalPDF_1S", "", mass, mean_1S, sigma_1S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar yield1S("yield1S", "N 1S", initYield, 0, nEntries);

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "massScaling_2S*mean_1S", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigma_2S("sigma_2S", "massScaling_2S*sigma_1S", RooArgSet(massScaling_2S, sigma_1S));

	RooCrystalBall signalPDF_2S("signalPDF_2S", "", mass, mean_2S, sigma_2S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar yield2S("yield2S", "N 2S", initYield / 4, 0, nEntries);
	// RooRealVar yield2S("yield2S", "N 2S", 0);

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "massScaling_3S*mean_1S", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigma_3S("sigma_3S", "massScaling_3S*sigma_1S", RooArgSet(massScaling_3S, sigma_1S));

	RooCrystalBall signalPDF_3S("signalPDF_3S", "", mass, mean_3S, sigma_3S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar yield3S("yield3S", "N 3S", initYield / 12, 0, nEntries);
	// RooRealVar yield3S("yield3S", "N 3S", 0);

	RooAddPdf signalModel("SymDSCBModel", "PDF of the sum of the three Y signal PDFs", {signalPDF_1S, signalPDF_2S, signalPDF_3S}, {yield1S, yield2S, yield3S});

	wspace.import(signalModel, RecycleConflictNodes()); // including all signal yield variables!!

	return signalModel;
}

RooArgList ChebychevCoefList(int order = 1) {
	RooArgList coefList;

	for (int i = 0; i <= order; i++) {
		RooRealVar* coef = new RooRealVar(Form("ChebychevCoef_%d", i), " ", 0.1, -10, 10);
		coefList.addOwned(*coef);
	}

	return coefList;
}

RooAddPdf BackgroundModel(RooWorkspace& wspace, const char* bkgShapeName, Long64_t nEntries = 1e6) {
	RooRealVar invMass = *wspace.var("mass");

	RooRealVar yieldBkg("yieldBkg", "N background events", 0, nEntries);

	// Chebychev Nth order polynomial
	if (strncmp(bkgShapeName, "ChebychevOrder", 14) == 0) {
		int charSize = strlen(bkgShapeName);
		const char* lastDigit = bkgShapeName + charSize - 1;
		std::stringstream converter(lastDigit);

		int order = 0;
		converter >> order;

		RooArgList coefList = ChebychevCoefList(order);
		RooChebychev bkgPDF("bkgPDF", " ", invMass, coefList);

		RooAddPdf bkgModel("bkgModel", "PDF of the backgroud", {bkgPDF}, {yieldBkg});

		wspace.import(bkgModel, RecycleConflictNodes());

		return bkgModel;
	}

	// exponential x err function
	else if (strcmp(bkgShapeName, "ExpTimesErr") == 0) {
<<<<<<< HEAD
		RooRealVar err_mu("err_mu", " ", 7, 5, 13);
		RooRealVar err_sigma("err_sigma", " ", 1, 0, 10);
		RooRealVar exp_lambda("exp_lambda", " ", 1.5, 0, 10);
=======
		RooRealVar err_mu("err_mu", " ", 0, 13);
		RooRealVar err_sigma("err_sigma", " ", 1, 0, 10);
		RooRealVar exp_lambda("exp_lambda", " ", 10, 0, 500);
>>>>>>> b3346481d0082db909096b565c02f49b6bee13d9

		ErrorFuncTimesExp bkgPDF("bkgPDF", " ", invMass, err_mu, err_sigma, exp_lambda);

		RooAddPdf bkgModel("bkgModel", "PDF of the backgroud", {bkgPDF}, {yieldBkg});

		wspace.import(bkgModel, RecycleConflictNodes());

		return bkgModel;
	}

	else if (strcmp(bkgShapeName, "Exponential") == 0) {
		RooRealVar exp_lambda("exp_lambda", " ", -10, 10);

		RooExponential bkgPDF("bkgPDF", " ", invMass, exp_lambda);

		RooAddPdf bkgModel("bkgModel", "PDF of the backgroud", {bkgPDF}, {yieldBkg});

		wspace.import(bkgModel, RecycleConflictNodes());

		return bkgModel;
	}

	else {
		std::cout << "No matching background model" << std::endl;
		return nullptr;
	}
}

// build the main PDF based on all PDFs above

RooAddPdf MassFitModel(RooWorkspace& wspace, const char* signalShapeName, const char* bkgShapeName, Int_t ptMin = 0, Int_t ptMax = 30, Long64_t nEntries = 1e6) {
	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	if (BeVerbose) std::cout << "\nBuilding invariant mass fit model with " << signalShapeName << " for the signal and " << bkgShapeName << " for the background\n";

	// tail parameters fixed to MC extracted values, and identical for the three resonances

	ImportAndFixMCSignalParameters(wspace, signalShapeName, ptMin, ptMax);

	auto signalModel = NominalSignalModel(wspace, nEntries);

	RooAbsPdf* signalPDF_1S = wspace.pdf("signalPDF_1S");
	RooAbsPdf* signalPDF_2S = wspace.pdf("signalPDF_2S");
	RooAbsPdf* signalPDF_3S = wspace.pdf("signalPDF_3S");

	RooRealVar* yield1S = wspace.var("yield1S");
	RooRealVar* yield2S = wspace.var("yield2S");
	RooRealVar* yield3S = wspace.var("yield3S");

	// background
	auto bkgModel = BackgroundModel(wspace, bkgShapeName, nEntries);

	RooAbsPdf* bkgPDF = wspace.pdf("bkgPDF");

	RooRealVar* yieldBkg = wspace.var("yieldBkg");

	RooAddPdf model("invMassModel", "", {*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, *bkgPDF}, {*yield1S, *yield2S, *yield3S, *yieldBkg});

	wspace.import(model, RecycleConflictNodes());

	if (BeVerbose) std::cout << "\nModel exported to " << wspace.GetName() << std::endl;

	return model;
}

#endif
