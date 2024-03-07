/*

All RooFit PDFs and models used in many macros compiled in a single header!

*/
using namespace RooFit;

#include "../Parameters/PhysicsConstants.h"

#include "ErrorFuncTimesExp.h"

/// PDFs for Y signal extraction

// one symmetric double-sided Crystal Ball PDF per Y resonance with PDG mass scaling for the mean and width of the excited states
// tail parameters fixed to MC extracted values, and identical for the three resonances

RooAddPdf NominalSignalModel(RooWorkspace& wspace, RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, Long64_t yieldMax = 1e6) {
	RooRealVar mass = *wspace.var("mass");

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.35, 9.55);
	RooRealVar sigma_1S("sigma_1S", "", .04, .13);

	RooCrystalBall signalPDF_1S("signalPDF_1S", "", mass, mean_1S, sigma_1S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar yield1S("yield1S", "N 1S", yieldMax / 5, 0, yieldMax);

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "massScaling_2S*mean_1S", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigma_2S("sigma_2S", "massScaling_2S*sigma_1S", RooArgSet(massScaling_2S, sigma_1S));

	RooCrystalBall signalPDF_2S("signalPDF_2S", "", mass, mean_2S, sigma_2S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar yield2S("yield2S", "N 2S", yieldMax / 10, 0, yieldMax / 2);

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "massScaling_3S*mean_1S", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigma_3S("sigma_3S", "massScaling_3S*sigma_1S", RooArgSet(massScaling_3S, sigma_1S));

	RooCrystalBall signalPDF_3S("signalPDF_3S", "", mass, mean_3S, sigma_3S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar yield3S("yield3S", "N 3S", yieldMax / 20, 0, yieldMax / 4);

	RooAddPdf signalModel("SymDSCBModel", "PDF of the sum of the three Y signal PDFs", {signalPDF_1S, signalPDF_2S, signalPDF_3S}, {yield1S, yield2S, yield3S});

	wspace.import(signalModel); // including all signal yield variables!!

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

RooAddPdf BackgroundModel(RooWorkspace& wspace, const char* bkgShapeName, Long64_t yieldMax = 1e6) {
	RooRealVar invMass = *wspace.var("mass");

	RooRealVar yieldBkg("yieldBkg", "N background events", 0, yieldMax);

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
		RooRealVar err_mu("err_mu", " ", 7, 0, 15);
		RooRealVar err_sigma("err_sigma", " ", 0, 10);
		RooRealVar exp_lambda("exp_lambda", " ", 3, 0, 20);

		ErrorFuncTimesExp bkgPDF("bkgPDF", " ", invMass, err_mu, err_sigma, exp_lambda);

		RooAddPdf bkgModel("bkgModel", "PDF of the backgroud", {bkgPDF}, {yieldBkg});

		wspace.import(bkgModel, RecycleConflictNodes());

		return bkgModel;
	}

	else {
		cout << "No matching background model" << endl;
		return nullptr;
	}
}

// build the main PDF based on all PDFs above

RooAddPdf* MassFitModel(RooWorkspace& wspace, const char* signalShapeName, const char* bkgShapeName, Int_t ptMin = 0, Int_t ptMax = 30, Long64_t yieldMax = 1e6) {
	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	cout << endl
	     << "Building invariant mass fit model with " << signalShapeName << " for the signal and " << bkgShapeName << " for the background" << endl;

	RooRealVar* alphaInf = new RooRealVar("alphaInf", "", 1);
	RooRealVar* orderInf = new RooRealVar("orderInf", "", 1);
	RooRealVar* alphaSup = new RooRealVar("alphaSup", "", 1);
	RooRealVar* orderSup = new RooRealVar("orderSup", "", 1);

	RooArgSet tailParams = GetMCSignalTailParameters(alphaInf, orderInf, alphaSup, orderSup, signalShapeName, ptMin, ptMax);

	auto signalModel = NominalSignalModel(wspace, alphaInf, orderInf, alphaSup, orderSup, yieldMax);

	RooAbsPdf* signalPDF_1S = wspace.pdf("signalPDF_1S");
	RooAbsPdf* signalPDF_2S = wspace.pdf("signalPDF_2S");
	RooAbsPdf* signalPDF_3S = wspace.pdf("signalPDF_3S");

	RooRealVar* yield1S = wspace.var("yield1S");
	RooRealVar* yield2S = wspace.var("yield2S");
	RooRealVar* yield3S = wspace.var("yield3S");

	// background
	auto bkgModel = BackgroundModel(wspace, bkgShapeName, yieldMax);

	RooAbsPdf* bkgPDF = wspace.pdf("bkgPDF");

	RooRealVar* yieldBkg = wspace.var("yieldBkg");

	RooAddPdf* model = new RooAddPdf("invMassModel", "", {*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, *bkgPDF}, {*yield1S, *yield2S, *yield3S, *yieldBkg});

	wspace.import(*model, RecycleConflictNodes());

	cout << endl
	     << "Model exported to " << wspace.GetName() << endl;

	return model;
}
