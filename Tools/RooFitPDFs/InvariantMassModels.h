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
