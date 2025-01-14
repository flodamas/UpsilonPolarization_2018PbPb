// when the header is included several times, to avoid the redefinition error
#ifndef Fit_Shortcuts_GUARD
#define Fit_Shortcuts_GUARD

#include "../Tools/BasicHeaders.h"

#include "../Tools/Parameters/PhysicsConstants.h"

#include "../AnalysisParameters.h"

using namespace RooFit;

// for data
RooFitResult* RawInvariantMassFit(RooWorkspace& wspace, RooDataSet data, RooArgSet varsForMinos = RooArgSet(), float massMin = MassBinMin, float massMax = MassBinMax) {
	auto model = wspace.pdf("invMassModel");

	RooRealVar* invMass = wspace.var("mass");
	invMass->setRange("MassFitRange", MassBinMin, MassBinMax);

	if (model == nullptr) {
		std::cerr << "\n[FitShortcuts] the workspace does not contain the proper invariant mass model for fitting!\n";
		return nullptr;
	}

	if (BeVerbose) cout << "\nFitting the raw invariant mass distribution...\n";

	//model.Print("v");
	cout << "varsForMinos: " << varsForMinos.getSize() << endl;

	auto* fitResult = model->fitTo(data, Save(), PrintLevel(-1), NumCPU(NCPUs), Range("MassFitRange"), Minos(varsForMinos), Extended(true));
	// only run Minos over the parsed variables

	if (BeVerbose) fitResult->Print("v");

	return fitResult;
}

RooFitResult* WeightedInvariantMassFit(RooDataSet data, RooAddPdf model, float massMin = MassBinMin, float massMax = MassBinMax) {
	if (BeVerbose) cout << "\nFitting the invariant mass distribution with weighted entries...\n\n";

	bool doWeightedError = true;

	auto* fitResult = model.fitTo(data, Save(), Extended(true) /*, PrintLevel(-1)*/, Minos(!doWeightedError), NumCPU(NCPUs), Range(massMin, massMax), AsymptoticError(doWeightedError));

	if (BeVerbose) fitResult->Print("v");

	return fitResult;
}

// for MC
RooFitResult* MCWeightedInvariantMassFit(RooDataSet data, RooAbsPdf& model, float massMin = MassBinMin, float massMax = MassBinMax) {
	if (BeVerbose) cout << "\nFitting the MC invariant mass distribution with weighted entries...\n\n";

	auto* fitResult = model.fitTo(data, Save(), PrintLevel(-1), Minos(!DoMCWeightedError), NumCPU(NCPUs), Range(massMin, massMax), AsymptoticError(DoMCWeightedError));
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit
	// stratiges
	if (BeVerbose) fitResult->Print("v");

	return fitResult;
}

RooFitResult* SymDSCBfit(RooWorkspace& wspace, RooDataSet data, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax) {
	// fit
	RooRealVar mean("meanSymDSCB", "", PDGmass_1S, 9., 10.);
	RooRealVar sigma("sigmaSymDSCB", "", 0.08, .03, .15);
	RooRealVar alphaInf("alphaInfSymDSCB", "", 1.14, 0.1, 10);
	RooRealVar orderInf("orderInfSymDSCB", "", 5.0, 0.1, 20);
	RooRealVar alphaSup("alphaSupSymDSCB", "", 1.3, 0.1, 10);
	RooRealVar orderSup("orderSupSymDSCB", "", 33, 0.1, 200);
	// RooRealVar orderSup("orderSupSymDSCB", "", 33.284);

	RooCrystalBall signal("SymDSCB", "SymDSCB", *wspace.var("mass"), mean, sigma, alphaInf, orderInf, alphaSup, orderSup);

	if (BeVerbose) cout << "\nFitting the MC signal shape (weighted entries!!) with a double-sided Crystal Ball PDF made of a symmetric Gaussian core and asymmetric tail distributions...\n";

	auto* fitResult = MCWeightedInvariantMassFit(data, signal, massMin, massMax);

	wspace.import(signal);

	return fitResult;
}

RooFitResult* AsymDSCBfit(RooRealVar* massVar, RooWorkspace* wspace, RooDataSet* data, Float_t massMin, Float_t massMax) {
	// fit
	RooRealVar mean("meanAsymDSCB", "", PDGmass_1S, 9., 10.);
	RooRealVar sigmaInf("sigmaInfAsymDSCB", "", 0.08, .05, .15);
	RooRealVar alphaInf("alphaInfAsymDSCB", "", 1.5, 0.1, 10);
	RooRealVar orderInf("orderInfAsymDSCB", "", 1.5, 0.1, 10);
	RooRealVar sigmaSup("sigmaSupAsymDSCB", "", 0.08, .05, .15);
	RooRealVar alphaSup("alphaSupAsymDSCB", "", 1.5, 0.1, 10);
	RooRealVar orderSup("orderSupAsymDSCB", "", 3, 0.1, 40);

	RooCrystalBall signal("AsymDSCB", "AsymDSCB", *massVar, mean, sigmaInf, sigmaSup, alphaInf, orderInf, alphaSup, orderSup);

	if (BeVerbose) cout << endl
		                  << "Fitting the MC signal shape (weighted entries!!) with a double-sided Crystal Ball PDF made of an asymmetric Gaussian core and asymmetric tail distributions..." << endl;

	bool doWeightedError = true;

	auto* fitResult = signal.fitTo(*data, Save(), Extended(true) /*, PrintLevel(-1)*/, Minos(!doWeightedError), NumCPU(3), Range(massMin, massMax), AsymptoticError(doWeightedError));
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit
	wspace->import(signal);

	if (BeVerbose) fitResult->Print("v");

	return fitResult;
}

RooFitResult* SymDSCBGaussfit(RooWorkspace& wspace, RooDataSet data, Float_t massMin, Float_t massMax) {
	/// fit
	/// (DSCB variables)
	RooRealVar mean("meanDSCBGauss", "", PDGmass_1S, 9., 10.);
	RooRealVar sigma("sigmaDSCBGauss", "", 0.08, .05, .15);
	RooRealVar alphaInf("alphaInfDSCBGauss", "", 1.5, 0.1, 10);
	RooRealVar orderInf("orderInfDSCBGauss", "", 1.5, 0.1, 10);
	RooRealVar alphaSup("alphaSupDSCBGauss", "", 1.5, 0.1, 10);
	RooRealVar orderSup("orderSupDSCBGauss", "", 3, 0.1, 40);

	/// (Gaussian variable (used the same mean as DSCB))
	RooRealVar sigma_gauss("sigma_gauss", "", 0.11, .05, .3);

	/// (fraction between two PDFs)
	RooRealVar normFraction("normFraction", "", 0.6, 0.01, 1);

	RooCrystalBall DSCB("DSCB", "DSCB", *wspace.var("mass"), mean, sigma, alphaInf, orderInf, alphaSup, orderSup);
	RooGaussian gauss("gauss", "gaussian", *wspace.var("mass"), mean, sigma_gauss);

	RooAddPdf signal("DSCBGauss", "sum of DSCB and CB PDF", RooArgList(DSCB, gauss), RooArgList(normFraction), kTRUE);

	if (BeVerbose) cout << "\nFitting the MC signal shape (weighted entries!!) with the sum of a double-sided Crystal Ball PDF made of a symmetric Gaussian core and asymmetric tail distributions, and a Gaussian PDF..." << endl;

	auto* fitResult = MCWeightedInvariantMassFit(data, signal, massMin, massMax);

	wspace.import(signal);

	return fitResult;
}

RooFitResult* HypatiaFit(RooWorkspace& wspace, RooDataSet data, Float_t massMin, Float_t massMax) {
	// fit
	RooRealVar mean("meanHypatia", "", PDGmass_1S, 9., 10.);
	RooRealVar lambda("lambdaHypatia", "lambda of hypatia PDF", -1.0, -10.0, 10);
	RooRealVar zeta("zetaHypatia", "zeta of hypatia PDF", 0.01, 0.0, 1.0);
	RooRealVar beta("betaHypatia", "beta of hypatia PDF", -0.01, -20.0, 0.0);
	RooRealVar sigma("sigmaHypatia", "sigma of hypatia PDF", 0.15, 0.03, 0.3);
	RooRealVar alphaInf("alphaInfHypatia", "al1s of hypatia PDF", 1.5, 0.1, 10.0);
	RooRealVar alphaSup("alphaSupHypatia", "ar1s of hypatia PDF", 1.5, 0.1, 10.0);
	RooRealVar orderInf("orderInfHypatia", "nl1s of hypatia PDF", 1.0, 0.2, 10.0);
	RooRealVar orderSup("orderSupHypatia", "nr1s of hypatia PDF", 5., 0.1, 100);

	//mean.setConstant();
	sigma.setConstant();
	zeta.setVal(0);
	zeta.setConstant();
	beta.setVal(0);
	beta.setConstant();
	/*
	lambda.setVal(-1.729);
	zeta.setVal(0);
	beta.setVal(-1.59);
	sigma.setVal(0.158);
	alphaInf.setVal(1.769);
	alphaSup.setVal(5.398);
	orderInf.setVal(1.422);
	orderSup.setVal(0.012);
*/
	//zeta.setConstant(1);

	RooHypatia2 signal("Hypatia", "Hypatia", *wspace.var("mass"), lambda, zeta, beta, sigma, mean, alphaInf, orderInf, alphaSup, orderSup);

	if (BeVerbose) cout << "\nFitting the MC signal shape (weighted entries!!) with a Hypatia PDF" << endl;

	auto* fitResult = MCWeightedInvariantMassFit(data, signal, massMin, massMax);

	wspace.import(signal);

	return fitResult;
}

RooFitResult* VoigtianFit(RooWorkspace& wspace, RooDataSet data, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax) {
	// fit
	RooRealVar mean("meanVoigtian", "", PDGmass_1S, 9., 10.);
	RooRealVar sigma("sigmaVoigtian", "", 0.08, .03, .15);
	RooRealVar width("widthVoigtian", "", 0.08, .03, .15);

	RooVoigtian signal("Voigtian", "Voigtian", *wspace.var("mass"), mean, sigma, width);

	if (BeVerbose) cout << "\nFitting the MC signal shape (weighted entries!!) with a Voigtian PDF...\n";

	auto* fitResult = MCWeightedInvariantMassFit(data, signal, massMin, massMax);

	wspace.import(signal);

	return fitResult;
}

RooFitResult* JohnsonFit(RooWorkspace& wspace, RooDataSet data, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax) {
	// fit
	RooRealVar mean("meanJohnson", "", PDGmass_1S, 9.3, 9.6);
	RooRealVar lambda("lambdaJohnson", "", 0.1, 0.001, 5);
	RooRealVar gamma("gammaJohnson", "", 0.4, -5, 5);
	RooRealVar delta("deltaJohnson", "", 1, 0.01, 10);

	RooJohnson signal("Johnson", "Johnson", *wspace.var("mass"), mean, lambda, gamma, delta);

	if (BeVerbose) cout << "\nFitting the MC signal shape (weighted entries!!) with a Johnson PDF...\n";

	auto* fitResult = MCWeightedInvariantMassFit(data, signal, massMin, massMax);

	wspace.import(signal);

	return fitResult;
}

/// Helpers to get the relevant signal shape parameters and to define an unique fit name

const char* GetSignalFitName(const char* signalShapeName = "SymDSCB", Int_t ptMin = 0, Int_t ptMax = 30) {
	return Form("%s_cent%dto%d_absy%dp%dto%dp%d_pt%dto%d", signalShapeName, gCentralityBinMin, gCentralityBinMax, (Int_t)gRapidityMin, (Int_t)(10 * (gRapidityMin - (Int_t)gRapidityMin)), (Int_t)gRapidityMax, (Int_t)(10 * (gRapidityMax - (Int_t)gRapidityMax)), ptMin, ptMax);
}

const char* GetFitModelName(const char* signalShapeName = "SymDSCB", Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "HX", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* signalFitName = GetSignalFitName(signalShapeName, ptMin, ptMax);

	// just need to append the specific (cos theta, phi) bin name

	return Form("%s_cosTheta%.2fto%.2f_phi%dto%d_%s", signalFitName, cosThetaMin, cosThetaMax, phiMin, phiMax, refFrameName);
}

const char* GetTotalFitModelName(const char* bkgShapeName = "ExpTimesErr", const char* signalShapeName = "SymDSCB", Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "HX", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	// total fit name = bkg + signal names (in this order)

	return Form("%s_%s", bkgShapeName, fitModelName);
}

const char** GetFitModelNames(const char* signalShapeName = "SymDSCB", Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Int_t nCosThetaBins = 10, const vector<Double_t>& cosThetaBinEdges = {}, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char** signalFitNames = new const char*[nCosThetaBins];

	// just need to append the specific (cos theta, phi) bin name
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		signalFitNames[iCosTheta] = Form("%s_cosTheta%.2fto%.2f_phi%dto%d_%s", GetSignalFitName(signalShapeName, ptMin, ptMax), cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], phiMin, phiMax, (isCSframe) ? "CS" : "HX");
	}

	return signalFitNames;
}

const char** GetFitModelNames(const char* signalShapeName = "SymDSCB", Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Double_t cosThetaMin = -1, Double_t cosThetaMax = 1, Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {}) {
	const char** signalFitNames = new const char*[nPhiBins];

	// just need to append the specific (cos theta, phi) bin name
	for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
		signalFitNames[iPhi] = Form("%s_cosTheta%.2fto%.2f_phi%dto%d_%s", GetSignalFitName(signalShapeName, ptMin, ptMax), cosThetaMin, cosThetaMax, (Int_t)phiBinEdges[iPhi], (Int_t)phiBinEdges[iPhi + 1], (isCSframe) ? "CS" : "HX");
	}

	return signalFitNames;
}

// small dummy function to avoid repetiting the same piece of code everywhere...
const char* GetMCFileName(const char* fitModelName = "SymDSCB_cent0to90_absy0p0to2p4_pt6to12") {
	return Form("../MonteCarlo/SignalParameters/absPhi/%s.txt", fitModelName);
}

// Define and read the signal parameters from MC extracted values, define Gaussian constraints on the parameters that are set to constants

RooArgList* ListOfSignalContraints(RooWorkspace& wspace, const char* signalShapeName = "SymDSCB", const char* fitModelName = "SymDSCB_cent0to90_absy0p0to2p4_pt6to12", bool fixSigmaToMC = false) {
	
	RooArgList* constraintsList = new RooArgList();

	if (strcmp(signalShapeName, "SymDSCB") == 0) {
		// Declare the parameter variables
		RooRealVar* sigma = new RooRealVar(Form("sigma%s", signalShapeName), "", 1);
		RooRealVar* alphaInf = new RooRealVar(Form("alphaInf%s", signalShapeName), "", 1);
		RooRealVar* orderInf = new RooRealVar(Form("orderInf%s", signalShapeName), "", 1);
		RooRealVar* alphaSup = new RooRealVar(Form("alphaSup%s", signalShapeName), "", 1);
		RooRealVar* orderSup = new RooRealVar(Form("orderSup%s", signalShapeName), "", 1);

		RooArgSet parametersList(*sigma, *alphaInf, *orderInf, *alphaSup, *orderSup);

		// if the .txt file for this specific fit model exists, just read the tail parameters from it

		const char* mcFileName = GetMCFileName(fitModelName);

		if (fopen(mcFileName, "r")) {
			if (BeVerbose) cout << "\nFound " << mcFileName << " file, will read the signal parameters from it\n";
			parametersList.readFromFile(mcFileName);
		} else {
			cout << endl
				<< mcFileName << " file does not seem to exist, you need to extract the signal paramaters from MC fit first!" << endl;
			exit(1);
		}

		// fix the tail parameters
		if (fixSigmaToMC) sigma->setConstant();
		alphaInf->setConstant();
		orderInf->setConstant();
		alphaSup->setConstant();
		orderSup->setConstant();

		if (BeVerbose) {
			cout << "\nSignal parameters fixed to the following MC values:\n";
			parametersList.Print("v");
		}
		wspace.import(parametersList);

		// build the Gaussian constraints on parameters to propagate the corresponding uncertainties
		RooGaussian* alphaInfConstraint = new RooGaussian("alphaInfConstraint", "Gaussian constraint on alphaInf", *alphaInf, alphaInf->getVal(), alphaInf->getError());
		RooGaussian* orderInfConstraint = new RooGaussian("orderInfConstraint", "Gaussian constraint on orderInf", *orderInf, orderInf->getVal(), orderInf->getError());
		RooGaussian* alphaSupConstraint = new RooGaussian("alphaSupConstraint", "Gaussian constraint on alphaSup", *alphaSup, alphaSup->getVal(), alphaSup->getError());
		RooGaussian* orderSupConstraint = new RooGaussian("orderSupConstraint", "Gaussian constraint on orderSup", *orderSup, orderSup->getVal(), orderSup->getError());

		if (fixSigmaToMC) {
			RooGaussian* sigmaConstraint = new RooGaussian("sigmaConstraint", "Gaussian constraint on width parameter", *sigma, sigma->getVal(), sigma->getError());
			constraintsList->add(*sigmaConstraint);
		}
		constraintsList->add(*alphaInfConstraint);
		constraintsList->add(*orderInfConstraint);
		constraintsList->add(*alphaSupConstraint);
		constraintsList->add(*orderSupConstraint);
	}

	else if (strcmp(signalShapeName, "Johnson") == 0) {
		// Declare the parameter variables
		RooRealVar* lambda = new RooRealVar(Form("lambda%s", signalShapeName), "", 1);
		RooRealVar* gamma = new RooRealVar(Form("gamma%s", signalShapeName), "", 1);
		RooRealVar* delta = new RooRealVar(Form("delta%s", signalShapeName), "", 1);

		RooArgSet parametersList(*lambda, *gamma, *delta);

		// if the .txt file for this specific fit model exists, just read the tail parameters from it

		const char* mcFileName = GetMCFileName(fitModelName);

		if (fopen(mcFileName, "r")) {
			if (BeVerbose) cout << "\nFound " << mcFileName << " file, will read the signal parameters from it\n";
			parametersList.readFromFile(mcFileName);
		} else {
			cout << endl
			     << mcFileName << " file does not seem to exist, you need to extract the signal paramaters from MC fit first!" << endl;
			exit(1);
		}

		// fix the tail parameters
		lambda->setConstant();
		gamma->setConstant();
		delta->setConstant();

		if (BeVerbose) {
			cout << "\nSignal parameters fixed to the following MC values:\n";
			parametersList.Print("v");
		}
		wspace.import(parametersList);

		// build the Gaussian constraints on parameters to propagate the corresponding uncertainties
		RooGaussian* lambdaConstraint = new RooGaussian("lambdaConstraint", "Gaussian constraint on lambda parameter", *lambda, lambda->getVal(), lambda->getError());
		RooGaussian* gammaConstraint = new RooGaussian("gammaConstraint", "Gaussian constraint on gamma parameter", *gamma, gamma->getVal(), gamma->getError());
		RooGaussian* deltaConstraint = new RooGaussian("deltaConstraint", "Gaussian constraint on delta parameter", *delta, delta->getVal(), delta->getError());

		constraintsList->add(*lambdaConstraint);
		constraintsList->add(*gammaConstraint);
		constraintsList->add(*deltaConstraint);
	}

	wspace.import(*constraintsList, RecycleConflictNodes());

	cout << "\nList of Gaussian constraints:\n";

	constraintsList->Print("v");

	//wspace.Print("v");

	return constraintsList;
}

/*
RooArgSet GetMCSignalTailParameters(RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, const char* signalShapeName = "SymDSCB", Int_t ptMin = 0, Int_t ptMax = 30) {
	RooArgSet parametersList(*alphaInf, *orderInf, *alphaSup, *orderSup);

	// if the .txt file for this specific fit model exists, just read the tail parameters from it

	const char* mcFileName = GetMCFileName(signalShapeName, ptMin, ptMax);

	if (fopen(mcFileName, "r")) {
		if (BeVerbose) cout << endl
			                  << "Found " << mcFileName << " file, will read the signal tail parameters from it" << endl;
		parametersList.readFromFile(mcFileName);
	} else {
		cout << endl
		     << mcFileName << " file does not seem to exist, you need to extract the signal tail paramaters from MC fit first!" << endl;
		exit(1);
	}

	// fix the tail parameters
	alphaInf->setConstant();
	orderInf->setConstant();
	alphaSup->setConstant();
	orderSup->setConstant();

	if (BeVerbose) {
		cout << endl
		     << "Tail parameters fixed to the following MC signal values:" << endl;
		parametersList.Print("v");
	}
	return parametersList;
}

RooArgSet GetMCSignalParameters(RooRealVar* sigma, RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, RooRealVar* sigma_gauss, RooRealVar* normFraction, const char* signalShapeName = "SymDSCB", Int_t ptMin = 0, Int_t ptMax = 30) {
	// RooArgSet GetMCSignalParameters(RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, RooRealVar* sigma_gauss, TString signalShapeName = "SymDSCB", Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
	// RooArgSet GetMCSignalParameters(RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, RooFormulaVar* ratio_sigma, TString signalShapeName = "SymDSCB", Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
	// RooArgSet GetMCSignalParameters(RooRealVar* sigma, RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, RooRealVar* sigma_gauss, TString signalShapeName = "SymDSCB", Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {

	RooArgSet Params(*sigma, *alphaInf, *orderInf, *alphaSup, *orderSup, *sigma_gauss, *normFraction);
	// RooArgSet Params(*alphaInf, *orderInf, *alphaSup, *orderSup, *sigma_gauss);
	// RooArgSet Params(*alphaInf, *orderInf, *alphaSup, *orderSup, *ratio_sigma);
	// RooArgSet Params(*sigma, *alphaInf, *orderInf, *alphaSup, *orderSup, *sigma_gauss);

	// if the .txt file for this specific fit model exists, just read the tail parameters from it

	const char* mcFileName = GetMCFileName(signalShapeName, ptMin, ptMax);

	if (fopen(mcFileName, "r")) {
		if (BeVerbose) cout << endl
			                  << "Found " << mcFileName << " file, will read the signal tail parameters from it" << endl;
		Params.readFromFile(mcFileName);
	} else {
		cout << endl
		     << mcFileName << " file does not seem to exist, you need to extract the signal tail parameters from MC fit first!" << endl;
	}
	// fix the parameters
	sigma->setConstant();
	alphaInf->setConstant();
	orderInf->setConstant();
	alphaSup->setConstant();
	orderSup->setConstant();
	sigma_gauss->setConstant();
	normFraction->setConstant();
	// ratio_sigma -> setConstant();

	if (BeVerbose) {
		cout << endl
		     << "Parameters fixed to the following MC signal values:" << endl;
		Params.Print("v");
	}
	return Params;
}
*/
Double_t ComputeSignalSignificance(RooWorkspace& wspace, Int_t iState = 1) {
	RooRealVar mass = *wspace.var("mass");

	auto* signalPDF = wspace.pdf(Form("signalPDF_%dS", iState));

	RooRealVar mean1S = *(wspace.var("mean_1S"));
	RooFormulaVar mean2S = *(RooFormulaVar*)(wspace.function("mean_2S"));
	double mean = (iState == 1) ? mean1S.getVal() : mean2S.getVal();

	RooRealVar sigma1S = *(wspace.var("sigmaSymDSCB")); // to be fixed
	RooFormulaVar sigma2S = *(RooFormulaVar*)(wspace.function("sigma_2S"));
	double width = (iState == 1) ? sigma1S.getVal() : sigma2S.getVal();

	// mass window for integral
	mass.setRange("integral", mean - 3 * width, mean + 3 * width);

	Double_t signalYield = (signalPDF->createIntegral(mass, RooFit::NormSet(mass), RooFit::Range("integral")))->getVal() * wspace.var(Form("yield%dS", iState))->getVal();

	Double_t bkgYield = (wspace.pdf("bkgPDF")->createIntegral(mass, RooFit::NormSet(mass), RooFit::Range("integral")))->getVal() * wspace.var("yieldBkg")->getVal();

	return signalYield / sqrt(signalYield + bkgYield);
}

void SaveSignalYields(RooArgSet* signalYields, const char* bkgShapeName, const char* fitModelName, const char* extraString = "") {
	gSystem->mkdir("../SignalExtraction/SignalYields/", kTRUE);
	signalYields->writeToFile(Form("../SignalExtraction/SignalYields/%s_%s%s.txt", bkgShapeName, fitModelName, extraString));
}

void SaveRawSignalYields(RooArgSet* signalYields, const char* fitModelName, const char* extraString = gMuonAccName) {
	gSystem->mkdir("../SignalExtraction/RawYields/", kTRUE);
	signalYields->writeToFile(Form("../SignalExtraction/RawYields/%s%s.txt", fitModelName, extraString));
}

void SaveRawFitResults(RooFitResult* fitResult, const char* fitModelName, const char* extraString = gMuonAccName) {
	gSystem->mkdir("../SignalExtraction/RawFitResults", kTRUE);
	TFile fitResultsFile(Form("RawFitResults/%s%s.root", fitModelName, gMuonAccName), "RECREATE");
	
	fitResult->Write();
	fitResultsFile.Close();
}

void SavePolarizationFitParameters(RooArgSet* parameters, const char* methodName, const char* modelName, const char* extraString = "") {
	gSystem->mkdir("../Polarization/ParametersResults/", kTRUE);

	const char* fileName = Form("../Polarization/ParametersResults/%s_%s_%s.txt", methodName, modelName, extraString);
	parameters->writeToFile(fileName);

	auto list = parameters->contentsString();

	cout << "\n[Polarization] fit results for parameters (" << list << ") saved in " << fileName << endl;
}

RooArgSet GetSignalYields(RooRealVar* yield1S, RooRealVar* yield2S, RooRealVar* yield3S, const char* bkgShapeName, const char* fitModelName, const char* extraString = "") {
	RooArgSet signalYields(*yield1S, *yield2S, *yield3S);

	char yieldsFileName[512];
	snprintf(yieldsFileName, sizeof(yieldsFileName), "../SignalExtraction/RawYields/%s_%s%s.txt", bkgShapeName, fitModelName, extraString);

	cout << yieldsFileName << endl;
	if (fopen(yieldsFileName, "r")) {
		cout << endl
		     << "Found" << yieldsFileName << " file, will read signal yields from it" << endl;
		signalYields.readFromFile(yieldsFileName);
	} else {
		cout << endl
		     << yieldsFileName << " file does not seem to exist, you need to perform the signal extraction first!" << endl;
		exit(1);
	}

	cout << endl
	     << "Signal yields values:" << endl;
	signalYields.Print("v");

	return signalYields;
}

RooFitResult* GetFitResults(const char* totalFitModelName, const char* extraString = "") {
    TFile *fitResultsFile = TFile::Open(Form("../SignalExtraction/RawFitResults/%s%s.root", totalFitModelName, extraString), "READ");  

	char fitResultsFileName[512];
	snprintf(fitResultsFileName, sizeof(fitResultsFileName), "../SignalExtraction/RawFitResults/%s%s.root", totalFitModelName, extraString);

	cout << fitResultsFileName << endl;
    if (!fitResultsFile) {
		cout << endl
       		 << fitResultsFileName << " file does not seem to exist, you need to perform the signal extraction first!" << endl;
		exit(1);
	}

    else {
		cout << endl
		     << "Found" << fitResultsFileName << " file, will read fit results from it" << endl;	
	}

    RooFitResult* fitResults = (RooFitResult*)fitResultsFile->Get("fitresult_invMassModel_dataset");

	cout << endl
         << "Fit result values:" << endl;
	fitResults->Print("v");

	return fitResults;
}

RooArgSet GetPolarParams(RooRealVar* lambdaTheta, RooRealVar* lambdaPhi, RooRealVar* lambdaThetaPhi, RooRealVar* lambdaTilde, const char* methodName, const char* modelName, const char* extraString = "") {
	RooArgSet polarParams(*lambdaTheta, *lambdaPhi, *lambdaThetaPhi, *lambdaTilde);

	char paramsFileName[512];
	snprintf(paramsFileName, sizeof(paramsFileName), "../Polarization/ParametersResults/%s_%s_%s.txt", methodName, modelName, extraString);

	cout << paramsFileName << endl;
	if (fopen(paramsFileName, "r")) {
		cout << endl
		     << "Found" << paramsFileName << " file, will read polarization parameters from it" << endl;
		polarParams.readFromFile(paramsFileName);
	} else {
		cout << endl
		     << paramsFileName << " file does not seem to exist, you need to perform the signal extraction first!" << endl;
		exit(1);
	}

	cout << endl
	     << "Polarization parameter values:" << endl;
	polarParams.Print("v");

	return polarParams;
}

void SaveCanvas(TCanvas* canvasName, const char* fitModelName) {
	gSystem->mkdir("InvMassFits/CorrectedData/", kTRUE);
	canvasName->SaveAs(Form("InvMassFits/CorrectedData/%s.png", fitModelName), "RECREATE");
}

void SaveRawDataFitCanvas(TCanvas* canvasName, const char* totalFitModelName, const char* extraString = gMuonAccName) {
	gSystem->mkdir("InvMassFits/RawData/", kTRUE);
	canvasName->SaveAs(Form("InvMassFits/RawData/%s%s.png", totalFitModelName, extraString), "RECREATE");
}

void calculateChi2(TH1D* standardCorrectedHist, TF1* PolarFunc, Int_t nCosThetaBins = 10) {
	double chiSqr = 0;

	for (int i = 1; i <= nCosThetaBins; i++) {
		Double_t x = standardCorrectedHist->GetBinCenter(i);
		Double_t res = (standardCorrectedHist->GetBinContent(i) - PolarFunc->Eval(x)) / (standardCorrectedHist->GetBinError(i));
		if (res == 0) continue;

		chiSqr += TMath::Power(res, 2);

		// cout << "x: " << x << endl;
		// cout << "res: " << res << endl;
	}

	cout << "chi2: " << chiSqr << endl;
	cout << "nDOF: " << (nCosThetaBins - PolarFunc->GetNpar()) << endl;
	cout << "reduced chi2: " << chiSqr / (nCosThetaBins - PolarFunc->GetNpar()) << endl;
}

#endif
