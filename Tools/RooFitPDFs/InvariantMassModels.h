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

// return the list of pdf names for internal usage
std::vector<std::string> BuildSignalPdfs(RooWorkspace& wspace, const char* signalShapeName = "SymDSCB") {
	std::vector<std::string> list;

	RooRealVar mass = *wspace.var("mass");
	mass.setRange("MassFitRange", MassBinMin, MassBinMax);

	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.35, 9.55);

	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);
	RooFormulaVar mean_2S("mean_2S", "@0*@1", RooArgSet(massScaling_2S, mean_1S));

	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);
	RooFormulaVar mean_3S("mean_3S", "@0*@1", RooArgSet(massScaling_3S, mean_1S));

	if (strcmp(signalShapeName, "SymDSCB") == 0) {
		// get the tail parameters, assuming that they have been imported to the workspace first!!
		// RooRealVar sigma_1S = *wspace.var("sigmaSymDSCB");
		RooRealVar sigma = *wspace.var("sigmaSymDSCB");

		RooRealVar alphaInf = *wspace.var("alphaInfSymDSCB");
		RooRealVar orderInf = *wspace.var("orderInfSymDSCB");
		RooRealVar alphaSup = *wspace.var("alphaSupSymDSCB");
		RooRealVar orderSup = *wspace.var("orderSupSymDSCB");

		// Y(1S) signal shape
		// RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.35, 9.55);
		RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.35, 9.60);
		// RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S);
		RooRealVar sigma_1S("sigma_1S", "", sigma.getVal(), .03, .15);

		RooCrystalBall signalPDF_1S("signalPDF_1S", "Symmetric DSCB pdf for Y(1S) mass peak", mass, mean_1S, sigma_1S, alphaInf, orderInf, alphaSup, orderSup);
		list.push_back(signalPDF_1S.GetName());

		// Y(2S) signal shape, mass scaling for mean and widths

		RooFormulaVar sigma_2S("sigma_2S", "@0*@1", RooArgSet(massScaling_2S, sigma_1S));

		RooCrystalBall signalPDF_2S("signalPDF_2S", "Symmetric DSCB pdf for Y(2S) mass peak", mass, mean_2S, sigma_2S, alphaInf, orderInf, alphaSup, orderSup);
		list.push_back(signalPDF_2S.GetName());

		// Y(3S) signal shape, mass scaling for mean and widths

		RooFormulaVar sigma_3S("sigma_3S", "@0*@1", RooArgSet(massScaling_3S, sigma_1S));

		RooCrystalBall signalPDF_3S("signalPDF_3S", "Symmetric DSCB pdf for Y(3S) mass peak", mass, mean_3S, sigma_3S, alphaInf, orderInf, alphaSup, orderSup);
		list.push_back(signalPDF_3S.GetName());

		wspace.import({signalPDF_1S, signalPDF_2S, signalPDF_3S}, RecycleConflictNodes());
	}

	else if (strcmp(signalShapeName, "Johnson") == 0) {
		// get the parameters, assuming that they have been imported to the workspace first!!
		RooRealVar lambda = *wspace.var("lambdaJohnson");
		RooRealVar gamma = *wspace.var("gammaJohnson");
		RooRealVar delta = *wspace.var("deltaJohnson");

		// Y(1S) signal shape

		RooJohnson signalPDF_1S("signalPDF_1S", "Johnson pdf for Y(1S) mass peak", mass, mean_1S, lambda, gamma, delta);
		list.push_back(signalPDF_1S.GetName());

		// Y(2S) signal shape, mass scaling for mean

		RooJohnson signalPDF_2S("signalPDF_2S", "Johnson pdf for Y(2S) mass peak", mass, mean_2S, lambda, gamma, delta);
		list.push_back(signalPDF_2S.GetName());

		// Y(3S) signal shape, mass scaling for mean

		RooJohnson signalPDF_3S("signalPDF_3S", "Johnson pdf for Y(3S) mass peak", mass, mean_3S, lambda, gamma, delta);
		list.push_back(signalPDF_3S.GetName());

		wspace.import({signalPDF_1S, signalPDF_2S, signalPDF_3S}, RecycleConflictNodes());
	}

	else if (strcmp(signalShapeName, "Voigtian") == 0) {
		// this function is not complete yet, need to add the mass scaling for the excited states
		// Voigtian distribution
		RooRealVar sigma("sigmaVoigtian", "", 0.08, .03, .15);
		RooRealVar width("widthVoigtian", "", 0.08, .03, .15);

		RooVoigtian signalPDF("signalPDF", "Voigtian distribution", mass, mean_1S, sigma, width);
		list.push_back(signalPDF.GetName());

		wspace.import(signalPDF, RecycleConflictNodes());
	}

	else {
		std::cerr << "No matching Signal shape name!" << std::endl;
		exit(1);
	}

	return list;
}

// Multiply the signal PDFs by constraints terms (if any)

void ApplyInternalConstraintsToSignal(RooWorkspace& wspace, std::vector<std::string>& signalPdfNameList, RooArgList* listOfConstraints) {
	// ideally we would like a simple function which can "guess" the constrained parameters (signal tail parameters, etc) and the corresponding terms (e.g. Gaussian terms which act as nuisance parameters)

	for (int i = 0; i < signalPdfNameList.size(); i++) {
		auto pdfName = signalPdfNameList.at(i);

		// get the original signal pdf

		if (wspace.pdf(pdfName) == nullptr) {
			std::cerr << "\n[InvariantMassModels] pdf " << pdfName << " not found in " << wspace.GetName() << ", check the calling order of the different internal functions.\nSkipping " << pdfName << " for internal constraints.\n";
			continue;
		}
		RooAbsPdf* originalPdf = wspace.pdf(pdfName.c_str());

		// list of pdfs to be multiplied
		RooArgList pdfList(*listOfConstraints);
		pdfList.add(*originalPdf);

		// cannot overwrite the signal pdf so use another name

		RooProdPdf* newPdf = new RooProdPdf(Form("constrained_%s", pdfName.c_str()), Form("(constrained) %s", originalPdf->GetTitle()), pdfList);
		newPdf->setNormRange("MassFitRange");

		wspace.import(*newPdf, RecycleConflictNodes());

		// replace the pdf name in the list
		signalPdfNameList[i] = newPdf->GetName();
	}
}

/// Background shape modeling

RooArgList* ChebychevCoefList(int order = 1) {
	RooArgList* coefList = new RooArgList();

	for (int i = 0; i <= order; i++) {
		RooRealVar* coef = new RooRealVar(Form("ChebychevCoef_%d", i), " ", -0.1, -10, 10);
		coefList->addOwned(*coef);
	}

	return coefList;
}

RooAbsPdf* BackgroundPDF(RooWorkspace& wspace, const char* bkgShapeName, const char* fitModelName) {
	RooRealVar* invMass = wspace.var("mass");
	invMass->setRange("MassFitRange", MassBinMin, MassBinMax);

	// Chebychev Nth order polynomial
	if (strncmp(bkgShapeName, "ChebychevOrder", 14) == 0) {
		int charSize = strlen(bkgShapeName);
		const char* lastDigit = bkgShapeName + charSize - 1;
		std::stringstream converter(lastDigit);

		int order = 0;
		converter >> order;

		RooArgList* coefList = ChebychevCoefList(order);
		RooChebychev* bkgPDF = new RooChebychev("bkgPDF", Form("Chebychev polynomial of order %d", order), *invMass, *coefList);

		bkgPDF->setNormRange("MassFitRange");

		wspace.import(*bkgPDF, RecycleConflictNodes());

		return bkgPDF;
	}

	// exponential x err function
	else if (strcmp(bkgShapeName, "ExpTimesErr") == 0) {
		RooRealVar* err_mu = new RooRealVar("err_mu", " ", 7, 4, 14);
		// RooRealVar* err_mu = new RooRealVar("err_mu", " ", 7);
		RooRealVar* err_sigma = new RooRealVar("err_sigma", " ", 1.5, 0.1, 5);
		// RooRealVar* err_sigma = new RooRealVar("err_sigma", " ", 0.9);
		RooRealVar* exp_lambda = new RooRealVar("exp_lambda", " ", 1.46, 0, 100);

		ErrorFuncTimesExp* bkgPDF = new ErrorFuncTimesExp("bkgPDF", "Product of an error function with an exponential", *invMass, *err_mu, *err_sigma, *exp_lambda);

		bkgPDF->setNormRange("MassFitRange");

		wspace.import(*bkgPDF, RecycleConflictNodes());

		return bkgPDF;
	}

	else if (strcmp(bkgShapeName, "Argus") == 0) {
		RooRealVar* m0 = new RooRealVar("m0_Argus", " ", 100, 0, 1000);
		RooRealVar* c = new RooRealVar("c_Argus", " ", 10, 0, 100);
		RooRealVar* p = new RooRealVar("p_Argus", " ", 1, 0, 20);

		RooArgusBG* bkgPDF = new RooArgusBG("bkgPDF", "ARGUS background shape", *invMass, *m0, *c);

		bkgPDF->setNormRange("MassFitRange");

		wspace.import(*bkgPDF, RecycleConflictNodes());

		return bkgPDF;
	}

	else if (strcmp(bkgShapeName, "Exponential") == 0) {
		RooRealVar* exp_lambda = new RooRealVar("exp_lambda", " ", -10, 10);

		RooExponential* bkgPDF = new RooExponential("bkgPDF", " ", *invMass, *exp_lambda);

		bkgPDF->setNormRange("MassFitRange");

		wspace.import(*bkgPDF, RecycleConflictNodes());

		return bkgPDF;
	}

	else {
		std::cout << "No matching background model name!" << std::endl;
		return nullptr;
	}
}

// build the main PDF based on all PDFs above

void BuildInvariantMassModel(RooWorkspace& wspace, const char* signalShapeName, const char* bkgShapeName, const char* fitModelName, Long64_t nEntries = 1e6, bool fixSigmaToMC = false, const char* muonAccName = "UpsilonTriggerThresholds") {
	RooRealVar* invMass = wspace.var("mass");
	invMass->setRange("MassFitRange", MassBinMin, MassBinMax);

	if (BeVerbose) std::cout << "\nBuilding invariant mass fit model with " << signalShapeName << " for the Y signal shapes and " << bkgShapeName << " for the background modeling\n";

	// 1. tail parameters fixed to MC extracted values, and identical for the three resonances

	RooArgList* constraintsList = ListOfSignalContraints(wspace, signalShapeName, fitModelName, fixSigmaToMC, muonAccName);

	// 2. signal model: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	auto signalPdfNameList = BuildSignalPdfs(wspace, signalShapeName);

	// 2.1 apply Gaussian penalty to internal constraint terms

	if (constraintsList->getSize()) ApplyInternalConstraintsToSignal(wspace, signalPdfNameList, constraintsList);

	auto signalPDF_1S = wspace.pdf(signalPdfNameList.at(0));
	auto signalPDF_2S = wspace.pdf(signalPdfNameList.at(1));
	auto signalPDF_3S = wspace.pdf(signalPdfNameList.at(2));
	signalPDF_1S->setNormRange("MassFitRange");
	signalPDF_2S->setNormRange("MassFitRange");
	signalPDF_3S->setNormRange("MassFitRange");

	// 3. background

	auto bkgPDF = BackgroundPDF(wspace, bkgShapeName, fitModelName);
	bkgPDF->setNormRange("MassFitRange");

	// complete invariant mass model

	Long64_t initYield = nEntries / 100;

	RooRealVar* yield1S = new RooRealVar("yield1S", "Number of Y(1S) signal candidates", initYield, 0, nEntries);
	RooRealVar* yield2S = new RooRealVar("yield2S", "Number of Y(2S) signal candidates", initYield / 4, 0, nEntries);
	RooRealVar* yield3S = new RooRealVar("yield3S", "Number of Y(3S) signal candidates", initYield / 12, 0, nEntries);

	RooRealVar* yieldBkg = new RooRealVar("yieldBkg", "Background yield", 0, nEntries);

	RooAddPdf model("invMassModel", "", {*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, *bkgPDF}, {*yield1S, *yield2S, *yield3S, *yieldBkg});
	model.setNormRange("MassFitRange");

	RooArgSet normSet(*invMass); // Define normalization set for the mass variable
	model.fixCoefNormalization(normSet);

	// wspace.import(model);
	wspace.import(model, RecycleConflictNodes());

	/// if the fit was already performed, get the seeds of the parameters from the fit results
	RooRealVar* err_mu = (RooRealVar*)wspace.var("err_mu");
	RooRealVar* err_sigma = (RooRealVar*)wspace.var("err_sigma");
	RooRealVar* exp_lambda = (RooRealVar*)wspace.var("exp_lambda");

	const char* totalFitModelName = Form("%s_%s", bkgShapeName, fitModelName);
	const char* savedFitResultFileName = Form("%s%s.root", totalFitModelName, gMuonAccName);

	// /// set the initial values of the parameters to the already fitted values
	// if (savedFitResultFileName) {
	// 	RooFitResult* fitResults = GetFitResults(totalFitModelName, gMuonAccName);

	// 	RooArgSet* fitParams = new RooArgSet(fitResults->floatParsFinal());

	// 	err_mu->setVal(((RooRealVar*)fitParams->find("err_mu"))->getVal());
	// 	err_sigma->setVal(((RooRealVar*)fitParams->find("err_sigma"))->getVal());
	// 	exp_lambda->setVal(((RooRealVar*)fitParams->find("exp_lambda"))->getVal());

	// 	cout << "err_mu: " << err_mu->getVal() << endl;
	// 	cout << "err_sigma: " << err_sigma->getVal() << endl;
	// 	cout << "exp_lambda: " << exp_lambda->getVal() << endl;
	// }

	// else {
	// 	std::cerr << savedFitResultFileName << "does not exist!!!" << std::endl;
	// }

	// /// set the initial values of the parameters to the fixed values
	// err_mu->setVal(6.7);
	// err_sigma->setVal(1.3);
	// exp_lambda->setVal(30);

	// cout << "err_mu: " << err_mu->getVal() << endl;
	// cout << "err_sigma: " << err_sigma->getVal() << endl;
	// cout << "exp_lambda: " << exp_lambda->getVal() << endl;

	if (BeVerbose) std::cout << "\nInvariant mass model exported to " << wspace.GetName() << std::endl;
	//wspace.Print("v");
}

#endif
