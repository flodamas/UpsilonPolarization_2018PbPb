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
std::vector<std::string> BuildSignalPdfs(RooWorkspace& wspace) {
	std::vector<std::string> list;

	RooRealVar mass = *wspace.var("mass");
	mass.setRange("fitRange", 7, 13);

	// get the tail parameters, assuming that they have been imported to the workspace first!!
	RooRealVar sigma_1S = *wspace.var("sigmaSymDSCB");

	RooRealVar alphaInf = *wspace.var("alphaInfSymDSCB");
	RooRealVar orderInf = *wspace.var("orderInfSymDSCB");
	RooRealVar alphaSup = *wspace.var("alphaSupSymDSCB");
	RooRealVar orderSup = *wspace.var("orderSupSymDSCB");

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.35, 9.55);
	//RooRealVar sigma_1S("sigma_1S", "", sigma.getVal(), .03, .15);

	RooCrystalBall signalPDF_1S("signalPDF_1S", "Symmetric DSCB pdf for Y(1S) mass peak", mass, mean_1S, sigma_1S, alphaInf, orderInf, alphaSup, orderSup);
	list.push_back(signalPDF_1S.GetName());

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "@0*@1", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigma_2S("sigma_2S", "@0*@1", RooArgSet(massScaling_2S, sigma_1S));

	RooCrystalBall signalPDF_2S("signalPDF_2S", "Symmetric DSCB pdf for Y(2S) mass peak", mass, mean_2S, sigma_2S, alphaInf, orderInf, alphaSup, orderSup);
	list.push_back(signalPDF_2S.GetName());

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "@0*@1", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigma_3S("sigma_3S", "@0*@1", RooArgSet(massScaling_3S, sigma_1S));

	RooCrystalBall signalPDF_3S("signalPDF_3S", "Symmetric DSCB pdf for Y(3S) mass peak", mass, mean_3S, sigma_3S, alphaInf, orderInf, alphaSup, orderSup);
	list.push_back(signalPDF_3S.GetName());

	wspace.import({signalPDF_1S, signalPDF_2S, signalPDF_3S}, RecycleConflictNodes());

	return list;
}

RooAddPdf* NominalSignalModel(RooWorkspace& wspace, Long64_t nEntries = 1e6) {
	RooRealVar mass = *wspace.var("mass");
	mass.setRange("fitRange", 7, 13);

	Long64_t initYield = nEntries / 100;

	// get the tail parameters, assuming that they have been imported to the workspace first!!
	RooRealVar sigma = *wspace.var("sigmaSymDSCB");
	//	sigma.Print();
	//	Double_t sigmaSeed = sigma.getVal();

	//RooRealVar sigma_1S = *wspace.var("sigma_1S");
	RooRealVar alphaInf = *wspace.var("alphaInfSymDSCB");
	RooRealVar orderInf = *wspace.var("orderInfSymDSCB");
	RooRealVar alphaSup = *wspace.var("alphaSupSymDSCB");
	RooRealVar orderSup = *wspace.var("orderSupSymDSCB");

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.35, 9.55);
	RooRealVar sigma_1S("sigma_1S", "", sigma.getVal(), .02, .15);

	RooCrystalBall* signalPDF_1S = new RooCrystalBall("signalPDF_1S", "", mass, mean_1S, sigma_1S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar* yield1S = new RooRealVar("yield1S", "N 1S", initYield, 0, nEntries);

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "@0*@1", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigma_2S("sigma_2S", "@0*@1", RooArgSet(massScaling_2S, sigma_1S));

	RooCrystalBall* signalPDF_2S = new RooCrystalBall("signalPDF_2S", "", mass, mean_2S, sigma_2S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar* yield2S = new RooRealVar("yield2S", "N 2S", initYield / 4, 0, nEntries);

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "@0*@1", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigma_3S("sigma_3S", "@0*@1", RooArgSet(massScaling_3S, sigma_1S));

	RooCrystalBall* signalPDF_3S = new RooCrystalBall("signalPDF_3S", "", mass, mean_3S, sigma_3S, alphaInf, orderInf, alphaSup, orderSup);
	RooRealVar* yield3S = new RooRealVar("yield3S", "N 3S", initYield / 12, 0, nEntries);

	RooAddPdf* SymDSCBModel = new RooAddPdf("signalPDF", "PDF of the sum of the three Y signal PDFs", {*signalPDF_1S, *signalPDF_2S, *signalPDF_3S}, {*yield1S, *yield2S, *yield3S});
	SymDSCBModel->setNormRange("fitRange");

	wspace.import(*SymDSCBModel, RecycleConflictNodes());

	return SymDSCBModel;
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

		wspace.import(*newPdf, RecycleConflictNodes());

		// replace the pdf name in the list
		signalPdfNameList[i] = newPdf->GetName();
	}
}

/// Background shape modeling

RooArgList* ChebychevCoefList(int order = 1) {
	RooArgList* coefList = new RooArgList();

	for (int i = 0; i <= order; i++) {
		RooRealVar* coef = new RooRealVar(Form("ChebychevCoef_%d", i), " ", 0.1, -3, 3);
		coefList->addOwned(*coef);
	}

	return coefList;
}

RooAbsPdf* BackgroundPDF(RooWorkspace& wspace, const char* bkgShapeName) {
	RooRealVar* invMass = wspace.var("mass");
	invMass->setRange("fitRange", 7, 13);

	// Chebychev Nth order polynomial
	if (strncmp(bkgShapeName, "ChebychevOrder", 14) == 0) {
		int charSize = strlen(bkgShapeName);
		const char* lastDigit = bkgShapeName + charSize - 1;
		std::stringstream converter(lastDigit);

		int order = 0;
		converter >> order;

		RooArgList* coefList = ChebychevCoefList(order);
		RooChebychev* bkgPDF = new RooChebychev("bkgPDF", Form("Chebychev polynomial of order %d", order), *invMass, *coefList);
		bkgPDF->setNormRange("fitRange");

		wspace.import(*bkgPDF, RecycleConflictNodes());

		return bkgPDF;
	}

	// exponential x err function
	else if (strcmp(bkgShapeName, "ExpTimesErr") == 0) {
		RooRealVar* err_mu = new RooRealVar("err_mu", " ", 6.9, 2, 13);
		// RooRealVar* err_mu = new RooRealVar("err_mu", " ", 7.2);
		RooRealVar* err_sigma = new RooRealVar("err_sigma", " ", 0.8, 0.0001, 10);
		RooRealVar* exp_lambda = new RooRealVar("exp_lambda", " ", 4.0, 0, 15);

		ErrorFuncTimesExp* bkgPDF = new ErrorFuncTimesExp("bkgPDF", "Product of an error function with an exponential", *invMass, *err_mu, *err_sigma, *exp_lambda);
		bkgPDF->setNormRange("fitRange");

		wspace.import(*bkgPDF, RecycleConflictNodes());

		return bkgPDF;
	}

	else if (strcmp(bkgShapeName, "Exponential") == 0) {
		RooRealVar* exp_lambda = new RooRealVar("exp_lambda", " ", -10, 10);

		RooExponential* bkgPDF = new RooExponential("bkgPDF", " ", *invMass, *exp_lambda);
		bkgPDF->setNormRange("fitRange");

		wspace.import(*bkgPDF, RecycleConflictNodes());

		return bkgPDF;
	}

	else {
		std::cout << "No matching background model name!" << std::endl;
		return nullptr;
	}
}

// build the main PDF based on all PDFs above

void BuildInvariantMassModel(RooWorkspace& wspace, const char* signalShapeName, const char* bkgShapeName, const char* fitModelName, Long64_t nEntries = 1e6, bool fixSigmaToMC = false) {
	if (BeVerbose) std::cout << "\nBuilding invariant mass fit model with " << signalShapeName << " for the Y signal shapes and " << bkgShapeName << " for the background modeling\n";

	// 1. tail parameters fixed to MC extracted values, and identical for the three resonances

	RooArgList* constraintsList = ListOfSignalContraints(wspace, signalShapeName, fitModelName, fixSigmaToMC);

	// 2. signal model: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	auto signalPdfNameList = BuildSignalPdfs(wspace);

	// 2.1 apply Gaussian penalty to internal constraint terms

	if (constraintsList->getSize()) ApplyInternalConstraintsToSignal(wspace, signalPdfNameList, constraintsList);

	auto signalPDF_1S = wspace.pdf(signalPdfNameList.at(0));
	auto signalPDF_2S = wspace.pdf(signalPdfNameList.at(1));
	auto signalPDF_3S = wspace.pdf(signalPdfNameList.at(2));

	// 3. background
	RooRealVar* invMass = wspace.var("mass");
	invMass->setRange("fitRange", 7, 13);

	auto bkgPDF = BackgroundPDF(wspace, bkgShapeName);
	bkgPDF->setNormRange("fitRange");

	// complete invariant mass model

	Long64_t initYield = nEntries / 100;

	RooRealVar* yield1S = new RooRealVar("yield1S", "Number of Y(1S) signal candidates", initYield, 0, nEntries);
	RooRealVar* yield2S = new RooRealVar("yield2S", "Number of Y(2S) signal candidates", initYield / 4, 0, nEntries);
	RooRealVar* yield3S = new RooRealVar("yield3S", "Number of Y(3S) signal candidates", initYield / 12, 0, nEntries);

	RooRealVar* yieldBkg = new RooRealVar("yieldBkg", "Background yield", 0, nEntries);

	RooAddPdf model("invMassModel", "", {*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, *bkgPDF}, {*yield1S, *yield2S, *yield3S, *yieldBkg});

	wspace.import(model, RecycleConflictNodes());

	if (BeVerbose) std::cout << "\nInvariant mass model exported to " << wspace.GetName() << std::endl;
	//wspace.Print("v");
}

RooAddPdf* MassFitModel(RooWorkspace& wspace, const char* signalShapeName, const char* bkgShapeName, const char* fitModelName, Long64_t nEntries = 1e6) {
	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	if (BeVerbose) std::cout << "\nBuilding invariant mass fit model with " << signalShapeName << " for the signal and " << bkgShapeName << " for the background\n";

	// tail parameters fixed to MC extracted values, and identical for the three resonances

	//ImportAndFixMCSignalParameters(wspace, signalShapeName, fitModelName);

	RooArgList* constraintsList = ListOfSignalContraints(wspace, signalShapeName, fitModelName, true);

	auto signalModel = NominalSignalModel(wspace, nEntries);

	// apply Gaussian penalty to  internal constraint terms
	constraintsList->add(*wspace.pdf("SymDSCBModel"));

	cout << "\nContent of constraints list:\n";

	constraintsList->Print("v");

	RooProdPdf constrainedSignalModel("constrainedSignalModel", " ", *constraintsList);

	wspace.import(constrainedSignalModel, RecycleConflictNodes()); // including all signal yield variables!!

	cout << "\nConstraints applied to the signal model\n";

	RooAbsPdf* signalPDF_1S = wspace.pdf("signalPDF_1S");
	RooAbsPdf* signalPDF_2S = wspace.pdf("signalPDF_2S");
	RooAbsPdf* signalPDF_3S = wspace.pdf("signalPDF_3S");

	RooRealVar* yield1S = wspace.var("yield1S");
	RooRealVar* yield2S = wspace.var("yield2S");
	RooRealVar* yield3S = wspace.var("yield3S");

	cout << "\nList of PDFs in constrainedSignalModel:\n";

	constrainedSignalModel.pdfList().Print("v");

	// background
	RooRealVar* invMass = wspace.var("mass");
	invMass->setRange("fitRange", 7, 13);

	auto bkgPDF = BackgroundPDF(wspace, bkgShapeName);
	bkgPDF->setNormRange("fitRange");

	RooRealVar* yieldBkg = new RooRealVar("yieldBkg", "Background yield", 0, nEntries);

	// complete invariant mass model

	cout << "\n[InvariantMassModels] Adding signal and background PDFs... ";

	//RooAddPdf* model = new RooAddPdf("invMassModel", "(constrained) signal + background PDF", RooArgList(constrainedSignalModel, bkgModel));

	cout << "Done!!\n";

	//model->Print("v");

	RooAddPdf* model = new RooAddPdf("invMassModel", "", {*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, *bkgPDF}, {*yield1S, *yield2S, *yield3S, *yieldBkg});
	
	wspace.import(*model, RecycleConflictNodes());

	if (BeVerbose) std::cout << "\nModel exported to " << wspace.GetName() << std::endl;
	//wspace.Print("v");

	return model;
}

#endif
