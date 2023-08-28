// reduce the input dataset (N dimensions) to the apply desired kinematic cuts
RooDataSet* ReducedDataset(RooDataSet* allDataset, RooWorkspace* wspace, Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	// not cutting in CS and HF variables at the same time!!! Either you're analyzing CS or HX frame, but not both at once

	Float_t minCosThetaCS, maxCosThetaCS;
	Int_t minPhiCS, maxPhiCS;
	Float_t minCosThetaHX, maxCosThetaHX;
	Int_t minPhiHX, maxPhiHX;
	Int_t minMass = 8, maxMass = 14;

	if (isCSframe) {
		minCosThetaCS = cosThetaMin;
		maxCosThetaCS = cosThetaMax;
		minPhiCS = phiMin;
		maxPhiCS = phiMax;

		minCosThetaHX = -1;
		maxCosThetaHX = 1;
		minPhiHX = -180;
		maxPhiHX = 180;
	} else {
		minCosThetaCS = -1;
		maxCosThetaCS = 1;
		minPhiCS = -180;
		maxPhiCS = 180;

		minCosThetaHX = cosThetaMin;
		maxCosThetaHX = cosThetaMax;
		minPhiHX = phiMin;
		maxPhiHX = phiMax;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (mass > %d && mass < %d) && (pt > %d && pt < %d) && (cosThetaCS > %f && cosThetaCS < %f) && (phiCS > %d && phiCS < %d)&& (cosThetaHX > %f && cosThetaHX < %f) && (phiHX > %d && phiHX < %d)", 2 * centMin, 2 * centMax, minMass, maxMass, ptMin, ptMax, minCosThetaCS, maxCosThetaCS, minPhiCS, maxPhiCS, minCosThetaHX, maxCosThetaHX, minPhiHX, maxPhiHX);
	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace->var("centrality")),*(wspace->var("mass")),*(wspace->var("rapidity")),*(wspace->var("pt")),*(wspace->var("cosThetaLab")),*(wspace->var("phiLab")),*(wspace->var("cosThetaCS")),*(wspace->var("phiCS")),*(wspace->var("cosThetaHX")),*(wspace->var("phiHX"))), kinematicCut);
	reducedDataset->SetName("reducedDataset");

	wspace->import(*reducedDataset);

	return reducedDataset;
}

// reduce the input dataset (N dimensions) to the mass dimension only dataset and apply desired kinematic cuts
RooDataSet* ReducedMassDataset(RooDataSet* allDataset, RooWorkspace* wspace, Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, double cosThetaMin = -1, double cosThetaMax = 1, double phiMin = -180, double phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	// not cutting in CS and HF variables at the same time!!! Either you're analyzing CS or HX frame, but not both at once

	double minCosThetaCS, maxCosThetaCS;
	double minPhiCS, maxPhiCS;
	double minCosThetaHX, maxCosThetaHX;
	double minPhiHX, maxPhiHX;

	if (isCSframe) {
		minCosThetaCS = cosThetaMin;
		maxCosThetaCS = cosThetaMax;
		minPhiCS = phiMin;
		maxPhiCS = phiMax;

		minCosThetaHX = -1;
		maxCosThetaHX = 1;
		minPhiHX = -180;
		maxPhiHX = 180;
	} else {
		minCosThetaCS = -1;
		maxCosThetaCS = 1;
		minPhiCS = -180;
		maxPhiCS = 180;

		minCosThetaHX = cosThetaMin;
		maxCosThetaHX = cosThetaMax;
		minPhiHX = phiMin;
		maxPhiHX = phiMax;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (pt > %d && pt < %d) && (cosThetaCS > %f && cosThetaCS < %f) && (phiCS > %f && phiCS < %f)&& (cosThetaHX > %f && cosThetaHX < %f) && (phiHX > %f && phiHX < %f)", 2 * centMin, 2 * centMax, ptMin, ptMax, minCosThetaCS, maxCosThetaCS, minPhiCS, maxPhiCS, minCosThetaHX, maxCosThetaHX, minPhiHX, maxPhiHX);

	RooDataSet* massDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace->var("mass"))), kinematicCut);
	massDataset->SetName(kinematicCut); // just to make it unique

	wspace->import(*massDataset);

	return massDataset;
}


RooFitResult* SymDSCBfit(RooRealVar* massVar, RooWorkspace* wspace, RooDataSet* massDataset, Float_t massMin, Float_t massMax){
	using namespace RooFit;

	// fit
	RooRealVar mean("meanSymDSCB", "", 9.457, 9., 10.);
	RooRealVar sigma("sigmaSymDSCB", "", 0.08, .05, .15);
	RooRealVar alphaInf("alphaInfSymDSCB", "", 1.5, 0.1, 10);
	RooRealVar orderInf("orderInfSymDSCB", "", 1.5, 0.1, 10);
	RooRealVar alphaSup("alphaSupSymDSCB", "", 1.5, 0.1, 10);
	RooRealVar orderSup("orderSupSymDSCB", "", 3, 0.1, 40);

	RooCrystalBall signal("SymDSCB", "SymDSCB", *massVar, mean, sigma, alphaInf, orderInf, alphaSup, orderSup);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with a double-sided Crystal Ball PDF made of a symmetric Gaussian core and asymmetric tail distributions..." << endl;

	bool doWeightedError = true;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true)/*, PrintLevel(-1)*/, Minos(!doWeightedError), NumCPU(3), Range(massMin, massMax), AsymptoticError(doWeightedError)); 
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit
	wspace -> import(signal);


	fitResult -> Print("v");

 	return fitResult;
}


RooFitResult* AsymDSCBfit(RooRealVar* massVar, RooWorkspace* wspace, RooDataSet* massDataset, Float_t massMin, Float_t massMax){
	using namespace RooFit;
	// fit
	RooRealVar mean("meanAsymDSCB", "", 9.457, 9., 10.);
	RooRealVar sigmaInf("sigmaInfAsymDSCB", "", 0.08, .05, .15);
	RooRealVar alphaInf("alphaInfAsymDSCB", "", 1.5, 0.1, 10);
	RooRealVar orderInf("orderInfAsymDSCB", "", 1.5, 0.1, 10);
	RooRealVar sigmaSup("sigmaSupAsymDSCB", "", 0.08, .05, .15);
	RooRealVar alphaSup("alphaSupAsymDSCB", "", 1.5, 0.1, 10);
	RooRealVar orderSup("orderSupAsymDSCB", "", 3, 0.1, 40);

	RooCrystalBall signal("AsymDSCB", "AsymDSCB", *massVar, mean, sigmaInf, sigmaSup, alphaInf, orderInf, alphaSup, orderSup);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with a double-sided Crystal Ball PDF made of an asymmetric Gaussian core and asymmetric tail distributions..." << endl;

	bool doWeightedError = true;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true)/*, PrintLevel(-1)*/, Minos(!doWeightedError), NumCPU(3), Range(massMin, massMax), AsymptoticError(doWeightedError)); 
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit
	wspace -> import(signal);

	fitResult -> Print("v");

 	return fitResult;
}


RooFitResult* SymDSCBGaussfit(RooRealVar* massVar, RooWorkspace* wspace, RooDataSet* massDataset, Float_t massMin, Float_t massMax){
	using namespace RooFit;

	/// fit
	/// (DSCB variables)
	RooRealVar mean("meanDSCBGauss", "", 9.457, 9., 10.);
	RooRealVar sigma("sigmaDSCBGauss", "", 0.08, .05, .15);
	RooRealVar alphaInf("alphaInfDSCBGauss", "", 1.5, 0.1, 10);
	RooRealVar orderInf("orderInfDSCBGauss", "", 1.5, 0.1, 10);
	RooRealVar alphaSup("alphaSupDSCBGauss", "", 1.5, 0.1, 10);
	RooRealVar orderSup("orderSupDSCBGauss", "", 3, 0.1, 40);
	
	/// (Gaussian variable (used the same mean as DSCB))
	RooRealVar sigma_gauss("sigma_gauss", "", 0.11, .05, .3);

	/// (fraction between two PDFs)
	RooRealVar normFraction("normFraction", "", 0.6, 0.01, 1);


	RooCrystalBall DSCB("DSCB", "DSCB", *massVar, mean, sigma, alphaInf, orderInf, alphaSup, orderSup);
	RooGaussian gauss("gauss","gaussian", *massVar, mean, sigma_gauss);

	RooAddPdf signal("DSCBGauss", "sum of DSCB and CB PDF", RooArgList(DSCB, gauss), RooArgList(normFraction), kTRUE);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with the sum of a double-sided Crystal Ball PDF made of a symmetric Gaussian core and asymmetric tail distributions, and a Gaussian PDF..." << endl;

	bool doWeightedError = true;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true)/*, PrintLevel(-1)*/, Minos(!doWeightedError), NumCPU(3), Range(massMin, massMax), AsymptoticError(doWeightedError)); 
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit
	wspace -> import(signal);
	
	fitResult->Print("v");

	return fitResult;
}



RooFitResult* Hypatiafit(RooRealVar* massVar, RooWorkspace* wspace, RooDataSet* massDataset, Float_t massMin, Float_t massMax){
	using namespace RooFit;

	// fit
	RooRealVar mean("meanHypatia", "", 9.451, 9., 10.);
	RooRealVar lambda("lambdaHypatia","lambda of hypatia PDF",-1.0, -20.0, -0.1);
	RooRealVar zeta("zetaHypatia","zeta of hypatia PDF",0.01, 0.0, 1.0);
	RooRealVar beta("betaHypatia","beta of hypatia PDF",-0.01, -20.0, 0.0);
	RooRealVar sigma("sigmaHypatia","sigma of hypatia PDF",0.15, 0.1, 0.3);
	RooRealVar alphaInf("alphaInfHypatia","al1s of hypatia PDF", 3.0 , 0.1, 5.0);
	RooRealVar alphaSup("alphaSupHypatia","ar1s of hypatia PDF", 3.0 , 0.1, 7.0);
	RooRealVar orderInf("orderInfHypatia","nl1s of hypatia PDF", 1.0 , 0.2, 14.718);
	RooRealVar orderSup("orderSupHypatia","nr1s of hypatia PDF", 1.0 , 0.0, 14.718);

	lambda.setVal(-1.729); 
	zeta.setVal(0); 
	beta.setVal(-1.59); 
	sigma.setVal(0.158); 
	alphaInf.setVal(1.769); 
	alphaSup.setVal(5.398); 
	orderInf.setVal(1.422); 
	orderSup.setVal(0.012); 

	zeta.setConstant(1);

	RooHypatia2 signal("Hypatia","Hypatia", *massVar, lambda, zeta, beta, sigma, mean, alphaInf, orderInf, alphaSup, orderSup);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with a Hypatia PDF" << endl;

	bool doWeightedError = true;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true)/*, PrintLevel(-1)*/, Minos(!doWeightedError), NumCPU(8), Range(massMin, massMax), AsymptoticError(doWeightedError)); 
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit
	wspace -> import(signal);

	fitResult -> Print("v");

 	return fitResult;
}

RooArgSet GetMCSignalTailParameters(RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, TString signalShapeName = "symCoreDSCB", Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
	RooArgSet tailParams(*alphaInf, *orderInf, *alphaSup, *orderSup);

	// if the .txt file for this specific fit model exists, just read the tail parameters from it
	const char* mcFileName = Form("../MonteCarlo/SignalParameters/%s_cent%dto%d_pt%dto%d.txt", signalShapeName.Data(), centMin, centMax, ptMin, ptMax);

	if (fopen(mcFileName, "r")) {
		cout << endl
		     << "Found " << mcFileName << " file, will read the signal tail parameters from it" << endl;
		tailParams.readFromFile(mcFileName);
	} else {
		cout << endl
		     << mcFileName << " file does not seem to exist, you need to extract the signal tail paramaters from MC fit first!" << endl;
	}

	// fix the tail parameters
	alphaInf->setConstant();
	orderInf->setConstant();
	alphaSup->setConstant();
	orderSup->setConstant();

	cout << endl
	     << "Tail parameters fixed to the following MC signal values:" << endl;
	tailParams.Print("v");

	return tailParams;
}

RooArgSet GetMCSignalParameters(RooRealVar* sigma, RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, RooRealVar* sigma_gauss, RooRealVar* normFraction, TString signalShapeName = "symCoreDSCB", Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
// RooArgSet GetMCSignalParameters(RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, RooRealVar* sigma_gauss, TString signalShapeName = "symCoreDSCB", Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
// RooArgSet GetMCSignalParameters(RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, RooFormulaVar* ratio_sigma, TString signalShapeName = "symCoreDSCB", Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
// RooArgSet GetMCSignalParameters(RooRealVar* sigma, RooRealVar* alphaInf, RooRealVar* orderInf, RooRealVar* alphaSup, RooRealVar* orderSup, RooRealVar* sigma_gauss, TString signalShapeName = "symCoreDSCB", Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {

	RooArgSet Params(*sigma, *alphaInf, *orderInf, *alphaSup, *orderSup, *sigma_gauss, *normFraction);
	// RooArgSet Params(*alphaInf, *orderInf, *alphaSup, *orderSup, *sigma_gauss);
	// RooArgSet Params(*alphaInf, *orderInf, *alphaSup, *orderSup, *ratio_sigma);
	// RooArgSet Params(*sigma, *alphaInf, *orderInf, *alphaSup, *orderSup, *sigma_gauss);

	// if the .txt file for this specific fit model exists, just read the tail parameters from it
	const char* mcFileName = Form("../MonteCarlo/SignalParameters/%s_cent%dto%d_pt%dto%d.txt", signalShapeName.Data(), centMin, centMax, ptMin, ptMax);

	if (fopen(mcFileName, "r")) {
		cout << endl
		     << "Found " << mcFileName << " file, will read the signal tail parameters from it" << endl;
		Params.readFromFile(mcFileName);
	} else {
		cout << endl
		     << mcFileName << " file does not seem to exist, you need to extract the signal tail paramaters from MC fit first!" << endl;
	}

	// fix the parameters
	sigma -> setConstant();
	alphaInf->setConstant();
	orderInf->setConstant();
	alphaSup->setConstant();
	orderSup->setConstant();
	sigma_gauss -> setConstant();
	normFraction -> setConstant();
	// ratio_sigma -> setConstant();

	cout << endl
	     << "Parameters fixed to the following MC signal values:" << endl;
	Params.Print("v");

	return Params;
}
