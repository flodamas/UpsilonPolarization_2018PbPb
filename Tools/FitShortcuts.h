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


void drawFitGraph(RooRealVar* massVar, RooWorkspace* wspace, RooDataSet* massDataset, TString pdfName, RooFitResult* fitResult, Float_t massMin, Float_t massMax, Int_t nBins){
	
	using namespace RooFit;
	
	/// draw the fit to see if the fit is reasonable
	string spdfName = pdfName.Data(); // convert TString to string
	auto* canvas = new TCanvas(Form("canvas%s",spdfName.c_str()), Form("canvas%s",spdfName.c_str()), 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	/// define frame for the invariant mass plot with a fit
	RooPlot* frame = massVar -> frame(Title(" "), Range(massMin, massMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	// frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset -> plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	RooAbsPdf* signal = nullptr;
	/// get the signal pdf depending on the pdf model
	if(pdfName=="SymDSCB") signal = wspace -> pdf("SymDSCB"); // symmetric core Double-sided Crystal Ball
	else if(pdfName=="SymDSCBGauss") signal = wspace -> pdf("DSCBGauss"); // Double-sided Crystal Ball + Gaussian
	else if(pdfName=="AsymDSCB") signal = wspace -> pdf("AsymDSCB"); // asymmetric core Double-sided Crystal Ball
	else( "Error! can't find the pdf model..");

	wspace -> Print();	

	/// put the signal pdf on the frame
	signal -> plotOn(frame, LineColor(kBlue));

	/// convert RooRealVals in wspace to Int_t
	RooRealVar* centMinVar = wspace->var("centMinVar");
	Int_t centMin = centMinVar -> getVal();
	RooRealVar* centMaxVar = wspace->var("centMaxVar");
	Int_t centMax = centMaxVar -> getVal();
	
	RooRealVar* ptMinVar = wspace->var("ptMinVar");
	Int_t ptMin = ptMinVar -> getVal();
	RooRealVar* ptMaxVar = wspace->var("ptMaxVar");
	Int_t ptMax = ptMaxVar -> getVal();
	
	RooRealVar* cosThetaMinVar = wspace->var("cosThetaMinVar");
	Int_t cosThetaMin = cosThetaMinVar -> getVal();
	RooRealVar* cosThetaMaxVar = wspace->var("cosThetaMaxVar");
	Int_t cosThetaMax = cosThetaMaxVar -> getVal();

	RooRealVar* isCSframeVar = wspace->var("isCSframeVar");
	Int_t isCSframe = isCSframeVar -> getVal();

	RooRealVar* phiMinVar = wspace->var("phiMinVar");
	Int_t phiMin = phiMinVar -> getVal();
	RooRealVar* phiMaxVar = wspace->var("phiMaxVar");
	Int_t phiMax = phiMaxVar -> getVal();

	/// Get variables in the signal
	RooArgSet* params = nullptr;
	params = signal -> getVariables() ;
	// cout << "params of " << pdfName << ": " << *params << endl;

	/// texts on the plot
	frame->addObject(KinematicsText(centMin, centMax, ptMin, ptMax));
	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	if(pdfName=="SymDSCB"){
		RooRealVar* meanVar = (RooRealVar*) params->find("meanSymDSCB");
		RooRealVar* sigmaVar = (RooRealVar*) params->find("sigmaSymDSCB");
		RooRealVar* alphaInfVar = (RooRealVar*) params->find("alphaInfSymDSCB");
		RooRealVar* orderInfVar = (RooRealVar*) params->find("orderInfSymDSCB");
		RooRealVar* alphaSupVar = (RooRealVar*) params->find("alphaSupSymDSCB");
		RooRealVar* orderSupVar = (RooRealVar*) params->find("orderSupSymDSCB");
		frame->addObject(SymCoreDoubleCBParamsText(*meanVar, *sigmaVar, *alphaInfVar, *orderInfVar, *alphaSupVar, *orderSupVar));
	}
	else if(pdfName=="SymDSCBGauss"){
		RooRealVar* meanVar = (RooRealVar*) params->find("meanDSCBGauss");
		RooRealVar* sigmaVar = (RooRealVar*) params->find("sigmaDSCBGauss");
		RooRealVar* alphaInfVar = (RooRealVar*) params->find("alphaInfDSCBGauss");
		RooRealVar* orderInfVar = (RooRealVar*) params->find("orderInfDSCBGauss");
		RooRealVar* alphaSupVar = (RooRealVar*) params->find("alphaSupDSCBGauss");
		RooRealVar* orderSupVar = (RooRealVar*) params->find("orderSupDSCBGauss");
		RooRealVar* sigma_gaussVar = (RooRealVar*) params->find("sigma_gauss");
		RooRealVar* normFractionVar = (RooRealVar*) params->find("normFraction");
		frame->addObject(SymCoreDoubleCBGaussParamsText(*meanVar, *sigmaVar, *alphaInfVar, *orderInfVar, *alphaSupVar, *orderSupVar, *sigma_gaussVar, *normFractionVar));
	}
	else if(pdfName=="AsymDSCB"){
		RooRealVar* meanVar = (RooRealVar*) params->find("meanAsymDSCB");
		RooRealVar* sigmaInfVar = (RooRealVar*) params->find("sigmaInfAsymDSCB");
		RooRealVar* alphaInfVar = (RooRealVar*) params->find("alphaInfAsymDSCB");
		RooRealVar* orderInfVar = (RooRealVar*) params->find("orderInfAsymDSCB");
		RooRealVar* sigmaSupVar = (RooRealVar*) params->find("sigmaSupAsymDSCB");
		RooRealVar* alphaSupVar = (RooRealVar*) params->find("alphaSupAsymDSCB");
		RooRealVar* orderSupVar = (RooRealVar*) params->find("orderSupAsymDSCB");
		frame->addObject(AsymDoubleCBParamsText(*meanVar, *sigmaInfVar, *alphaInfVar, *orderInfVar, *sigmaSupVar, *alphaSupVar, *orderSupVar));	
	}
	else("Error! can't find the pdf model..");

	frame->Draw();

	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->cd();
	pad1->Draw();
	pad2->Draw();
	canvas->SaveAs(Form("SignalParameters/MCfit_%s_%dto%d.png", spdfName.c_str(), ptMin, ptMax), "RECREATE");
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
