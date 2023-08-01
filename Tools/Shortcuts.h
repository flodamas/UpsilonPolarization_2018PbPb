
// reduce the input dataset (N dimensions) to the mass dimension only dataset and apply desired kinematic cuts
RooDataSet* ReducedMassDataset(RooDataSet* allDataset, RooWorkspace* wspace, Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	// not cutting in CS and HF variables at the same time!!! Either you're analyzing CS or HX frame, but not both at once

	Float_t minCosThetaCS, maxCosThetaCS;
	Int_t minPhiCS, maxPhiCS;
	Float_t minCosThetaHX, maxCosThetaHX;
	Int_t minPhiHX, maxPhiHX;

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

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (pt > %d && pt < %d) && (cosThetaCS > %f && cosThetaCS < %f) && (phiCS > %d && phiCS < %d)&& (cosThetaHX > %f && cosThetaHX < %f) && (phiHX > %d && phiHX < %d)", 2 * centMin, 2 * centMax, ptMin, ptMax, minCosThetaCS, maxCosThetaCS, minPhiCS, maxPhiCS, minCosThetaHX, maxCosThetaHX, minPhiHX, maxPhiHX);

	RooDataSet* massDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace->var("mass"))), kinematicCut);
	massDataset->SetName("massDataset");

	wspace->import(*massDataset);

	return massDataset;
}


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


// return the signal shape tail parameters (alphaL, orderL, alphaR, orderR)
// based on the name of the fit model, if there is no .txt file corresponding, it will perform the MC fit in order to extract the parameters

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
