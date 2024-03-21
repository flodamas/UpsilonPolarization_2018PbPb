#include "../../Tools/BasicHeaders.h"

#include "../../AnalysisParameters.h"

using namespace RooFit;

// open the data file, create a workspace, and import the dataset
RooWorkspace SetUpWorkspace(const char* filename, const char* refFrameName = "CS") {
	cout << "\nSetting up the workspace...\n\n";

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return RooWorkspace("emptyWorkspace");
	}

	cout << "File " << filename << " opened" << endl;

	const char* datasetName = RawDatasetName(refFrameName);
	RooDataSet* data = (RooDataSet*)f->Get(datasetName);

	RooWorkspace wspace(Form("workspace%s", refFrameName));

	wspace.import(*data);

	cout << endl
	     << datasetName << " imported into " << wspace.GetName() << endl;

	f->Close();

	return wspace;
}

// reduce the whole dataset (N dimensions) to (invariant mass, cos theta, phi), allowing cutting off along phi
RooDataSet* InvMassCosThetaPhiDataset(RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180) {
	RooDataSet* allDataset = (RooDataSet*)wspace.data(RawDatasetName(refFrameName));

	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (%s > %d && %s < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass")), *(wspace.var(CosThetaVarName(refFrameName))), *(wspace.var(PhiVarName(refFrameName)))), kinematicCut);

	wspace.import(*reducedDataset, RooFit::Rename(Form("(inv mass, cos theta, phi) %s reduced dataset", refFrameName)));

	return reducedDataset;
}

// reduce the input dataset (N dimensions) to the apply desired kinematic cuts
RooDataSet* ReducedDataset(RooDataSet* allDataset, RooWorkspace* wspace, Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30, Double_t massMin = 8, Double_t massMax = 14, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
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

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (mass > %f && mass < %f) && (pt > %d && pt < %d) && (cosThetaCS > %f && cosThetaCS < %f) && (phiCS > %d && phiCS < %d)&& (cosThetaHX > %f && cosThetaHX < %f) && (phiHX > %d && phiHX < %d)", 2 * centMin, 2 * centMax, massMin, massMax, ptMin, ptMax, minCosThetaCS, maxCosThetaCS, minPhiCS, maxPhiCS, minCosThetaHX, maxCosThetaHX, minPhiHX, maxPhiHX);
	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace->var("centrality")), *(wspace->var("mass")), *(wspace->var("rapidity")), *(wspace->var("pt")), *(wspace->var("cosThetaLab")), *(wspace->var("phiLab")), *(wspace->var("etaLabMupl")), *(wspace->var("etaLabMumi")), *(wspace->var("cosThetaCS")), *(wspace->var("phiCS")), *(wspace->var("cosThetaHX")), *(wspace->var("phiHX"))), kinematicCut);

	reducedDataset->SetName("reducedDataset");

	wspace->import(*reducedDataset);

	return reducedDataset;
}

// reduce the input dataset (N dimensions) to the mass dimension only dataset and apply desired kinematic cuts
RooDataSet* ReducedMassDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, double cosThetaMin = -1, double cosThetaMax = 1, double phiMin = -180, double phiMax = 180) {
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

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (cosThetaCS > %f && cosThetaCS < %f) && (phiCS > %f && phiCS < %f)&& (cosThetaHX > %f && cosThetaHX < %f) && (phiHX > %f && phiHX < %f)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, minCosThetaCS, maxCosThetaCS, minPhiCS, maxPhiCS, minCosThetaHX, maxCosThetaHX, minPhiHX, maxPhiHX);

	RooDataSet* massDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass"))), kinematicCut);
	massDataset->SetName(kinematicCut); // just to make it unique

	wspace.import(*massDataset);

	return massDataset;
}
