#ifndef ROODATASETHELPERS_H
#define ROODATASETHELPERS_H

#include "../../Tools/BasicHeaders.h"

#include "../../AnalysisParameters.h"

using namespace RooFit;

// open the data file, create a workspace, and import the dataset
RooWorkspace SetUpWorkspace(const char* filename, const char* refFrameName = "") {
	if (BeVerbose) std::cout << "\nSetting up the workspace...\n\n";

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		std::cout << "File " << filename << " not found. Check the directory of the file.\n";
		return RooWorkspace("emptyWorkspace");
	}

	if (BeVerbose) std::cout << "File " << filename << " opened\n";

	const char* datasetName = RawDatasetName(refFrameName);
	RooDataSet* data = (RooDataSet*)f->Get(datasetName);

	RooWorkspace wspace(Form("workspace%s", refFrameName));

	wspace.import(*data);

	if (BeVerbose) std::cout << datasetName << " imported into " << wspace.GetName() << std::endl;

	f->Close();

	return wspace;
}

// reduce the whole raw dataset (N dimensions) to (invariant mass, cos theta, phi)
RooDataSet InvMassCosThetaPhiDataset(RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "") {
	RooDataSet* allDataset = (RooDataSet*)wspace.data(RawDatasetName(refFrameName));

	if (allDataset == nullptr) {
		std::cerr << "Null RooDataSet provided to the reducer method!!" << std::endl;
		return RooDataSet();
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass")), *(wspace.var(CosThetaVarName("CS"))), *(wspace.var(PhiVarName("CS"))), *(wspace.var(CosThetaVarName("HX"))), *(wspace.var(PhiVarName("HX")))), kinematicCut);

	wspace.import(reducedDataset, RooFit::Rename(Form("reducedDataset%s", refFrameName)));

	return reducedDataset;
}

// reduce the whole dataset (N dimensions) to (cos theta, phi) for a given reference frame
RooDataSet ReducedRecoMCDataset(RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS") {
	RooDataSet* allDataset = (RooDataSet*)wspace.data("MCdataset");

	if (allDataset == nullptr) {
		std::cerr << "Null RooDataSet provided to the reducer method!!" << std::endl;
		return RooDataSet();
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) ", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var(CosThetaVarName(refFrameName))), *(wspace.var(PhiVarName(refFrameName)))), kinematicCut);

	wspace.import(reducedDataset, RooFit::Rename(Form("RecoMCDataset%s", refFrameName)));

	return reducedDataset;
}

// reduce the input dataset (N dimensions) to the mass dimension only dataset and apply desired kinematic cuts
RooDataSet ReducedMassDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, double cosThetaMin = -1, double cosThetaMax = 1, double phiMin = -180, double phiMax = 180) {
	if (allDataset == nullptr) {
		std::cerr << "Null RooDataSet provided to the reducer method!!" << std::endl;
		return RooDataSet();
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

	RooDataSet massDataset = *(RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass"))), kinematicCut);
	massDataset.SetName(kinematicCut); // just to make it unique

	wspace.import(massDataset);

	return massDataset;
}

#endif
