#ifndef SPlotHelpers_h
#define SPlotHelpers_h

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/FitDistributions.h"
#include "../Tools/RooFitPDFs/InvariantMassModels.h"

#include "RooStats/SPlot.h"

using namespace RooStats;

/// sPlot class documentation https://root.cern/doc/master/classRooStats_1_1SPlot.html

const char* sWeightedDatasetsFileName = "../sPlot/Datasets/sWeightedDatasets.root";

const char* sWeightedDatasetName(const char* fitModelName) {
	const char* name = Form("sWeightedDataset_%s", fitModelName);

	return name;
}

RooDataSet* CreateSWeights(RooWorkspace& wspace, Int_t ptMin, Int_t ptMax, const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ExpTimesErr") {
	// first, check that the dataset is present in the workspace

	if (!wspace.data("reducedDataset")) {
		std::cout << "\n[sPlot] ERROR: cannot find the reduced dataset in the workspace! Returning nothing now!!!\n\n";
		return nullptr;
	}

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax);

	auto data = *(RooDataSet*)wspace.data("reducedDataset");
	Long64_t nEntries = data.sumEntries();

	// then, check if the invariant mass model is present in the workspace. If not, build it!

	BuildInvariantMassModel(wspace, signalShapeName, bkgShapeName, fitModelName, nEntries);

	// need a "standard fit" first

	auto* fitResult = RawInvariantMassFit(wspace, data, RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S")));

	// draw and save the invariant mass distribution, to check the fit
	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	const char* totalFitModelName = GetTotalFitModelName(bkgShapeName, signalShapeName, ptMin, ptMax);

	gSystem->mkdir("../sPlot/InvMassFits", kTRUE);
	massCanvas->SaveAs(Form("../sPlot/InvMassFits/%s.png", totalFitModelName), "RECREATE");

	delete massCanvas;

	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments

	if (BeVerbose) {
		std::cout << "\n---------------------sWeighted data based on yields-------------------\n\nDataset content before creating sWeights:\n";

		data.Print();
	}

	RooRealVar yield1S = *wspace.var("yield1S");
	RooRealVar yield2S = *wspace.var("yield2S");
	RooRealVar yield3S = *wspace.var("yield3S");
	RooRealVar yieldBkg = *wspace.var("yieldBkg");

	RooStats::SPlot mySPlot{"mySPlot", "", data, wspace.pdf("invMassModel"), RooArgList(yield1S, yield2S, yield3S, yieldBkg), RooArgSet(), true, false, "", Range("MassFitRange"), NumCPU(NCPUs)};

	// copy the newly created dataset

	const char* fitName = GetTotalFitModelName(bkgShapeName, signalShapeName, ptMin, ptMax);

	auto* sData = new RooDataSet(data, sWeightedDatasetName(fitName));

	if (BeVerbose) {
		std::cout << "\nAfter creating sWeights:\n";
		sData->Print();

		// check that our weights have the desired properties

		std::cout << "\nChecking sWeights:" << std::endl;

		std::cout << "\nYield of Y(1S) is\t" << yield1S.getVal() << ".  From sWeights it is "
		          << mySPlot.GetYieldFromSWeight("yield1S") << std::endl;

		std::cout << "Yield of background is\t" << yieldBkg.getVal() << ".  From sWeights it is "
		          << mySPlot.GetYieldFromSWeight("yieldBkg") << std::endl
		          << std::endl;
	}

	// import the sWeighted dataset in the workspace
	wspace.import(*sData);

	// save it in a ROOT file for later usage
	gSystem->mkdir("../sPlot/Datasets", kTRUE);
	TFile file(sWeightedDatasetsFileName, "UPDATE");

	sData->Write();

	file.Close();

	if (BeVerbose) std::cout << "[sPlot] INFO: " << sData->GetName() << " imported into " << wspace.GetName() << " and saved in " << sWeightedDatasetsFileName << std::endl;

	return sData;
}

RooDataSet* SWeightedDataset(RooWorkspace& wspace, Int_t ptMin, Int_t ptMax, const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ExpTimesErr", bool update = false) {
	RooDataSet* sData = new RooDataSet();

	const char* fitName = GetTotalFitModelName(bkgShapeName, signalShapeName, ptMin, ptMax);

	const char* datasetName = sWeightedDatasetName(fitName);

	// recreate it upon request
	if (update) {
		if (BeVerbose) std::cout << "\n[sWeights] updating " << datasetName << " as requested\n";

		sData = CreateSWeights(wspace, ptMin, ptMax, signalShapeName, bkgShapeName);
		return sData;
	}

	// check if the sWeighted dataset has been already been imported to the workspace
	else if (wspace.data(datasetName)) {
		if (BeVerbose) std::cout << "\n[sWeights] found " << datasetName << " in workspace\n";

		sData = (RooDataSet*)wspace.data(datasetName);
	}

	// then, if it exists in the ROOT file
	else if (TFile::Open(sWeightedDatasetsFileName, "READ")) {
		TFile* file = TFile::Open(sWeightedDatasetsFileName, "READ");

		if (file->Get(datasetName)) {
			if (BeVerbose) std::cout << "\n[sWeights] found " << datasetName << " in " << sWeightedDatasetsFileName << std::endl;

			sData = (RooDataSet*)file->Get(datasetName);
			wspace.import(*sData);
			file->Close();
		}

		else {
			if (BeVerbose) std::cout << "\n[sWeights] could not find " << datasetName << ", will create it now\n";

			sData = CreateSWeights(wspace, ptMin, ptMax, signalShapeName, bkgShapeName);
		}

	}

	// otherwise, (re)create the sWeights assuming that the original dataset is in the workspace!
	else {
		if (BeVerbose) std::cout << "\n[sWeights] could not find " << datasetName << ", will create it now\n";

		sData = CreateSWeights(wspace, ptMin, ptMax, signalShapeName, bkgShapeName);
	}

	return sData;
}

RooDataSet GetSpeciesSWeightedDataset(RooDataSet* sData, const char* species = "1S") {
	return RooDataSet(sData->GetName(), sData->GetTitle(), sData, *sData->get(), nullptr, Form("yield%s_sw", species));
}

// for polarization likelihood fit
RooConstVar GetSPlotScaleFactor(RooDataSet* sData, Int_t iState = gUpsilonState) {
	double scaleFactor = 0, sumOfSWeightsSquared = 0;

	for (int i = 0; i < sData->numEntries(); i++) {
		const RooArgSet* values = sData->get(i);

		double sWeight = dynamic_cast<const RooRealVar*>(values->find(Form("yield%dS_sw", iState)))->getVal();

		scaleFactor += sWeight;
		sumOfSWeightsSquared += sWeight * sWeight;
	}

	scaleFactor /= sumOfSWeightsSquared;

	return RooConstVar("sPlotScaleFactor", "sPlot scale factor", scaleFactor);
}

#endif
