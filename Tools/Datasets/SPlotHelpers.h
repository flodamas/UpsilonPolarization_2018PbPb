#include "../../Tools/BasicHeaders.h"

#include "../../AnalysisParameters.h"

#include "../../Tools/FitShortcuts.h"
#include "../../Tools/Style/FitDistributions.h"

#include "RooStats/SPlot.h"

using namespace RooStats;

/// sPlot class documentation https://root.cern/doc/master/classRooStats_1_1SPlot.html

SPlot CreateSWeights(RooWorkspace& wspace, RooDataSet& data) {
	// check first that the invariant mass fit model is present in the workspace

	if (!wspace.pdf("invMassModel")) {
		cout << "\n ERROR: cannot find the invariant mass fit model in the provided workspace for sPlot!\n";
		return SPlot();
	}

	if (BeVerbose) {
		std::cout << "\n---------------------sWeighted data based on yields-------------------\n\nDataset content before creating sWeights:\n";

		data.Print();
	}

	auto invMassModel = *(RooAddPdf*)wspace.pdf("invMassModel");

	// need a "standard fit" first
	auto* fitResult = RawInvariantMassFit(data, invMassModel, RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S")));
	/*
	// draw and save the invariant mass distribution, to check the fit
	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	gSystem->mkdir("InvMassFits", kTRUE);
	massCanvas->SaveAs(Form("InvMassFits/rawInvMassFit_%s_cent%dto%d_pt%dto%dGeV.png", bkgShapeName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
*/
	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments

	RooRealVar yield1S = *wspace.var("yield1S");
	RooRealVar yield2S = *wspace.var("yield2S");
	RooRealVar yield3S = *wspace.var("yield3S");
	RooRealVar yieldBkg = *wspace.var("yieldBkg");

	RooStats::SPlot mySPlot{"mySPlot", "", data, &invMassModel, RooArgList(yield1S, yield2S, yield3S, yieldBkg), RooArgSet(), true, false, "", Range(MassBinMin, MassBinMax), NumCPU(NCPUs)};

	if (BeVerbose) {
		std::cout << "\nAfter creating sWeights:\n";
		data.Print();

		// check that our weights have the desired properties

		std::cout << "\nChecking sWeights:" << std::endl;

		std::cout << "\nYield of Y(1S) is\t" << yield1S.getVal() << ".  From sWeights it is "
		          << mySPlot.GetYieldFromSWeight("yield1S") << std::endl;

		std::cout << "Yield of background is\t" << yieldBkg.getVal() << ".  From sWeights it is "
		          << mySPlot.GetYieldFromSWeight("yieldBkg") << std::endl
		          << std::endl;
	}

	// import the sWeighted dataset in the workspace
	wspace.import(data, RooFit::Rename(Form("sWeighted%s", data.GetName())));

	// return the sPlot instance, in case we need to access some info from outside
	return mySPlot;
}

RooDataSet GetSWeightedDataset(RooDataSet* data, const char* species = "1S") {
	// TO-DO: modify the code such that this function calls 'CreateSWeights' if the sWeighted RooDataSet is not found in the workspace (to be parsed as input)
	return RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), nullptr, Form("yield%s_sw", species));
}
