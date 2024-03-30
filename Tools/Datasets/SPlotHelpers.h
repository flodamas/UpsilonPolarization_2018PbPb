#include "../../Tools/BasicHeaders.h"

#include "../../AnalysisParameters.h"

#include "RooStats/SPlot.h"

using namespace RooStats;

/// sPlot class documentation https://root.cern/doc/master/classRooStats_1_1SPlot.html

SPlot CreateSPlot(RooWorkspace& wspace, RooDataSet& data, RooAddPdf fitModel) {
	if (BeVerbose) {
		std::cout << "\n--------------------------sWeighted data based on yields------------------------\n\nDataset content before creating sWeights:\n";

		data.Print();
	}

	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments

	RooRealVar yield1S = *wspace.var("yield1S");
	RooRealVar yield2S = *wspace.var("yield2S");
	RooRealVar yield3S = *wspace.var("yield3S");
	RooRealVar yieldBkg = *wspace.var("yieldBkg");

	RooStats::SPlot mySPlot{"mySPlot", "", data, &fitModel, RooArgList(yield1S, yield2S, yield3S, yieldBkg), RooArgSet(), true, false, "", Range(MassBinMin, MassBinMax), NumCPU(NCPUs)};

	if (BeVerbose) {
		std::cout << "\nAfter creating sWeights:\n";
		data.Print();

		// check that our weights have the desired properties

		std::cout << "\nChecking sWeights:" << std::endl;

		std::cout << std::endl
		          << "Yield of Y(1S) is\t" << yield1S.getVal() << ".  From sWeights it is "
		          << mySPlot.GetYieldFromSWeight("yield1S") << std::endl;

		std::cout << "Yield of background is\t" << yieldBkg.getVal() << ".  From sWeights it is "
		          << mySPlot.GetYieldFromSWeight("yieldBkg") << std::endl
		          << std::endl;
	}

	// return the sPlot instance, in case we need to access some info from outside
	return mySPlot;
}

RooDataSet GetSWeightedDataset(RooDataSet* data, const char* species = "1S") {
	return RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), nullptr, Form("yield%s_sw", species));
}
