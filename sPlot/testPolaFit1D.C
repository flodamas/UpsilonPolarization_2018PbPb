#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "RooStats/SPlot.h"

// reduce the whole dataset (N dimensions)
RooDataSet* InvMassCosThetaWeightedDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS") {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass")), *(wspace.var(Form("cosTheta%s", refFrameName)))), kinematicCut);

	wspace.import(*reducedDataset, RooFit::Rename("(inv mass, cos theta) dataset"));

	return reducedDataset;
}

// based on the tutorial https://github.com/root-project/root/blob/master/tutorials/roostats/rs301_splot.C
void testPolaFit1D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t nCosThetaBins = 10, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root") {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Float_t phiMin = 0, phiMax = 180;

	/// Set up the data
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	// Read skimmed dataset (contains angular distributions in CS and HX after kinematic cuts due to acceptance)
	RooDataSet* allDatasetCS = (RooDataSet*)f->Get("datasetCS");

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDatasetCS);

	auto* data = InvMassCosThetaWeightedDataset(allDatasetCS, wspace, ptMin, ptMax, refFrameName);

	std::cout << "\n------------------------------------------\nThe dataset before creating sWeights:\n";
	data->Print();

	// read variables in the reduced dataset in the workspace
	RooRealVar* invMass = wspace.var("mass");

	RooRealVar cosTheta = *wspace.var(Form("cosTheta%s", refFrameName));

	Long64_t nEntries = data->sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	const char* signalShapeName = "SymDSCB";

	// get the tail parameters of the signal shape first in case the MC fit is needed
	RooRealVar* alphaInf = new RooRealVar("alphaInf", "", 1);
	RooRealVar* orderInf = new RooRealVar("orderInf", "", 1);
	RooRealVar* alphaSup = new RooRealVar("alphaSup", "", 1);
	RooRealVar* orderSup = new RooRealVar("orderSup", "", 1);

	RooArgSet tailParams = GetMCSignalTailParameters(alphaInf, orderInf, alphaSup, orderSup, signalShapeName, ptMin, ptMax);

	auto signalModel = NominalSignalModel(wspace, alphaInf, orderInf, alphaSup, orderSup, nEntries);

	RooAbsPdf* signalMassPDF_1S = wspace.pdf("signalPDF_1S");
	RooAbsPdf* signalMassPDF_2S = wspace.pdf("signalPDF_2S");
	RooAbsPdf* signalMassPDF_3S = wspace.pdf("signalPDF_3S");

	RooRealVar yield1S = *wspace.var("yield1S");
	RooRealVar yield2S = *wspace.var("yield2S");
	RooRealVar yield3S = *wspace.var("yield3S");

	// background: Chebychev polynomial

	int order = 3;

	RooArgList coefList = ChebychevCoefList(order);

	RooChebychev bkgMassPDF("bkgPDF", " ", *invMass, coefList);

	RooRealVar yieldBkg("yieldBkg", "N background events", 0, nEntries);

	RooAddPdf* invMassModel = new RooAddPdf("fitModel", "", RooArgList(*signalMassPDF_1S, *signalMassPDF_2S, *signalMassPDF_3S, bkgMassPDF), {yield1S, yield2S, yield3S, yieldBkg});

	/// SPlot time!

	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments
	RooStats::SPlot sData{"sData", "", *data, invMassModel, RooArgList(yield1S, yield2S, yield3S, yieldBkg), RooArgSet(), true, false, "dataWithSWeights", Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError), NumCPU(NCPUs)};

	std::cout << "\n\nThe dataset after creating sWeights:\n";
	data->Print();

	// check that our weights have the desired properties

	std::cout << "\n------------------------------------------\n\nCheck SWeights:" << std::endl;

	std::cout << std::endl
	          << "Yield of Y(1S) is\t" << yield1S.getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yield1S") << std::endl;

	std::cout << "Yield of background is\t" << yieldBkg.getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yieldBkg") << std::endl
	          << std::endl;

	// import this new dataset with sWeights
	std::cout << "import new dataset with sWeights" << std::endl;
	wspace.import(*data, Rename("dataWithSWeights"));

	/// Draw the cos theta distribution with and without sWeights

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta.frame(Title(" "), Range(cosThetaMin, cosThetaMax));
	frame->SetXTitle(Form("cos #theta_{%s}", refFrameName));
	data->plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), Name("data"), DataError(RooAbsData::SumW2));

	// create weighted data sets
	RooDataSet data_weight1S{data->GetName(), data->GetTitle(), data, *data->get(), nullptr, "yield1S_sw"};

	data_weight1S.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(kRed), DataError(RooAbsData::SumW2), Name("data1S"));

	frame->GetYaxis()->SetMaxDigits(3);

	//CMS_lumi(canvas, gCMSLumiText);

	// polarization fit, just to give it a try but weird results obtained so far, WIP
	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -2., 2.);
	RooRealVar normCosTheta("normCosTheta", "normCosTheta", 0, nEntries);
	RooGenericPdf cosThetaPDF_1S("cosThetaPDF_1S", "cosThetaPDF_1S", "(@0/(3 + @1)) * (1 + @1*@2*@2)", {normCosTheta, lambdaTheta, cosTheta});

	auto* polarizationFitResult = cosThetaPDF_1S.fitTo(data_weight1S, Save(), Extended(kTRUE), PrintLevel(+1), NumCPU(NCPUs), AsymptoticError(DoAsymptoticError), Range(cosThetaMin, cosThetaMax));

	polarizationFitResult->Print("v");

	cosThetaPDF_1S.plotOn(frame, LineColor(kBlue));

	frame->Draw();

	gPad->RedrawAxis();

	frame->SetMaximum(nEntries / 30);

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.05);
	text.DrawLatexNDC(.55, .85, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	//text.DrawLatexNDC(.48, .8, Form("%.0f < m_{#mu#mu} < %.0f GeV/c^{2}", massMin, massMax));

	TLegend legend(.3, .8, .5, .65);
	legend.SetTextSize(.05);

	legend.AddEntry(frame->findObject("data"), "corrected dimuon events", "lp");
	legend.AddEntry(frame->findObject("data1S"), "with #varUpsilon(1S) yield sWeights", "lp");
	legend.DrawClone();

	gPad->Update();
}
