#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "RooStats/SPlot.h"

// compare the sPlot and raw yield extraction methods
void compareRawCosThetaDistrib(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180) { //possible refFrame names: CS or HX

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 10;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	/// Set up the data
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/UpsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	const char* datasetName = Form("dataset%s", refFrameName);
	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	// import the dataset to a workspace
	RooWorkspace wspace(Form("workspace_%s", refFrameName));
	wspace.import(*allDataset);

	auto* data = InvMassCosThetaPhiDataset(allDataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	std::cout << "\n------------------------------------------\nThe dataset before creating sWeights:\n";
	data->Print();

	// read variables in the reduced dataset in the workspace
	RooRealVar* invMass = wspace.var("mass");

	RooRealVar* cosTheta = wspace.var(Form("cosTheta%s", refFrameName));

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

	RooAbsPdf* signalPDF_1S = wspace.pdf("signalPDF_1S");
	RooAbsPdf* signalPDF_2S = wspace.pdf("signalPDF_2S");
	RooAbsPdf* signalPDF_3S = wspace.pdf("signalPDF_3S");

	RooRealVar* yield1S = wspace.var("yield1S");
	RooRealVar* yield2S = wspace.var("yield2S");
	RooRealVar* yield3S = wspace.var("yield3S");

	// background: error function x exponential
	const char* bkgShapeName = "ExpTimesErr";
	auto bkgModel = NominalBkgModel(wspace, bkgShapeName, nEntries);

	RooAbsPdf* bkgPDF = wspace.pdf("bkgPDF");

	RooRealVar* yieldBkg = wspace.var("yieldBkg");

	RooAddPdf* invMassModel = new RooAddPdf("fitModel", "", RooArgList(*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, *bkgPDF), {*yield1S, *yield2S, *yield3S, *yieldBkg});

	auto* fitResult = invMassModel->fitTo(*data, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), Minos(kTRUE));

	fitResult->Print("v");

	wspace.import(*invMassModel, RecycleConflictNodes());

	/// Draw the invariant mass distribution, to check the fit
	TCanvas* massCanvas = DrawMassFitDistributions(wspace, data, fitResult->floatParsFinal().getSize(), ptMin, ptMax);

	massCanvas->SaveAs(Form("InvMassFits/rawInvMassFit_%s_cent%dto%d_pt%dto%dGeV.png", bkgShapeName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	/// SPlot time!

	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments
	RooStats::SPlot sData{"sData", "", *data, invMassModel, RooArgList(*yield1S, *yield2S, *yield3S, *yieldBkg), RooArgSet(), true, false, "dataWithSWeights", Range(MassBinMin, MassBinMax), NumCPU(NCPUs)};

	std::cout << "\n\nThe dataset after creating sWeights:\n";
	data->Print();

	// check that our weights have the desired properties

	std::cout << "\n------------------------------------------\n\nCheck SWeights:" << std::endl;

	std::cout << std::endl
	          << "Yield of Y(1S) is\t" << yield1S->getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yield1S") << std::endl;

	std::cout << "Yield of background is\t" << yieldBkg->getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yieldBkg") << std::endl
	          << std::endl;

	// import this new dataset with sWeights
	std::cout << "import new dataset with sWeights" << std::endl;
	wspace.import(*data, Rename("dataWithSWeights"));

	// create weighted data sets
	RooDataSet data_weight1S{data->GetName(), data->GetTitle(), data, *data->get(), nullptr, "yield1S_sw"};

	/// Standard extraction of the raw yield per bin of cos theta

	// loop and fit

	TH1D yieldHist("yieldHist", " ", nCosThetaBins, cosThetaMin, cosThetaMax);

	Float_t cosThetaStep = ((cosThetaMax - cosThetaMin) / nCosThetaBins);

	Float_t maxYield = 0;

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		Float_t cosThetaVal = cosThetaMin + iCosTheta * cosThetaStep;

		cout << "Invariant mass fit for cos theta = [" << cosThetaVal << ", " << cosThetaVal + cosThetaStep << "]" << endl;

		std::unique_ptr<RooAbsData> reducedDataset{data->reduce({*invMass}, Form("cosTheta%s > %f && cosTheta%s < %f", refFrameName, cosThetaVal, refFrameName, cosThetaVal + cosThetaStep))};

		auto* massFitResult = invMassModel->fitTo(*reducedDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), Minos(kTRUE));

		massFitResult->Print("v");

		double rawYield = yield1S->getVal();
		double errorRawYield = yield1S->getError();

		yieldHist.SetBinContent(iCosTheta + 1, rawYield);
		yieldHist.SetBinError(iCosTheta + 1, errorRawYield);

		if (yield1S->getVal() > maxYield) maxYield = rawYield;
	}

	RooDataHist yieldRooDataHist("yieldRooDataHist", " ", *cosTheta, Import(yieldHist));

	/// Draw the cos theta distributions

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosTheta->frame(Title(" "), Range(cosThetaMin, cosThetaMax));

	data_weight1S.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(kRed), Name("sData"));

	yieldRooDataHist.plotOn(frame, DrawOption("P0Z"), MarkerColor(kAzure + 2), Name("rawYield"));

	frame->GetYaxis()->SetMaxDigits(3);
	frame->SetMaximum(2 * maxYield);

	frame->Draw();

	gPad->RedrawAxis();

	//frame->SetMaximum(nEntries / 15);

	TLegend legend(.25, .85, .5, .65);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("%d < p_{T}^{#mu#mu} < %d GeV/c, %d < #varphi_{%s} < %d", ptMin, ptMax, phiMin, refFrameName, phiMax));
	legend.AddEntry(frame->findObject("sData"), "sWeighted data", "lp");
	legend.AddEntry(frame->findObject("rawYield"), "#varUpsilon(1S) raw yield", "lp");

	legend.DrawClone();

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	gSystem->mkdir("1D", kTRUE);
	canvas->SaveAs(Form("1D/compareRawCosTheta%s_cent%dto%d_pt%dto%dGeV_phi_%dto%d.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax), "RECREATE");
}
