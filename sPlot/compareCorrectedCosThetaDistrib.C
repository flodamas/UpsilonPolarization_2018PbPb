#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "RooStats/SPlot.h"

RooDataSet* InvMassCosThetaWeightedDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass")), *(wspace.var("cosThetaCS"))), kinematicCut);

	wspace.import(*reducedDataset, RooFit::Rename("(inv mass, cos theta) dataset"));

	return reducedDataset;
}

RooDataSet* InvMassDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Float_t cosThetaMin = -0.1, Float_t cosThetaMax = 0.1) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (cosThetaCS > %f && cosThetaCS < %f)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, cosThetaMin, cosThetaMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass"))), kinematicCut);

	wspace.import(*reducedDataset, Rename(Form("dataset_cosTheta_%.1fto%.1f", cosThetaMin, cosThetaMax)));

	return reducedDataset;
}

// compare the resulting distributions before and after acc x eff correction
void compareCorrectedCosThetaDistrib(Int_t ptMin = 0, Int_t ptMax = 30, const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root", const char* datasetName = "datasetCS") {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 10;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	/// Set up the data
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}

	cout << "File " << filename << " opened" << endl;

	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	// import the dataset to a workspace
	RooWorkspace wspace(Form("workspace_%s", datasetName));
	wspace.import(*allDataset);

	RooRealVar invMass = *wspace.var("mass");

	RooRealVar cosThetaCS = *wspace.var("cosThetaCS");

	auto* data = InvMassCosThetaWeightedDataset(allDataset, wspace, ptMin, ptMax);

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

	// background: Chebychev polynomial

	int order = 2;

	RooArgList coefList = ChebychevCoefList(order);

	RooChebychev bkgPDF("bkgPDF", " ", invMass, coefList);

	RooRealVar yieldBkg("yieldBkg", "N background events", 0, nEntries);

	RooAddPdf* invMassModel = new RooAddPdf("fitModel", "", RooArgList(*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, bkgPDF), {*yield1S, *yield2S, *yield3S, yieldBkg});
	wspace.import(*invMassModel, RecycleConflictNodes());

	/// "Standard" procedure: extract the yields per bin

	TH1D standardCorrectedHist("standardCorrectedHist", " ", nCosThetaBins, cosThetaMin, cosThetaMax);

	Float_t cosThetaStep = ((cosThetaMax - cosThetaMin) / nCosThetaBins);

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		Float_t cosThetaVal = cosThetaMin + iCosTheta * cosThetaStep;

		cout << "Invariant mass fit for cos theta = [" << cosThetaVal << ", " << cosThetaVal + cosThetaStep << "]" << endl;

		RooDataSet* reducedDataset = InvMassDataset(allDataset, wspace, ptMin, ptMax, cosThetaVal, cosThetaVal + cosThetaStep);

		auto* fitResult = invMassModel->fitTo(*reducedDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

		fitResult->Print("v");

		// save the invariant mass distribution fit for further checks
		// one pad for the invariant mass data distribution with fit components, one for the pull distribution
		TCanvas* massCanvas = new TCanvas("massCanvas", "", 600, 600);
		TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
		pad1->SetBottomMargin(0.03);
		pad1->Draw();
		pad1->cd();

		RooPlot* frame = InvariantMassRooPlot(wspace, reducedDataset);

		frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

		//frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

		frame->addObject(FitResultText(*wspace.var("yield1S"), ComputeSignalSignificance(wspace, 1), *wspace.var("yield2S"), ComputeSignalSignificance(wspace, 2)));

		frame->Draw();

		gPad->RedrawAxis();

		// pull distribution
		massCanvas->cd();

		TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

		//canvas->Modified();
		//canvas->Update();
		massCanvas->cd();
		pad1->Draw();
		pad2->Draw();

		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, true, cosThetaVal, cosThetaVal + cosThetaStep, 0, 180);

		massCanvas->SaveAs(Form("InvMassFits/CorrectedData_ChebychevOrder%d_%s.png", order, fitModelName), "RECREATE");

		standardCorrectedHist.SetBinContent(iCosTheta + 1, yield1S->getVal());
		standardCorrectedHist.SetBinError(iCosTheta + 1, yield1S->getError());
	}

	RooDataHist correctedHist("correctedHist", " ", cosThetaCS, Import(standardCorrectedHist));

	/// SPlot time!

	// copy for sPlot
	RooDataSet dataForsPlot = RooDataSet(*data, "dataForsPlot");

	RooWorkspace wsPlot("wsPlot");
	wsPlot.import(dataForsPlot);

	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments
	RooStats::SPlot sData{"sData", "", dataForsPlot, invMassModel, RooArgList(*yield1S, *yield2S, *yield3S, yieldBkg), RooArgSet(), true, false, "dataWithSWeights", Range(MassBinMin, MassBinMax), NumCPU(NCPUs), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError)};

	std::cout << "\n\nThe dataset after creating sWeights:\n";
	dataForsPlot.Print();

	// check that our weights have the desired properties

	std::cout << "\n------------------------------------------\n\nCheck SWeights:" << std::endl;

	std::cout << std::endl
	          << "Yield of Y(1S) is\t" << yield1S->getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yield1S") << std::endl;

	std::cout << "Yield of background is\t" << yieldBkg.getVal() << ".  From sWeights it is "
	          << sData.GetYieldFromSWeight("yieldBkg") << std::endl
	          << std::endl;

	// import this new dataset with sWeights
	//std::cout << "import new dataset with sWeights" << std::endl;
	//wspace.import(*reducedDataset, Rename("dataWithSWeights"));

	// create weighted data sets
	RooDataSet data_weight1S{dataForsPlot.GetName(), dataForsPlot.GetTitle(), &dataForsPlot, *dataForsPlot.get(), nullptr, "yield1S_sw"};

	/// Draw the cos theta distributions

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosThetaCS.frame(Title(" "), Range(cosThetaMin, cosThetaMax));
	frame->SetXTitle("cos #theta_{CS}");

	data_weight1S.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(kRed), DataError(RooAbsData::SumW2), Name("fromsPlot"));

	correctedHist.plotOn(frame, DrawOption("P0Z"), MarkerColor(kAzure + 2), Name("standard"));

	frame->GetYaxis()->SetMaxDigits(3);

	frame->Draw();

	gPad->RedrawAxis();

	//frame->SetMaximum(nEntries / 15);

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.05);
	text.DrawLatexNDC(.55, .85, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	//text.DrawLatexNDC(.48, .8, Form("%.0f < m_{#mu#mu} < %.0f GeV/c^{2}", massMin, massMax));

	TLegend legend(.25, .8, .5, .63);
	legend.SetTextSize(.05);
	legend.SetHeader("dimuon events corrected for acc x eff");
	legend.AddEntry(frame->findObject("fromsPlot"), "sWeighted with #varUpsilon(1S) yield", "lp");
	legend.AddEntry(frame->findObject("standard"), "standard yield extraction", "lp");

	legend.DrawClone();

	gPad->Update();

	//CMS_lumi(canvas, gCMSLumiText);
	canvas->SaveAs(Form("1D/compareCorrectedCosThetaCS_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
}