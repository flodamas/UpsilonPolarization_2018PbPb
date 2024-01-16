#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "RooStats/SPlot.h"

// reduce the whole dataset (N dimensions)
RooDataSet sWeightedCosThetaDataset(Int_t ptMin = 0, Int_t ptMax = 30, const char* filename = "../Files/UpsilonSkimmedDataset.root", const char* datasetName = "datasetCS") {
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}

	cout << "File " << filename << " opened" << endl;

	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	// import the dataset to a workspace
	RooWorkspace wspace(Form("workspace_%s", datasetName));
	wspace.import(*allDataset);

	// read variables in the reduced dataset in the workspace
	RooRealVar invMass = *wspace.var("mass");

	RooRealVar cosThetaCS = *wspace.var("cosThetaCS");

	RooRealVar phiCS = *wspace.var("phiCS");

	RooRealVar pt = *wspace.var("pt");

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(invMass, cosThetaCS, phiCS, pt), kinematicCut);

	//	wspace.import(*reducedDataset);

	Long64_t nEntries = reducedDataset->sumEntries();

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

	auto* fitResult = invMassModel->fitTo(*reducedDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

	fitResult->Print("v");

	/// SPlot time!

	// yields from invariant mass distribution fit as sWeights
	// see https://root.cern/doc/master/classRooStats_1_1SPlot.html#a5b30f5b1b2a3723bbebef17ffb6507b2 constructor for the arguments
	RooStats::SPlot sData{"sData", "", *reducedDataset, invMassModel, RooArgList(*yield1S, *yield2S, *yield3S, yieldBkg), RooArgSet(), true, false, "dataWithSWeights", Range(MassBinMin, MassBinMax), NumCPU(NCPUs)};

	std::cout << "\n\nThe dataset after creating sWeights:\n";
	reducedDataset->Print();

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
	wspace.import(*reducedDataset, Rename("dataWithSWeights"));

	// create weighted data sets
	RooDataSet data_weight1S{reducedDataset->GetName(), reducedDataset->GetTitle(), reducedDataset, *reducedDataset->get(), nullptr, "yield1S_sw"};

	//delete wspace;

	return data_weight1S;
}

RooDataSet* InvMassRawDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Float_t cosThetaMin = -0.1, Float_t cosThetaMax = 0.1) {
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
void compareCosThetaDistrib(Int_t ptMin = 0, Int_t ptMax = 30) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	cout << endl
	     << "************* sPlot of corrected data ************" << endl
	     << endl;

	RooDataSet correctedBeforeDataset = sWeightedCosThetaDataset(ptMin, ptMax, "../Files/WeightedUpsilonSkimmedDataset.root", "datasetCS");

	cout << endl
	     << "************* sPlot of raw data ************" << endl
	     << endl;

	// to be corrected for acceptance x efficiency
	RooDataSet sWeightedRawDataset = sWeightedCosThetaDataset(ptMin, ptMax, "../Files/UpsilonSkimmedDataset.root", "dataset");

	Int_t nCosThetaBins = 10;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	RooRealVar cosThetaCS("cosThetaCS", "cos theta in the Collins-Soper frame", cosThetaMin, cosThetaMax);

	RooRealVar correctionCSVar("correctionCS", "1 / (acc x eff) weight in the CS frame", -1000000, 1000000);

	RooDataSet correctedAfterDataset("correctedAfterDataset", " ", RooArgSet(cosThetaCS, correctionCSVar), RooFit::WeightVar("correctionCS"), RooFit::StoreAsymError(RooArgSet(correctionCSVar)));

	// acceptance maps
	TFile* acceptanceFile = TFile::Open("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root", "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return;
	}

	auto* accMapCS = (TEfficiency*)acceptanceFile->Get("AccMatrixCS");

	// efficiency maps
	TFile* efficiencyFile = TFile::Open("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root", "READ");
	if (!efficiencyFile) {
		cout << "Efficiency file not found. Check the directory of the file." << endl;
		return;
	}

	auto* effMapCS = (TEfficiency*)efficiencyFile->Get("NominalEff_CS");
	auto* systEffCS = (TH3D*)efficiencyFile->Get("RelatSystEff_CS");
	auto* effMapHX = (TEfficiency*)efficiencyFile->Get("NominalEff_HX");
	auto* systEffHX = (TH3D*)efficiencyFile->Get("RelatSystEff_HX");

	// run over the sWeighted dataset and fill a new one with corrected data
	Double_t weightCS = 0;
	Double_t errorWeightLowCS = 0, errorWeightHighCS = 0;

	// Retrieve and parse event buffer before loop
	RooArgSet* obs = (RooArgSet*)sWeightedRawDataset.get();
	RooRealVar* tempCosThetaCS = (RooRealVar*)obs->find("cosThetaCS");
	RooRealVar* tempPhiCS = (RooRealVar*)obs->find("phiCS");
	RooRealVar* tempPt = (RooRealVar*)obs->find("pt");

	for (int i = 0; i < sWeightedRawDataset.numEntries(); i++) {
		sWeightedRawDataset.get(i);

		//cout << tempCosThetaCS->getVal() << ": weight = " << sWeightedRawDataset.weight() << endl;

		// get the corresponding weights
		int accBinCS = accMapCS->FindFixBin(tempCosThetaCS->getVal(), tempPhiCS->getVal(), tempPt->getVal());
		double acceptanceCS = accMapCS->GetEfficiency(accBinCS);

		int effBinCS = effMapCS->FindFixBin(tempCosThetaCS->getVal(), tempPhiCS->getVal(), tempPt->getVal());
		double efficiencyCS = effMapCS->GetEfficiency(effBinCS);

		weightCS = ((acceptanceCS == 0) || (efficiencyCS == 0)) ? 0 : sWeightedRawDataset.weight() / (acceptanceCS * efficiencyCS); // IMPORTANT!

		//		weightCS = (efficiencyCS == 0) ? 0 : sWeightedRawDataset.weight() / (efficiencyCS); // IMPORTANT!

		// propagate both scale factor uncertainties and efficiency stat errors to the weight
		errorWeightLowCS = weightCS * TMath::Hypot(systEffCS->GetBinContent(effBinCS), effMapCS->GetEfficiencyErrorUp(effBinCS) / efficiencyCS);

		errorWeightHighCS = weightCS * TMath::Hypot(systEffCS->GetBinContent(effBinCS), effMapCS->GetEfficiencyErrorLow(effBinCS) / efficiencyCS);

		correctionCSVar = weightCS;

		correctionCSVar.setAsymError(errorWeightLowCS, errorWeightHighCS);

		correctedAfterDataset.add(RooArgSet(*tempCosThetaCS, correctionCSVar), weightCS, errorWeightLowCS, errorWeightHighCS);
	}

	/// Last check: extract the raw yields for each cos theta bins, correct for acceptance times efficiency, and plot the resulting distribution

	const char* filename = "../Files/UpsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}

	cout << "File " << filename << " opened" << endl;

	RooDataSet* allDataset = (RooDataSet*)f->Get("dataset");

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	// read variables in the reduced dataset in the workspace
	RooRealVar invMass = *wspace.var("mass");

	Long64_t nEntries = 3e4;

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

	int order = 1;

	RooArgList coefList = ChebychevCoefList(order);

	RooChebychev bkgPDF("bkgPDF", " ", invMass, coefList);

	RooRealVar yieldBkg("yieldBkg", "N background events", 0, nEntries);

	RooAddPdf* invMassModel = new RooAddPdf("fitModel", "", RooArgList(*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, bkgPDF), {*yield1S, *yield2S, *yield3S, yieldBkg});

	// get the correction histograms (1D!)
	auto* efficiency1D = (TEfficiency*)efficiencyFile->Get(Form("hEffCosThetaCS_pt%dto%d", ptMin, ptMax));

	auto* acceptance1D = (TEfficiency*)acceptanceFile->Get(Form("AccCosThetaCS_pt%dto%d", ptMin, ptMax));

	// loop and fit

	TH1D standardCorrectedHist("standardCorrectedHist", " ", nCosThetaBins, cosThetaMin, cosThetaMax);

	Float_t cosThetaStep = ((cosThetaMax - cosThetaMin) / nCosThetaBins);

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		Float_t cosThetaVal = cosThetaMin + iCosTheta * cosThetaStep;

		cout << "Invariant mass fit for cos theta = [" << cosThetaVal << ", " << cosThetaVal + cosThetaStep << "]" << endl;

		RooDataSet* reducedDataset = InvMassRawDataset(allDataset, wspace, ptMin, ptMax, cosThetaVal, cosThetaVal + cosThetaStep);

		auto* fitResult = invMassModel->fitTo(*reducedDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), Minos(kTRUE));

		fitResult->Print("v");

		double rawYield = yield1S->getVal();
		double errorRawYield = yield1S->getError();

		int globalBin = efficiency1D->GetGlobalBin(iCosTheta + 1);

		double effValue = efficiency1D->GetEfficiency(globalBin);

		double accValue = acceptance1D->GetEfficiency(globalBin);

		double normYield = (effValue == 0 || accValue == 0) ? 0 : rawYield / (effValue * accValue);

		standardCorrectedHist.SetBinContent(iCosTheta + 1, normYield);
		standardCorrectedHist.SetBinError(iCosTheta + 1, errorRawYield / (effValue * accValue));
	}

	RooDataHist correctedHist("correctedHist", " ", cosThetaCS, Import(standardCorrectedHist));

	/// Draw the cos theta distributions

	TCanvas* canvas = new TCanvas("canvas", "canvas", 650, 600);

	RooPlot* frame = cosThetaCS.frame(Title(" "), Range(cosThetaMin, cosThetaMax));
	frame->SetXTitle("cos #theta_{CS}");

	correctedBeforeDataset.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(kRed), DataError(RooAbsData::SumW2), Name("dataBefore"));

	correctedAfterDataset.plotOn(frame, Binning(nCosThetaBins), DrawOption("P0Z"), MarkerColor(kAzure + 2), DataError(RooAbsData::SumW2), Name("dataAfter"));

	correctedHist.plotOn(frame, DrawOption("P0Z"), MarkerColor(kGreen + 2), Name("corrected"));

	frame->GetYaxis()->SetMaxDigits(3);

	frame->Draw();

	gPad->RedrawAxis();

	//frame->SetMaximum(nEntries / 15);

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.05);
	text.DrawLatexNDC(.55, .85, Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	//text.DrawLatexNDC(.48, .8, Form("%.0f < m_{#mu#mu} < %.0f GeV/c^{2}", massMin, massMax));

	TLegend legend(.3, .8, .55, .65);
	legend.SetTextSize(.05);
	//legend.SetHeader("sWeighted data, efficiency-corrected");
	legend.AddEntry(frame->findObject("dataBefore"), "2D correction + sPlot", "lp");
	legend.AddEntry(frame->findObject("dataAfter"), "sPlot + 2D correction", "lp");
	legend.AddEntry(frame->findObject("corrected"), "raw yield + 1D correction", "lp");

	legend.DrawClone();

	gPad->Update();

	//CMS_lumi(canvas, gCMSLumiText);
	canvas->SaveAs(Form("1D/compareCosThetaCS_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
}
