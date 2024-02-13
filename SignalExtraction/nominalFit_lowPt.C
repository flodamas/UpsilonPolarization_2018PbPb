#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/RooFitPDFs/ErrorFuncTimesExp.h"
#include "../Tools/Style/FitDistributions.h"

#include "RooStats/SPlot.h"

RooDataSet* InvMassCosThetaWeightedDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t phiMin = 0, Int_t phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (phi%s > %d && phi%s < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, phiMin, refFrameName, phiMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass")), *(wspace.var(Form("cosTheta%s", refFrameName)))), kinematicCut);

	wspace.import(*reducedDataset, RooFit::Rename("(inv mass, cos theta) dataset"));

	return reducedDataset;
}

RooDataSet* InvMassDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Float_t cosThetaMin = -0.1, Float_t cosThetaMax = 0.1, const char* refFrameName = "CS", Int_t phiMin = 0, Int_t phiMax = 180) {
	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return nullptr;
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (cosTheta%s > %f && cosTheta%s < %f) && (phi%s > %d && phi%s < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, refFrameName, phiMin, refFrameName, phiMax);

	RooDataSet* reducedDataset = (RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass"))), kinematicCut);

	wspace.import(*reducedDataset, Rename(Form("dataset_cosTheta_%.1fto%.1f_phi_%dto%d", cosThetaMin, cosThetaMax, phiMin, phiMax)));

	return reducedDataset;
}


void nominalFit_lowPt(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Float_t binMin = 7, binMax = 13;
	Int_t nBins = 80;

	/// Set up the data
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	const char* refFrameName = isCSframe ? "CS":"HX" ;

	const char* datasetName = Form("dataset%s", refFrameName);
	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	// import the dataset to a workspace
	RooWorkspace wspace(Form("workspace_%s", datasetName));
	wspace.import(*allDataset);

	RooRealVar invMass = *wspace.var("mass");

	RooRealVar cosTheta = *wspace.var(Form("cosTheta%s", refFrameName));

	auto* data = InvMassCosThetaWeightedDataset(allDataset, wspace, ptMin, ptMax, refFrameName);

	Long64_t nEntries = data->sumEntries();
	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	const char* signalShapeName = "symDSCB";

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

	// background: Choose "ChebychevOrderN" or "ExpTimesErr"

	const char* bkgShapeName = "ChebychevOrder2";

	auto bkgModel = NominalBkgModel(wspace, bkgShapeName, nEntries);

	RooAbsPdf* bkgPDF = wspace.pdf("bkgPDF");

	RooRealVar* yieldBkg = wspace.var("yieldBkg");

	// sig + bkg model

	RooAddPdf* invMassModel = new RooAddPdf("fitModel", "", RooArgList(*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, *bkgPDF), {*yield1S, *yield2S, *yield3S, *yieldBkg});

	RooDataSet* reducedDataset = InvMassDataset(allDataset, wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, refFrameName, phiMin, phiMax);

	auto* fitResult = invMassModel->fitTo(*reducedDataset, Save(), Extended(kTRUE)/*, PrintLevel(-1)*/, NumCPU(NCPUs), Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

	fitResult->Print("v");

	auto* massCanvas = new TCanvas("massCanvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	wspace.import(*invMassModel, RecycleConflictNodes());
	// RooPlot* frame = InvariantMassRooPlot(wspace, reducedDataset);

	RooPlot* frame = (*wspace.var("mass")).frame(Title(" "), Range(MassBinMin, MassBinMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	reducedDataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"), DataError(RooAbsData::SumW2));

	invMassModel->plotOn(frame, Components(*wspace.pdf("bkgPDF")), LineColor(kGray + 2), LineStyle(kDashed));
	invMassModel->plotOn(frame, Components(*wspace.pdf("signalPDF_1S")), LineColor(kRed));
	invMassModel->plotOn(frame, Components(*wspace.pdf("signalPDF_2S")), LineColor(kRed));
	invMassModel->plotOn(frame, Components(*wspace.pdf("signalPDF_3S")), LineColor(kRed));
	invMassModel->plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

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

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	gSystem->mkdir("InvMassFits", kTRUE);
	massCanvas->SaveAs(Form("InvMassFits/CorrectedData_%s_%s.png", bkgShapeName, fitModelName), "RECREATE");
}

void scanNominalFit_lowPt(){

	Int_t ptEdges[9] = {0, 2, 4, 6, 8, 12, 16, 20, 30};
	Double_t cosThetaEdges[11] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
	Int_t phiEdges[7] = {-180, -120, -60, 0, 60, 120, 180};
	Int_t numPhiEle = sizeof(phiEdges)/sizeof(Int_t);
	for(Int_t idx =0; idx < numPhiEle-1; idx++){
		nominalFit_lowPt(ptEdges[1], ptEdges[2], kFALSE, cosThetaEdges[3], cosThetaEdges[4], phiEdges[idx], phiEdges[idx+1]);
	}
}
