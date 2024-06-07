#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../Tools/RooFitPDFs/InvariantMassModels.h"

RooDataSet InvMassDataset(RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Float_t cosThetaMin = -0.1, Float_t cosThetaMax = 0.1, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180) {
	RooDataSet* allDataset = (RooDataSet*)wspace.data(RawDatasetName(""));

	if (allDataset == nullptr) {
		cerr << "Null RooDataSet provided to the reducer method!!" << endl;
		return RooDataSet();
	}

	const char* kinematicCut = Form("centrality >= %d && centrality < %d && rapidity > %f && rapidity < %f && pt > %d && pt < %d && cosTheta%s > %f && cosTheta%s < %f && %s > %d && %s < %d", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);

	RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass"))), kinematicCut);

	wspace.import(reducedDataset, Rename(Form("dataset_cosTheta_%.1fto%.1f_phi_%dto%d", cosThetaMin, cosThetaMax, phiMin, phiMax)));

	return reducedDataset;
}

void nominalFit_2D_RawDataset(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax, Bool_t isBkgExpTimesErr = kTRUE, Int_t ChebychevOrder = 1) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/UpsilonSkimmedDataset.root";

	const char* refFrameName = isCSframe ? "CS" : "HX";

	RooWorkspace wspace = SetUpWorkspace(filename, "");

	RooRealVar invMass = *wspace.var("mass");

	RooRealVar cosTheta = *wspace.var(Form("cosTheta%s", refFrameName));

	RooRealVar phi = *wspace.var(Form("phi%s", refFrameName));

	wspace.Print();

	auto allDataset = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax, "");

	Long64_t nEntries = allDataset.sumEntries();

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	const char* signalShapeName = "SymDSCB";

	// background
	
	const char* bkgShapeName = nullptr;

	if (isBkgExpTimesErr) bkgShapeName = "ExpTimesErr";		
	else bkgShapeName = Form("ChebychevOrder%d", ChebychevOrder);

	auto invMassModel = MassFitModel(wspace, signalShapeName, bkgShapeName, ptMin, ptMax, nEntries);

	RooDataSet reducedDataset = InvMassDataset(wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, refFrameName, phiMin, phiMax);

	auto* fitResult = invMassModel.fitTo(reducedDataset, Save(), Extended(kTRUE) /*, PrintLevel(-1)*/, NumCPU(NCPUs), Range(massMin, massMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

	fitResult->Print("v");

	// save the invariant mass distribution fit for further checks
	// one pad for the invariant mass data distribution with fit components, one for the pull distribution
	auto* massCanvas = new TCanvas("massCanvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);

	pad1->Draw();
	pad1->cd();

	// RooPlot* frame = InvariantMassRooPlot(wspace, reducedDataset);
	RooPlot* frame = invMass.frame(Title(" "), Range(MassBinMin, MassBinMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	reducedDataset.plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"));

	invMassModel.plotOn(frame, Components(*wspace.pdf("bkgPDF")), LineColor(kGray + 2), LineStyle(kDashed));
	invMassModel.plotOn(frame, Components(*wspace.pdf("signalPDF_1S")), LineColor(kRed));
	invMassModel.plotOn(frame, Components(*wspace.pdf("signalPDF_2S")), LineColor(kRed));
	invMassModel.plotOn(frame, Components(*wspace.pdf("signalPDF_3S")), LineColor(kRed));
	invMassModel.plotOn(frame, LineColor(kBlue));

	frame->GetYaxis()->SetMaxDigits(3);

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

	RooArgSet* signalYields = new RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S"), *wspace.var("yield3S"));

	if (!isBkgExpTimesErr) bkgShapeName = Form("ChebychevOrder%d", ChebychevOrder);
	SaveRawDataSignalYields(signalYields, bkgShapeName, fitModelName);
	SaveRawDataCanvas(massCanvas, bkgShapeName, fitModelName);
}

void scanNominalFit_lowPt_RawDataset(Int_t ptMin = 0, Int_t ptMax = 2) {
	Int_t ptEdges[8] = {0, 2, 4, 6, 8, 12, 16, 30};
	// Float_t cosThetaEdges[21] = {-1.,-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
	Float_t cosThetaEdges[11] = {-1., -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.};
	Int_t phiEdges[7] = {-180, -120, -60, 0, 60, 120, 180};
	// Int_t phiEdges[7] = {-180, -120, -60, 0, 60, 120, 180};
	// Int_t phiEdges[2] = {-180, 180};

	Int_t numPtEle = sizeof(ptEdges) / sizeof(Int_t);
	Int_t numCosThetaEle = sizeof(cosThetaEdges) / sizeof(Float_t);
	Int_t numPhiEle = sizeof(phiEdges) / sizeof(Int_t);

	for (Int_t cosThetaIdx = 0; cosThetaIdx < numCosThetaEle - 1; cosThetaIdx++) {
		for (Int_t idx = 0; idx < numPhiEle - 1; idx++) {
			nominalFit_2D_RawDataset(ptMin, ptMax, kFALSE, cosThetaEdges[cosThetaIdx], cosThetaEdges[cosThetaIdx + 1], phiEdges[idx], phiEdges[idx + 1]);
		}
	}
}
