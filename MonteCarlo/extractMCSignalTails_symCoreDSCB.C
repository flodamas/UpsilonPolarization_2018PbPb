#include "../Tools/BasicHeaders.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"

#include "../Tools/FitShortcuts.h"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

// crystal ball shape with symmetric Gaussian core and asymmetric tails (just like RooDSCBShape)

void extractMCSignalTails_symCoreDSCB(Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30, bool saveParams = true) {
	/// open the MC skimmed file
	const char* filename = "../Files/Y1SSelectedMCWeightedDataset.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

	const char* signalShapeName = "SymDSCB";

	// we only extract the tail parameters for a specific pt bin, they do not vary significantly with cos theta or phi within this pt bin
	// Bool_t isCSframe = kTRUE;
	Bool_t isCSframe = kFALSE;
	Float_t cosThetaMin = -1, cosThetaMax = 1;
	Float_t phiMin = -180, phiMax = 180;

	Float_t massMin = 8.5, massMax = 10.5;
	Int_t nBins = 80;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)file->Get("MCdataset");

	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	const char* refFrameName = (isCSframe) ? "CS" : "HX";

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar* massVar = wspace.var("mass");

	// fit
	auto* fitResult = SymDSCBfit(wspace, *massDataset, massMin, massMax);

	RooRealVar mean = *wspace.var(Form("mean%s", signalShapeName));
	RooRealVar sigma = *wspace.var(Form("sigma%s", signalShapeName));
	RooRealVar alphaInf = *wspace.var(Form("alphaInf%s", signalShapeName));
	RooRealVar orderInf = *wspace.var(Form("orderInf%s", signalShapeName));
	RooRealVar alphaSup = *wspace.var(Form("alphaSup%s", signalShapeName));
	RooRealVar orderSup = *wspace.var(Form("orderSup%s", signalShapeName));

	/// draw the fit to see if the fit is reasonable (we can comment it (lines 79-105) out if drawing is not necessary)
	auto* canvas = new TCanvas("canvas", "", 600, 600);
	
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = massVar->frame(Title(" "), Range(massMin, massMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	// frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset->plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	wspace.pdf(signalShapeName)->plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(centMin, centMax, ptMin, ptMax));
	//frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	frame->addObject(SymCoreDoubleCBParamsText(mean, sigma, alphaInf, orderInf, alphaSup, orderSup));
	frame->Draw();

	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->Modified();
	canvas->Update();

	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	const char* outputName = GetSignalFitName(signalShapeName, ptMin, ptMax);

	gSystem->mkdir("SignalShapeFits", kTRUE);
	canvas->SaveAs(Form("SignalShapeFits/%s.png", outputName), "RECREATE");

	if (saveParams) {
		// save signal shape parameters in a txt file to be read for data fit
		RooArgSet* tailParams = new RooArgSet(alphaInf, orderInf, alphaSup, orderSup);

		SaveMCSignalParameters(tailParams, outputName); // so that we don't have to refit later
	}

	// file->Close();
}

void scanExtractMCSignalTails_symCoreDSCB() {
	Int_t PtEdges[8] = {0, 2, 4, 6, 8, 12, 16, 30};
	Int_t NumPtEle = sizeof(PtEdges) / sizeof(Int_t);
	for (Int_t idx = 0; idx < NumPtEle - 1; idx++) {
		extractMCSignalTails_symCoreDSCB(gCentralityBinMin, gCentralityBinMax, PtEdges[idx], PtEdges[idx + 1]);
	}
}
