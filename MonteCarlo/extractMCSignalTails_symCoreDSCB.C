#include "../Tools/BasicHeaders.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"

#include "../Tools/FitShortcuts.h"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

// crystal ball shape with symmetric Gaussian core and asymmetric tails (just like RooDSCBShape)

void extractMCSignalTails_symCoreDSCB(Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30, const char* filename = "../Files/Y1SReconstructedMCWeightedDataset_TriggerAcc_Lambda_Theta0.00_Phi0.00_ThetaPhi0.00.root", Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Float_t phiMin = -180, Float_t phiMax = 180, bool saveParams = true) {
	/// open the MC skimmed file

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}

	cout << "File " << filename << " opened\n";

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

	const char* signalShapeName = "SymDSCB";

	// we only extract the tail parameters for a specific pt bin, they do not vary significantly with cos theta or phi within this pt bin
	// Bool_t isCSframe = kTRUE;
	// Float_t cosThetaMin = -1, cosThetaMax = 1;
	// Float_t phiMin = -180, phiMax = 180;

	Float_t massMin = 8.5, massMax = 10.5;
	Int_t nBins = 80;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* refFrameName = (isCSframe) ? "CS" : "HX";

	RooDataSet* allDataset = (RooDataSet*)file->Get(Form("MCdataset%s", refFrameName));

	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooDataSet massDataset = ReducedMassDataset(allDataset, wspace, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar massVar = *wspace.var("mass");

	// fit
	auto* fitResult = SymDSCBfit(wspace, massDataset, massMin, massMax);

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

	RooPlot* frame = massVar.frame(Title(" "), Range(massMin, massMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	// frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset.plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	wspace.pdf(signalShapeName)->plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(centMin, centMax, ptMin, ptMax));
	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	frame->addObject(SymCoreDoubleCBParamsText(mean, sigma, alphaInf, orderInf, alphaSup, orderSup));
	frame->Draw();

	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->Modified();
	canvas->Update();

	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	// const char* outputName = GetSignalFitName(signalShapeName, ptMin, ptMax);
	const char* outputName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);
	//cout << outputName << endl;
	gSystem->mkdir("SignalShapeFits", kTRUE);
	canvas->SaveAs(Form("SignalShapeFits/%s.png", outputName), "RECREATE");

	if (saveParams) {
		// save signal shape parameters in a txt file to be read for data fit
		RooArgSet* Params = new RooArgSet(sigma, alphaInf, orderInf, alphaSup, orderSup);

		SaveMCSignalParameters(Params, outputName); // so that we don't have to refit later
	}

	// file->Close();
}

void scanExtractMCSignalTails_symCoreDSCB(Bool_t isCSframe = kFALSE) {
	Float_t cosThetaEdges[6] = {-0.7, -0.42, -0.14, 0.14, 0.42, 0.7};
	Float_t phiEdges[7] = {-180, -120, -60, 0, 60, 120, 180};

	Int_t NumCosThetaEle = sizeof(cosThetaEdges) / sizeof(Float_t);
	Int_t NumPhiEle = sizeof(phiEdges) / sizeof(Float_t);

	for (Int_t ptIdx = 0; ptIdx < NPtBins; ptIdx++) {
		for (Int_t cosThetaIdx = 0; cosThetaIdx < NumCosThetaEle - 1; cosThetaIdx++) {
			for (Int_t phiIdx = 0; phiIdx < NumPhiEle - 1; phiIdx++) {
				extractMCSignalTails_symCoreDSCB(gCentralityBinMin, gCentralityBinMax, gPtBinning[ptIdx], gPtBinning[ptIdx + 1], "../Files/Y1SReconstructedMCWeightedDataset_TriggerAcc_Lambda_Theta0.00_Phi0.00_ThetaPhi0.00.root", isCSframe, cosThetaEdges[cosThetaIdx], cosThetaEdges[cosThetaIdx + 1], phiEdges[phiIdx], phiEdges[phiIdx + 1], true);
			}
		}
	}
}
