#include "../Tools/BasicHeaders.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"

#include "../Tools/FitShortcuts.h"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../MonteCarlo/AccEffHelpers.h"

// crystal ball shape with symmetric Gaussian core and asymmetric tails (just like RooDSCBShape)

void extractMCSignalTails_symCoreDSCB(Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30, const char* muonAccName = "UpsilonTriggerThresholds", Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Float_t phiMin = -180, Float_t phiMax = 180, bool saveParams = true, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0) {
	/// open the MC skimmed file

	const char* filename = Form("../Files/Y1SReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", muonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}

	cout << "File " << filename << " opened\n";

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

	const char* signalShapeName = "SymDSCB";

	Float_t massMin = 8.8, massMax = 10.2;
	// Float_t massMin = 8.65, massMax = 10.5;
	Int_t nBins = 40;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* refFrameName = (isCSframe) ? "CS" : "HX";

	RooDataSet* allDataset = (RooDataSet*)file->Get(Form("MCdataset%s", refFrameName));

	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooDataSet massDataset = ReducedMassDataset(allDataset, wspace, ptMin, ptMax, refFrameName, massMin, massMax, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar massVar = *wspace.var("mass");
	massVar.setRange("MassFitRange", massMin, massMax);

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
	pad1->SetRightMargin(0.035);
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = massVar.frame(Title(" "), Range(massMin, massMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	// frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset.plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	wspace.pdf(signalShapeName)->plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(centMin, centMax, ptMin, ptMax));
	frame->addObject(RefFrameTextPhiFolded(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	frame->addObject(SymCoreDoubleCBParamsText(mean, sigma, alphaInf, orderInf, alphaSup, orderSup));
	frame->GetYaxis()->SetMaxDigits(3);

	frame->Draw();

	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->Modified();
	canvas->Update();

	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	const char* path = Form("SignalShapeFits/%s/", muonAccName);
	gSystem->mkdir(path, kTRUE);

	// const char* outputName = GetSignalFitName(signalShapeName, ptMin, ptMax);
	const char* outputName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);
	//cout << outputName << endl;
	canvas->SaveAs(Form("%s/%s.png", path, outputName), "RECREATE");

	if (saveParams) {
		// save signal shape parameters in a txt file to be read for data fit
		RooArgSet* Params = new RooArgSet(sigma, alphaInf, orderInf, alphaSup, orderSup);

		const char* fullPath = Form("SignalParameters/%s/%s.txt", muonAccName, outputName);

		gSystem->mkdir(Form("SignalParameters/%s/", muonAccName), kTRUE);

		SaveMCSignalParameters(Params, fullPath); // so that we don't have to refit later
	}

	// file->Close();
}

void scanExtractMCSignalTails_symCoreDSCB(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180, const char* muonAccName = "UpsilonTriggerThresholds") {
	std::vector<Double_t> cosThetaEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	std::vector<Double_t> phiEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	for (Int_t cosThetaIdx = 0; cosThetaIdx < nCosThetaBins; cosThetaIdx++) {
		for (Int_t phiIdx = 0; phiIdx < nPhiBins; phiIdx++) {
			extractMCSignalTails_symCoreDSCB(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, muonAccName, isCSframe, cosThetaEdges[cosThetaIdx], cosThetaEdges[cosThetaIdx + 1], phiEdges[phiIdx], phiEdges[phiIdx + 1], true);
		}
	}
}
