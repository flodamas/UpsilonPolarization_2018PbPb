// #include "../Tools/Style/tdrStyle.C"
// #include "../Tools/Style/CMS_lumi.C"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/FitShortcuts.h"

// #include "../Tools/Parameters/PhysicsConstants.h"

// crystal ball shape with symmetric Gaussian core and asymmetric tails (just like RooDSCBShape)

RooArgSet* extractMCSignalTails_asymCoreDSCB(Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
	/// open the MC skimmed file
	const char* filename = "../Files/MCUpsilonSkimmedWeightedDataset.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

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

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	RooDataSet* massDataset = ReducedMassDataset(allDataset, *wspace, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar* massVar = wspace->var("mass");

	using namespace RooFit;
	// fit
	RooRealVar mean("mean", "", 9., 10.);
	RooRealVar sigmaInf("sigmaInf", "", .05, .15);
	RooRealVar alphaInf("alphaInf", "", 0.1, 10);
	RooRealVar orderInf("orderInf", "", 0.1, 10);
	RooRealVar sigmaSup("sigmaSup", "", .05, .15);
	RooRealVar alphaSup("alphaSup", "", 0.1, 10);
	RooRealVar orderSup("orderSup", "", 0.1, 40);

	RooCrystalBall signal("CB", "", *massVar, mean, sigmaInf, sigmaSup, alphaInf, orderInf, alphaSup, orderSup);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with a double-sided Crystal Ball PDF made of an asymmetric Gaussian core and asymmetric tail distributions..." << endl;

	bool doWeightedError = true;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true)/*, PrintLevel(-1)*/, Minos(!doWeightedError), NumCPU(3), Range(massMin, massMax), AsymptoticError(doWeightedError)); 
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit

	fitResult -> Print("v");

	TString outputName = Form("asymCoreDSCB_cent%dto%d_pt%dto%d", centMin, centMax, ptMin, ptMax);

	// save signal shape parameters in a txt file to be read for data fit
	RooArgSet* tailParams = new RooArgSet(alphaInf, orderInf, alphaSup, orderSup);

	SaveMCSignalParameters(tailParams, outputName.Data()); // so that we don't have to refit later

 	/// draw the fit to see if the fit is reasonable
	auto* canvas = new TCanvas("canvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = massVar -> frame(Title(" "), Range(massMin, massMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	// frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset -> plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	signal.plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(centMin, centMax, ptMin, ptMax));
	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	frame->addObject(AsymDoubleCBParamsText(mean, sigmaInf, alphaInf, orderInf, sigmaSup, alphaSup, orderSup));
	frame->Draw();

	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->cd();
	pad1->Draw();
	pad2->Draw();
	canvas->SaveAs(Form("SignalParameters/MCfit_%s_%dto%d.png", "asymCoreDSCB", ptMin, ptMax), "RECREATE");

	// file->Close();
	return tailParams;
}
