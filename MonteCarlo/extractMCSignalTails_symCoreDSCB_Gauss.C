//#include "../Tools/Style/tdrStyle.C"
//#include "../Tools/Style/CMS_lumi.C"
#include "../Tools/FitShortcuts.h"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

// #include "../Tools/Parameters/PhysicsConstants.h"

// crystal ball shape with symmetric Gaussian core and asymmetric tails (just like RooDSCBShape)

RooArgSet* extractMCSignalTails_symCoreDSCB_Gauss(Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {
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

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar* massVar = wspace->var("mass");

	/// fit
	/// (DSCB variables)
	RooRealVar mean("mean", "", 9.457, 9., 10.);
	RooRealVar sigma("sigma", "", 0.08, .05, .15);
	RooRealVar alphaInf("alphaInf", "", 1.5, 0.1, 10);
	RooRealVar orderInf("orderInf", "", 1.5, 0.1, 10);
	RooRealVar alphaSup("alphaSup", "", 1.5, 0.1, 10);
	RooRealVar orderSup("orderSup", "", 3, 0.1, 40);

	/// (Gaussian variable (used the same mean as DSCB))
	RooRealVar sigma_gauss("sigma_gauss", "", 0.11, .05, .3);

	RooRealVar normFraction("normFraction", "", 0.6, 0.01, 1);

	RooCrystalBall DSCB("DSCB", "", *massVar, mean, sigma, alphaInf, orderInf, alphaSup, orderSup);
	RooGaussian gauss("gauss", "gauss", *massVar, mean, sigma_gauss);

	// RooCBShape CB_1("CB_1", "", *massVar, mean, sigma_1, alphaPar, orderPar);
	RooAddPdf signal("signal", "sum of DSCB and Gaussian PDF", RooArgList(DSCB, gauss), RooArgList(normFraction), kTRUE);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with a double-sided Crystal Ball PDF made of a symmetric Gaussian core and asymmetric tail distributions..." << endl;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(false), PrintLevel(-1), Minos(!DoMCWeightedError), NumCPU(DoMCWeightedError), Range(massMin, massMax), AsymptoticError(DoMCWeightedError));
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit

	fitResult->Print("v");

	//TString outputName = Form("symCoreDSCB_Gaus_cent%dto%d_pt%dto%d", centMin, centMax, ptMin, ptMax);
	const char* outputName = GetSignalFitName("DSCB_Gauss", ptMin, ptMax);

	// save signal shape parameters in a txt file to be read for data fit
	RooArgSet* Params = new RooArgSet(sigma, alphaInf, orderInf, alphaSup, orderSup, sigma_gauss, normFraction);

	SaveMCSignalParameters(Params, outputName); // so that we don't have to refit later

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

	signal.plotOn(frame, Components("gauss"), LineColor(kRed));

	signal.plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(centMin, centMax, ptMin, ptMax));
	//frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	frame->addObject(SymCoreDoubleCBGaussParamsText(mean, sigma, alphaInf, orderInf, alphaSup, orderSup, sigma_gauss, normFraction));
	frame->Draw();

	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	const char* outputName2 = GetSignalFitName("DSCB_Gauss", ptMin, ptMax);
	canvas->SaveAs(Form("SignalShapeFits/%s.png", outputName2), "RECREATE");

	// file->Close();
	return Params;
}

void scanExtractMCSignalTails_symCoreDSCB_Gauss(){

	Int_t PtEdges[9] = {0, 2, 4, 6, 8, 12, 16, 20, 30};
	Int_t NumPtEle = sizeof(PtEdges)/sizeof(Int_t);
	for(Int_t idx =0; idx < NumPtEle-1; idx++){
		extractMCSignalTails_symCoreDSCB_Gauss(0, 90, PtEdges[idx], PtEdges[idx+1]);
	}
}