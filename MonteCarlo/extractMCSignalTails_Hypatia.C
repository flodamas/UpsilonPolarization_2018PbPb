#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/FitShortcuts.h"

// #include "../Tools/Parameters/PhysicsConstants.h"

// crystal ball shape with symmetric Gaussian core and asymmetric tails (just like RooDSCBShape)

RooArgSet* extractMCSignalTails_Hypatia(Int_t centMin = 0, Int_t centMax = 90, Int_t ptMin = 0, Int_t ptMax = 30) {

	/// Start measuring time
	clock_t start, end, cpu_time;
	start = clock();

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

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, centMin, centMax, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar* massVar = wspace->var("mass");

	/// fit
	/// Hypatia variables
	RooRealVar mean("meanHypatia", "", 9.47, 9., 10.);
	RooRealVar lambda("lambdaHypatia","lambda of hypatia PDF",-1.0, -20.0, -0.1);
	RooRealVar zeta("zetaHypatia","zeta of hypatia PDF",0.01, 0.0, 1.0);
	RooRealVar beta("betaHypatia","beta of hypatia PDF",-0.01, -20.0, 0.0);
	RooRealVar sigma("sigmaHypatia","sigma of hypatia PDF",0.15, 0.1, 0.3);
	RooRealVar alphaInf("alphaInfHypatia","al1s of hypatia PDF", 3.0 , 0.1, 5.0);
	RooRealVar alphaSup("alphaSupHypatia","ar1s of hypatia PDF", 3.0 , 0.1, 7.0);
	RooRealVar orderInf("orderInfHypatia","nl1s of hypatia PDF", 1.0 , 0.2, 14.718);
	RooRealVar orderSup("orderSupHypatia","nr1s of hypatia PDF", 1.0 , 0.0, 14.718);

	lambda.setVal(-1.754); 
	zeta.setVal(0); 
	beta.setVal(-1.545); 
	sigma.setVal(0.160); 
	alphaInf.setVal(1.66); 
	alphaSup.setVal(5.360); 
	orderInf.setVal(1.611); 
	orderSup.setVal(0.003); 

	// zeta.setConstant(1);

	RooHypatia2 signal("Hypatia","Hypatia", *massVar, lambda, zeta, beta, sigma, mean, alphaInf, orderInf, alphaSup, orderSup);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with a Hypatia PDF" << endl;

	bool doWeightedError = true;

	/// fit
	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true)/*, PrintLevel(-1)*/, Minos(!doWeightedError), NumCPU(8), Range(massMin, massMax), AsymptoticError(doWeightedError)); 
	// quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit

	fitResult->Print("v");

	TString outputName = Form("Hypatia_cent%dto%d_pt%dto%d", centMin, centMax, ptMin, ptMax);

	// save signal shape parameters in a txt file to be read for data fit
	RooArgSet* Params = new RooArgSet(mean, lambda, zeta, beta, sigma, alphaInf, orderInf, alphaSup, orderSup);

	SaveMCSignalParameters(Params, outputName.Data()); // so that we don't have to refit later

	/// draw the fit to see if the fit is reasonable (we can comment it (lines 79-105) out if drawing is not necessary)
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
	frame->addObject(HypatiaParamsText(mean, lambda, zeta, beta, sigma, alphaInf, orderInf, alphaSup, orderSup));
	frame->Draw();

	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->cd();
	pad1->Draw();
	pad2->Draw();
	canvas->SaveAs(Form("SignalParameters/plots/MCfit_%s_%dto%d.png", "Hypatia", ptMin, ptMax), "RECREATE");

	// file->Close();

	/// End measuring time
	end = clock();
	cpu_time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time / 60. << "minutes" << endl;

	return Params;
}


void extractMCSignalPara_symCoreDSCB_Gauss_scan(){

	Int_t ptCuts[9] = {0, 2, 4, 6, 8, 12, 16, 20, 30};
	Int_t NumPtArr = sizeof(ptCuts)/sizeof(Int_t);
	for(int idx=0; idx< NumPtArr-1; idx++){
		
	extractMCSignalTails_Hypatia(0, 90, ptCuts[idx], ptCuts[idx+1]);
	}

}