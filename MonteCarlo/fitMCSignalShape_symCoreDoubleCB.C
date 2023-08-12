#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Shortcuts.h"

#include "../AnalysisParameters.h"

// crystal ball shape with symmetric Gaussian core and asymmetric tails (just like RooDSCBShape)

void fitMCSignalShape_symCoreDoubleCB(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* filename = "../Files/MCUpsilonSkimmedWeightedDataset.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

	Float_t massMin = 8.5, massMax = 10.5;
	Int_t nBins = 75;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)file->Get("MCdataset");

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar* massVar = wspace->var("mass");

	Long64_t nEntries = massDataset->sumEntries();

	auto* canvas = new TCanvas("canvas", "", 600, 650);
	//canvas->Divide(2);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.01);
	pad1->SetLogy();
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = massVar->frame(Title(" "), Range(massMin, massMax));
	//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset->plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	// fit
	RooRealVar mean("mean", "", 9., 10.);
	RooRealVar sigma("sigma", "", .05, .15);
	RooRealVar alphaInf("alphaInf", "", 0.1, 10);
	RooRealVar orderInf("orderInf", "", 0.1, 10);
	RooRealVar alphaSup("alphaSup", "", 0.1, 10);
	RooRealVar orderSup("orderSup", "", 0.1, 40);

	RooCrystalBall signal("CB", "", *massVar, mean, sigma, alphaInf, orderInf, alphaSup, orderSup);

	cout << endl
	     << "Fitting the MC signal shape (weighted entries!!) with a double-sided Crystal Ball PDF with a symmetric Gaussian core and asymmetric tail distributions..." << endl;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true), PrintLevel(-1), Minos(!DoMCWeightedError), NumCPU(NCPUs), Range(massMin, massMax), AsymptoticError(DoMCWeightedError)); // quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit

	fitResult->Print("v");

	signal.plotOn(frame, LineColor(kBlue));

	// add legends
	frame->addObject(KinematicsText(centMin, centMax));

	frame->addObject(RefFrameText(ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->addObject(SymCoreDoubleCBParamsText(mean, sigma, alphaInf, orderInf, alphaSup, orderSup));

	frame->Draw();

	frame->SetMaximum(3 * nEntries);
	frame->SetMinimum(20);

	CMS_lumi(pad1, "#varUpsilon(1S) Hydjet-embedded PbPb MC");

	// pull distribution
	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	canvas->Modified();
	canvas->Update();
	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	TString outputName = Form("symCoreDSCB_cent%dto%d_pt%dto%d_cosTheta%.1fto%.1f_phi%dto%d_%s", centMin, centMax, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax, (isCSframe) ? "CS" : "HX");

	canvas->SaveAs(Form("../MonteCarlo/SignalShapeFits/%s.png", outputName.Data()), "RECREATE");

	// save signal shape parameters in a txt file to be read for data fit
	RooArgSet tailParams(alphaInf, orderInf, alphaSup, orderSup);

	SaveMCSignalTailParameters(tailParams, outputName.Data());

	file->Close();
}
