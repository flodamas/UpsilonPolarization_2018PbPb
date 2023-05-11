#include "../Tools/Style/tdrStyle.C"
#include "../Tools/Style/CMS_lumi.C"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Shortcuts.h"

// sum of two crystal ball PDF with common pole and tails, but "scaling" Gaussian widths

void fitMCSignalShape_twoCB(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* filename = "../Files/MCUpsilonSkimmedWeightedDataset.root";

	TFile* file = TFile::Open(filename, "READ");
	if (!file) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "      Simulation Internal";

	Int_t centMin = 0, centMax = 90;
	Float_t massMin = 9, massMax = 10.;
	Int_t nBins = 100;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)file->Get("MCdataset");

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, centMin, centMax, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

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

	// first CB
	RooRealVar mean("mean", "", 9., 10.);
	RooRealVar sigma_1("sigma_1", "", .05, .3);

	RooRealVar alphaPar("alphaPar", "", 0.1, 10);
	RooRealVar orderPar("orderPar", "", 0.1, 10);

	RooCBShape CB_1("CB_1", "", *massVar, mean, sigma_1, alphaPar, orderPar);

	// second CB
	RooRealVar resFraction("resFraction", "", 0.01, 1);
	RooFormulaVar sigma_2("sigma_2", "resFraction*sigma_1", RooArgSet(resFraction, sigma_1));
	//RooRealVar sigma_2("sigma_2", "", .05, .3);

	RooCBShape CB_2("CB_2", "", *massVar, mean, sigma_2, alphaPar, orderPar);

	// add the two
	RooRealVar normFraction("normFraction", "", 0.01, 1);
	RooAddPdf signal("signal", "sum of two CBs PDF", RooArgList(CB_1, CB_2), RooArgList(normFraction), kTRUE);

	// fit
	bool doWeightedError = true;

	auto* fitResult = signal.fitTo(*massDataset, Save(), Extended(true), PrintLevel(-1), Minos(!doWeightedError), NumCPU(4), Range(massMin, massMax), AsymptoticError(doWeightedError)); // quoting RooFit: "sum-of-weights and asymptotic error correction do not work with MINOS errors", so let's turn off Minos, no need to estimate asymmetric errors with MC fit

	fitResult->Print("v");

	signal.plotOn(frame, Components("CB_1"), LineColor(kRed));
	signal.plotOn(frame, Components("CB_2"), LineColor(kOrange + 1));

	signal.plotOn(frame, LineColor(kBlue));

	// add legends
	frame->addObject(KinematicsText(centMin, centMax, ptMin, ptMax));

	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->addObject(TwoCBParamsText(mean, sigma_1, alphaPar, orderPar, resFraction, normFraction));

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

	TString outputName = Form("twoCB_cent%dto%d_pt%dto%d_cosTheta%.1fto%.1f_phi%dto%d_%s", centMin, centMax, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax, (isCSframe) ? "CS" : "HX");

	canvas->SaveAs(Form("../MonteCarlo/SignalShapeFits/%s.png", outputName.Data()), "RECREATE");

	// save signal shape parameters in a txt file to be read for data fit
	SaveMCSignalTailParameters(RooArgSet(alphaPar, orderPar, resFraction, normFraction), outputName.Data());

	file->Close();
}
