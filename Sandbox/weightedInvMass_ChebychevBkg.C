//#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"

void weightedInvMass_ChebychevBkg(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* signalShapeName = "SymDSCB";

	// get the tail parameters of the signal shape first in case the MC fit is needed
	RooRealVar* alphaInf = new RooRealVar("alphaInf", "", 1);
	RooRealVar* orderInf = new RooRealVar("orderInf", "", 1);
	RooRealVar* alphaSup = new RooRealVar("alphaSup", "", 1);
	RooRealVar* orderSup = new RooRealVar("orderSup", "", 1);

	RooArgSet tailParams = GetMCSignalTailParameters(alphaInf, orderInf, alphaSup, orderSup, signalShapeName, ptMin, ptMax);

	const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)f->Get(isCSframe ? "datasetCS" : "datasetHX");

	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooDataSet* massDataset = (isCSframe) ? ReducedWeightedMassDatasetCS(allDataset, wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax) : ReducedWeightedMassDatasetHX(allDataset, wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar massVar = *wspace.var("mass");

	Long64_t nEntries = massDataset->sumEntries();

	/// fitting model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	auto signalModel = NominalSignalModel(wspace, alphaInf, orderInf, alphaSup, orderSup, nEntries);

	// get variables
	RooAbsPdf* signalPDF_1S = wspace.pdf("signalPDF_1S");
	RooAbsPdf* signalPDF_2S = wspace.pdf("signalPDF_2S");
	RooAbsPdf* signalPDF_3S = wspace.pdf("signalPDF_3S");

	RooRealVar yield1S = *(wspace.var("yield1S"));
	RooRealVar yield2S = *(wspace.var("yield2S"));
	RooRealVar yield3S = *(wspace.var("yield3S"));

	// background: Chebychev polynomial

	int order = 3;

	RooArgList coefList = ChebychevCoefList(order);

	RooChebychev bkgPDF("bkgPDF", " ", massVar, coefList);

	RooRealVar yieldBkg("yieldBkg", "N background events", 0, nEntries);

	RooAddPdf fitModel("fitModel", "", {*signalPDF_1S, *signalPDF_2S, *signalPDF_3S, bkgPDF}, RooArgList(yield1S, yield2S, yield3S, yieldBkg));

	auto* fitResult = fitModel.fitTo(*massDataset, Save(), Extended(kTRUE), PrintLevel(-1), NumCPU(NCPUs), Range(MassBinMin, MassBinMax), AsymptoticError(DoAsymptoticError), SumW2Error(!DoAsymptoticError));

	fitResult->Print("v");

	wspace.import(fitModel);

	// draw
	auto* canvas = new TCanvas("canvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = massVar.frame(Title(" "), Range(MassBinMin, MassBinMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	massDataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"));

	fitModel.plotOn(frame, Components(bkgPDF), LineColor(kGray + 2), LineStyle(kDashed));
	fitModel.plotOn(frame, Components(*signalPDF_1S), LineColor(kRed));
	fitModel.plotOn(frame, Components(*signalPDF_2S), LineColor(kRed));
	fitModel.plotOn(frame, Components(*signalPDF_3S), LineColor(kRed));
	fitModel.plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	//frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->addObject(FitResultText(yield1S, ComputeSignalSignificance(wspace, 1), yield2S, ComputeSignalSignificance(wspace, 2)));

	frame->Draw();
	frame->GetYaxis()->SetMaxDigits(3);

	gPad->RedrawAxis();

	// pull distribution
	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	//canvas->Modified();
	//canvas->Update();
	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	canvas->SaveAs(Form("mass_distrib/weightedEvents_ChebychevBkg_%s.png", fitModelName), "RECREATE");
}
