//#include "../Tools/Style/tdrStyle.C"
//#include "../Tools/Style/CMS_lumi.C"
#include "../Tools/FitShortcuts.h"

#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Parameters/PhysicsConstants.h"

#include "../Tools/CustomRoofitPDFs/ErrorFuncTimesExp.h"

void SymDSCB_ErrxExp(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* signalShapeName = "SymDSCB";

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	// get the tail parameters of the signal shape first in case the MC fit is needed
	RooRealVar* alphaInf = new RooRealVar("alphaInf", "", 1);
	RooRealVar* orderInf = new RooRealVar("orderInf", "", 1);
	RooRealVar* alphaSup = new RooRealVar("alphaSup", "", 1);
	RooRealVar* orderSup = new RooRealVar("orderSup", "", 1);

	RooArgSet tailParams = GetMCSignalTailParameters(alphaInf, orderInf, alphaSup, orderSup, signalShapeName, ptMin, ptMax);

	const char* filename = "../Files/UpsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)f->Get("dataset");

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar* massVar = wspace->var("mass");

	Long64_t nEntries = massDataset->sumEntries();

	auto* canvas = new TCanvas("canvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = massVar->frame(Title(" "), Range(6, MassBinMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"));

	/// fitting model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.3, 9.6);
	RooRealVar sigma_1S("sigma_1S", "", .13, .01, .2);

	RooCrystalBall signal_1S("signalPDF_1S", "", *massVar, mean_1S, sigma_1S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar nSignal_1S("nSignal_1S", "N 1S", nEntries / 5, 0, nEntries);

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "massScaling_2S*mean_1S", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigma_2S("sigma_2S", "massScaling_2S*sigma_1S", RooArgSet(massScaling_2S, sigma_1S));

	RooCrystalBall signal_2S("signalPDF_2S", "", *massVar, mean_2S, sigma_2S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar nSignal_2S("nSignal_2S", "N 2S", nEntries / 10, 0, nEntries / 2);

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "massScaling_3S*mean_1S", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigma_3S("sigma_3S", "massScaling_3S*sigma_1S", RooArgSet(massScaling_3S, sigma_1S));

	RooCrystalBall signal_3S("signalPDF_3S", "", *massVar, mean_3S, sigma_3S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar nSignal_3S("nSignal_3S", "N 3S", nEntries / 20, 0, nEntries / 2);

	// background: error function x exponential
	RooRealVar err_mu("err_mu", "err_mu", 8, 6, 12);
	RooRealVar err_sigma("err_sigma", "err_sigma", 1, 0, 10);
	RooRealVar exp_lambda("exp_lambda", "m_lambda", 5, 0, 40);

	ErrorFuncTimesExp bkgPDF("bkgPDF", "", *massVar, err_mu, err_sigma, exp_lambda);
	RooRealVar nBkg("nBkg", "N background events", 0, nEntries);

	RooAddPdf fitModel("fitModel", "", RooArgList(signal_1S, signal_2S, signal_3S, bkgPDF), RooArgList(nSignal_1S, nSignal_2S, nSignal_3S, nBkg));

	auto* fitResult = fitModel.fitTo(*massDataset, Save(), Extended(kTRUE), PrintLevel(-1), Minos(kTRUE), NumCPU(NCPUs), Range(MassBinMin, MassBinMax));

	fitResult->Print("v");

	wspace->import(fitModel);

	// compute significance
	Double_t significance1S = ComputeSignalSignificance(wspace, mean_1S.getVal(), sigma_1S.getVal(), 1);

	Double_t significance2S = ComputeSignalSignificance(wspace, mean_2S.getVal(), sigma_2S.getVal(), 2);

	fitModel.plotOn(frame, Components(bkgPDF), LineColor(kGray + 2), LineStyle(kDashed));
	fitModel.plotOn(frame, Components(signal_1S), LineColor(kRed));
	fitModel.plotOn(frame, Components(signal_2S), LineColor(kRed));
	fitModel.plotOn(frame, Components(signal_3S), LineColor(kRed));
	fitModel.plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	//frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->addObject(FitResultText(nSignal_1S, significance1S, nSignal_2S, significance2S));

	frame->Draw();
	frame->GetYaxis()->SetMaxDigits(3);

	gPad->RedrawAxis();

	//frame->SetMaximum(nEntries / 15);
	//frame->SetMinimum(0.8);

	//	CMS_lumi(pad1, "2018 PbPb miniAOD, DoubleMuon PD");

	// pull distribution
	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	//canvas->Modified();
	//canvas->Update();
	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	canvas->SaveAs(Form("FitPlots/%s.png", fitModelName), "RECREATE");
}
