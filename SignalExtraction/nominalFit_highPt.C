#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/RooFitPDFs/ErrorFuncTimesExp.h"
#include "../Tools/Style/FitDistributions.h"

void nominalFit_highPt(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	TString fitModelName = Form("symCoreDSCB_expo_cent%dto%d_pt%dto%d_cosTheta%.1fto%.1f_phi%dto%d_%s", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax, (isCSframe) ? "CS" : "HX");

	// get the tail parameters of the signal shape first in case the MC fit is needed
	RooRealVar* alphaInf = new RooRealVar("alphaInf", "", 1);
	RooRealVar* orderInf = new RooRealVar("orderInf", "", 1);
	RooRealVar* alphaSup = new RooRealVar("alphaSup", "", 1);
	RooRealVar* orderSup = new RooRealVar("orderSup", "", 1);

	RooArgSet tailParams = GetMCSignalTailParameters(alphaInf, orderInf, alphaSup, orderSup, "symCoreDSCB", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);

	const char* filename = "../Files/upsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* datasetName = Form("dataset%s", isCSframe? "CS":"HX");
	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

	RooRealVar* massVar = wspace->var("mass");

	Long64_t nEntries = massDataset->sumEntries();

	auto* canvas = new TCanvas("canvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = massVar->frame(Title(" "), Range(MassBinMin, MassBinMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad

	//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (binMax - binMin) / nBins));
	massDataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"));

	/// fitting model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.3, 9.6);
	RooRealVar sigma_1S("sigma_1S", "", .03, .15);

	RooCrystalBall signal_1S("signal_1S", "", *massVar, mean_1S, sigma_1S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar nSignal_1S("nSignal_1S", "N 1S", 0, nEntries);

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "massScaling_2S*mean_1S", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigma_2S("sigma_2S", "massScaling_2S*sigma_1S", RooArgSet(massScaling_2S, sigma_1S));

	RooCrystalBall signal_2S("signal_2S", "", *massVar, mean_2S, sigma_2S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar nSignal_2S("nSignal_2S", "N 2S", 0, nEntries / 2);

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "massScaling_3S*mean_1S", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigma_3S("sigma_3S", "massScaling_3S*sigma_1S", RooArgSet(massScaling_3S, sigma_1S));

	RooCrystalBall signal_3S("signal_3S", "", *massVar, mean_3S, sigma_3S, *alphaInf, *orderInf, *alphaSup, *orderSup);
	RooRealVar nSignal_3S("nSignal_3S", "N 3S", 0, nEntries / 4);

	// background: exponential
	RooRealVar exp_lambda("exp_lambda", "m_lambda", -10, 10);

	RooExponential bkgPDF("bkgPDF", "", *massVar, exp_lambda);
	RooRealVar nBkg("nBkg", "N background events", 0, nEntries);

	RooAddPdf fitModel("fitModel", "", RooArgList(signal_1S, signal_2S, signal_3S, bkgPDF), RooArgList(nSignal_1S, nSignal_2S, nSignal_3S, nBkg));

	auto* fitResult = fitModel.fitTo(*massDataset, Save(), Extended(kTRUE), PrintLevel(-1), Minos(kTRUE), NumCPU(NCPUs), Range(MassBinMin, MassBinMax));

	fitResult->Print("v");

	// compute significance
	massVar->setRange("integral", mean_1S.getVal() - 3 * sigma_1S.getVal(), mean_1S.getVal() + 3 * sigma_1S.getVal());

	Double_t signal = signal_1S.createIntegral(*massVar, NormSet(*massVar), Range("integral"))->getVal() * nSignal_1S.getVal();

	Double_t bkg = bkgPDF.createIntegral(*massVar, NormSet(*massVar), Range("integral"))->getVal() * nBkg.getVal();

	Double_t significance = signal / sqrt(signal + bkg);

	cout << endl
	     << "Y(1S) yield significance = " << significance << endl;

	// for Y(2S)
	massVar->setRange("integral2S", mean_2S.getVal() - 3 * sigma_2S.getVal(), mean_2S.getVal() + 3 * sigma_2S.getVal());

	Double_t signal2S = signal_2S.createIntegral(*massVar, NormSet(*massVar), Range("integral2S"))->getVal() * nSignal_2S.getVal();

	Double_t bkg2S = bkgPDF.createIntegral(*massVar, NormSet(*massVar), Range("integral2S"))->getVal() * nBkg.getVal();

	Double_t significance2S = signal2S / sqrt(signal2S + bkg2S);

	cout << endl
	     << "Y(2S) yield significance = " << significance2S << endl;

	fitModel.plotOn(frame, Components(bkgPDF), LineColor(kGray + 2), LineStyle(kDashed));
	fitModel.plotOn(frame, Components(signal_1S), LineColor(kRed));
	fitModel.plotOn(frame, Components(signal_2S), LineColor(kRed));
	fitModel.plotOn(frame, Components(signal_3S), LineColor(kRed));
	fitModel.plotOn(frame, LineColor(kBlue));

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax));

	frame->addObject(RefFrameText(ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->addObject(FitResultText(nSignal_1S, significance, nSignal_2S, significance2S));

	frame->Draw();
	gPad->RedrawAxis();

	frame->SetMaximum(nEntries / 15);
	//frame->SetMinimum(0.8);

	CMS_lumi(pad1, lumi_PbPb);

	// pull distribution
	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	//canvas->Modified();
	//canvas->Update();
	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	canvas->SaveAs(Form("FitPlots/%s.png", fitModelName.Data()), "RECREATE");
}

void scanNominalFit_highPt(){

	Int_t PtEdges[9] = {0, 2, 4, 6, 8, 12, 16, 20, 30};
	Int_t NumPtEle = sizeof(PtEdges)/sizeof(Int_t);
	for(Int_t idx =0; idx < NumPtEle-1; idx++){
		nominalFit_highPt(0, 90, PtEdges[idx], PtEdges[idx+1]);
	}
}


