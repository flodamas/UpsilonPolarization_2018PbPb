
#include "../Tools/Style/FitDistributions.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Shortcuts.h"

#include "../Tools/Parameters/PhysicsConstants.h"

#include "../AnalysisParameters.h"

#include "../Tools/CustomRoofitPDFs/ErrorFuncTimesExp.h"

void nominalFit_lowPt_2Danalysis(Bool_t isCSframe = kTRUE) {
	// get the skimmed dataset

	const char* filename = "../Files/upsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	string color = "\033[1;31m";

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)f->Get("dataset");

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	/// definition of the fitting model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
	// tail parameters fixed to MC extracted values, and identical for the three resonances
	RooRealVar* alphaInf = new RooRealVar("alphaInf", "", 1);
	RooRealVar* orderInf = new RooRealVar("orderInf", "", 1);
	RooRealVar* alphaSup = new RooRealVar("alphaSup", "", 1);
	RooRealVar* orderSup = new RooRealVar("orderSup", "", 1);

	RooArgSet tailParams = GetMCSignalTailParameters(alphaInf, orderInf, alphaSup, orderSup, "symCoreDSCB", gCentralityBinMin, gCentralityBinMax, gPtMin, gPtMax);

	// Y(1S) signal shape
	RooRealVar mean_1S("mean_1S", "mean 1S", PDGmass_1S, 9.3, 9.6);
	RooRealVar sigma_1S("sigma_1S", "", .03, .15);

	// Y(2S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_2S("massScaling_2S", "", PDGmass_2S / PDGmass_1S);

	RooFormulaVar mean_2S("mean_2S", "massScaling_2S*mean_1S", RooArgSet(massScaling_2S, mean_1S));
	RooFormulaVar sigma_2S("sigma_2S", "massScaling_2S*sigma_1S", RooArgSet(massScaling_2S, sigma_1S));

	// Y(3S) signal shape, mass scaling for mean and widths
	RooConstVar massScaling_3S("massScaling_3S", "", PDGmass_3S / PDGmass_1S);

	RooFormulaVar mean_3S("mean_3S", "massScaling_3S*mean_1S", RooArgSet(massScaling_3S, mean_1S));
	RooFormulaVar sigma_3S("sigma_3S", "massScaling_3S*sigma_1S", RooArgSet(massScaling_3S, sigma_1S));

	// background: error function x exponential
	RooRealVar err_mu("err_mu", "err_mu", 2, 12);
	RooRealVar err_sigma("err_sigma", "err_sigma", 0, 5);
	RooRealVar exp_lambda("exp_lambda", "m_lambda", 0, 10);

	// store the upsilon's raw yields in 2D histograms for the polarization extraction procedure
	TH2D* hYield1S = new TH2D(Form("RawYield1S_%s_pt%dto%d", (isCSframe) ? "CS" : "HX", gPtMin, gPtMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);corrected yield", NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);

	TH2D* hYield2S = new TH2D(Form("RawYield2S_%s_pt%dto%d", (isCSframe) ? "CS" : "HX", gPtMin, gPtMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);corrected yield", NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);

	// store the upsilon's yield significance too, quite useful
	TH2D* hSignificance1S = new TH2D(Form("Significance1S_%s_pt%dto%d", (isCSframe) ? "CS" : "HX", gPtMin, gPtMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);corrected yield", NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);

	TH2D* hSignificance2S = new TH2D(Form("Significance2S_%s_pt%dto%d", (isCSframe) ? "CS" : "HX", gPtMin, gPtMax), ";cos #theta_{CS}; #varphi_{CS} (#circ);corrected yield", NCosThetaBinsCS, CosThetaBinningCS, NPhiBinsCS, PhiBinningCS);

	// perform the signal extraction for every (cos theta, phi) bins defined in AnalysisParameters.h
	for (int iCosTheta = 0; iCosTheta < NCosThetaBinsCS; iCosTheta++) {
		for (int iPhi = 0; iPhi < NPhiBinsCS; iPhi++) {
			// get the boundaries corresponding to those bins
			double cosThetaMin = CosThetaBinningCS[iCosTheta], cosThetaMax = CosThetaBinningCS[iCosTheta + 1];
			double phiMin = PhiBinningCS[iPhi], phiMax = PhiBinningCS[iPhi + 1];

			TString fitModelName = Form("symCoreDSCB_cent%dto%d_pt%dto%d_cosTheta%.2fto%.2f_phi%.0fto%.0f_%s", gCentralityBinMin, gCentralityBinMax, gPtMin, gPtMax, cosThetaMin, cosThetaMax, phiMin, phiMax, (isCSframe) ? "CS" : "HX");

			cout << endl
			     << color
			     << cosThetaMin << " < cos theta < " << cosThetaMax << "; " << phiMin << " < phi < " << phiMax << "\033[0m" << endl;

			// reduce the dataset to the corresponding kinematic region
			RooDataSet* massDataset = ReducedMassDataset(allDataset, wspace, gCentralityBinMin, gCentralityBinMax, gPtMin, gPtMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);

			RooRealVar* massVar = wspace->var("mass");

			Long64_t nEntries = massDataset->sumEntries();

			auto* canvas = new TCanvas("canvas", "", 600, 600);
			TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
			pad1->SetBottomMargin(0.03);
			pad1->Draw();
			pad1->cd();

			RooPlot* frame = massVar->frame(Title(" "), Range(MassBinMin, MassBinMax));
			frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
			//frame->SetYTitle(Form("Candidates / (%d MeV)", (int)1000 * (MassBinMax - MassBinMin) / NMassBins));
			massDataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"));

			/// fitting model

			// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance
			// tail parameters fixed to MC extracted values, and identical for the three resonances

			// Y(1S) signal shape
			RooCrystalBall signal_1S("signal_1S", "", *massVar, mean_1S, sigma_1S, *alphaInf, *orderInf, *alphaSup, *orderSup);
			RooRealVar nSignal_1S("nSignal_1S", "N 1S", nEntries / 5, 0, nEntries);

			// Y(2S) signal shape
			RooCrystalBall signal_2S("signal_2S", "", *massVar, mean_2S, sigma_2S, *alphaInf, *orderInf, *alphaSup, *orderSup);
			RooRealVar nSignal_2S("nSignal_2S", "N 2S", nEntries / 10, 0, nEntries / 2);

			// Y(3S) signal shape
			RooCrystalBall signal_3S("signal_3S", "", *massVar, mean_3S, sigma_3S, *alphaInf, *orderInf, *alphaSup, *orderSup);
			RooRealVar nSignal_3S("nSignal_3S", "N 3S", nEntries / 20, 0, nEntries / 2);

			// background: error function x exponential
			ErrorFuncTimesExp bkgPDF("bkgPDF", "", *massVar, err_mu, err_sigma, exp_lambda);
			RooRealVar nBkg("nBkg", "N background events", 0, nEntries);

			// add all the fit components
			RooAddPdf fitModel("fitModel", "", RooArgList(signal_1S, signal_2S, signal_3S, bkgPDF), RooArgList(nSignal_1S, nSignal_2S, nSignal_3S, nBkg));

			// fit!!!
			auto* fitResult = fitModel.fitTo(*massDataset, Save(), Extended(kTRUE), PrintLevel(-1), Minos(kTRUE), NumCPU(NCPUs), Range(MassBinMin, MassBinMax));

			fitResult->Print("v");

			// compute significance (TO DO: move this part to a dedicated macro, would be much more convenient)
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

			// store the yield results
			// WATCH OUT: histogram binning starts from 1
			int globalBin = hYield1S->GetBin(iCosTheta + 1, iPhi + 1);
			hYield1S->SetBinContent(globalBin, nSignal_1S.getVal());
			hYield1S->SetBinError(globalBin, nSignal_1S.getError()); // symmetric errors for now, will need to switch to TGraph2DAsymmErrors to store asymmetric ones...

			hYield2S->SetBinContent(globalBin, nSignal_2S.getVal());
			hYield2S->SetBinError(globalBin, nSignal_2S.getError()); // symmetric errors for now, will need to switch to TGraph2DAsymmErrors to store asymmetric ones...

			hSignificance1S->SetBinContent(globalBin, significance);
			hSignificance2S->SetBinContent(globalBin, significance2S);

			// draw stuff
			fitModel.plotOn(frame, Components(bkgPDF), LineColor(kGray + 2), LineStyle(kDashed));
			fitModel.plotOn(frame, Components(signal_1S), LineColor(kRed));
			fitModel.plotOn(frame, Components(signal_2S), LineColor(kRed));
			fitModel.plotOn(frame, Components(signal_3S), LineColor(kRed));
			fitModel.plotOn(frame, LineColor(kBlue));

			frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax));

			frame->addObject(RefFrameText(gPtMin, gPtMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

			frame->addObject(FitResultText(nSignal_1S, significance, nSignal_2S, significance2S));

			frame->Draw();
			gPad->RedrawAxis();

			frame->SetMaximum(nEntries / 15);

			CMS_lumi(pad1, "2018 PbPb miniAOD, DoubleMuon PD");

			// pull distribution
			canvas->cd();

			TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

			canvas->cd();
			pad1->Draw();
			pad2->Draw();

			canvas->SaveAs(Form("FitPlots/%s.png", fitModelName.Data()), "RECREATE");
			canvas->Close();
		} // end of the loop over the phi dimension
	}   // end of the loop over the cos theta dimension

	// save the fit results in a TFile
	const char* outputFileName = "YieldResults/NominalFit_2DAnalysis.root";

	TFile outputFile(outputFileName, "RECREATE");

	hYield1S->Write();
	hYield2S->Write();
	hSignificance1S->Write();
	hSignificance2S->Write();

	cout << endl
	     << "Distributions of the signal extraction results saved in " << outputFileName << endl;
}
