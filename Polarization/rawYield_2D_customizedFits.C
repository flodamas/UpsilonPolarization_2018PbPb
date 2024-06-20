#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../sPlot/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "PolarFitHelpers.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"

void rawYield_2D_customizedFits(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 6, Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = gUpsilonState) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Bin edges and width
	// Set the bin edges along the cosTheta/phi axis depending on the number of bins, min and max values
	// (If want to use non-uniform bin width, the bin edges should be pre-defined in PolarFitHelpers.h)
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	/// Set up the variables
	RooRealVar cosTheta("cosTheta", "", cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	cosTheta.setRange(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	cosTheta.setBins(nCosThetaBins);

	RooRealVar* yield1S = new RooRealVar("yield1S", "", 1000);
	RooRealVar* yield2S = new RooRealVar("yield2S", "", 100);
	RooRealVar* yield3S = new RooRealVar("yield3S", "", 10);

	/// Assign signal and background shape name to read the file for the yield extraction results
	const char* signalShapeName = "SymDSCB";

	// background shape array: ChebychevOrderN or ExpTimesErr

	const Int_t nCosThetaBinsMax = 20;
	const Int_t nPhiBinsMax = 10;

	std::string bkgShapeName[nCosThetaBinsMax][nPhiBinsMax];

	// fill the background shape array with ChebychevOrder2
	// std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ChebychevOrder2");

	// exceptions
	// bkgShapeName[1][1] = "ChebychevOrder1";
	// bkgShapeName[1][3] = "ChebychevOrder1";
	// bkgShapeName[1][4] = "ChebychevOrder1";
	// bkgShapeName[1][5] = "ChebychevOrder1";

	// bkgShapeName[2][2] = "ChebychevOrder1";

	// bkgShapeName[3][0] = "ChebychevOrder1";
	// bkgShapeName[3][5] = "ChebychevOrder1";

	// fill the background shape array with ExpTimesErr
	std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ExpTimesErr");

	// // exceptions
	// bkgShapeName[3][0] = "ChebychevOrder1";

	/// "Standard" procedure: extract the yields per bin
	TH2D* yieldMap = new TH2D("yieldMap", " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	TH2D* standardCorrectedMap = new TH2D("standardCorrectedMap", " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get acceptance maps
	TFile* acceptanceFile = openFile("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root");
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);

	// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// get efficiency maps
	TFile* efficiencyFile = openFile("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root");
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);

	// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// get relative systematic uncertainty of efficiency
	auto* systEff = (TH3D*)efficiencyFile->Get(RelativeSystTEfficiency3DName(refFrameName));

	// rebin uncertainty map based on costheta, phi, and pT selection
	TH1D* systEffCosThetaPhi = rebinRel3DUncCosTheta(systEff, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], ptMin, ptMax, nCosThetaBins, cosThetaBinEdges);

	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	// draw acc and eff histograms to check if the rebinning works well
	DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kTRUE);
	DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE);

	TH2D* hTotalCosThetaPhi = (TH2D*)accMapCosThetaPhi->GetTotalHistogram();

	// define histograms to draw uncertainty plots
	TH1D* statHighEffCosTheta = new TH1D("statHighEffCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* statLowEffCosTheta = new TH1D("statLowEffCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* statHighAccCosTheta = new TH1D("statHighAccCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* statLowAccCosTheta = new TH1D("statLowAccCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* yield1SUncCosTheta = new TH1D("yield1SUncCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* totalRelUncCosTheta = new TH1D("totalRelUncCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());

	/// apply weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			Double_t weight = 1;

			// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotalCosThetaPhi->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotalCosThetaPhi->GetYaxis()->GetBinCenter(iPhi + 1);

			Int_t iGlobalBin = hTotalCosThetaPhi->FindFixBin(binCenterCosTheta, binCenterPhi);

			// get the corresponding weights
			double acceptance = accMapCosThetaPhi->GetEfficiency(iGlobalBin);
			double efficiency = effMapCosThetaPhi->GetEfficiency(iGlobalBin);

			// calculate weight
			weight = 1. / (acceptance * efficiency);

			// // propagate both scale factor uncertainties and efficiency stat errors to the weight
			// double relSystUnc = systEffCosTheta->GetBinContent(iCosTheta + 1);

			// double relEffUncHigh = effMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) / efficiency;
			// double relEffUncLow = effMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) / efficiency;

			// double relAccUncHigh = accMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) / acceptance;
			// double relAccUncLow = accMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) / acceptance;

			// totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

			// totalUncHigh = totalRelUncHigh * efficiency * acceptance;
			// totalUncLow = totalRelUncLow * efficiency * acceptance;

			// get yields and their uncertainties
			const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[iPhi], (Int_t)phiBinEdges[iPhi + 1]);

			RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("RawData_%s", bkgShapeName[iCosTheta][iPhi].c_str()), fitModelName);

			double yield1SVal = yield1S->getVal();

			double yield1SUnc = yield1S->getError();

			// set the bin contents reflecting weights
			// only raw yield itself before correction
			yieldMap->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal);

			// yield with acceptance x efficiency correction
			standardCorrectedMap->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal * weight);

			// yieldMap->SetBinError(iCosTheta + 1, yield1SErr * weight);
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh));
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));

			yieldMap->SetBinError(iCosTheta + 1, iPhi + 1, yield1SUnc);

			standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, yield1SUnc);

			// yieldMap->SetBinError(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

			// yieldMap->SetBinError(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			// // fill uncertainty histograms
			// statHighEffCosTheta->SetBinContent(iCosTheta + 1, relEffUncHigh);
			// statLowEffCosTheta->SetBinContent(iCosTheta + 1, relEffUncLow);
			// statHighAccCosTheta->SetBinContent(iCosTheta + 1, relAccUncHigh);
			// statLowAccCosTheta->SetBinContent(iCosTheta + 1, relAccUncLow);
			// yield1SUncCosTheta->SetBinContent(iCosTheta + 1, yield1SUnc / yield1SVal);
			// totalRelUncCosTheta->SetBinContent(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

			if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
		}
	}

	/// Polarization fit
	// with Root Fit function

	// TVirtualFitter::SetDefaultFitter("Minuit");

	TCanvas* yieldCanvas = drawYieldMap(yieldMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	yieldCanvas->SetLogz();

	displayYieldUncertainties(yieldMap, nCosThetaBins, nPhiBins);

	TPaveText* kinematicsText = new TPaveText(0.14, 0.84, 0.81, 0.93, "NDCNB");
	kinematicsText->SetFillColor(4000);
	kinematicsText->SetBorderSize(0);
	kinematicsText->AddText(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	// kinematicsText->AddText(DimuonPtRangeText(ptMin, ptMax));
	kinematicsText->SetAllWith("", "align", 12);
	kinematicsText->Draw("SAME");

	yieldCanvas->Modified();
	yieldCanvas->Update();

	TCanvas* correctedMapCanvas = drawYieldMap(standardCorrectedMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	standardCorrectedMap->GetZaxis()->SetTitle("Corrected #varUpsilon(1S) Yields");
	// standardCorrectedMap->GetZaxis()->SetTitleOffset(1.1);

	standardCorrectedMap->GetZaxis()->SetRangeUser(1e-6, maxYield * 2);

	standardCorrectedMap->SetMinimum(1e-6);

	TF2* polarFunc2D = generalPolarFunc(maxYield);

	TFitResultPtr fitResults = standardCorrectedMap->Fit("polarFunc2D", "ESV");

	// Fit results

	double chi2 = fitResults->Chi2();
	double nDOF = nCosThetaBins * nPhiBins - polarFunc2D->GetNpar();

	double normVal = fitResults->Parameter(0);
	double normErr = fitResults->ParError(0);

	double lambdaThetaVal = fitResults->Parameter(1);
	double lambdaThetaErr = fitResults->ParError(1);

	double lambdaPhiVal = fitResults->Parameter(2);
	double lambdaPhiErr = fitResults->ParError(2);

	double lambdaThetaPhiVal = fitResults->Parameter(3);
	double lambdaThetaPhiErr = fitResults->ParError(3);

	kinematicsText->Draw("SAME");

	// correctedMapCanvas->Modified();
	// correctedMapCanvas->Update();

	TLegend legend2(.17, .60, .28, .84);
	legend2.SetTextSize(.05);
	// legend2.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend2.AddEntry(standardCorrectedMap, "#varUpsilon(1S) corrected yield", "lp");
	legend2.AddEntry(polarFunc2D, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaThetaVal, lambdaThetaErr), "l");
	legend2.AddEntry((TObject*)0, Form("                       #lambda_{#varphi} = %.2f #pm %.2f", lambdaPhiVal, lambdaPhiErr), "");
	legend2.AddEntry((TObject*)0, Form("                       #lambda_{#theta#varphi} = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr), "");
	legend2.AddEntry((TObject*)0, Form("                       n  = %.2f #pm %.2f", normVal, normErr), "");

	legend2.DrawClone();

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.045);
	textChi2.DrawLatexNDC(0.74, 0.044, Form("#chi^{2} / n_{dof} = %.2f", chi2 / nDOF));

	gPad->Update();

	// draw uncertainties

	// TCanvas* errCanvas = drawUncertaintyPlot(refFrameName, systEffCosTheta, statHighEffCosTheta, statLowEffCosTheta, statHighAccCosTheta, statLowAccCosTheta, yield1SUncCosTheta, totalRelUncCosTheta);

	// TLegend legend3(.22, .9, .5, .61);
	// legend3.SetTextSize(.04);
	// legend3.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	// legend3.AddEntry(systEffCosTheta, "syst uncer of Efficiency (muon SF)", "lp");
	// legend3.AddEntry(statHighEffCosTheta, "stat uncer high of Efficiency", "lp");
	// legend3.AddEntry(statLowEffCosTheta, "stat uncer low of Efficiency", "lp");
	// legend3.AddEntry(statHighAccCosTheta, "stat uncer high of Acceptance", "lp");
	// legend3.AddEntry(statLowAccCosTheta, "stat uncer low of Acceptance", "lp");
	// legend3.AddEntry(yield1SUncCosTheta, "raw yield extraction uncertainty", "lp");
	// legend3.AddEntry(totalRelUncCosTheta, "total uncertainty", "lp");

	// legend3.DrawClone();

	// gPad->Update();

	/// calculate chi2 / nDOF by hand for cross-check

	// calculateChi2(standardCorrectedHist, PolarFunc, nCosThetaBins);

	/// contour plot
	// (ref: https://root-forum.cern.ch/t/roofit-minos-errors-for-2-parameters-of-interest/16157)

	// // set the confidence level
	// gMinuit->SetErrorDef(2.30); // 1 sigma corresponds to delchi2 = 2.30
	// TGraph* contourPlot1 = (TGraph*)gMinuit->Contour(1000, 1, 0); // Contour(number of points, lambda_theta, normalization factor)

	// gMinuit->SetErrorDef(6.18); // 2 sigma corresponds to delchi2 = 6.18
	// TGraph* contourPlot2 = (TGraph*)gMinuit->Contour(1000, 1, 0);

	// TCanvas* contourCanvas = drawContourPlots(ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], refFrameName, contourPlot1, contourPlot2);

	//

	// save canvas

	gSystem->mkdir("YieldMap/2D", kTRUE);
	yieldCanvas->SaveAs(Form("YieldMap/2D/YieldMapCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");

	gSystem->mkdir("YieldMap/2D", kTRUE);
	correctedMapCanvas->SaveAs(Form("YieldMap/2D/CorrectedMapCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");

	// gSystem->mkdir("UncertaintyPlots/2D", kTRUE);
	// errCanvas->SaveAs(Form("UncertaintyPlots/2D/uncertainty%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiBinEdges[0], phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");

	// gSystem->mkdir("ContourPlots/2D", kTRUE);
	// contourCanvas->SaveAs(Form("ContourPlots/2D/contour%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiBinEdges[0], phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
}
