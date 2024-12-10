#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
//#include "../sPlot/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "PolarFitHelpers.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"

void rawYield_2D_customizedFits(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 6, Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = gUpsilonState, Bool_t LEGOplot = kTRUE) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	//using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Bin edges and width
	// Set the bin edges along the cosTheta/phi axis depending on the number of bins, min and max values
	// (If want to use non-uniform bin width, the bin edges should be pre-defined in PolarFitHelpers.h)
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	RooRealVar* yield1S = new RooRealVar("yield1S", "", 1000);
	RooRealVar* yield2S = new RooRealVar("yield2S", "", 100);
	RooRealVar* yield3S = new RooRealVar("yield3S", "", 10);

	/// Assign signal and background shape name to read the file for the yield extraction results
	const char* signalShapeName = "SymDSCB";

	// background shape array: ChebychevOrderN or ExpTimesErr

	const Int_t nCosThetaBinsMax = 20;
	const Int_t nPhiBinsMax = 10;

	std::string bkgShapeName[nCosThetaBinsMax][nPhiBinsMax];

	// // fill the background shape array with ChebychevOrder2
	// std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ChebychevOrder2");

	// // exceptions
	// // bkgShapeName[1][1] = "ChebychevOrder1";
	// // bkgShapeName[1][3] = "ChebychevOrder1";
	// // bkgShapeName[1][4] = "ChebychevOrder1";
	// // bkgShapeName[1][5] = "ChebychevOrder1";

	// // bkgShapeName[2][2] = "ChebychevOrder1";

	// // bkgShapeName[3][0] = "ChebychevOrder1";
	// // bkgShapeName[3][5] = "ChebychevOrder1";

	// // HX, 2 < pT < 6 GeV/c, |phi|
	// bkgShapeName[4][0] = "ExpTimesErr";
	// bkgShapeName[4][1] = "ExpTimesErr";
	// bkgShapeName[4][2] = "ExpTimesErr";

	// // HX, 12 < pT < 20
	// bkgShapeName[1][1] = "ChebychevOrder1";
	// bkgShapeName[2][2] = "ChebychevOrder1";
	// bkgShapeName[3][2] = "ChebychevOrder1";

	// // CS, 6 < pT < 12 GeV/c, |phi|
	// bkgShapeName[4][0] = "ChebychevOrder1";

	// // CS, 12 < pT < 20
	// bkgShapeName[1][2] = "ChebychevOrder1";
	// bkgShapeName[3][2] = "ChebychevOrder1";

	// fill the background shape array with ExpTimesErr
	std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ExpTimesErr");

	// exceptions

	// // HX, 2 < pT < 6 GeV/c (-180 to 180)
	// bkgShapeName[0][2] = "ChebychevOrder1";
	// bkgShapeName[0][3] = "ChebychevOrder1";

	// bkgShapeName[3][2] = "ChebychevOrder2";
	// bkgShapeName[3][3] = "ChebychevOrder2";

	// bkgShapeName[3][0] = "ChebychevOrder2";
	// bkgShapeName[3][5] = "ChebychevOrder2";

	// // HX, 2 < pT < 6 GeV/c (abs(phi))
	// bkgShapeName[0][0] = "ChebychevOrder1";

	// bkgShapeName[3][0] = "ChebychevOrder2";

	// bkgShapeName[4][2] = "ChebychevOrder1";

	/// "Standard" procedure: extract the yields per bin
	TH2D* yieldMap = new TH2D("yieldMap", " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	TH2D* standardCorrectedMap = new TH2D("standardCorrectedMap", " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	Double_t lambdaTheta = 0., lambdaPhi = 0., lambdaThetaPhi = 0.;

	// get acceptance maps
	TFile* acceptanceFile = openFile(Form("../MonteCarlo/AcceptanceMaps/%dS/AcceptanceResults%s.root", iState, gMuonAccName));
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);

	// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// get efficiency maps
	TFile* efficiencyFile = openFile(Form("../MonteCarlo/EfficiencyMaps/%dS/EfficiencyResults%s.root", iState, gMuonAccName));
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);

	// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// get relative systematic uncertainty of efficiency
	auto* systEff = (TH3D*)efficiencyFile->Get(RelativeSystTEfficiency3DName(refFrameName));

	// rebin uncertainty map based on costheta, phi, and pT selection
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	//	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	// draw acc and eff histograms to check if the rebinning works well
	DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kTRUE, kTRUE);

	DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);

	TH2D* hTotalCosThetaPhi = (TH2D*)accMapCosThetaPhi->GetTotalHistogram();

	// define a histogram to draw weight map
	TH2D* weightMap = new TH2D("weightMap", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	// define histograms to draw uncertainty plots
	TH2D* relSystEffCosThetaPhi = new TH2D("relSystEffCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* statHighEffCosThetaPhi = new TH2D("statHighEffCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* statLowEffCosThetaPhi = new TH2D("statLowEffCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* statHighAccCosThetaPhi = new TH2D("statHighAccCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* statLowAccCosThetaPhi = new TH2D("statLowAccCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* yield1SUncCosThetaPhi = new TH2D("yield1SUncCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* totalRelUncCosThetaPhi = new TH2D("totalRelUncCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

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
			if (acceptance == 0 || efficiency == 0)
				weight = 0.;
			else
				weight = 1. / (acceptance * efficiency);

			weightMap->SetBinContent(iCosTheta + 1, iPhi + 1, weight);

			// propagate both scale factor uncertainties and efficiency stat errors to the weight
			double relSystUnc = systEffCosThetaPhi->GetBinContent(iCosTheta + 1, iPhi + 1) / efficiency;

			double relEffUncHigh = effMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
			double relEffUncLow = effMapCosThetaPhi->GetEfficiencyErrorLow(iGlobalBin) / efficiency;

			double relAccUncHigh = accMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
			double relAccUncLow = accMapCosThetaPhi->GetEfficiencyErrorLow(iGlobalBin) / acceptance;

			// totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

			totalRelUncHigh = TMath::Hypot(relEffUncHigh, relAccUncHigh);
			totalRelUncLow = TMath::Hypot(relEffUncLow, relAccUncLow);

			totalUncHigh = totalRelUncHigh * efficiency * acceptance;
			totalUncLow = totalRelUncLow * efficiency * acceptance;

			// get yields and their uncertainties
			// const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[iPhi], (Int_t)phiBinEdges[iPhi + 1]);

			// RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("RawData_%s", bkgShapeName[iCosTheta][iPhi].c_str()), fitModelName, Form("%s%s", gMuonAccName, "_absphi"));

			// double yield1SVal = (yield1S->getVal());

			// double yield1SUnc = (yield1S->getError());

			Int_t absiPhi = yieldMap->GetYaxis()->FindBin(fabs(binCenterPhi)) - 1;
			cout << "absiPhi: " << absiPhi << endl;

			const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[absiPhi], (Int_t)phiBinEdges[absiPhi + 1]);

			RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("%s", bkgShapeName[iCosTheta][absiPhi].c_str()), fitModelName, Form("%s", gMuonAccName));

			cout << "signalYields: " << signalYields << endl;

			// double yield1SVal = (yield1S->getVal())/2.;

			// double yield1SUnc = (yield1S->getError())/2.;
			double yield1SVal = (yield1S->getVal());

			double yield1SUnc = (yield1S->getError());

			// set the bin contents reflecting weights
			// only raw yield itself before correction
			yieldMap->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal);

			// yield with acceptance x efficiency correction
			standardCorrectedMap->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal * weight);

			// yieldMap->SetBinError(iCosTheta + 1, yield1SErr * weight);
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh));
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));

			yieldMap->SetBinError(iCosTheta + 1, iPhi + 1, yield1SUnc);

			// yieldMap->SetBinError(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

			// yieldMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			// fill uncertainty histograms
			relSystEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relSystUnc);
			statHighEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relEffUncHigh);
			statLowEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relEffUncLow);
			statHighAccCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relAccUncHigh);
			statLowAccCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relAccUncLow);
			yield1SUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SUnc / yield1SVal);
			totalRelUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

			// if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight / 2.;
			if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
		}
	}

	/// Polarization fit
	// with Root Fit function (E: minos on, S: save results, V: Verbose, I: integral, M: imporve algorithm, R: specificed range)

	// TVirtualFitter::SetDefaultFitter("Minuit");

	// draw 1/(acceptance x efficiency) map

	TCanvas* weightCanvas = draw2DMap(weightMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE);

	weightMap->GetZaxis()->SetTitle("1 / (acc x eff)");
	weightMap->GetZaxis()->SetTitleOffset(1.);

	display2DMapContents(weightMap, nCosThetaBins, nPhiBins, kFALSE);

	TPaveText* kinematicsText = new TPaveText(0.14, 0.84, 0.81, 0.93, "NDCNB");
	kinematicsText->SetFillColor(4000);
	kinematicsText->SetBorderSize(0);
	kinematicsText->AddText(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	kinematicsText->SetAllWith("", "align", 12);
	kinematicsText->Draw("SAME");

	weightCanvas->Modified();
	weightCanvas->Update();

	// draw yield map before applying corrections

	TCanvas* yieldCanvas = draw2DMap(yieldMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, LEGOplot);

	yieldMap->GetZaxis()->SetTitle(Form("#varUpsilon(%dS) yields", iState));
	yieldMap->GetZaxis()->SetTitleOffset(1.);

	yieldMap->GetZaxis()->SetRangeUser(0, maxYield * 2);
	// yieldMap->SetMinimum(0);

	// yieldCanvas->SetLogz(); // useful when the one of the bins has an exceptionally high value :')

	if (!LEGOplot) display2DMapContents(yieldMap, nCosThetaBins, nPhiBins, kTRUE);

	kinematicsText->Draw("SAME");

	yieldCanvas->Modified();
	yieldCanvas->Update();

	// draw yield map corrected by acceptance and efficiency

	TCanvas* correctedMapCanvas = draw2DMap(standardCorrectedMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, LEGOplot);

	standardCorrectedMap->GetZaxis()->SetTitle("Corrected #varUpsilon(1S) Yields");
	if (!LEGOplot) standardCorrectedMap->GetZaxis()->SetTitleOffset(1.1);

	standardCorrectedMap->GetZaxis()->SetRangeUser(0, maxYield * 2);

	// standardCorrectedMap->SetMinimum(0);

	if (!LEGOplot) display2DMapContents(standardCorrectedMap, nCosThetaBins, nPhiBins, kTRUE);

	kinematicsText->Draw("SAME");

	/// fit!!!

	TF2* polarFunc2D = getGeneralPolarFunc(maxYield);

	if (LEGOplot) {
		TFitResultPtr fitResults = standardCorrectedMap->Fit("generalPolarFunc", "ESVIMR");

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

		double lambdaTildeVal = (lambdaThetaVal + 3. * lambdaPhiVal) / (1. - lambdaPhiVal);
		double lambdaTildeErr = TMath::Hypot(1. / (1. - lambdaPhiVal) * lambdaPhiErr, (3. - lambdaThetaVal - 6. * lambdaPhiVal) / TMath::Power((1. - lambdaPhiVal), 2) * lambdaPhiErr);

		RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1.2, 1.2);
		RooRealVar lambdaPhi("lambdaPhi", "lambdaPhi", -1.2, 1.2);
		RooRealVar lambdaThetaPhi("lambdaThetaPhi", "lambdaThetaPhi", -1.2, 1.2);
		RooRealVar lambdaTilde("lambdaTilde", "lambdaTilde", -1.2, 1.2);

		lambdaTheta.setVal(lambdaThetaVal);
		lambdaTheta.setError(lambdaThetaErr);

		lambdaPhi.setVal(lambdaPhiVal);
		lambdaPhi.setError(lambdaPhiErr);

		lambdaThetaPhi.setVal(lambdaThetaPhiVal);
		lambdaThetaPhi.setError(lambdaThetaPhiErr);

		lambdaTilde.setVal(lambdaTildeVal);
		lambdaTilde.setError(lambdaTildeErr);

		RooArgSet* savedParams = new RooArgSet(lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde);

		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

		SavePolarizationFitParameters(savedParams, "rootFit", fitModelName);

		TLegend legend2(.17, .55, .28, .84);
		legend2.SetTextSize(.05);
		legend2.AddEntry(standardCorrectedMap, "#varUpsilon(1S) corrected yield", "lp");
		legend2.AddEntry(polarFunc2D, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaThetaVal, lambdaThetaErr), "l");
		legend2.AddEntry((TObject*)0, Form("                       #lambda_{#varphi} = %.2f #pm %.2f", lambdaPhiVal, lambdaPhiErr), "");
		legend2.AddEntry((TObject*)0, Form("                       #lambda_{#theta#varphi} = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr), "");
		legend2.AddEntry((TObject*)0, Form("                       #tilde{#lambda} = %.2f #pm %.2f", lambdaTildeVal, lambdaTildeErr), "");
		legend2.AddEntry((TObject*)0, Form("                       n  = %.2f #pm %.2f", normVal, normErr), "");

		legend2.DrawClone();

		TLatex textChi2;
		textChi2.SetTextAlign(12);
		textChi2.SetTextSize(0.045);
		textChi2.DrawLatexNDC(0.74, 0.044, Form("#chi^{2} / n_{dof} = %.2f", chi2 / nDOF));
	}

	gPad->Update();

	/// draw uncertainty 2D plots

	// statistical uncertainty of acceptance up
	TCanvas* statHighAccCanvas = draw2DMap(statHighAccCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	statHighAccCosThetaPhi->GetZaxis()->SetTitle("stat uncer high of acceptance");
	statHighAccCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText->Draw("SAME");

	display2DMapContents(statHighAccCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	statHighAccCanvas->Modified();
	statHighAccCanvas->Update();

	// statistical uncertainty of acceptance down
	TCanvas* statLowAccCanvas = draw2DMap(statLowAccCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	statLowAccCosThetaPhi->GetZaxis()->SetTitle("stat uncer low of acceptance");
	statLowAccCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText->Draw("SAME");

	display2DMapContents(statLowAccCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	statLowAccCanvas->Modified();
	statLowAccCanvas->Update();

	// statistical uncertainty of efficiency up
	TCanvas* statHighEffCanvas = draw2DMap(statHighEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	statHighEffCosThetaPhi->GetZaxis()->SetTitle("stat uncer high of efficiency");
	statHighEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText->Draw("SAME");

	display2DMapContents(statHighEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	statHighEffCanvas->Modified();
	statHighEffCanvas->Update();

	// statistical uncertainty of efficiency down
	TCanvas* statLowEffCanvas = draw2DMap(statLowEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	statLowEffCosThetaPhi->GetZaxis()->SetTitle("stat uncer low of efficiency");
	statLowEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText->Draw("SAME");

	display2DMapContents(statLowEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	statLowEffCanvas->Modified();
	statLowEffCanvas->Update();

	// systematic uncertainty of efficiency (muon scale factor)
	TCanvas* systEffCanvas = draw2DMap(relSystEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	relSystEffCosThetaPhi->GetZaxis()->SetTitle("syst uncer of efficiency (muon SF)");
	relSystEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText->Draw("SAME");

	display2DMapContents(relSystEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	systEffCanvas->Modified();
	systEffCanvas->Update();

	// yield extraction uncertainty
	TCanvas* yieldUncCanvas = draw2DMap(yield1SUncCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	yield1SUncCosThetaPhi->GetZaxis()->SetTitle("raw yield statistical uncertainty");
	//yield1SUncCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText->Draw("SAME");

	display2DMapContents(yield1SUncCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	yieldUncCanvas->Modified();
	yieldUncCanvas->Update();

	// total uncertatinty

	TCanvas* totalUncCanvas = draw2DMap(totalRelUncCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	totalRelUncCosThetaPhi->GetZaxis()->SetTitle("total uncertainty");
	totalRelUncCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText->Draw("SAME");

	display2DMapContents(totalRelUncCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	totalUncCanvas->Modified();
	totalUncCanvas->Update();

	/// save canvas
	gSystem->mkdir(Form("EfficiencyMaps/%dS", iState), kTRUE);
	weightCanvas->SaveAs(Form("EfficiencyMaps/%dS/WeightsMapCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");

	gSystem->mkdir("YieldMap/2D", kTRUE);
	yieldCanvas->SaveAs(Form("YieldMap/2D/YieldMapCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");

	gSystem->mkdir("YieldMap/2D", kTRUE);
	correctedMapCanvas->SaveAs(Form("YieldMap/2D/CorrectedMapCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");

	gSystem->mkdir("UncertaintyPlots/2D", kTRUE);
	statHighAccCanvas->SaveAs(Form("UncertaintyPlots/2D/statHighAcc%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	statLowAccCanvas->SaveAs(Form("UncertaintyPlots/2D/statLowAcc%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	statHighEffCanvas->SaveAs(Form("UncertaintyPlots/2D/statHighEff%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	statLowEffCanvas->SaveAs(Form("UncertaintyPlots/2D/statLowEff%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	systEffCanvas->SaveAs(Form("UncertaintyPlots/2D/sysEff%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	yieldUncCanvas->SaveAs(Form("UncertaintyPlots/2D/yieldUnc%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	totalUncCanvas->SaveAs(Form("UncertaintyPlots/2D/totalUnc%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");

	// save the histograms to the root files to see them with the root viewer
	TFile* EfficiencyOutFile = new TFile(Form("EfficiencyMaps/%dS/efficiencyHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.root", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	accMapCosThetaPhi->Write();
	effMapCosThetaPhi->Write();
	systEffCosThetaPhi->Write();
	weightMap->Write();
	EfficiencyOutFile->Close();

	TFile* ResultOutFile = new TFile(Form("YieldMap/2D/resultsHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.root", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	yieldMap->Write();
	standardCorrectedMap->Write();
	ResultOutFile->Close();

	TFile* uncOutFile = new TFile(Form("UncertaintyPlots/2D/uncertaintyHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.root", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	statLowAccCosThetaPhi->Write();
	statLowEffCosThetaPhi->Write();
	statHighAccCosThetaPhi->Write();
	statHighEffCosThetaPhi->Write();
	relSystEffCosThetaPhi->Write();
	yield1SUncCosThetaPhi->Write();
	totalRelUncCosThetaPhi->Write();
	uncOutFile->Close();
}
