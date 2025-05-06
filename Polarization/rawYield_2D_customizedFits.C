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

void rawYield_2D_customizedFits(Int_t ptMin = 0, Int_t ptMax = 2, const char* muonAccName = "UpsilonTriggerThresholds", const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180, Int_t iState = gUpsilonState, Bool_t LEGOplot = kTRUE, const char* defaultBkgShapeName = "ExpTimesErr") { //Chebychev
	writeExtraText = true;                                                                                                                                                                                                                                                                                                                                                                                                   // if extra text
	extraText = "       Preliminary";

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
	if (strcmp(defaultBkgShapeName, "Chebychev") == 0) {
		std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ChebychevOrder2");

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
	}

	if (strcmp(defaultBkgShapeName, "ExpTimesErr") == 0) {
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
	}

	/// define a TH2D for the yield map before any correction
	TH2D* yieldMap = new TH2D("yieldMap", " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	/// define a TH2D for the yield map after correction
	/// "Standard" procedure means extract the yields per bin
	TH2D* standardCorrectedMap = new TH2D("standardCorrectedMap", " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	/// Read acceptance and efficiency files for correction
	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	/// polarization weight of the acceptance and efficiency maps (for now all 0)
	Double_t lambdaTheta = 0., lambdaPhi = 0., lambdaThetaPhi = 0.;

	Bool_t isPhiFolded = kTRUE;

	/// get acceptance maps
	const char* accMapPath = AcceptanceResultsPath(muonAccName);

	TFile* acceptanceFile = openFile(Form("%s/AcceptanceResults%s.root", accMapPath, isPhiFolded ? "" : "_fullPhi"));
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);

	/// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	/// get efficiency maps
	const char* effMapPath = EfficiencyResultsPath(muonAccName);

	TFile* efficiencyFile = openFile(Form("%s/EfficiencyResults%s.root", effMapPath, isPhiFolded ? "" : "_fullPhi"));
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);

	/// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	/// get relative systematic uncertainty of efficiency
	auto* systEff = (TH3D*)efficiencyFile->Get(SystTEfficiency3DName(refFrameName));

	/// rebin uncertainty map based on costheta, phi, and pT selection
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

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

	/// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	/// apply weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			Double_t weight = 1;

			/// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotalCosThetaPhi->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotalCosThetaPhi->GetYaxis()->GetBinCenter(iPhi + 1);

			Int_t iGlobalBin = hTotalCosThetaPhi->FindFixBin(binCenterCosTheta, binCenterPhi);

			/// get the corresponding weights
			double acceptance = accMapCosThetaPhi->GetEfficiency(iGlobalBin);
			double efficiency = effMapCosThetaPhi->GetEfficiency(iGlobalBin);

			/// calculate weight
			if (acceptance == 0 || efficiency == 0)
				/// if acceptance or efficiency is 0, set the weight to 0
				weight = 0.;
			else
				weight = 1. / (acceptance * efficiency);

			/// fill the weight map
			weightMap->SetBinContent(iCosTheta + 1, iPhi + 1, weight);

			/// Error propagation
			// propagate both scale factor uncertainties and efficiency stat errors to the weight
			double relSystUnc = systEffCosThetaPhi->GetBinContent(iCosTheta + 1, iPhi + 1) / efficiency;

			double relEffUncHigh = effMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
			double relEffUncLow = effMapCosThetaPhi->GetEfficiencyErrorLow(iGlobalBin) / efficiency;

			double relAccUncHigh = accMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
			double relAccUncLow = accMapCosThetaPhi->GetEfficiencyErrorLow(iGlobalBin) / acceptance;

			// totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

			/// apply only the statistical uncertainties of the efficiency and acceptance
			totalRelUncHigh = TMath::Hypot(relEffUncHigh, relAccUncHigh);
			totalRelUncLow = TMath::Hypot(relEffUncLow, relAccUncLow);

			totalUncHigh = totalRelUncHigh * efficiency * acceptance;
			totalUncLow = totalRelUncLow * efficiency * acceptance;

			/// Get yields and their uncertainties
			Int_t absiPhi = yieldMap->GetYaxis()->FindBin(fabs(binCenterPhi)) - 1;
			cout << "absiPhi: " << absiPhi << endl;

			const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[absiPhi], (Int_t)phiBinEdges[absiPhi + 1]);

			RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("%s", bkgShapeName[iCosTheta][absiPhi].c_str()), fitModelName, Form("%s", muonAccName));

			cout << "signalYields: " << signalYields << endl;

			/// get the yield value
			double yield1SVal = (yield1S->getVal());

			double yield1SUnc = (yield1S->getError());

			/// Set the bin contents reflecting weights
			/// only raw yield itself before correction
			yieldMap->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal);

			/// Set the bin error of the raw yield
			// yieldMap->SetBinError(iCosTheta + 1, yield1SErr * weight);
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh));
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));

			yieldMap->SetBinError(iCosTheta + 1, iPhi + 1, yield1SUnc);
			// yieldMap->SetBinError(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));
			// yieldMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			/// yield with acceptance x efficiency correction
			standardCorrectedMap->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal * weight);

			/// set the bin error of the corrected yield
			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);
			standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, yield1SUnc * weight);

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

	/// draw 1/(acceptance x efficiency) map

	TCanvas* weightCanvas = draw2DMap(weightMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE);

	weightMap->GetZaxis()->SetTitle("1 / (acc x eff)");
	weightMap->GetZaxis()->SetTitleOffset(0.8);

	display2DMapContents(weightMap, nCosThetaBins, nPhiBins, kFALSE);

	TPaveText* kinematicsText = new TPaveText(0.15, 0.79, 0.77, 0.99, "NDCNB");
	kinematicsText->SetFillColor(4000);
	kinematicsText->SetBorderSize(0);
	kinematicsText->AddText(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	kinematicsText->SetAllWith("", "align", 31);
	kinematicsText->Draw("SAME");

	weightCanvas->Modified();
	weightCanvas->Update();

	/// draw yield map before applying corrections

	TCanvas* yieldCanvas = draw2DMap(yieldMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, LEGOplot);

	yieldMap->GetZaxis()->SetTitle(Form("#varUpsilon(%dS) yields", iState));
	yieldMap->GetZaxis()->SetTitleOffset(1.);

	yieldMap->GetZaxis()->SetRangeUser(0, maxYield * 2);
	// yieldMap->SetMinimum(0);

	// yieldCanvas->SetLogz(); // useful when the one of the bins has an exceptionally high value :')

	TPaveText* kinematicsText_2D = new TPaveText(0.17, 0.83, 0.79, 0.93, "NDCNB");
	kinematicsText_2D->SetFillColor(4000);
	kinematicsText_2D->SetBorderSize(0);
	kinematicsText_2D->AddText(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonRapidityRangeText(gRapidityMin, gRapidityMax)));
	kinematicsText_2D->AddText(Form("%s, %s frame", DimuonPtRangeText(ptMin, ptMax), (strstr(refFrameName, "CS")) ? "Collins-Soper" : "Helicity"));

	kinematicsText_2D->SetAllWith("", "align", 12);

	if (LEGOplot)
		kinematicsText->Draw("SAME");

	else if (!LEGOplot) {
		// display2DMapContents(yieldMap, nCosThetaBins, nPhiBins, kTRUE);
		// kinematicsText_2D ->Draw("SAME");
		TLatex legend;
		legend.SetTextAlign(22);
		legend.SetTextSize(0.04);
		legend.DrawLatexNDC(.48, .86, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));
		legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
		// if (strcmp(extraString, "_TriggerAcc") == 0)
		// 	legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
		// else if (strcmp(extraString, "_SimpleAcc") == 0)
		// 	legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", iState));
		// else if (strcmp(extraString, "_2018PbPbAcc") == 0)
		// 	legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", iState));
	}

	yieldCanvas->Modified();
	yieldCanvas->Update();

	/// draw yield map corrected by acceptance and efficiency

	TCanvas* correctedMapCanvas = draw2DMap(standardCorrectedMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, LEGOplot);

	standardCorrectedMap->GetZaxis()->SetTitle("#varUpsilon(1S) corrected yield");
	standardCorrectedMap->GetZaxis()->SetTitleSize(0.055);
	standardCorrectedMap->GetZaxis()->SetTitleOffset(1.4);
	standardCorrectedMap->GetZaxis()->SetLabelSize(0.044);

	standardCorrectedMap->GetYaxis()->SetTitleSize(0.055);
	standardCorrectedMap->GetYaxis()->SetTitleOffset(1.3);
	standardCorrectedMap->GetYaxis()->SetLabelSize(0.044);
	standardCorrectedMap->GetYaxis()->SetLabelOffset(0);

	standardCorrectedMap->GetXaxis()->SetTitleSize(0.055);
	standardCorrectedMap->GetXaxis()->SetTitleOffset(1.2);
	standardCorrectedMap->GetXaxis()->SetLabelSize(0.044);

	if (!LEGOplot) standardCorrectedMap->GetZaxis()->SetTitleOffset(1.2);

	standardCorrectedMap->GetZaxis()->SetRangeUser(0, maxYield * 1.6);

	// standardCorrectedMap->SetMinimum(0);

	if (LEGOplot)
		kinematicsText->Draw("SAME");

	else if (!LEGOplot) {
		// display2DMapContents(standardCorrectedMap, nCosThetaBins, nPhiBins, kTRUE);
		// kinematicsText_2D->Draw("SAME");
		TLatex legend;
		legend.SetTextAlign(22);
		legend.SetTextSize(0.04);
		legend.DrawLatexNDC(.48, .86, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));
		legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
		// if (strcmp(extraString, "_TriggerAcc") == 0)
		// 	legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
		// else if (strcmp(extraString, "_SimpleAcc") == 0)
		// 	legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", iState));
		// else if (strcmp(extraString, "_2018PbPbAcc") == 0)
		// 	legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", iState));
	}

	/// fit!!!

	TF2* polarFunc2D = getGeneralPolarFunc(maxYield);

	if (LEGOplot) {
		TFitResultPtr fitResults = standardCorrectedMap->Fit("generalPolarFunc", "ESVIMR" /*, "ESVIR"*/);

		// Fit results

		double chi2 = fitResults->Chi2();
		double nDOF = nCosThetaBins * nPhiBins - polarFunc2D->GetNpar();

		double pValue = TMath::Prob(chi2, nDOF); // https://root.cern.ch/doc/master/namespaceTMath.html#a3aba6abaf2bc605680c7be8fc5fd98aa

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

		RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -2., 2.);
		RooRealVar lambdaPhi("lambdaPhi", "lambdaPhi", -2., 2.);
		RooRealVar lambdaThetaPhi("lambdaThetaPhi", "lambdaThetaPhi", -2., 2.);
		RooRealVar lambdaTilde("lambdaTilde", "lambdaTilde", -2., 2.);

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

		SavePolarizationFitParameters(savedParams, "rootFit", fitModelName, muonAccName);
		// TLegend legend2(.17, .7, .23, .87);
		TLegend legend2(.18, .75, .25, .8);
		legend2.SetTextSize(.045);
		legend2.SetFillColor(0);
		legend2.SetFillStyle(1001);

		//legend2.AddEntry(standardCorrectedMap, "#varUpsilon(1S) corrected yield", "lp");
		// legend2.AddEntry(polarFunc2D, Form("fit: #lambda_{#theta}  = %.2f #pm %.2f    #lambda_{#varphi} = %.2f #pm %.2f", lambdaThetaVal, lambdaThetaErr, lambdaPhiVal, lambdaPhiErr), "l");
		// legend2.AddEntry((TObject*)0, Form("     #lambda_{#theta#varphi} = %.2f #pm %.2f  #tilde{#lambda}  = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr, lambdaTildeVal, lambdaTildeErr), "");
		legend2.AddEntry(polarFunc2D, "fit: ", "l");

		legend2.DrawClone();

		TPaveText* resultTextRight = new TPaveText(0.25, 0.70, 0.60, 0.81, "NDC"); // Adjust coordinates
		resultTextRight->SetFillColor(0);                                          // White background
		resultTextRight->SetFillStyle(1001);                                       // Solid fill
		resultTextRight->SetBorderSize(0);                                         // Optional: Thin border
		resultTextRight->SetTextSize(0.045);
		resultTextRight->SetTextAlign(12); // Align text left
		resultTextRight->SetMargin(0.03);
		resultTextRight->AddText(Form("#lambda_{#theta}  = %.2f #pm %.2f   #lambda_{#varphi} = %.2f #pm %.2f", lambdaThetaVal, lambdaThetaErr, lambdaPhiVal, lambdaPhiErr));
		resultTextRight->AddText(Form("#lambda_{#theta#varphi} = %.2f #pm %.2f   #tilde{#lambda}  = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr, lambdaTildeVal, lambdaTildeErr));
		resultTextRight->Draw();

		TPaveText* resultTextLeft = new TPaveText(0.52, 0.72, 0.90, 0.805, "NDC"); // Adjust coordinates
		resultTextLeft->SetFillColor(0);                                           // White background
		resultTextLeft->SetFillStyle(1001);                                        // Solid fill
		resultTextLeft->SetBorderSize(0);                                          // Optional: Thin border
		resultTextLeft->SetTextSize(0.045);
		resultTextLeft->SetTextAlign(12); // Align text left
		resultTextLeft->SetMargin(0.03);
		resultTextLeft->AddText(Form("#lambda_{#varphi} = %.2f #pm %.2f", lambdaPhiVal, lambdaPhiErr));
		resultTextLeft->AddText(Form("#tilde{#lambda}  = %.2f #pm %.2f", lambdaTildeVal, lambdaTildeErr));
		//		resultTextLeft->Draw();

		TLatex textChi2;
		textChi2.SetTextAlign(22);
		// textChi2.SetTextAlign(12);
		textChi2.SetTextSize(0.04);
		textChi2.DrawLatexNDC(0.72, 0.044, Form("#chi^{2} / n_{dof} = %.2f, p-value = %.2f", chi2 / nDOF, pValue));

		correctedMapCanvas->SetTopMargin(0.2);

		gPad->Modified();
		gPad->Update();
	}

	gPad->Update();

	/// draw uncertainty 2D plots

	/// statistical uncertainty of acceptance up
	TCanvas* statHighAccCanvas = draw2DMap(statHighAccCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	statHighAccCosThetaPhi->GetZaxis()->SetTitle("stat uncer high of acceptance");
	statHighAccCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText_2D->Draw("SAME");

	display2DMapContents(statHighAccCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	statHighAccCanvas->Modified();
	statHighAccCanvas->Update();

	/// statistical uncertainty of acceptance down
	TCanvas* statLowAccCanvas = draw2DMap(statLowAccCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	statLowAccCosThetaPhi->GetZaxis()->SetTitle("stat uncer low of acceptance");
	statLowAccCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText_2D->Draw("SAME");

	display2DMapContents(statLowAccCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	statLowAccCanvas->Modified();
	statLowAccCanvas->Update();

	/// statistical uncertainty of efficiency up
	TCanvas* statHighEffCanvas = draw2DMap(statHighEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	statHighEffCosThetaPhi->GetZaxis()->SetTitle("stat uncer high of efficiency");
	statHighEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText_2D->Draw("SAME");

	display2DMapContents(statHighEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	statHighEffCanvas->Modified();
	statHighEffCanvas->Update();

	/// statistical uncertainty of efficiency down
	TCanvas* statLowEffCanvas = draw2DMap(statLowEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	statLowEffCosThetaPhi->GetZaxis()->SetTitle("stat uncer low of efficiency");
	statLowEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText_2D->Draw("SAME");

	display2DMapContents(statLowEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	statLowEffCanvas->Modified();
	statLowEffCanvas->Update();

	/// systematic uncertainty of efficiency (muon scale factor)
	TCanvas* systEffCanvas = draw2DMap(relSystEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	relSystEffCosThetaPhi->GetZaxis()->SetTitle("syst uncer of efficiency (muon SF)");
	relSystEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText_2D->Draw("SAME");

	display2DMapContents(relSystEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	systEffCanvas->Modified();
	systEffCanvas->Update();

	/// yield extraction uncertainty
	TCanvas* yieldUncCanvas = draw2DMap(yield1SUncCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	yield1SUncCosThetaPhi->GetZaxis()->SetTitle("raw yield statistical uncertainty");
	//yield1SUncCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText_2D->Draw("SAME");

	display2DMapContents(yield1SUncCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	yieldUncCanvas->Modified();
	yieldUncCanvas->Update();

	/// total uncertatinty

	TCanvas* totalUncCanvas = draw2DMap(totalRelUncCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	totalRelUncCosThetaPhi->GetZaxis()->SetTitle("total uncertainty");
	totalRelUncCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	kinematicsText_2D->Draw("SAME");

	display2DMapContents(totalRelUncCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	totalUncCanvas->Modified();
	totalUncCanvas->Update();

	/// save canvas

	const char* commonOutputName = Form("cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax, cosThetaMin, cosThetaMax);

	gSystem->mkdir(Form("EfficiencyMaps/%s", muonAccName), kTRUE);
	weightCanvas->SaveAs(Form("EfficiencyMaps/%s/WeightsMapCosTheta%s_%s", muonAccName, refFrameName, commonOutputName), "RECREATE");

	const char* yieldMapsPath = Form("YieldMaps/%s", muonAccName);
	gSystem->mkdir(yieldMapsPath, kTRUE);
	if (!LEGOplot)
		yieldCanvas->SaveAs(Form("%s/YieldMapCosTheta%s_%s_%s", yieldMapsPath, refFrameName, defaultBkgShapeName, commonOutputName), "RECREATE");
	else
		yieldCanvas->SaveAs(Form("%s/YieldMapCosTheta%s_3D_%s_%s", yieldMapsPath, refFrameName, defaultBkgShapeName, commonOutputName), "RECREATE");

	if (!LEGOplot)
		correctedMapCanvas->SaveAs(Form("%s/CorrectedMapCosTheta%s_%s_%s", yieldMapsPath, refFrameName, defaultBkgShapeName, commonOutputName), "RECREATE");
	else
		correctedMapCanvas->SaveAs(Form("%s/CorrectedMapCosTheta%s_fit_%s_%s", yieldMapsPath, refFrameName, defaultBkgShapeName, commonOutputName), "RECREATE");

	/// Uncertainty plots
	const char* uncertaintyPath = Form("UncertaintyPlots/%s", muonAccName);

	gSystem->mkdir(uncertaintyPath, kTRUE);
	statHighAccCanvas->SaveAs(Form("%s/statHighAcc%s_%s", uncertaintyPath, refFrameName, commonOutputName), "RECREATE");
	statLowAccCanvas->SaveAs(Form("%s/statLowAcc%s_%s", uncertaintyPath, refFrameName, commonOutputName), "RECREATE");
	statHighEffCanvas->SaveAs(Form("%s/statHighEff%s_%s", uncertaintyPath, refFrameName, commonOutputName), "RECREATE");
	statLowEffCanvas->SaveAs(Form("%s/statLowEff%s_%s", uncertaintyPath, refFrameName, commonOutputName), "RECREATE");
	systEffCanvas->SaveAs(Form("%s/sysEff%s_%s", uncertaintyPath, refFrameName, commonOutputName), "RECREATE");
	yieldUncCanvas->SaveAs(Form("%s/yieldUnc%s_%s", uncertaintyPath, refFrameName, commonOutputName), "RECREATE");
	totalUncCanvas->SaveAs(Form("%s/totalUnc%s_%s", uncertaintyPath, refFrameName, commonOutputName), "RECREATE");

	// save the histograms to the root files to see them with the root viewer
	TFile* EfficiencyOutFile = new TFile(Form("EfficiencyMaps/%s/efficiencyHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.root", muonAccName, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
	accMapCosThetaPhi->Write();
	effMapCosThetaPhi->Write();
	systEffCosThetaPhi->Write();
	weightMap->Write();
	EfficiencyOutFile->Close();

	TFile* ResultOutFile = new TFile(Form("%s/resultsHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.root", yieldMapsPath, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
	yieldMap->Write();
	standardCorrectedMap->Write();
	ResultOutFile->Close();

	TFile* uncOutFile = new TFile(Form("%s/uncertaintyHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.root", uncertaintyPath, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
	statLowAccCosThetaPhi->Write();
	statLowEffCosThetaPhi->Write();
	statHighAccCosThetaPhi->Write();
	statHighEffCosThetaPhi->Write();
	relSystEffCosThetaPhi->Write();
	yield1SUncCosThetaPhi->Write();
	totalRelUncCosThetaPhi->Write();
	uncOutFile->Close();
}

void scanRawYield_2D_customizedFits(const char* refFrameName = "CS") {
	/// loop over
	for (Int_t ptIdx = 0; ptIdx < NPtBins; ptIdx++) {
		for (Int_t idx = 0; idx < 2; idx++) {
			if (idx == 0)
				rawYield_2D_customizedFits(gPtBinning[ptIdx], gPtBinning[ptIdx + 1], gMuonAccName, refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kFALSE, "ExpTimesErr"); // 2D plot
			else
				rawYield_2D_customizedFits(gPtBinning[ptIdx], gPtBinning[ptIdx + 1], gMuonAccName, refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kTRUE, "ExpTimesErr"); // LEGO plot + fit
		}
	}
}
