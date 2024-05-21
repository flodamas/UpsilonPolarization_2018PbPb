#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/cosThetaPolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"

TEfficiency* rebinTEff3DMap(TEfficiency* TEff3DMap, Int_t phiMin = -180, Int_t phiMax = 180, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 10, const vector<Double_t>& cosThetaBinEdges = {}) {
	/// rebin efficiency maps based on costheta, phi, and pT selection

	// extract the numerator and the denominator from the 3D TEfficiency Map
	TH3D* hPassed = (TH3D*)TEff3DMap->GetPassedHistogram();
	TH3D* hTotal = (TH3D*)TEff3DMap->GetTotalHistogram();

	// obtain the bin numbers of the boundaries on phi and pt
	Int_t iPhiMin = hPassed->GetYaxis()->FindBin(phiMin);
	Int_t iPhiMax = hPassed->GetYaxis()->FindBin(phiMax);

	Int_t iPtMin = hPassed->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = hPassed->GetZaxis()->FindBin(ptMax);

	// obtain the projection histogram along the costheta axis within boundaries of phi and pt
	// (option e: calculate errors, o: only bins inside the selected range will be filled)
	TH1D* hPassedCosTheta = (TH1D*)hPassed->ProjectionX("hPassedCosTheta", iPhiMin, iPhiMax - 1, iPtMin, iPtMax - 1, "eo");
	TH1D* hTotalCosTheta = (TH1D*)hTotal->ProjectionX("hTotalCosTheta", iPhiMin, iPhiMax - 1, iPtMin, iPtMax - 1, "eo");

	// rebin the projection histogram (default: 20 bins from -1 to 1)
	TH1D* hPassedCosTheta_Rebin = (TH1D*)hPassedCosTheta->Rebin(nCosThetaBins, "hPassedCosTheta_Rebin", cosThetaBinEdges.data());
	TH1D* hTotalCosTheta_Rebin = (TH1D*)hTotalCosTheta->Rebin(nCosThetaBins, "hTotalCosTheta_Rebin", cosThetaBinEdges.data());

	// define TEfficiency using the final numerator and denominator
	TEfficiency* TEffCosTheta = new TEfficiency("TEffCosTheta", "cos #theta_{CS}; efficiency", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);

	TEffCosTheta->SetPassedHistogram(*hPassedCosTheta_Rebin, "f");
	TEffCosTheta->SetTotalHistogram(*hTotalCosTheta_Rebin, "f");

	return TEffCosTheta;
}

TH1D* rebinRel3DUnc(TH3D* systEff, Int_t phiMin = -180, Int_t phiMax = 180, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 10, const vector<Double_t>& cosThetaBinEdges = {}) {
	/// rebin efficiency maps based on costheta, phi, and pT selection
	// uncertainty addition is sqrt(pow(unc1, 2) + pow(unc2, 2)), so fold it manually

	TH1D* h1DSystEff = new TH1D("h1DSystEff", "", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);

	// obtain the bin numbers of the boundaries on phi and pt
	Int_t iPhiMin = systEff->GetYaxis()->FindBin(phiMin);
	Int_t iPhiMax = systEff->GetYaxis()->FindBin(phiMax) - 1;

	Int_t iPtMin = systEff->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = systEff->GetZaxis()->FindBin(ptMax) - 1;

	Int_t ibin = 1;

	// calculate the systematic uncertainties merging phi and pt bins
	// and set the bin content of 1D costheta hist using the calculated value
	for (int ibin = 1; ibin <= nCosThetaBins; ibin++){
		Double_t cosThetaSumSystEff = 0;

		Int_t iCosThetaMin = systEff->GetXaxis()->FindBin(cosThetaBinEdges[ibin-1]);
		Int_t iCosThetaMax = systEff->GetXaxis()->FindBin(cosThetaBinEdges[ibin])-1;

		// merge bins along the cosTheta axis
		for (int iCosTheta = iCosThetaMin; iCosTheta <= iCosThetaMax; iCosTheta++) {
			Double_t phiSumSystEff = 0;

			// sum uncertainties along the phi axis
			for (int iPhi = iPhiMin; iPhi <= iPhiMax; iPhi++) {
				Double_t ptSumSystEff = 0;

				// sum uncertainties along the pt axis
				for (int iPt = iPtMin; iPt <= iPtMax; iPt++) {
					ptSumSystEff = TMath::Hypot(ptSumSystEff, systEff->GetBinContent(iCosTheta, iPhi, iPt));
				}

				phiSumSystEff = TMath::Hypot(phiSumSystEff, ptSumSystEff);
			}

			cosThetaSumSystEff = TMath::Hypot(cosThetaSumSystEff, phiSumSystEff);
		}

		h1DSystEff->SetBinContent(ibin, cosThetaSumSystEff);

	}

	return h1DSystEff;
}

void DrawEfficiency1DHist(TEfficiency* effHist, Int_t ptMin, Int_t ptMax, Int_t iState = gUpsilonState, Bool_t	isAcc = kTRUE) {
	TCanvas* canvas = new TCanvas(effHist->GetName(), "", 600, 600);
	canvas->SetRightMargin(0.05);

	// empty frame for the axes
	TH1D* frameHist = new TH1D("frameHist", "", NCosThetaBinsHX, CosThetaBinningHX);

	frameHist->Draw(); 

	// draw efficiency plot on top of the histogram frame
	effHist->SetLineWidth(3);
	effHist->Draw("PL E0 SAME");

	// cosmetics of the histogram
	CMS_lumi(canvas, Form("Unpolarized #varUpsilon(%dS) Pythia 8 MC", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.55, .88, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));
	legend.DrawLatexNDC(.55, .8, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));

	if (strstr(effHist->GetName(), "CS")) frameHist->SetXTitle(CosThetaVarTitle("CS"));
	else frameHist->SetXTitle(CosThetaVarTitle("HX"));

	frameHist->SetYTitle(TEfficiencyMainTitle(iState));
	
	frameHist->GetXaxis()->CenterTitle();
	frameHist->GetYaxis()->CenterTitle();

	frameHist->GetXaxis()->SetRangeUser(-1, 1);
	frameHist->GetYaxis()->SetRangeUser(0, 1);

	frameHist->GetXaxis()->SetNdivisions(510, kTRUE);

	// save the plot
	gSystem->mkdir(Form("EfficiencyMaps/%dS", iState), kTRUE);
	if (isAcc) canvas->SaveAs(Form("EfficiencyMaps/%dS/acc_%s_pt%dto%d.png", iState, effHist->GetName(), ptMin, ptMax), "RECREATE");
	else canvas->SaveAs(Form("EfficiencyMaps/%dS/eff_%s_pt%dto%d.png", iState, effHist->GetName(), ptMin, ptMax), "RECREATE");
}

vector<Double_t> setCosThetaBinEdges(Int_t nCosThetaBins){

	vector<Double_t> cosThetaBinEdges;

	// define the bin edges along the cosTheta axis depending on the number of bins
	if (nCosThetaBins == 5) cosThetaBinEdges = {-0.5, -0.3, -0.1, 0.1, 0.3, 0.5};
	else if (nCosThetaBins == 6) cosThetaBinEdges = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
	else if (nCosThetaBins == 7) cosThetaBinEdges = {-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7};
	else if (nCosThetaBins == 8) cosThetaBinEdges = {-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8};
	else if (nCosThetaBins == 9) cosThetaBinEdges = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9};
	else if (nCosThetaBins == 10) cosThetaBinEdges = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
	else {
		cout << "Pre-defined binning not found. Check the input nCosThetaBins" << endl;
		exit(1);
	}

	return cosThetaBinEdges;
}

TCanvas* drawUncertainties(const char* refFrameName, TH1D* uncPlot1, TH1D* uncPlot2, TH1D* uncPlot3, TH1D* uncPlot4, TH1D* uncPlot5, TH1D* uncPlot6, TH1D* uncPlot7){

	TCanvas *errCanvas = new TCanvas("errCanvas", "errCanvas", 650, 600);

	uncPlot1->GetYaxis()->SetRangeUser(0, 1);

	uncPlot1->SetXTitle(Form("cos #theta_{%s}", refFrameName));
	uncPlot1->SetYTitle("Relative Uncertainty");

	uncPlot1->GetXaxis()->CenterTitle();
	uncPlot1->GetYaxis()->CenterTitle();

	Int_t lineWidth = 6;

	uncPlot1->SetLineWidth(lineWidth);
	uncPlot1->SetLineColor(kRed+1);

	uncPlot1->Draw();

	uncPlot2->SetLineWidth(lineWidth);
	uncPlot2->SetLineColor(kRed-7);

	uncPlot2->Draw("SAME");

	uncPlot3->SetLineWidth(lineWidth);
	uncPlot3->SetLineColor(kOrange);

	uncPlot3->Draw("SAME");

	uncPlot4->SetLineWidth(lineWidth);
	uncPlot4->SetLineColor(kAzure-9);

	uncPlot4->Draw("SAME");	

	uncPlot5->SetLineWidth(lineWidth);
	uncPlot5->SetLineColor(kBlue-3);

	uncPlot5->Draw("SAME");

	uncPlot6->SetLineWidth(lineWidth);
	uncPlot6->SetLineColor(kViolet-8);

	uncPlot6->Draw("SAME");

	uncPlot7->SetLineWidth(lineWidth);
	uncPlot7->SetLineColor(kBlack);

	uncPlot7->Draw("SAME");

	return errCanvas;
}

void rawYield_1D_customizedFits(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 10, Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = gUpsilonState) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Set up the variables
	RooRealVar* yield1S = new RooRealVar("yield1S", "", 1000);
	RooRealVar* yield2S = new RooRealVar("yield2S", "", 100);
	RooRealVar* yield3S = new RooRealVar("yield3S", "", 10);

	/// Bin edges and width 
	// Set the bin edges along the cosTheta axis depending on the number of bins 
	// (The bin edges are pre-defined, so need to modify them if different bin edges are required)
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	/// Set up the variables
	RooRealVar cosTheta("cosTheta", "", cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	cosTheta.setBins(nCosThetaBins);

	/// Assign signal and background shape name to read the file for the yield extraction results
	const char* signalShapeName = "SymDSCB";

	// background shape array: ChebychevOrderN or ExpTimesErr	
	// const char* bkgShapeName[] = {
	//   // "ChebychevOrder1",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   // "ChebychevOrder1"
	// };

	const char* bkgShapeName[] = {
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr"
	};

	/// "Standard" procedure: extract the yields per bin
	TH1D* standardCorrectedHist = new TH1D("standardCorrectedHist", " ", nCosThetaBins, cosThetaBinEdges.data());

	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get acceptance maps
	TFile* acceptanceFile = openFile("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root");
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);

	// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosTheta = rebinTEff3DMap(accMap, phiMin, phiMax, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges);

	// get efficiency maps
	TFile* efficiencyFile = openFile("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root");
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);

	// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effMapCosTheta = rebinTEff3DMap(effMap, phiMin, phiMax, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges);

	// get relative systematic uncertainty of efficiency
	auto* systEff = (TH3D*)efficiencyFile->Get(RelativeSystTEfficiency3DName(refFrameName));

	// rebin uncertainty map based on costheta, phi, and pT selection
	TH1D* systEffCosTheta = rebinRel3DUnc(systEff, phiMin, phiMax, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges);

	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;

	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	// // define arrays for TGraph (to draw AsymmError)
	// double finalDataPoints[nCosThetaBins];
	// double finalErrHigh[nCosThetaBins];
	// double finalErrLow[nCosThetaBins];
	// double cosThetaBinCenter[nCosThetaBins];

	// draw acc and eff histograms to check if the rebinning works well
	DrawEfficiency1DHist(accMapCosTheta, ptMin, ptMax, iState, kTRUE);
	DrawEfficiency1DHist(effMapCosTheta, ptMin, ptMax, iState, kFALSE);

	TH1D* statHighEffCosTheta = new TH1D("statHighEffCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* statLowEffCosTheta = new TH1D("statLowEffCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* statHighAccCosTheta = new TH1D("statHighAccCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* statLowAccCosTheta = new TH1D("statLowAccCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* yield1SUncCosTheta = new TH1D("yield1SUncCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	TH1D* totalRelUncCosTheta = new TH1D("totalRelUncCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());

	/// apply weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		Double_t weight = 0;

		// get the corresponding weights
		double acceptance = accMapCosTheta->GetEfficiency(iCosTheta + 1);
		double efficiency = effMapCosTheta->GetEfficiency(iCosTheta + 1);

		// calculate weight
		weight = 1. / (acceptance * efficiency);

		// print acc, eff, and weight values
		// cout << "bin " << iCosTheta << endl;
		// cout << "acceptance: " << acceptance << endl;
		// cout << "efficiency: " << efficiency << endl;
		// cout << "weight: " << weight << endl;
		// cout << endl;

		// propagate both scale factor uncertainties and efficiency stat errors to the weight
		double relSystUnc = systEffCosTheta->GetBinContent(iCosTheta + 1);

		double relEffUncHigh = effMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) / efficiency;
		double relEffUncLow = effMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) / efficiency;

		double relAccUncHigh = accMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) / acceptance;
		double relAccUncLow = accMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) / acceptance;

		// print statistical uncertainties of nominal efficiency and acceptance
		// cout << "eff error up: " << effMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) << endl;
		// cout << "eff error low: " << effMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) << endl;
		// cout << "acc error up: " << accMapCosTheta->GetEfficiencyErrorUp(iCosTheta + 1) << endl;
		// cout << "acc error low: " << accMapCosTheta->GetEfficiencyErrorLow(iCosTheta + 1) << endl;

		totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
		totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

		totalUncHigh = totalRelUncHigh * efficiency * acceptance;
		totalUncLow = totalRelUncLow * efficiency * acceptance;

		// get yields and their uncertainties
		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, isCSframe, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], phiMin, phiMax);

		RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("RawData_%s", bkgShapeName[iCosTheta]), fitModelName);

		double yield1SVal = yield1S->getVal();

		double yield1SUnc = yield1S->getError();

		// set the bin contents reflecting weights
		standardCorrectedHist->SetBinContent(iCosTheta + 1, yield1SVal * weight);

		// standardCorrectedHist->SetBinError(iCosTheta + 1, yield1SErr * weight);
		// standardCorrectedHist->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh));
		// standardCorrectedHist->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));
		
		// standardCorrectedHist->SetBinError(iCosTheta + 1, yield1SUnc * weight);

		// standardCorrectedHist->SetBinError(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

		standardCorrectedHist->SetBinError(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

		/// fill arrays for TGraphAsymmError
		// cosThetaBinCenter[iCosTheta] = (cosThetaBinEdges[iCosTheta] + cosThetaBinEdges[iCosTheta + 1]) / 2.;

		// finalDataPoints[iCosTheta] = yield1SVal * weight;

		// finalErrHigh[iCosTheta] = yield1SVal * TMath::Hypot(weight * yield1SErr / yield1SVal, errorWeightHigh);
		// finalErrLow[iCosTheta] = yield1SVal * TMath::Hypot(weight * yield1SErr / yield1SVal, errorWeightLow);

		// standardCorrectedHist->SetBinError(iCosTheta + 1, finalErrHigh[iCosTheta]);

		// fill uncertainty histograms
		statHighEffCosTheta->SetBinContent(iCosTheta + 1, relEffUncHigh);
		statLowEffCosTheta->SetBinContent(iCosTheta + 1, relEffUncLow);
		statHighAccCosTheta->SetBinContent(iCosTheta + 1, relAccUncHigh);
		statLowAccCosTheta->SetBinContent(iCosTheta + 1, relAccUncLow);
		yield1SUncCosTheta->SetBinContent(iCosTheta + 1, yield1SUnc / yield1SVal);
		totalRelUncCosTheta->SetBinContent(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

		// cout << "bin: " << iCosTheta+1 << endl;
		// cout << "sys: " << relSystUnc << endl;
		// cout << "effHigh: " << relEffUncHigh << endl;
		// cout << "effLow: " << relEffUncLow << endl;
		// cout << "accHigh: " << relAccUncHigh << endl;
		// cout << "accLow: " << relAccUncLow << endl;
		// cout << "yield: " << yield1SUnc/yield1SVal << endl;
		// cout << "total: " << TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) << endl;

		if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
		// cout << "Uncertainty: " << finalErrHigh[iCosTheta] << endl;
	}
	cout << "Uncertainty: " << yield1S->getError() << endl;

	/// TGraphAsymmErrors - comment out for now

	// TGraphAsymmErrors *correctedGraph = new TGraphAsymmErrors(nCosThetaBins, cosThetaBinCenter, finalDataPoints, 0, 0, finalErrLow, finalErrHigh);

	/// Polarization fit
	// with Root Fit function

	TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 650, 600);

	TF1* PolarFunc = cosThetaPolarFunc(maxYield);

	TFitResultPtr fitResults = standardCorrectedHist->Fit("PolarFunc", "ESV", "", cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]); //L:log likelihood fit (default: chi2 method), E: NINOS

	cout << "Error of hist (bin1): " << standardCorrectedHist->GetBinError(1) << endl;
	
	// Fit results

	double chi2 = fitResults->Chi2();
	double nDOF = nCosThetaBins - PolarFunc->GetNpar();

	double lambdaVal = fitResults->Parameter(0);
	double lambdaErr = fitResults->ParError(0);

	gStyle->SetOptFit(1011);

	// cosmetics

	standardCorrectedHist->SetMarkerStyle(20);
	standardCorrectedHist->SetMarkerSize(1);
	standardCorrectedHist->SetMarkerColor(kAzure + 2);

	standardCorrectedHist->GetYaxis()->SetRangeUser(0, 2 * maxYield);
	standardCorrectedHist->GetYaxis()->SetMaxDigits(3);

	standardCorrectedHist->SetXTitle(Form("cos #theta_{%s}", refFrameName));
	standardCorrectedHist->SetYTitle(Form("Events / ( %0.1f )", cosThetaStep));

	standardCorrectedHist->GetXaxis()->CenterTitle();
	standardCorrectedHist->GetYaxis()->CenterTitle();

	TLegend legend2(.22, .88, .5, .68);
	legend2.SetTextSize(.05);
	legend2.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend2.AddEntry(standardCorrectedHist, "#varUpsilon(1S) corrected yield", "lp");
	legend2.AddEntry(PolarFunc, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaVal, lambdaErr), "l");

	legend2.DrawClone();

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.045);
	textChi2.DrawLatexNDC(0.74, 0.044, Form("#chi^{2} / n_{dof} = %.2f", chi2 / nDOF));
	
	gPad->Update();

	// draw uncertainties
	
	TCanvas* errCanvas = drawUncertainties(refFrameName, systEffCosTheta, statHighEffCosTheta, statLowEffCosTheta, statHighAccCosTheta, statLowAccCosTheta, yield1SUncCosTheta, totalRelUncCosTheta);

	TLegend legend3(.22, .9, .5, .61);
	legend3.SetTextSize(.04);
	legend3.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legend3.AddEntry(systEffCosTheta, "syst uncer of Efficiency (muon SF)", "lp");
	legend3.AddEntry(statHighEffCosTheta, "stat uncer high of Efficiency", "lp");
	legend3.AddEntry(statLowEffCosTheta, "stat uncer low of Efficiency", "lp");
	legend3.AddEntry(statHighAccCosTheta, "stat uncer high of Acceptance", "lp");
	legend3.AddEntry(statLowAccCosTheta, "stat uncer low of Acceptance", "lp");
	legend3.AddEntry(yield1SUncCosTheta, "raw yield extraction uncertainty", "lp");	
	legend3.AddEntry(totalRelUncCosTheta, "total uncertainty", "lp");
	
	legend3.DrawClone();

	gPad->Update();

	// calculate chi2 / nDOF by hand for cross-check

	calculateChi2(standardCorrectedHist, PolarFunc, nCosThetaBins);	

	gSystem->mkdir("DistributionFits/1D", kTRUE);
	canvas2->SaveAs(Form("DistributionFits/1D/ROOTFIT_compareCorrectedCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
	errCanvas->SaveAs(Form("DistributionFits/1D/uncertainty%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
}
