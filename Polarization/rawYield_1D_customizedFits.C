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
// #include "../Tools/RooFitPDFs/cosThetaPolarFunc.h"
#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"

void rawYield_1D_customizedFits(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 10, Int_t phiMin = -180, Int_t phiMax = 180, const Int_t nPhiBins = 6, Int_t cosThetaMin = -1, Int_t cosThetaMax = 1, Int_t iState = gUpsilonState) {
	
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Bin edges and width 
	// Set the bin edges along the cosTheta/phi axis depending on the number of bins 
	// (The bin edges are pre-defined, so need to modify them if different bin edges are required)
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);
	
	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	/// Set up the variables
	RooRealVar cosTheta("cosTheta", "", cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	cosTheta.setRange(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	cosTheta.setBins(nCosThetaBins);

	/// Assign signal and background shape name to read the file for the yield extraction results
	const char* signalShapeName = "SymDSCB";

	// background shape array: ChebychevOrderN or ExpTimesErr	
	// const char* bkgShapeNamesCosTheta[] = {
	//   // "ChebychevOrder1",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",

	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",

	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",

	//   "ChebychevOrder2"
	//   // "ChebychevOrder1"
	// };
	// const char* bkgShapeNamesPhi[] = {
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2",
	//   "ChebychevOrder2"
	// };
	const char* bkgShapeNamesCosTheta[] = {
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
	  "ExpTimesErr",
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

	const char* bkgShapeNamesPhi[] = {
		"ExpTimesErr",
		"ExpTimesErr",
		"ExpTimesErr",
		"ExpTimesErr",
		"ExpTimesErr",
		"ExpTimesErr",
	};

	/// "Standard" procedure: extract the yields per bin
	TH1D* standardCorrectedCosThetaHist = new TH1D("standardCorrectedCosThetaHist", " ", nCosThetaBins, cosThetaBinEdges.data());

	TH1D* standardCorrectedPhiHist = new TH1D("standardCorrectedPhiHist", " ", nPhiBins, phiBinEdges.data());

	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get acceptance maps
	TFile* acceptanceFile = openFile(Form("%s/AcceptanceResults%s.root", AcceptanceResultsPath(gMuonAccName.Data()), "_fullPhi"));
	
	if (!acceptanceFile) {
		std::cerr << "Error: acceptanceFile is null." << std::endl;
		acceptanceFile->Close();
		delete acceptanceFile;

		// acceptanceMap_noGenFilter(0, 30, gUpsilonState, lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded, "MuonUpsilonTriggerAcc");

		return;
	}
	
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);

	// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accHistCosTheta = rebinTEff3DMapCosTheta(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);

	TEfficiency* accHistPhi = rebinTEff3DMapPhi(accMap, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	// get efficiency maps
	TFile* efficiencyFile = openFile(Form("%s/EfficiencyResults%s.root", EfficiencyResultsPath(gMuonAccName.Data()), "_fullPhi"));
	
	if (!efficiencyFile) {
		std::cerr << "Error: efficiencyFile is null." << std::endl;

		return;
	}	
	
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);

	// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effHistCosTheta = rebinTEff3DMapCosTheta(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);

	TEfficiency* effHistPhi = rebinTEff3DMapPhi(effMap, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	// get relative systematic uncertainty of efficiency
	auto* systEff = (TH3D*)efficiencyFile->Get(SystTEfficiency3DName(refFrameName));

	// rebin uncertainty map based on costheta, phi, and pT selection
	TH1D* systEffCosTheta = rebinRel3DUncCosTheta(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);

	TH1D* systEffPhi = rebinRel3DUncPhi(effMap, systEff, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	Bool_t isCSframe = (strcmp(refFrameName, "CS") == 0) ? kTRUE : kFALSE;


	// // draw acc and eff histograms to check if the rebinning works well
	DrawEfficiency1DHist(accHistCosTheta, ptMin, ptMax, iState, kTRUE, kTRUE);
	DrawEfficiency1DHist(effHistCosTheta, ptMin, ptMax, iState, kFALSE, kTRUE);

	DrawEfficiency1DHist(accHistPhi, ptMin, ptMax, iState, kTRUE, kFALSE);
	DrawEfficiency1DHist(effHistPhi, ptMin, ptMax, iState, kFALSE, kFALSE);

	// // define histograms to draw uncertainty plots
	// TH1D* statHighEffCosTheta = new TH1D("statHighEffCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	// TH1D* statLowEffCosTheta = new TH1D("statLowEffCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	// TH1D* statHighAccCosTheta = new TH1D("statHighAccCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	// TH1D* statLowAccCosTheta = new TH1D("statLowAccCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	// TH1D* yield1SUncCosTheta = new TH1D("yield1SUncCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());
	// TH1D* totalRelUncCosTheta = new TH1D("totalRelUncCosTheta", "", nCosThetaBins, cosThetaBinEdges.data());

	// TH1D* statHighEffPhi = new TH1D("statHighEffPhi", "", nPhiBins, phiBinEdges.data());
	// TH1D* statLowEffPhi = new TH1D("statLowEffPhi", "", nPhiBins, phiBinEdges.data());
	// TH1D* statHighAccPhi = new TH1D("statHighAccPhi", "", nPhiBins, phiBinEdges.data());
	// TH1D* statLowAccPhi = new TH1D("statLowAccPhi", "", nPhiBins, phiBinEdges.data());
	// TH1D* yield1SUncPhi = new TH1D("yield1SUncPhi", "", nPhiBins, phiBinEdges.data());
	// TH1D* totalRelUncPhi = new TH1D("totalRelUncPhi", "", nPhiBins, phiBinEdges.data());

	const char** fitModelNamesCosTheta = GetFitModelNames(signalShapeName, ptMin, ptMax, isCSframe, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);

	const char** fitModelNamesPhi = GetFitModelNames(signalShapeName, ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	/// apply weights and errors to each costheta bin

	Float_t maxYieldCosTheta = correctRawYield1DHist(standardCorrectedCosThetaHist, accHistCosTheta, effHistCosTheta, systEffCosTheta, nCosThetaBins, bkgShapeNamesCosTheta, fitModelNamesCosTheta);
	
	Float_t maxYieldPhi = correctRawYield1DHist(standardCorrectedPhiHist, accHistPhi, effHistPhi, systEffPhi, nPhiBins, bkgShapeNamesPhi, fitModelNamesPhi);

	/// Polarization fit
	// with Root Fit function

	TVirtualFitter::SetDefaultFitter("Minuit"); 

	TCanvas* polarCanvas = new TCanvas(standardCorrectedCosThetaHist->GetName(), "", 1200, 600);

	polarCanvas->Divide(2, 1);

	polarCanvas->cd(1);

	TF1* PolarFunc = getCosThetaPolarFunc(maxYieldCosTheta);

	TFitResultPtr fitResults = standardCorrectedCosThetaHist->Fit("PolarFunc", "ESV", "", cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]); //L:log likelihood fit (default: chi2 method), E: NINOS

	// Fit results

	double chi2 = fitResults->Chi2();
	double nDOF = nCosThetaBins - PolarFunc->GetNpar();

	double normVal = fitResults->Parameter(0);
	double normErr = fitResults->ParError(0);

	double lambdaVal = fitResults->Parameter(1);
	double lambdaErr = fitResults->ParError(1);

	// 1 Sigma band
    TH1D* errorBand = new TH1D("errorBand", "", 1000, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
   	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(errorBand, 0.68);

	gStyle->SetOptFit(1011);

	// cosmetics

	errorBand->SetFillColorAlpha(kRed-9, 0.2);
	errorBand->SetMarkerSize(0);

	errorBand->Draw("E3 SAME");

	standardCorrectedCosThetaHist->SetMarkerStyle(20);
	standardCorrectedCosThetaHist->SetMarkerSize(1);
	standardCorrectedCosThetaHist->SetMarkerColor(kAzure + 2);

	standardCorrectedCosThetaHist->GetYaxis()->SetRangeUser(0, 2 * maxYieldCosTheta);
	standardCorrectedCosThetaHist->GetYaxis()->SetMaxDigits(3);

	standardCorrectedCosThetaHist->SetXTitle(Form("cos #theta_{%s}", refFrameName));
	standardCorrectedCosThetaHist->SetYTitle(Form("Events / ( %0.1f )", cosThetaStep));

	standardCorrectedCosThetaHist->GetXaxis()->CenterTitle();
	standardCorrectedCosThetaHist->GetYaxis()->CenterTitle();

	TLegend legendCosTheta(.2, .91, .48, .67);
	legendCosTheta.SetTextSize(.05);
	legendCosTheta.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legendCosTheta.AddEntry(standardCorrectedCosThetaHist, "#varUpsilon(1S) corrected yield", "lp");
	legendCosTheta.AddEntry(PolarFunc, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaVal, lambdaErr), "l");
	legendCosTheta.AddEntry((TObject*)0, Form("                       n  = %.2f #pm %.2f", normVal, normErr), "");

	legendCosTheta.DrawClone();

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.045);
	textChi2.DrawLatexNDC(0.74, 0.044, Form("#chi^{2} / n_{dof} = %.2f", chi2 / nDOF));

	gPad->Update();

	polarCanvas->cd(2);

	standardCorrectedPhiHist->Draw();

	standardCorrectedPhiHist->SetMarkerStyle(20);
	standardCorrectedPhiHist->SetMarkerSize(1);
	standardCorrectedPhiHist->SetMarkerColor(kAzure + 2);

	standardCorrectedPhiHist->GetYaxis()->SetRangeUser(0, 2 * maxYieldPhi);
	standardCorrectedPhiHist->GetYaxis()->SetMaxDigits(3);

	standardCorrectedPhiHist->SetXTitle(Form("#varphi_{%s}", refFrameName));
	standardCorrectedPhiHist->SetYTitle(Form("Events / ( %0.0f )", phiStep));

	standardCorrectedPhiHist->GetXaxis()->CenterTitle();
	standardCorrectedPhiHist->GetYaxis()->CenterTitle();

	TLegend* legendPhi = new TLegend(.2, .91, .48, .67);
	legendPhi->SetTextSize(.05);
	legendPhi->SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	legendPhi->AddEntry(standardCorrectedPhiHist, "#varUpsilon(1S) corrected yield", "lp");
	// legendPhi.AddEntry(PolarFunc, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaVal, lambdaErr), "l");
	// legendPhi.AddEntry((TObject*)0, Form("                       n  = %.2f #pm %.2f", normVal, normErr), "");

	legendPhi->Draw();

	polarCanvas->Update();

	// // draw uncertainties
	
	// TCanvas* errCanvas = drawUncertaintyPlot1D(refFrameName, systEffCosTheta, statHighEffCosTheta, statLowEffCosTheta, statHighAccCosTheta, statLowAccCosTheta, yield1SUncCosTheta, totalRelUncCosTheta);

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
	
	// calculateChi2(standardCorrectedCosThetaHist, PolarFunc, nCosThetaBins);	
	
	/// contour plot
	// (ref: https://root-forum.cern.ch/t/roofit-minos-errors-for-2-parameters-of-interest/16157)
	
	// set the confidence level
   	gMinuit->SetErrorDef(2.30); // 1 sigma corresponds to delchi2 = 2.30 
   	TGraph* contourPlot1 = (TGraph*)gMinuit->Contour(1000, 1, 0); // Contour(number of points, lambda_theta, normalization factor)

   	gMinuit->SetErrorDef(6.18); // 2 sigma corresponds to delchi2 = 6.18
   	TGraph* contourPlot2 = (TGraph*)gMinuit->Contour(1000, 1, 0);	

   	TCanvas* contourCanvas = drawContourPlots(ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], refFrameName, contourPlot1, contourPlot2);

	// save canvas

	gSystem->mkdir("DistributionFits/1D", kTRUE);
	polarCanvas->SaveAs(Form("DistributionFits/1D/ROOTFIT_compareCorrectedCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
	
	// gSystem->mkdir("UncertaintyPlots/1D", kTRUE);
	// errCanvas->SaveAs(Form("UncertaintyPlots/1D/uncertainty%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");

	gSystem->mkdir("ContourPlots/1D", kTRUE);
	contourCanvas->SaveAs(Form("ContourPlots/1D/contour%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, phiMin, phiMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]), "RECREATE");
}
