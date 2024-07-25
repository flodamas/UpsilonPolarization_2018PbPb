#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../sPlot/SPlotHelpers.h"

#include "../Tools/Style/Legends.h"
#include "../Tools/Style/Figures.h"

#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "PolarFitHelpers.h"

// extraction of the polarization parameters based on the LHCb method https://arxiv.org/abs/1709.01301

// MLL fit of the angular distibution of sWeighted signal events

void detailed(bool updatesWeights = false, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Float_t phiMin = -180, Float_t phiMax = 180, Int_t iState = 1) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 10, nPhiBins = 10;

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/UpsilonSkimmedDataset.root";

	RooWorkspace wspace = SetUpWorkspace(filename);

	auto data = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax);

	// read variables in the reduced dataset in the workspace

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));

	RooRealVar phi = *wspace.var(PhiVarName(refFrameName));

	// definition of the fit ranges (important!)
	cosTheta.setRange("PolaFitRange", cosThetaMin, cosThetaMax);
	phi.setRange("PolaFitRange", phiMin, phiMax);

	/// Invariant mass model

	// signal: one double-sided Crystal Ball PDF (symmetric Gaussian core) per Y resonance

	const char* signalShapeName = "SymDSCB";

	// background
	//int order = 2;
	//const char* bkgShapeName = Form("ChebychevOrder%d", order);
	const char* bkgShapeName = "ExpTimesErr";

	/// sWeighted data

	auto* sData = SWeightedDataset(wspace, ptMin, ptMax, signalShapeName, bkgShapeName, updatesWeights);

	RooDataSet sWeightedData = GetSpeciesSWeightedDataset(sData, Form("%dS", iState));

	//	RooRealVar sWeightVar = *(RooRealVar*)wspace.var(Form("yield%dS_sw", iState));

	/// Set up the likelihood function

	// 1. the polarization POIs and the PDF
	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1.2, 1.2);
	RooRealVar lambdaPhi("lambdaPhi", "lambdaPhi", -0.8, 0.8);
	RooRealVar lambdaThetaPhi("lambdaThetaPhi", "lambdaThetaPhi", -0.6, 0.6);

	RooFormulaVar lambdaTilde("lambdaTilde", "(@0 + 3*@1)/ (1-@1)", {lambdaTheta, lambdaPhi}); // frame-invariant parameter

	auto polarizationPDF = GeneralPolarizationPDF("polarizationPDF", " ", cosTheta, phi, lambdaTheta, lambdaPhi, lambdaThetaPhi);
	polarizationPDF.setNormRange("PolaFitRange");

	// 2. the sWeights scale factor (see LHCb measurement)
	// inflates the uncertainties to account for the statistical fluctuations in the background subtraction

	RooConstVar sPlotScaleFactor = GetSPlotScaleFactor(sData, iState);

	// 3. the acceptance x efficiency PDF

	const char* mapName = CosThetaPhiTEfficiency2DName(ptMin, ptMax, refFrameName);

	// acceptance maps
	TFile* acceptanceFile = TFile::Open(Form("../MonteCarlo/AcceptanceMaps/%dS/AcceptanceResults.root", iState), "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return;
	}

	auto* accMap = (TEfficiency*)acceptanceFile->Get(mapName);

	// efficiency maps
	TFile* efficiencyFile = TFile::Open(Form("../MonteCarlo/EfficiencyMaps/%dS/EfficiencyResults.root", iState), "READ");
	if (!efficiencyFile) {
		cout << "Efficiency file not found. Check the directory of the file." << endl;
		return;
	}

	auto* effMap = (TEfficiency*)efficiencyFile->Get(mapName);

	/// 2. do the product

	TH2* accTH2 = accMap->CreateHistogram();
	TH2* effTH2 = effMap->CreateHistogram();

	effTH2->Multiply(accTH2);

	/// 3. transform into a RooDataHist, then into a RooHistPdf
	RooDataHist effDataHist("effDataHist", "", {cosTheta, phi}, effTH2);

	RooHistPdf* accEffPDF = new RooHistPdf("effPDF", "", {cosTheta, phi}, effDataHist, 3);
	accEffPDF->setNormRange("PolaFitRange");

	// 3. the normalization factor

	auto* productPDF = new RooProdPdf("productPDF", "acc x eff x polarization PDF", polarizationPDF, *accEffPDF);
	productPDF->setNormRange("PolaFitRange");

	//auto* normFactor = productPDF->createIntegral(RooArgSet(cosTheta, phi), NormSet(RooArgSet(cosTheta, phi)), Range("PolaFitRange"));

	// 4. bind all components into a single variable

	RooGenericPdf totalPDF("totalPDF", " total polarization PDF", "@0*@1", {*productPDF, sPlotScaleFactor});
	totalPDF.setNormRange("PolaFitRange");

	auto totalNLL = totalPDF.createNLL(sWeightedData, Range("PolaFitRange"));

	/// Fit the angular distibution

	cout << "\n[Polarization] fitting the angular distribution, hold on...\n";

	// 1. minimizer settings

	RooMinimizer minimizer(*totalNLL);

	minimizer.setStrategy(2); // 1 for faster minimization, 2 for better convergence
	minimizer.setPrintLevel(0);
	minimizer.setRecoverFromNaNStrength(10.);
	minimizer.optimizeConst(true);
	minimizer.setOffsetting(true);
	//minimizer.setVerbose(true);

	// 2. perform the minimization
	minimizer.migrad();
	minimizer.hesse();

	auto* polarizationFitResult = minimizer.save();

	polarizationFitResult->Print("v");

	cout << "Frame-invariant parameter (\"lambda tilde\") = " << lambdaTilde.getVal() << " +/- " << lambdaTilde.getPropagatedError(*polarizationFitResult) << endl;

	RooRealVar savedLambdaTilde("lambdaTilde", "frame-invariant parameter", -5, 5); // can't directly dump the result for the lambdaTilde variable as RooFormulaVar, will write the analytical formula...
	savedLambdaTilde.setVal(lambdaTilde.getVal());
	savedLambdaTilde.setError(lambdaTilde.getPropagatedError(*polarizationFitResult));
	RooArgSet* savedParams = new RooArgSet(lambdaTheta, lambdaPhi, lambdaThetaPhi, savedLambdaTilde);

	const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	SavePolarizationFitParameters(savedParams, "detailedViaMinimizer", fitModelName);

	cout << "\nSECOND FIT METHOD: directly call fitTo() to the total PDF to enable AsymptoticError\n";

	auto* testResult = totalPDF.fitTo(sWeightedData, Save(), Range("PolaFitRange"), SumW2Error(true), PrintLevel(0), Offset("initial"), RecoverFromUndefinedRegions(10.), Strategy(1));

	testResult->Print("v");

	// plot likelihood scans
	auto canvasLambdaTheta = new TCanvas("canvasLambdaTheta", " ", 600, 600);
	RooPlot* frameLambdaTheta = lambdaTheta.frame(Title(" "));
	totalNLL->plotOn(frameLambdaTheta, ShiftToZero());
	frameLambdaTheta->Draw();

	auto canvasLambdaPhi = new TCanvas("canvasLambdaPhi", " ", 600, 600);

	RooPlot* frameLambdaPhi = lambdaPhi.frame(Title(" "));
	totalNLL->plotOn(frameLambdaPhi, ShiftToZero());
	frameLambdaPhi->Draw();
}
