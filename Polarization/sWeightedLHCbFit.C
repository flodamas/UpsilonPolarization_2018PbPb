#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../sPlot/SPlotHelpers.h"

#include "../Tools/Style/Legends.h"
#include "../Tools/Style/Figures.h"

#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.h"

#include "../MonteCarlo/prodAccEffPDF.C"

void drawAndSaveDistribution(TH2* histo, const char* name, Int_t ptMin, Int_t ptMax, const char* legend) {
	TCanvas* canvas = new TCanvas("canvas", "canvas", 700, 600);

	histo->Draw("COLZ");

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.045);
	text.SetTextColor(kWhite);

	text.DrawLatexNDC(.48, .8, Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	text.DrawLatexNDC(.48, .72, legend);

	gPad->Update();

	CMS_lumi(canvas, gCMSLumiText);
	gSystem->mkdir("2D", kTRUE);
	canvas->SaveAs(Form("2D/%s.png", name), "RECREATE");

	delete canvas;
}

// extraction of the polarization parameters based on the LHCb method https://arxiv.org/abs/1709.01301

// MLL fit of the angular distibution of sWeighted signal events

void sWeightedLHCbFit(bool updatesWeights = false, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Float_t phiMin = -180, Float_t phiMax = 180, Int_t iState = 1) {
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

	/// Set up the likelihood function

	// 1. the polarization POIs and the PDF
	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1., 1.);
	RooRealVar lambdaPhi("lambdaPhi", "lambdaPhi", -1., 1.);
	RooRealVar lambdaThetaPhi("lambdaThetaPhi", "lambdaThetaPhi", -1., 1.);

	RooFormulaVar lambdaTilde("lambdaTilde", "(@0 + 3*@1)/ (1-@1)", {lambdaTheta, lambdaPhi}); // invariant

	auto polarizationPDF = GeneralPolarizationPDF("polarizationPDF", " ", cosTheta, phi, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	auto dataNLL = polarizationPDF.createNLL(sWeightedData, Range("PolaFitRange"));

	//dataNLL->Print("v");

	// 2. the sWeights scale factor (see LHCb measurement)
	// inflates the uncertainties to account for the statistical fluctuations in the background subtraction

	RooConstVar sPlotScaleFactor = GetSPlotScaleFactor(sData, iState);

	// 2. the acceptance x efficiency PDF

	//auto* accEffPDF = prodAccEffPDF(ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, iState);

	auto* accEffPDF = prodAccEffPDF(cosTheta, phi, ptMin, ptMax, refFrameName, iState);

	// 3. the normalization factor

	auto* productPDF = new RooProdPdf("productPDF", "acc x eff x polarization PDF", polarizationPDF, *accEffPDF);

	auto* normFactor = productPDF->createIntegral(RooArgSet(cosTheta, phi), NormSet(RooArgSet(cosTheta, phi)), Range("PolaFitRange"));

	//normFactor->Print("v");

	// 4. bind all components into a single variable

	RooGenericPdf totalPDF("totalPDF", " total polarization PDF", "(@0/@1)*@2", {*productPDF, *normFactor, sPlotScaleFactor});

	RooFormulaVar totalNLL("totalNLL", "totalNLL", "(@0 + TMath::Log(@1))*@2", {*dataNLL, *normFactor, sPlotScaleFactor});

	//auto* totalNLL = totalPDF.createNLL(sWeightedData, Range("PolaFitRange"));

	//cout << "\n------------ Total NLL info ------------\n";

	//totalNLL.Print("v");

	//	RooFormulaVar totalNLL("totalNLL", "totalNLL", "@0", {dataNLL});

	/// Fit the angular distibution

	cout << "\n[Polarization] fitting the angular distribution, hold on...\n";

	// 1. minimizer settings

	RooMinimizer minimizer(totalNLL);

	minimizer.setStrategy(2); // better convergence
	minimizer.setPrintLevel(-1);
	minimizer.setRecoverFromNaNStrength(1.);
	minimizer.optimizeConst(true);
	minimizer.setOffsetting(true);
	//minimizer.setVerbose(true);

	// 2. perform the minimization
	minimizer.migrad();
	minimizer.hesse();

	auto* polarizationFitResult = minimizer.save();

	polarizationFitResult->Print("v");

	cout << "Frame-invariant polarization parameter (\"lambda tilde\") = " << lambdaTilde.getVal() << " +/- " << lambdaTilde.getPropagatedError(*polarizationFitResult) << endl;

	cout << "\nSECOND FIT METHOD: directly call fitTo() to the total PDF to enable AsymptoticError\n";

	auto* testResult = totalPDF.fitTo(sWeightedData, Save(), Range("PolaFitRange"), Extended(true), AsymptoticError(true), NumCPU(NCPUs), PrintLevel(-1), RecoverFromUndefinedRegions(1.), Offset());

	testResult->Print("v");

	cout << "Frame-invariant polarization parameter (\"lambda tilde\") = " << lambdaTilde.getVal() << " +/- " << lambdaTilde.getPropagatedError(*polarizationFitResult) << endl;

	/*
	/// Draw the (cos theta, phi) distributions with and without sWeights
	gStyle->SetPadLeftMargin(.14);
	gStyle->SetTitleYOffset(1.);
	gStyle->SetPadRightMargin(0.17);
	gStyle->SetTitleOffset(1., "z");
	SetColorPalette("TamDragon");

	const char* histoName = Form("rawCosThetaPhi%s_cent%dto%d_pt%dto%dGeV", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);

	// Y(1S)

	const char* name1S = Form("%s_1S", histoName);

	auto histo1S = TH2fromRooDataSet(data_weight1S, name1S, cosTheta, nCosThetaBins, cosThetaMin, cosThetaMax, phi, nPhiBins, phiMin, phiMax);

	drawAndSaveDistribution(histo1S, name1S, ptMin, ptMax, "#varUpsilon(1S) sWeights");
	*/
}
