#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"

#include "../Tools/Style/Legends.h"
#include "../Tools/Style/Figures.h"

#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "PolarFitHelpers.h"

/*

Test the LHCb's extraction method of the polarization parameters (https://arxiv.org/abs/1709.01301) with Monte Carlo, using the gen dimuon events as signal

Expectation: zero polarization

*/

RooDataSet ReducedGenMCDataset(RooDataSet* allDataset, RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS") {
	const char* kinematicCut = Form("(rapidity > %f && rapidity < %f) && (pt > %d && pt < %d) && (eventCat == 1)", gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var(Form("cosTheta%s", refFrameName))), *(wspace.var(Form("phi%s", refFrameName)))), kinematicCut);

	wspace.import(reducedDataset, RooFit::Rename(Form("GenMCDataset%s", refFrameName)));

	return reducedDataset;
}

void checkClosureGenMC(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Float_t phiMin = -180, Float_t phiMax = 180, Int_t iState = 1) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 20, nPhiBins = 18;

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = Form("../Files/Y%dSGeneratedMCDataset.root", iState);

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	RooDataSet* allDataset = (RooDataSet*)f->Get("MCdataset");

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	auto data = ReducedGenMCDataset(allDataset, wspace, ptMin, ptMax, refFrameName);

	// read variables in the reduced dataset in the workspace

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));

	RooRealVar phi = *wspace.var(PhiVarName(refFrameName));

	// draw

	//gStyle->SetPadLeftMargin(.15);
	//gStyle->SetTitleYOffset(1.2);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetTitleOffset(1.2, "z");
	SetColorPalette("TamDragon");

	TCanvas* canvas2D = new TCanvas("canvas2D", "canvas2D", 700, 600);

	TH1* histo2D = data.createHistogram("histo2D", cosTheta, Binning(nCosThetaBins, -1, 1), YVar(phi, Binning(nPhiBins, -180, 180)));

	histo2D->SetTitle(" ");
	histo2D->GetXaxis()->CenterTitle();
	//histo2D->GetYaxis()->SetRangeUser(-180, 300);
	histo2D->GetYaxis()->CenterTitle();
	histo2D->GetZaxis()->SetMaxDigits(3);
	histo2D->Draw("COLZ");

	// cosmetics
	/*
	TLegend legend(.18, .95, .45, .85);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("%d < p_{T}^{#mu#mu} < %d GeV/c", ptMin, ptMax));
	//legend.AddEntry(histoPDF, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	legend.DrawClone();
*/
	gPad->Update();

	CMS_lumi(canvas2D, Form("Unpolarized #varUpsilon(%dS) gen MC", iState));

	gSystem->mkdir("MCstudies", kTRUE);
	canvas2D->SaveAs(Form("MCstudies/DistribGenY%dSMC_%s_cent%dto%d_pt%dto%dGeV.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	// definition of the fit ranges (important!)
	cosTheta.setRange("PolaFitRange", cosThetaMin, cosThetaMax);
	phi.setRange("PolaFitRange", phiMin, phiMax);

	/// Set up the likelihood function

	// 1. the polarization POIs and the PDF
	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", 0.3, -1.1, 1.1);
	RooRealVar lambdaPhi("lambdaPhi", "lambdaPhi", -1.1, 1.1);
	RooRealVar lambdaThetaPhi("lambdaThetaPhi", "lambdaThetaPhi", -1., 1.);

	RooFormulaVar lambdaTilde("lambdaTilde", "(@0 + 3*@1)/ (1-@1)", {lambdaTheta, lambdaPhi}); // frame-invariant parameter

	auto polarizationPDF2D = GeneralPolarizationPDF("polarizationPDF2D", " ", cosTheta, phi, lambdaTheta, lambdaPhi, lambdaThetaPhi);
	polarizationPDF2D.setNormRange("PolaFitRange");

	//dataNLL->Print("v");

	// 3. the acceptance PDF

	const char* mapName = CosThetaPhiTEfficiency2DName(ptMin, ptMax, refFrameName);

	// acceptance maps
	TFile* acceptanceFile = TFile::Open(Form("../MonteCarlo/AcceptanceMaps/%dS/AcceptanceResults.root", iState), "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return;
	}

	auto* accMap = (TEfficiency*)acceptanceFile->Get(mapName);

	/// 2. do the product

	TH2* accTH2 = accMap->CreateHistogram();

	/// 3. transform into a RooDataHist, then into a RooHistPdf

	RooDataHist accDataHist2D("accDataHist2D", "", {cosTheta, phi}, accTH2);

	RooHistPdf accPDF2D("accPDF2D", "", {cosTheta, phi}, accDataHist2D, 3);
	accPDF2D.setNormRange("PolaFitRange");

	TCanvas* canvasAcc = new TCanvas("canvasAcc", " ", 700, 600);
	TH1* histoAccPDF2D = accPDF2D.createHistogram("histoAccPDF2D", cosTheta, RooFit::Binning(nCosThetaBins, -1, 1), RooFit::YVar(phi, RooFit::Binning(nPhiBins, -180, 180)));

	histoAccPDF2D->Draw("COLZ");

	RooProdPdf productPDF2D("productPDF2D", "acc x polarization PDF", polarizationPDF2D, accPDF2D);
	productPDF2D.setNormRange("PolaFitRange");

	/// Fit the angular distibution

	cout << "\n[Polarization] fitting the angular distribution, hold on...\n";

	// 1. minimizer settings

	auto nll2D = productPDF2D.createNLL(data, Range("PolaFitRange"));

	RooMinimizer minimizer(*nll2D);

	minimizer.setStrategy(2); // 1 for faster minimization, 2 for better convergence
	minimizer.setPrintLevel(0);
	minimizer.setRecoverFromNaNStrength(100.);
	minimizer.optimizeConst(true);
	minimizer.setOffsetting(true);
	//minimizer.setVerbose(true);

	// 2. perform the minimization
	minimizer.migrad();
	minimizer.hesse();

	auto* compactFitResult = minimizer.save();

	compactFitResult->Print("v");

	cout << "Frame-invariant parameter (\"lambda tilde\") = " << lambdaTilde.getVal() << " +/- " << lambdaTilde.getPropagatedError(*compactFitResult) << endl;

	//	SavePolarizationFitParameters(savedParams, "detailedViaMinimizer", fitModelName);

	cout << "\nSECOND FIT METHOD: directly call fitTo() to the total PDF\n";
	lambdaTheta.setVal(0.3);
	lambdaPhi.setVal(0);
	lambdaThetaPhi.setVal(0);

	auto* testResult2D = productPDF2D.fitTo(data, Save(), Range("PolaFitRange"), PrintLevel(0), Offset("bin"), RecoverFromUndefinedRegions(10.), NumCPU(3), Strategy(2));

	testResult2D->Print("v");

	// plot likelihood (and profile) scans

	RooPlot* frame1 = lambdaTheta.frame(Title(" "));
	productPDF2D.plotOn(frame1, ShiftToZero());

	std::unique_ptr<RooAbsReal> pll_lambdaTheta{productPDF2D.createProfile(lambdaTheta)};

	pll_lambdaTheta->plotOn(frame1, LineColor(kRed));
	//frame1->Draw();
}
