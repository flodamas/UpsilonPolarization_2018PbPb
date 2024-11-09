#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"

#include "../Tools/Style/Legends.h"
#include "../Tools/Style/Figures.h"

#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.h"

#include "../MonteCarlo/AccEffHelpers.h"

/*

Test the LHCb's extraction method of the polarization parameters (https://arxiv.org/abs/1709.01301) with Monte Carlo, using the reconstructed dimuon events as signal

This allows to both validate the code implementation AND the acceptance x efficiency closure (1D and 2D)

Expectation: zero polarization

*/

void testFullMCclosure(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Float_t phiMin = -180, Float_t phiMax = 180, Int_t iState = 1) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	Int_t nCosThetaBins = 20, nPhiBins = 20;

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = Form("../Files/Y%dSSelectedMCWeightedDataset%s.root", iState, gMuonAccName);

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

	auto data = ReducedRecoMCDataset(wspace, ptMin, ptMax, refFrameName);

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

	TH1* histo2D = data.createHistogram("histo2D", cosTheta, Binning(nCosThetaBins, -1, 1), YVar(phi, Binning(nPhiBins, -200, 200)));

	histo2D->SetTitle(" ");
	histo2D->GetXaxis()->CenterTitle();
	histo2D->GetYaxis()->SetRangeUser(-180, 300);
	histo2D->GetYaxis()->CenterTitle();
	histo2D->GetZaxis()->SetMaxDigits(3);
	histo2D->Draw("COLZ");

	// cosmetics

	TLegend legend(.18, .95, .45, .85);
	legend.SetTextSize(.05);
	legend.SetHeader(Form("centrality %d-%d%%, %d < p_{T}^{#mu#mu} < %d GeV/c", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));
	//legend.AddEntry(histoPDF, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaTheta.getVal(), lambdaTheta.getError()), "l");

	legend.DrawClone();

	gPad->Update();

	CMS_lumi(canvas2D, Form("#varUpsilon(%dS) Hydjet-embedded MC", iState));

	gSystem->mkdir("MCstudies", kTRUE);
	canvas2D->SaveAs(Form("MCstudies/DistribRecoY%dSMC_%s_cent%dto%d_pt%dto%dGeV.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	// definition of the fit ranges (important!)
	cosTheta.setRange("PolaFitRange", cosThetaMin, cosThetaMax);
	phi.setRange("PolaFitRange", phiMin, phiMax);

	/// Set up the likelihood function

	// 1. the polarization POIs and the PDF
	RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -1., 1.);
	RooRealVar lambdaPhi("lambdaPhi", "lambdaPhi", -1., 1.);
	RooRealVar lambdaThetaPhi("lambdaThetaPhi", "lambdaThetaPhi", -0.5, 0.5);

	RooFormulaVar lambdaTilde("lambdaTilde", "(@0 + 3*@1)/ (1-@1)", {lambdaTheta, lambdaPhi}); // frame-invariant parameter

	auto polarizationPDF2D = GeneralPolarizationPDF("polarizationPDF2D", " ", cosTheta, phi, lambdaTheta, lambdaPhi, lambdaThetaPhi);
	polarizationPDF2D.setNormRange("PolaFitRange");

	//auto dataNLL = polarizationPDF2D.createNLL(data, Range("PolaFitRange"));

	//dataNLL->Print("v");

	// 3. the acceptance x efficiency PDF

	const char* mapName = CosThetaPhiTEfficiency2DName(ptMin, ptMax, refFrameName);

	// acceptance maps
	TFile* acceptanceFile = TFile::Open(Form("../MonteCarlo/AcceptanceMaps/%dS/AcceptanceResults%s.root", iState, gMuonAccName), "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return;
	}

	auto* accMap = (TEfficiency*)acceptanceFile->Get(mapName);

	// efficiency maps
	TFile* efficiencyFile = TFile::Open(Form("../MonteCarlo/EfficiencyMaps/%dS/EfficiencyResults%s.root", iState, gMuonAccName), "READ");
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

	RooDataHist effDataHist2D("effDataHist2D", "", {cosTheta, phi}, effTH2);

	RooHistPdf accEffPDF2D("effPDF2D", "", {cosTheta, phi}, effDataHist2D, 3);
	accEffPDF2D.setNormRange("PolaFitRange");

	// 3. the normalization factor

	RooProdPdf productPDF2D("productPDF2D", "acc x eff x polarization PDF", polarizationPDF2D, accEffPDF2D);
	productPDF2D.setNormRange("PolaFitRange");
	//productPDF2D.forceNumInt(true);
	//RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-5);
	//RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-5);

	/// Fit the angular distibution

	cout << "\n[Polarization] fitting the angular distribution, hold on...\n";

	// 1. minimizer settings

	auto nll2D = productPDF2D.createNLL(data, Range("PolaFitRange"));

	RooMinimizer minimizer(*nll2D);

	minimizer.setStrategy(2); // 1 for faster minimization, 2 for better convergence
	minimizer.setPrintLevel(0);
	minimizer.setRecoverFromNaNStrength(1000.);
	minimizer.optimizeConst(true);
	minimizer.setOffsetting(true);
	//minimizer.setVerbose(true);

	// 2. perform the minimization
	minimizer.migrad();
	minimizer.hesse();

	auto* compactFitResult = minimizer.save();

	compactFitResult->Print("v");

	cout << "Frame-invariant parameter (\"lambda tilde\") = " << lambdaTilde.getVal() << " +/- " << lambdaTilde.getPropagatedError(*compactFitResult) << endl;

	cout << "\nSECOND FIT METHOD: directly call fitTo() to the total PDF to enable AsymptoticError\n";
	lambdaTheta.setVal(0.2);
	lambdaPhi.setVal(0);
	lambdaThetaPhi.setVal(0);

	auto* testResult2D = productPDF2D.fitTo(data, Save(), Range("PolaFitRange"), AsymptoticError(true), PrintLevel(0), RecoverFromUndefinedRegions(10.), Offset("bin"), Optimize(true), Strategy(2));

	testResult2D->Print("v");
}
