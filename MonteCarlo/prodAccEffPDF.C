#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "AccEffHelpers.h"

#include "../Tools/Style/Legends.h"

// transform the product of the acceptance and efficiency maps into a RooHistPdf

// see tutorial https://root.cern/doc/master/rf706__histpdf_8C.html

//RooHistPdf* prodAccEffPDF(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1., Float_t phiMin = -180, Float_t phiMax = 180, Int_t iState = gUpsilonState) { //possible refFrame names: CS or HX

RooHistPdf* prodAccEffPDF(RooRealVar cosTheta, RooRealVar phi, Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t iState = gUpsilonState) { //possible refFrame names: CS or HX

	writeExtraText = true;
	extraText = "       Internal";

	/// Set up the data
	//	using namespace RooFit;
	//	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	Int_t nCosThetaBins = 20;

	//	RooRealVar cosTheta(CosThetaVarName(refFrameName), CosThetaVarTitle(refFrameName), cosThetaMin, cosThetaMax);

	Int_t nPhiBins = 18;

	//	RooRealVar phi(PhiVarName(refFrameName), PhiVarTitle(refFrameName), phiMin, phiMax, gPhiUnit);

	/// 1. retrieve the 2D maps
	const char* mapName = CosThetaPhiTEfficiency2DName(refFrameName, ptMin, ptMax);

	// acceptance maps
	TFile* acceptanceFile = TFile::Open(Form("../MonteCarlo/AcceptanceMaps/%dS/AcceptanceResults%s.root", iState, gMuonAccName), "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return nullptr;
	}

	auto* accMap = (TEfficiency*)acceptanceFile->Get(mapName);

	// efficiency maps
	TFile* efficiencyFile = TFile::Open(Form("../MonteCarlo/EfficiencyMaps/%dS/EfficiencyResults%s.root", iState, gMuonAccName), "READ");
	if (!efficiencyFile) {
		cout << "Efficiency file not found. Check the directory of the file." << endl;
		return nullptr;
	}

	auto* effMap = (TEfficiency*)efficiencyFile->Get(mapName);

	/// 2. do the product

	TH2* accTH2 = accMap->CreateHistogram();
	TH2* effTH2 = effMap->CreateHistogram();

	effTH2->Multiply(accTH2);

	/// 3. transform into a RooDataHist, then into a RooHistPdf
	RooDataHist effDataHist("effDataHist", "", {cosTheta, phi}, effTH2);

	RooHistPdf* effPDF = new RooHistPdf("effPDF", "", {cosTheta, phi}, effDataHist, 3);

	/// Draw the distributions
	//gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(0.18);
	SetColorPalette(gPreferredColorPaletteName);

	// the acc x eff map
	TCanvas* canvas = new TCanvas("canvas", "canvas", 700, 600);
	effTH2->GetZaxis()->SetTitle(Form("#varUpsilon(%dS) acceptance x efficiency", iState));

	//effTH2->GetYaxis()->SetRangeUser(-180, 240);
	effTH2->GetXaxis()->CenterTitle();

	effTH2->GetYaxis()->CenterTitle();

	//effTH2->GetZaxis()->SetRangeUser(0, 0.8);
	effTH2->SetMinimum(0);

	effTH2->Draw("COLZ");

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.045);
	legend.DrawLatexNDC(.48, .8, Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	legend.DrawLatexNDC(.48, .72, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));

	CMS_lumi(canvas, Form("Unpolarized #varUpsilon(%dS) MC", iState));

	gSystem->mkdir(Form("../MonteCarlo/EfficiencyMaps/%dS", iState), kTRUE);
	canvas->SaveAs(Form("../MonteCarlo/EfficiencyMaps/%dS/AccEff_CosThetaPhi%s_cent%dto%d_pt%dto%dGeV.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	// draw the sampling of the PDF in 3D
	TCanvas* canvas2 = new TCanvas("canvas2", "canvas", 700, 600);
	TH1* histoPDF = effPDF->createHistogram("histoPDF", cosTheta, RooFit::Binning(nCosThetaBins, -1, 1), RooFit::YVar(phi, RooFit::Binning(nPhiBins, -180, 180)));

	histoPDF->SetTitle(" ");

	histoPDF->GetZaxis()->SetTitle(Form("#varUpsilon(%dS) acceptance #times efficiency PDF", iState));
	histoPDF->GetXaxis()->CenterTitle();

	histoPDF->GetYaxis()->CenterTitle();

	histoPDF->GetZaxis()->SetMaxDigits(3);

	histoPDF->Draw("COLZ");

	CMS_lumi(canvas2, Form("Unpolarized #varUpsilon(%dS) MC", iState));

	legend.DrawLatexNDC(.48, .8, Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	legend.DrawLatexNDC(.48, .72, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));

	canvas2->SaveAs(Form("../MonteCarlo/EfficiencyMaps/%dS/AccEffPDF_CosThetaPhi%s_cent%dto%d_pt%dto%dGeV.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	acceptanceFile->Close();
	efficiencyFile->Close();

	delete canvas;
	delete canvas2;

	return effPDF;
}
