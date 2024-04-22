#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "AccEffHelpers.h"

#include "../Tools/Style/Legends.h"

// transform the product of the acceptance and efficiency maps into a RooHistPdf

// see tutorial https://root.cern/doc/master/rf706__histpdf_8C.html

void prodAccEffPDF(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Int_t iState = 1) { //possible refFrame names: CS or HX
	writeExtraText = true;
	extraText = "       Internal";

	/// Set up the data
	//	using namespace RooFit;
	//	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	Int_t nCosThetaBins = 20;
	Float_t cosThetaMin = -1, cosThetaMax = 1;

	RooRealVar cosTheta(CosThetaVarName(refFrameName), CosThetaVarTitle(refFrameName), cosThetaMin, cosThetaMax);

	Int_t nPhiBins = 18;
	Int_t phiMin = -180, phiMax = 180;

	RooRealVar phi(PhiVarName(refFrameName), PhiVarTitle(refFrameName), phiMin, phiMax);

	/// 1. retrieve the 2D maps
	const char* mapName = CosThetaPhiTEfficiency2DName(ptMin, ptMax, refFrameName);

	// acceptance maps
	TFile* acceptanceFile = TFile::Open(Form("AcceptanceMaps/%dS/AcceptanceResults.root", iState), "READ");
	if (!acceptanceFile) {
		cout << "Acceptance file not found. Check the directory of the file." << endl;
		return;
	}

	auto* accMap = (TEfficiency*)acceptanceFile->Get(mapName);

	// efficiency maps
	TFile* efficiencyFile = TFile::Open(Form("EfficiencyMaps/%dS/EfficiencyResults.root", iState), "READ");
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

	RooHistPdf effPDF("effPDF", "", {cosTheta, phi}, effDataHist, 3);

	/// Draw the distributions
	gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPalette(gPreferredColorPalette);
	gStyle->SetNumberContours(256);

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
	legend.SetHeader(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));

	CMS_lumi(canvas, Form("Unpolarized #varUpsilon(%dS) MC", iState));

	gSystem->mkdir(Form("EfficiencyMaps/%dS", iState), kTRUE);
	canvas->SaveAs(Form("EfficiencyMaps/%dS/AccEff_CosThetaPhi%s_cent%dto%d_pt%dto%dGeV.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");

	// draw the sampling of the PDF in 3D
	TCanvas* canvas2 = new TCanvas("canvas2", "canvas", 700, 600);
	TH1* histoPDF = effPDF.createHistogram("histoPDF", cosTheta, RooFit::Binning(nCosThetaBins, cosThetaMin, cosThetaMax), RooFit::YVar(phi, RooFit::Binning(nPhiBins, phiMin, phiMax)));

	histoPDF->SetTitle(" ");

	histoPDF->GetZaxis()->SetTitle("acceptance #times efficiency PDF value");
	histoPDF->GetXaxis()->CenterTitle();

	histoPDF->GetYaxis()->CenterTitle();

	histoPDF->GetZaxis()->SetMaxDigits(3);

	histoPDF->Draw("COLZ");

	CMS_lumi(canvas2, Form("Unpolarized #varUpsilon(%dS) MC", iState));

	legend.DrawLatexNDC(.5, .75, Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));

	canvas2->SaveAs(Form("EfficiencyMaps/%dS/AccEffPDF_CosThetaPhi%s_cent%dto%d_pt%dto%dGeV.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
}
