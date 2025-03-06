#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "AccEffHelpers.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
// #include "../sPlot/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Polarization/PolarFitHelpers.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"

std::vector<TEfficiency*> accEffplots_3Dto1D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 5, Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = gUpsilonState, Bool_t isPhiFolded = kFALSE, TString accName = "MuonUpsilonTriggerAcc") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc"
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace std;
	using namespace RooFit;
	using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Bin edges and width
	// Set the bin edges along the cosTheta/phi axis depending on the number of bins, min and max values
	// (If want to use non-uniform bin width, the bin edges should be pre-defined in PolarFitHelpers.h)
	std::vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	std::vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get acceptance maps
	Double_t lambdaTheta = 0., lambdaPhi = 0., lambdaThetaPhi = 0.;	

	TString MuonAccName = "";
	TString accFileName = "";

	if (accName == TString("MuonUpsilonTriggerAcc")) MuonAccName = "_TriggerAcc";
	else if (accName == TString("MuonWithin2018PbPbAcc")) MuonAccName = "_2018PbPbAcc";
	else if (accName == TString("MuonSimpleAcc")) MuonAccName = "_SimpleAcc";
	else {
		cout << "Invalid acceptance name. Please choose from 'MuonUpsilonTriggerAcc', 'MuonWithin2018PbPbAcc', or 'MuonSimpleAcc'." << endl;
		return {};
	}

	if (isPhiFolded == kTRUE) accFileName = Form("./AcceptanceMaps/1S/AcceptanceResults%s.root", MuonAccName.Data());
	else accFileName = Form("./AcceptanceMaps/1S/AcceptanceResults%s_fullPhi.root", MuonAccName.Data());

	TFile* acceptanceFile = openFile(accFileName);
	cout << accFileName << " opened" << endl;
	
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);
	// cout << "Nominal map name: " << nominalMapName << endl;

	if (!accMap) {
    	std::cerr << "Error: accMap is null." << std::endl;
    	return {};
	}

	// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosTheta = rebinTEff3DMapCosTheta(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);
	TEfficiency* accMapPhi = rebinTEff3DMapPhi(accMap, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	// get efficiency maps
	// TFile* efficiencyFile = openFile(Form("./EfficiencyMaps/1S/EfficiencyResults%s.root", gMuonAccName));
	TString effFileName = "";

	if (isPhiFolded == kTRUE) effFileName = Form("./EfficiencyMaps/1S/EfficiencyResults%s.root", MuonAccName.Data());
	else effFileName = Form("./EfficiencyMaps/1S/EfficiencyResults%s_fullPhi.root", MuonAccName.Data());

	TFile* efficiencyFile = openFile(effFileName);
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);
	// cout << "Nominal map name: " << nominalMapName << endl;

	if (!effMap) {
    	std::cerr << "Error: effMap is null." << std::endl;
    	return {};
	}

	// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effMapCosTheta = rebinTEff3DMapCosTheta(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);
	TEfficiency* effMapPhi = rebinTEff3DMapPhi(effMap, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	/// draw acceptance 1D histograms
	DrawEfficiency1DHist(accMapCosTheta, ptMin, ptMax, iState, kTRUE, kTRUE, refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi, MuonAccName.Data(), isPhiFolded);
	DrawEfficiency1DHist(accMapPhi, ptMin, ptMax, iState, kTRUE, kFALSE, refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi, MuonAccName.Data(), isPhiFolded);

	/// draw efficiency 1D histograms
	DrawEfficiency1DHist(effMapCosTheta, ptMin, ptMax, iState, kFALSE, kTRUE, refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi, MuonAccName.Data(), isPhiFolded);
	DrawEfficiency1DHist(effMapPhi, ptMin, ptMax, iState, kFALSE, kFALSE, refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi, MuonAccName.Data(), isPhiFolded);

	/// draw acceptance x efficiency 1D histograms
	TEfficiency* phiHist = DrawEffxAcc1DHist(accMapCosTheta, effMapCosTheta, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, (double)phiMin, (double)phiMax, kTRUE, refFrameName, iState, MuonAccName.Data(), isPhiFolded);
	TEfficiency* cosThetaHist = DrawEffxAcc1DHist(accMapPhi, effMapPhi, ptMin, ptMax, nPhiBins, phiBinEdges, cosThetaMin, cosThetaMax, kFALSE, refFrameName, iState, MuonAccName.Data(), isPhiFolded);

	std::vector<TEfficiency*> hists = {phiHist, cosThetaHist};
	return hists;
}

void accEffplots_3Dto1D_comparison(Int_t ptMin = 2, Int_t ptMax = 6, const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 1, Int_t phiMin = 0, Int_t phiMax = 60, Int_t iState = gUpsilonState, Bool_t isPhiFolded = kFALSE, Bool_t isCosTheta = kTRUE, TString accName = "MuonUpsilonTriggerAcc") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc"
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Bin edges and width
	// Set the bin edges along the cosTheta/phi axis depending on the number of bins, min and max values
	// (If want to use non-uniform bin width, the bin edges should be pre-defined in PolarFitHelpers.h)
	std::vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	std::vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get acceptance maps
	Double_t lambdaTheta = 0., lambdaPhi = 0., lambdaThetaPhi = 0.;	

	TString MuonAccName = "";
	TString accFileName = "";

	if (accName == TString("MuonUpsilonTriggerAcc")) MuonAccName = "_TriggerAcc";
	else if (accName == TString("MuonWithin2018PbPbAcc")) MuonAccName = "_2018PbPbAcc";
	else if (accName == TString("MuonSimpleAcc")) MuonAccName = "_SimpleAcc";
	else {
		cout << "Invalid acceptance name. Please choose from 'MuonUpsilonTriggerAcc', 'MuonWithin2018PbPbAcc', or 'MuonSimpleAcc'." << endl;
		return;
	}

	cout << "Drawing acceptance x efficiency 1D histograms for |phi| range " << phiMin << " to " << phiMax << " degrees..." << endl;

	std::vector<TEfficiency*> positiveHists;
	std::vector<TEfficiency*> negativeHists;
	TEfficiency* positiveHist;
	TEfficiency* negativeHist;

	if (isCosTheta == kTRUE) {
		positiveHists = accEffplots_3Dto1D(ptMin, ptMax, refFrameName, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, iState, isPhiFolded, accName);
		negativeHists = accEffplots_3Dto1D(ptMin, ptMax, refFrameName, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, -phiMax, -phiMin, iState, isPhiFolded, accName);
		positiveHist = positiveHists[0];
		negativeHist = positiveHists[0];
	}
	else {
		positiveHists = accEffplots_3Dto1D(ptMin, ptMax, refFrameName, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, iState, isPhiFolded, accName);
		negativeHists = accEffplots_3Dto1D(ptMin, ptMax, refFrameName, nCosThetaBins, -cosThetaMax, -cosThetaMin, nPhiBins, phiMin, phiMax, iState, isPhiFolded, accName);
		positiveHist = positiveHists[1];
		negativeHist = negativeHists[1];
	}

	// empty frame for the axes
	TH1D* frameHist1D;
	
	if (isCosTheta == kTRUE) frameHist1D = new TH1D("frameHist1D", "", nCosThetaBins, cosThetaBinEdges.data());
	else frameHist1D = new TH1D("frameHist1D", "", nPhiBins, phiBinEdges.data());

	/// draw acceptance x efficiency 1D histograms for positive and negative phi
	TCanvas* comparisonCanvas = new TCanvas("comparisonCanvas", "comparisonCanvas", 600, 600);
	
	comparisonCanvas->SetRightMargin(0.05);

	frameHist1D->Draw();
	frameHist1D->GetYaxis()->SetRangeUser(0, 1);

	if (isCosTheta == kTRUE) {
		frameHist1D->GetXaxis()->SetTitle(CosThetaVarTitle(refFrameName));
		frameHist1D->GetXaxis()->SetNdivisions(-505);
	}
	else {
		frameHist1D->GetXaxis()->SetTitle(PhiVarTitle(refFrameName));
		frameHist1D->GetXaxis()->SetNdivisions(-6);
	}

	frameHist1D->GetXaxis()->CenterTitle();
	frameHist1D->GetYaxis()->SetTitle(TEfficiencyMainTitle(iState, "acceptance x efficiency"));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.55, .87, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));

	if (strcmp(MuonAccName.Data(), "_TriggerAcc") == 0) legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	else if (strcmp(MuonAccName.Data(), "_SimpleAcc") == 0) legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", iState));
	else if (strcmp(MuonAccName.Data(), "_2018PbPbAcc") == 0) legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", iState));

	positiveHist->Draw("EP same");
	positiveHist->SetMarkerColor(TColor::GetColor("#A559AA"));
	positiveHist->SetMarkerStyle(20);
	positiveHist->SetLineColor(TColor::GetColor("#A559AA"));
	
	negativeHist->Draw("PE same");
	negativeHist->SetMarkerColor(TColor::GetColor("#009ADE"));
	negativeHist->SetMarkerStyle(24);
	negativeHist->SetLineColor(TColor::GetColor("#009ADE"));

	TLegend* legend1D = new TLegend(0.17, 0.6, 0.4, 0.76);
	legend1D->SetBorderSize(0);
	legend1D->SetFillStyle(0);
	if (isCosTheta == kTRUE) {
		legend1D->AddEntry(positiveHist, Form("%d < #varphi < %d#circ", phiMin, phiMax), "lep");
		legend1D->AddEntry(negativeHist, Form("%d < #varphi < %d#circ", -phiMax, -phiMin), "lep");
	}
	else {
		legend1D->AddEntry(positiveHist, Form("%.2f < cos#theta < %.2f", cosThetaMin, cosThetaMax), "lep");
		legend1D->AddEntry(negativeHist, Form("%.2f < cos#theta < %.2f", -cosThetaMax, -cosThetaMin), "lep");
	}
	legend1D->Draw();
	
	CMS_lumi(comparisonCanvas, Form("#varUpsilon(%dS) Pythia 8 (5.02 TeV)", iState));

	gSystem->mkdir(Form("AccxEffMaps/%dS/analysisBin", iState), kTRUE);
	if (isCosTheta == kTRUE) comparisonCanvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/1DAccxEffComp_%s%s%s_pt%dto%d_|phi|%dto%d_%s.png", iState, positiveHist->GetName(), refFrameName, MuonAccName.Data(), ptMin, ptMax, phiMin, phiMax, isPhiFolded ? "folded" : "fullPhi"), "RECREATE");
	else comparisonCanvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/1DAccxEffComp_%s%s%s_pt%dto%d_|cosTheta|%.2fto%.2f_%s.png", iState, positiveHist->GetName(), refFrameName, MuonAccName.Data(), ptMin, ptMax, cosThetaMin, cosThetaMax, isPhiFolded ? "folded" : "fullPhi"), "RECREATE");

	delete comparisonCanvas;
	delete frameHist1D;
	delete positiveHist;
	delete negativeHist;

	return;
}

void scan_accEffplots_3Dto1D_comparison(){

	// // for (int ptBin = 0; ptBin < NPtBins; ++ptBin) {
	// for (int ptBin = 1; ptBin < NPtBins; ++ptBin) {
	// 	for (int phiBin = 0; phiBin < NPhiBins; ++phiBin) {
	// 		for (int refFrame = 0; refFrame < 2; ++refFrame) {
	// 			accEffplots_3Dto1D_comparison(gPtBinning[ptBin], gPtBinning[ptBin + 1], refFrame == 0 ? "CS" : "HX", 5, -0.7, 0.7, 1, gPhiBinning[phiBin], gPhiBinning[phiBin + 1], gUpsilonState, kFALSE, kTRUE, "MuonUpsilonTriggerAcc");
	// 		}
	// 	}
	// }	

	// for (int ptBin = 0; ptBin < NPtBins; ++ptBin) {
	for (int ptBin = 1; ptBin < NPtBins; ++ptBin) {
		for (int cosThetaBin = 3; cosThetaBin < NCosThetaBins; ++cosThetaBin) {
			for (int refFrame = 0; refFrame < 2; ++refFrame) {
				accEffplots_3Dto1D_comparison(gPtBinning[ptBin], gPtBinning[ptBin + 1], refFrame == 0 ? "CS" : "HX", 1, gCosThetaBinning[cosThetaBin], gCosThetaBinning[cosThetaBin + 1], 6, -180, 180, gUpsilonState, kFALSE, kFALSE, "MuonUpsilonTriggerAcc");
			}
		}
	}
}