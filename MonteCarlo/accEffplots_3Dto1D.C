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

TLine* drawLine(double x1, double y1, double x2, double y2) {
    TLine* line = new TLine(x1, y1, x2, y2);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);

    line->Draw("SAME");
    return line;
}

TPad* drawPulls(TEfficiency* positiveHist, TEfficiency* negativeHist, const char* xTitle, bool isCosTheta = kTRUE) {

	/// Create a new pad
	TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0, 0, 1, .25);
	bottomPad->SetTopMargin(0.01);
	bottomPad->SetBottomMargin(0.4);
	bottomPad->SetRightMargin(0.03);
	
	bottomPad->SetTicks(1, 1);
	
	bottomPad->Draw();
	bottomPad->cd();

	/// Get the number of bins and the range of the histogram
	int nBins = positiveHist->GetTotalHistogram()->GetNbinsX();
	double xmin = positiveHist->GetTotalHistogram()->GetXaxis()->GetBinLowEdge(1);
	double xmax = positiveHist->GetTotalHistogram()->GetXaxis()->GetBinUpEdge(nBins);

	// cout << nBins << " " << xmin << " " << xmax << endl;

	/// Create a new frame
	TH1D* frameHist = new TH1D("frameHist", "", nBins, xmin, xmax);

	frameHist->GetXaxis()->SetTitle(xTitle);
	frameHist->GetXaxis()->SetTitleSize(0.2);
	frameHist->GetXaxis()->SetLabelSize(0.15);
	frameHist->GetXaxis()->SetTitleOffset(0.8);
	frameHist->GetXaxis()->CenterTitle();

	frameHist->GetYaxis()->SetTitle("Pull");
	frameHist->GetYaxis()->SetTitleOffset(0.37);
	frameHist->GetYaxis()->SetTitleSize(0.18);
	frameHist->GetYaxis()->SetLabelSize(0.15);
	frameHist->GetYaxis()->CenterTitle();

	if (isCosTheta) frameHist->GetXaxis()->SetNdivisions(-5);
	else frameHist->GetXaxis()->SetNdivisions(-6);

	frameHist->GetYaxis()->SetRangeUser(-4.5, 4.5);
	frameHist->GetYaxis()->SetNdivisions(405);

	/// Create a new graph
	TGraphErrors* pullGraph = new TGraphErrors(nBins);

	/// Styling the graph
	pullGraph->SetMarkerStyle(20);
	pullGraph->SetMarkerSize(1.2);
	pullGraph->SetMarkerColor(TColor::GetColor("#333333"));
	pullGraph->SetLineColor(TColor::GetColor("#333333"));

	/// Fill the graph
	for (int iBin = 0; iBin < nBins; iBin++) {

		double positiveEff = positiveHist->GetEfficiency(positiveHist->GetGlobalBin(iBin + 1));
		double negativeEff = negativeHist->GetEfficiency(negativeHist->GetGlobalBin(iBin + 1));

		double positiveEffErrUp = positiveHist->GetEfficiencyErrorUp(positiveHist->GetGlobalBin(iBin + 1));
		double negativeEffErrUp = negativeHist->GetEfficiencyErrorUp(negativeHist->GetGlobalBin(iBin + 1));

		double positiveEffErrDown = positiveHist->GetEfficiencyErrorLow(positiveHist->GetGlobalBin(iBin + 1));
		double negativeEffErrDown = negativeHist->GetEfficiencyErrorLow(negativeHist->GetGlobalBin(iBin + 1));

		double xpoint = positiveHist->GetTotalHistogram()->GetXaxis()->GetBinCenter(iBin + 1);
		// cout << xpoint << " " << positiveEff << " " << negativeEff << " " << positiveEffErrUp << " " << negativeEffErrUp << " " << positiveEffErrDown << " " << negativeEffErrDown << endl;

		double pull = 0;
		double pullerr = 0;

		if (positiveEff > negativeEff) {
			pull = (positiveEff - negativeEff) / sqrt(positiveEffErrUp * positiveEffErrUp + negativeEffErrDown * negativeEffErrDown);
			pullerr = sqrt(pow(positiveEffErrUp / sqrt(pow(positiveEffErrUp, 2) + pow(negativeEffErrDown, 2)), 2) + pow(negativeEffErrDown / sqrt(pow(positiveEffErrUp, 2) + pow(negativeEffErrDown, 2)), 2));
		}

		else {
			pull = (positiveEff - negativeEff) / sqrt(positiveEffErrDown * positiveEffErrDown + negativeEffErrUp * negativeEffErrUp);
			pullerr = sqrt(pow(positiveEffErrDown / sqrt(pow(positiveEffErrDown, 2) + pow(negativeEffErrUp, 2)), 2) + pow(negativeEffErrUp / sqrt(pow(positiveEffErrDown, 2) + pow(negativeEffErrUp, 2)), 2));
		}
		// cout << "iBin: " << iBin << " xpoint: " << xpoint << " pull: " << pull << " pullerr: " << pullerr << endl;
		pullGraph->SetPoint(iBin, xpoint, pull);
		pullGraph->SetPointError(iBin, 0, pullerr);
	}

	/// Draw the frame
	frameHist->Draw();

	/// Draw white line at y = 0 to hide the x-axis
	TLine* whiteline = new TLine(xmin, 0, xmax, 0);  // Horizontal line at y = 0
	whiteline->SetLineColor(kWhite);  // Set to white or transparent
	whiteline->SetLineWidth(4);
	whiteline->Draw("SAME");

	/// Draw a horizontal line at y = 0
	TLine* line0 = drawLine(xmin, 0, xmax, 0);

	// TLine* line1 = drawLine(xmin, 1, xmax, 1);
	// TLine* line_1 = drawLine(xmin, -1, xmax, -1);

	TLine* line2 = drawLine(xmin, 2, xmax, 2);
	TLine* line_2 = drawLine(xmin, -2, xmax, -2);

	// TLine* line3 = drawLine(xmin, 3, xmax, 3);
	// TLine* line_3 = drawLine(xmin, -3, xmax, -3);

	TLine* line4 = drawLine(xmin, 4, xmax, 4);
	TLine* line_4 = drawLine(xmin, -4, xmax, -4);

	/// Draw the graph
	pullGraph->Draw("P");	

	/// Redraw the axis
	bottomPad->RedrawAxis();

	/// Return the pad
	return bottomPad;
}

std::vector<std::vector<TEfficiency*>> accEffplots_3Dto1D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 5, Int_t phiMin = -180, Int_t phiMax = 180, Int_t iState = gUpsilonState, Bool_t isPhiFolded = kFALSE, TString accName = "MuonUpsilonTriggerAcc") { // accName = "MuonSimpleAcc", "MuonWithin2018PbPbAcc", or "MuonUpsilonTriggerAcc"
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

	/// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	/// get acceptance maps
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

	/// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosTheta = rebinTEff3DMapCosTheta(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);
	TEfficiency* accMapPhi = rebinTEff3DMapPhi(accMap, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	/// get efficiency maps
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

	/// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effMapCosTheta = rebinTEff3DMapCosTheta(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);
	TEfficiency* effMapPhi = rebinTEff3DMapPhi(effMap, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	/// draw acceptance 1D histograms
	DrawEfficiency1DHist(accMapCosTheta, ptMin, ptMax, iState, kTRUE, kTRUE, refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi, MuonAccName.Data(), isPhiFolded);
	DrawEfficiency1DHist(accMapPhi, ptMin, ptMax, iState, kTRUE, kFALSE, refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi, MuonAccName.Data(), isPhiFolded);

	/// draw efficiency 1D histograms
	DrawEfficiency1DHist(effMapCosTheta, ptMin, ptMax, iState, kFALSE, kTRUE, refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi, MuonAccName.Data(), isPhiFolded);
	DrawEfficiency1DHist(effMapPhi, ptMin, ptMax, iState, kFALSE, kFALSE, refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi, MuonAccName.Data(), isPhiFolded);

	/// draw acceptance x efficiency 1D histograms
	TEfficiency* cosThetaHist = DrawEffxAcc1DHist(accMapCosTheta, effMapCosTheta, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, (double)phiMin, (double)phiMax, kTRUE, refFrameName, iState, MuonAccName.Data(), isPhiFolded);
	TEfficiency* phiHist = DrawEffxAcc1DHist(accMapPhi, effMapPhi, ptMin, ptMax, nPhiBins, phiBinEdges, cosThetaMin, cosThetaMax, kFALSE, refFrameName, iState, MuonAccName.Data(), isPhiFolded);

	std::vector<std::vector<TEfficiency*>> hists = {{accMapCosTheta, accMapPhi},
									                {effMapCosTheta, effMapPhi},
									                {cosThetaHist, phiHist}};
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

	/// get acceptance maps 
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

	/// get acceptance, efficiency, and, acceptance x efficiency 1D histograms: positive and negative sides
	std::vector<std::vector<TEfficiency*>> positiveHists;
	std::vector<std::vector<TEfficiency*>> negativeHists;
	std::vector<TEfficiency*> positiveHist;
	std::vector<TEfficiency*> negativeHist;

	if (isCosTheta == kTRUE) {
		positiveHists = accEffplots_3Dto1D(ptMin, ptMax, refFrameName, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, iState, isPhiFolded, accName);
		negativeHists = accEffplots_3Dto1D(ptMin, ptMax, refFrameName, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, -phiMax, -phiMin, iState, isPhiFolded, accName);
		positiveHist = {positiveHists[0][0], positiveHists[1][0], positiveHists[2][0]};
		negativeHist = {negativeHists[0][0], negativeHists[1][0], negativeHists[2][0]};
	}
	else {
		positiveHists = accEffplots_3Dto1D(ptMin, ptMax, refFrameName, nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax, iState, isPhiFolded, accName);
		negativeHists = accEffplots_3Dto1D(ptMin, ptMax, refFrameName, nCosThetaBins, -cosThetaMax, -cosThetaMin, nPhiBins, phiMin, phiMax, iState, isPhiFolded, accName);
		positiveHist = {positiveHists[0][1], positiveHists[1][1], positiveHists[2][1]};
		negativeHist = {negativeHists[0][1], negativeHists[1][1], negativeHists[2][1]};
	}

	/// draw plots
	/// empty frame for the axes
	std::vector<TH1D*> frameHist1D = {nullptr, nullptr, nullptr};
	std::vector<TCanvas*> comparisonCanvas = {nullptr, nullptr, nullptr};
	
	/// loop over for acceptance, efficiency, and acceptance x efficiency 1D histograms each
	for (int plotTypeBin = 0; plotTypeBin < 3; plotTypeBin++) {

		if (isCosTheta == kTRUE) {
			if (plotTypeBin == 0) frameHist1D[plotTypeBin] = new TH1D(Form("frameHist1D%s", "_acceptance"), "", nCosThetaBins, cosThetaBinEdges.data());
			else if (plotTypeBin == 1) frameHist1D[plotTypeBin] = new TH1D(Form("frameHist1D%s", "_efficiency"), "", nCosThetaBins, cosThetaBinEdges.data());
			else frameHist1D[plotTypeBin] = new TH1D(Form("frameHist1D%s", "_accxeff"), "", nCosThetaBins, cosThetaBinEdges.data());
		}
		else {
			if (plotTypeBin == 0) frameHist1D[plotTypeBin] = new TH1D(Form("frameHist1D%s", "_acceptance"), "", nPhiBins, phiBinEdges.data());
			else if (plotTypeBin == 1) frameHist1D[plotTypeBin] = new TH1D(Form("frameHist1D%s", "_efficiency"), "", nPhiBins, phiBinEdges.data());
			else frameHist1D[plotTypeBin] = new TH1D(Form("frameHist1D%s", "_accxeff"), "", nPhiBins, phiBinEdges.data());
		}

		/// draw acceptance x efficiency 1D histograms for positive and negative phi
		comparisonCanvas[plotTypeBin] = new TCanvas(Form("comparisonCanvas%d", plotTypeBin), Form("comparisonCanvas%d", plotTypeBin), 600, 600);
		
		comparisonCanvas[plotTypeBin]->SetRightMargin(0.05);

		TPad *pad1 = new TPad("pad1", "First Pad", 0.0, 0.25, 1.0, 1.0);
		pad1->Draw();
		pad1->cd(); 

		/// Set the margins
		pad1->SetTopMargin(0.07);
		pad1->SetBottomMargin(0.03);
		pad1->SetRightMargin(0.03);

		frameHist1D[plotTypeBin]->Draw("same");
		frameHist1D[plotTypeBin]->GetYaxis()->SetRangeUser(0, 1.4);

		/// Set the axis titles
		TString xTitle = "";

		if (isCosTheta == kTRUE) {
			xTitle = CosThetaVarTitle(refFrameName);
			frameHist1D[plotTypeBin]->GetXaxis()->SetTitle(CosThetaVarTitle(refFrameName));
			frameHist1D[plotTypeBin]->GetXaxis()->SetNdivisions(-5);
		}
		else {
			xTitle = PhiAxisTitle(refFrameName);
			frameHist1D[plotTypeBin]->GetXaxis()->SetTitle(PhiAxisTitle(refFrameName));
			frameHist1D[plotTypeBin]->GetXaxis()->SetNdivisions(-6);
		}

		frameHist1D[plotTypeBin]->GetXaxis()->CenterTitle();
		frameHist1D[plotTypeBin]->GetXaxis()->SetTitleSize(0);
		frameHist1D[plotTypeBin]->GetXaxis()->SetLabelSize(0);
		
		if (plotTypeBin == 0) frameHist1D[plotTypeBin]->GetYaxis()->SetTitle("Acceptance");
		else if (plotTypeBin == 1) frameHist1D[plotTypeBin]->GetYaxis()->SetTitle("Efficiency");
		else frameHist1D[plotTypeBin]->GetYaxis()->SetTitle(TEfficiencyMainTitle(iState, "Acceptance x Efficiency"));

		/// draw legend
		TLatex legend;
		legend.SetTextAlign(22);
		legend.SetTextSize(0.05);
		legend.DrawLatexNDC(.55, .85, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));

		if (strcmp(MuonAccName.Data(), "_TriggerAcc") == 0) legend.DrawLatexNDC(.55, .78, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
		else if (strcmp(MuonAccName.Data(), "_SimpleAcc") == 0) legend.DrawLatexNDC(.55, .78, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", iState));
		else if (strcmp(MuonAccName.Data(), "_2018PbPbAcc") == 0) legend.DrawLatexNDC(.55, .78, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", iState));

		/// draw acceptance x efficiency 1D histograms for positive and negative sides
		positiveHist[plotTypeBin]->Draw("EP same");
		positiveHist[plotTypeBin]->SetMarkerColor(TColor::GetColor("#A559AA"));
		positiveHist[plotTypeBin]->SetMarkerStyle(20);
		positiveHist[plotTypeBin]->SetLineColor(TColor::GetColor("#A559AA"));
		
		negativeHist[plotTypeBin]->Draw("PE same");
		negativeHist[plotTypeBin]->SetMarkerColor(TColor::GetColor("#009ADE"));
		negativeHist[plotTypeBin]->SetMarkerStyle(24);
		negativeHist[plotTypeBin]->SetLineColor(TColor::GetColor("#009ADE"));

		TLegend* legend1D = new TLegend(0.17, 0.58, 0.4, 0.73);
		legend1D->SetBorderSize(0);
		legend1D->SetFillStyle(0);
		if (isCosTheta == kTRUE) {
			legend1D->AddEntry(positiveHist[plotTypeBin], Form("%d < #varphi < %d#circ", phiMin, phiMax), "lep");
			legend1D->AddEntry(negativeHist[plotTypeBin], Form("%d < #varphi < %d#circ", -phiMax, -phiMin), "lep");
		}
		else {
			legend1D->AddEntry(positiveHist[plotTypeBin], Form("%.2f < cos#theta < %.2f", cosThetaMin, cosThetaMax), "lep");
			legend1D->AddEntry(negativeHist[plotTypeBin], Form("%.2f < cos#theta < %.2f", -cosThetaMax, -cosThetaMin), "lep");
		}
		legend1D->Draw();
		
		CMS_lumi(pad1, Form("#varUpsilon(%dS) Pythia 8 (5.02 TeV)", iState));

		gPad->Update();
		comparisonCanvas[plotTypeBin]->Update();

		/// Draw the pulls at the bottom
		comparisonCanvas[plotTypeBin]->cd();
		
		TPad* bottomPad = drawPulls(positiveHist[plotTypeBin], negativeHist[plotTypeBin], xTitle.Data(), isCosTheta);

		bottomPad->Draw();

		gPad->Update();
		comparisonCanvas[plotTypeBin]->Update();

		/// Save the canvas
		gSystem->mkdir(Form("AcceptanceMaps/%dS/analysisBin", iState), kTRUE);
		gSystem->mkdir(Form("EfficiencyMaps/%dS/analysisBin", iState), kTRUE);
		gSystem->mkdir(Form("AccxEffMaps/%dS/analysisBin", iState), kTRUE);

		if (isCosTheta == kTRUE) { 
			if (plotTypeBin == 0) comparisonCanvas[plotTypeBin]->SaveAs(Form("AcceptanceMaps/%dS/analysisBin/1DAccComp_%s%s%s_pt%dto%d_absPhi%dto%d_%s.png", iState, positiveHist[plotTypeBin]->GetName(), refFrameName, MuonAccName.Data(), ptMin, ptMax, phiMin, phiMax, isPhiFolded ? "folded" : "fullPhi"), "RECREATE");
			else if (plotTypeBin == 1) comparisonCanvas[plotTypeBin]->SaveAs(Form("EfficiencyMaps/%dS/analysisBin/1DEffComp_%s%s%s_pt%dto%d_absPhi%dto%d_%s.png", iState, positiveHist[plotTypeBin]->GetName(), refFrameName, MuonAccName.Data(), ptMin, ptMax, phiMin, phiMax, isPhiFolded ? "folded" : "fullPhi"), "RECREATE");
			else if (plotTypeBin == 2) comparisonCanvas[plotTypeBin]->SaveAs(Form("AccxEffMaps/%dS/analysisBin/1DAccxEffComp_%s%s%s_pt%dto%d_absPhi%dto%d_%s.png", iState, positiveHist[plotTypeBin]->GetName(), refFrameName, MuonAccName.Data(), ptMin, ptMax, phiMin, phiMax, isPhiFolded ? "folded" : "fullPhi"), "RECREATE");
		}
		
		else {
			if (plotTypeBin == 0) comparisonCanvas[plotTypeBin]->SaveAs(Form("AcceptanceMaps/%dS/analysisBin/1DAccComp_%s%s%s_pt%dto%d_absCosTheta%.2fto%.2f_%s.png", iState, positiveHist[plotTypeBin]->GetName(), refFrameName, MuonAccName.Data(), ptMin, ptMax, cosThetaMin, cosThetaMax, isPhiFolded ? "folded" : "fullPhi"), "RECREATE");
			else if (plotTypeBin == 1) comparisonCanvas[plotTypeBin]->SaveAs(Form("EfficiencyMaps/%dS/analysisBin/1DEffComp_%s%s%s_pt%dto%d_absCosTheta%.2fto%.2f_%s.png", iState, positiveHist[plotTypeBin]->GetName(), refFrameName, MuonAccName.Data(), ptMin, ptMax, cosThetaMin, cosThetaMax, isPhiFolded ? "folded" : "fullPhi"), "RECREATE");
			else if (plotTypeBin == 2) comparisonCanvas[plotTypeBin]->SaveAs(Form("AccxEffMaps/%dS/analysisBin/1DAccxEffComp_%s%s%s_pt%dto%d_absCosTheta%.2fto%.2f_%s.png", iState, positiveHist[plotTypeBin]->GetName(), refFrameName, MuonAccName.Data(), ptMin, ptMax, cosThetaMin, cosThetaMax, isPhiFolded ? "folded" : "fullPhi"), "RECREATE");
		}
		/// delete objects when running the scan function to avoid memory leak (comment out to check the plots)
		delete comparisonCanvas[plotTypeBin];
		delete frameHist1D[plotTypeBin];
		delete positiveHist[plotTypeBin];
		delete negativeHist[plotTypeBin];
	}
	return;
}

void scan_accEffplots_3Dto1D_comparison(){

	// for (int ptBin = 0; ptBin < NPtBins; ++ptBin) {
	for (int ptBin = 1; ptBin < NPtBins; ++ptBin) {
		for (int phiBin = 0; phiBin < NPhiBins; ++phiBin) {
			for (int refFrame = 0; refFrame < 2; ++refFrame) {
				accEffplots_3Dto1D_comparison(gPtBinning[ptBin], gPtBinning[ptBin + 1], refFrame == 0 ? "CS" : "HX", 5, -0.7, 0.7, 1, gPhiBinning[phiBin], gPhiBinning[phiBin + 1], gUpsilonState, kFALSE, kTRUE, "MuonUpsilonTriggerAcc");
			}
		}
	}	

	// for (int ptBin = 0; ptBin < NPtBins; ++ptBin) {
	for (int ptBin = 1; ptBin < NPtBins; ++ptBin) {
		for (int cosThetaBin = 3; cosThetaBin < NCosThetaBins; ++cosThetaBin) {
			for (int refFrame = 0; refFrame < 2; ++refFrame) {
				accEffplots_3Dto1D_comparison(gPtBinning[ptBin], gPtBinning[ptBin + 1], refFrame == 0 ? "CS" : "HX", 1, gCosThetaBinning[cosThetaBin], gCosThetaBinning[cosThetaBin + 1], 6, -180, 180, gUpsilonState, kFALSE, kFALSE, "MuonUpsilonTriggerAcc");
			}
		}
	}
}