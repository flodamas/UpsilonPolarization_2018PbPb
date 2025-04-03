#include "../Tools/FitShortcuts.h"

/// apply weights and errors to each costheta bin
Float_t correctRawYield1DHist(TH1D* standardCorrectedHist, TEfficiency* accMap, TEfficiency* effMap, TH1D* systEff, Int_t nBins = 10, const char** bkgShapeNames = nullptr, const char** fitModelNames = nullptr) {
	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	TCanvas* massCanvas = 0;

	Float_t maxYield = 0;

	for (Int_t iBin = 0; iBin < nBins; iBin++) {
		Double_t weight = 0;

		RooRealVar* yield1S = new RooRealVar("yield1S", "", 1000);
		RooRealVar* yield2S = new RooRealVar("yield2S", "", 100);
		RooRealVar* yield3S = new RooRealVar("yield3S", "", 10);

		// get the corresponding weights
		double acceptance = accMap->GetEfficiency(iBin + 1);
		double efficiency = effMap->GetEfficiency(iBin + 1);

		// calculate weight
		weight = 1. / (acceptance * efficiency);

		// propagate both scale factor uncertainties and efficiency stat errors to the weight
		double relSystUnc = systEff->GetBinContent(iBin + 1);

		double relEffUncHigh = effMap->GetEfficiencyErrorUp(iBin + 1) / efficiency;
		double relEffUncLow = effMap->GetEfficiencyErrorLow(iBin + 1) / efficiency;

		double relAccUncHigh = accMap->GetEfficiencyErrorUp(iBin + 1) / acceptance;
		double relAccUncLow = accMap->GetEfficiencyErrorLow(iBin + 1) / acceptance;

		totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
		totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

		totalUncHigh = totalRelUncHigh * efficiency * acceptance;
		totalUncLow = totalRelUncLow * efficiency * acceptance;

		// get yields and their uncertainties
		RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("RawData_%s", bkgShapeNames[iBin]), fitModelNames[iBin]);

		double yield1SVal = yield1S->getVal();

		double yield1SUnc = yield1S->getError();

		// set the bin contents reflecting weights
		standardCorrectedHist->SetBinContent(iBin + 1, yield1SVal * weight);

		// standardCorrectedHist->SetBinError(iBin + 1, yield1SErr * weight);
		// standardCorrectedHist->SetBinError(iBin + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh));
		// standardCorrectedHist->SetBinError(iBin + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));

		// standardCorrectedHist->SetBinError(iBin + 1, yield1SUnc * weight);

		// standardCorrectedHist->SetBinError(iBin + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

		standardCorrectedHist->SetBinError(iBin + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

		// // fill uncertainty histograms
		// statHighEff->SetBinContent(iBin + 1, relEffUncHigh);
		// statLowEff->SetBinContent(iBin + 1, relEffUncLow);
		// statHighAcc->SetBinContent(iBin + 1, relAccUncHigh);
		// statLowAcc->SetBinContent(iBin + 1, relAccUncLow);
		// yield1SUnc->SetBinContent(iBin + 1, yield1SUnc / yield1SVal);
		// totalRelUnc->SetBinContent(iBin + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

		if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;

		delete yield1S;
	}

	return maxYield;
}

TCanvas* drawUncertaintyPlot1D(const char* refFrameName, TH1D* uncPlot1, TH1D* uncPlot2, TH1D* uncPlot3, TH1D* uncPlot4, TH1D* uncPlot5, TH1D* uncPlot6, TH1D* uncPlot7) {
	TCanvas* errCanvas = new TCanvas("errCanvas", "errCanvas", 650, 600);

	uncPlot1->GetYaxis()->SetRangeUser(0, 1);

	uncPlot1->SetXTitle(Form("cos #theta_{%s}", refFrameName));
	uncPlot1->SetYTitle("Relative uncertainty");

	uncPlot1->GetXaxis()->CenterTitle();
	uncPlot1->GetYaxis()->CenterTitle();

	Int_t lineWidth = 6;

	uncPlot1->SetLineWidth(lineWidth);
	uncPlot1->SetLineColor(kRed + 1);

	uncPlot1->Draw();

	uncPlot2->SetLineWidth(lineWidth);
	uncPlot2->SetLineColor(kRed - 7);

	uncPlot2->Draw("SAME");

	uncPlot3->SetLineWidth(lineWidth);
	uncPlot3->SetLineColor(kOrange);

	uncPlot3->Draw("SAME");

	uncPlot4->SetLineWidth(lineWidth);
	uncPlot4->SetLineColor(kAzure - 9);

	uncPlot4->Draw("SAME");

	uncPlot5->SetLineWidth(lineWidth);
	uncPlot5->SetLineColor(kBlue - 3);

	uncPlot5->Draw("SAME");

	uncPlot6->SetLineWidth(lineWidth);
	uncPlot6->SetLineColor(kViolet - 8);

	uncPlot6->Draw("SAME");

	uncPlot7->SetLineWidth(lineWidth);
	uncPlot7->SetLineColor(kBlack);

	uncPlot7->Draw("SAME");

	return errCanvas;
}

// get maximum y value of the contour plot
double getMaxYValue(TGraph* graph) {
	double maxY = -1e20; // Initialize with a very small value

	// Get the number of points in the graph
	int nPoints = graph->GetN();

	// Iterate through each point in the graph
	for (int iPoint = 0; iPoint < nPoints; iPoint++) {
		double x, y;
		graph->GetPoint(iPoint, x, y); // Get x and y coordinates of the point

		// Update the maximum y-value if the current y-value is greater
		if (y > maxY) {
			maxY = y;
		}
	}

	return maxY;
}

// get minimum y value of the contour plot
double getMinYValue(TGraph* graph) {
	double minY = 1e20; // Initialize with a very large value

	// Get the number of points in the graph
	int nPoints = graph->GetN();

	// Iterate through each point in the graph
	for (int iPoint = 0; iPoint < nPoints; iPoint++) {
		double x, y;
		graph->GetPoint(iPoint, x, y); // Get x and y coordinates of the point

		// Update the maximum y-value if the current y-value is greater
		if (y < minY) {
			minY = y;
		}
	}

	return minY;
}

// draw contour plots
TCanvas* drawContourPlots(Int_t ptMin = 0, Int_t ptMax = 30, Double_t cosThetaMin = -1, Double_t cosThetaMax = 1, const char* refFrameName = "CS", TGraph* contour1 = nullptr, TGraph* contour2 = nullptr, TGraph* contour3 = nullptr) {
	TCanvas* contourCanvas = new TCanvas(contour1->GetName(), "", 650, 600);

	contourCanvas->SetLeftMargin(0.17);

	TH2D* contourPlotFrame = new TH2D("contourPlotFrame", ";#lambda_{#theta};Normalization Factor", 20, -2, 2, 100, 0, 7000);

	// contour1->SetTitle(";#lambda#theta;Normalization Factor");

	contourPlotFrame->GetXaxis()->CenterTitle();
	contourPlotFrame->GetYaxis()->CenterTitle();

	contourPlotFrame->GetYaxis()->SetTitleOffset(1.4);

	contourPlotFrame->GetXaxis()->SetRangeUser(-2, 2);
	contourPlotFrame->GetYaxis()->SetRangeUser(0, 7000);

	contourPlotFrame->Draw();

	if (contour3) {
		contour3->SetFillColorAlpha(kGreen - 8, 0.2);
		contour3->SetLineColor(kGreen - 8);
		contour3->SetLineWidth(2);
		contour3->Draw("FL SAME");
	}

	contour2->SetFillColorAlpha(kAzure - 9, 0.7);
	contour2->SetLineColor(kAzure - 9);
	contour2->SetLineWidth(2);
	contour2->Draw("FL SAME");

	contour1->SetFillColorAlpha(kOrange, 0.9);
	contour1->SetLineColor(kOrange + 1);
	contour1->SetLineWidth(2);
	contour1->Draw("FL SAME");

	TPaveText* text = new TPaveText(0.20, 0.70, 0.57, 0.90, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(CentralityRangeText(gCentralityBinMin, gCentralityBinMax));
	text->AddText(DimuonPtRangeText(ptMin, ptMax));
	text->AddText(CosThetaRangeText(refFrameName, cosThetaMin, cosThetaMax));

	text->SetAllWith("", "align", 12);
	text->Draw("SAME");

	TLegend contourLegend(.21, .54, .46, .70, NULL, "brNDC");
	contourLegend.SetTextSize(.05);

	contourLegend.AddEntry(contour1, "1#sigma", "F");
	contourLegend.AddEntry(contour2, "2#sigma", "F");

	contourLegend.DrawClone();

	gPad->Update();

	return contourCanvas;
}

TCanvas* draw2DMap(TH2D* mapCosThetaPhi, const char* refFrameName = "CS", Int_t nCosThetaBins = 5, const std::vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 5, const std::vector<Double_t>& phiBinEdges = {}, Bool_t LEGO = kFALSE, Bool_t isRange0to1 = kFALSE, Int_t iState = 1, Bool_t isPhiFolded = kTRUE) {
	TCanvas* map2DCanvas = new TCanvas(mapCosThetaPhi->GetName(), "", 680, 600);

	// SetColorPalette(gPreferredColorPaletteName);
	SetColorPalette("TamDragon");

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	TH2D* frameHist = new TH2D(Form("%sframeHist", mapCosThetaPhi->GetName()), " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins] + phiStep);

	if (LEGO) {
		map2DCanvas->SetLeftMargin(0.18);
		map2DCanvas->SetRightMargin(0.09);

		mapCosThetaPhi->Draw("LEGO E");

		mapCosThetaPhi->SetXTitle(Form("cos #theta_{%s}", refFrameName));
		if (isPhiFolded) mapCosThetaPhi->SetYTitle(Form("|#varphi_{%s}| (#circ)", refFrameName));
		else mapCosThetaPhi->SetYTitle(Form("#varphi_{%s} (#circ)", refFrameName));
		mapCosThetaPhi->SetZTitle(Form("#varUpsilon(%dS) yields", iState));

		mapCosThetaPhi->GetXaxis()->SetNdivisions(-500 - (nCosThetaBins));
		// mapCosThetaPhi->GetXaxis()->SetNdivisions(-200 - (nCosThetaBins / 2.));
		if (isPhiFolded) mapCosThetaPhi->GetYaxis()->SetNdivisions(-500 - (nPhiBins));
		else mapCosThetaPhi->GetYaxis()->SetNdivisions(-300 - (nPhiBins - 1));
		// else mapCosThetaPhi->GetYaxis()->SetNdivisions(-300 - (nPhiBins / 3.));

		mapCosThetaPhi->GetXaxis()->CenterTitle();
		mapCosThetaPhi->GetYaxis()->CenterTitle();

		// Set title offsets
		mapCosThetaPhi->GetXaxis()->SetTitleOffset(1.3);
		mapCosThetaPhi->GetYaxis()->SetTitleOffset(1.5);
		mapCosThetaPhi->GetZaxis()->SetTitleOffset(1.8);

		if (isPhiFolded) mapCosThetaPhi->GetYaxis()->SetRangeUser(phiBinEdges[0], phiBinEdges[nPhiBins]);
		// else mapCosThetaPhi->GetYaxis()->SetRangeUser(phiBinEdges[0], phiBinEdges[nPhiBins]);
		else mapCosThetaPhi->GetYaxis()->SetRangeUser(phiBinEdges[0], phiBinEdges[nPhiBins - 1]);

		mapCosThetaPhi->SetStats(0);
	}

	else {
		map2DCanvas->SetRightMargin(0.18);

		gStyle->SetPadRightMargin(0.2);

		gPad->Modified();
		gPad->Update();

		frameHist->Draw("COLZ");

		mapCosThetaPhi->Draw("SAME COLZ");

		frameHist->SetXTitle(Form("cos #theta_{%s}", refFrameName));
		if (isPhiFolded) frameHist->SetYTitle(Form("|#varphi_{%s}| (#circ)", refFrameName));
		else frameHist->SetYTitle(Form("#varphi_{%s} (#circ)", refFrameName));

		// frameHist->GetXaxis()->SetNdivisions(-500 - (nCosThetaBins));
		// frameHist->GetYaxis()->SetNdivisions(-500 - (nPhiBins + 1));

		/// Ndivision setting for the finer binning (costheta: 20 bins from -1 to 1, phi: 18 bins from -180 to 180)
		frameHist->GetXaxis()->SetNdivisions(-200 - (nCosThetaBins / 2.));
		frameHist->GetYaxis()->SetNdivisions(-300 - (nPhiBins / 3.));

		frameHist->GetXaxis()->CenterTitle();
		frameHist->GetYaxis()->CenterTitle();
		frameHist->GetZaxis()->CenterTitle();

		mapCosThetaPhi->GetZaxis()->SetMaxDigits(3);

		frameHist->SetStats(0);

		if (!isRange0to1)
			frameHist->GetZaxis()->SetRangeUser(0, mapCosThetaPhi->GetMaximum());
		else
			frameHist->GetZaxis()->SetRangeUser(0, 1);
	}

	// CMS_lumi(map2DCanvas, gCMSLumiText);
	CMS_lumi(map2DCanvas, "#varUpsilon(1S) Pythia 8 (5.02 TeV)");

	map2DCanvas->Modified();
	map2DCanvas->Update();

	return map2DCanvas;
}

// display the uncertainties signal extraction yield on each bin of 2D yield map
void display2DMapContents(TH2D* mapCosThetaPhi, Int_t nCosThetaBins = 10, Int_t nPhiBins = 6, Bool_t displayError = kFALSE) {
	if (!mapCosThetaPhi) {
		std::cout << "no 2D map found!!!" << std::endl;
		exit(1);
	}

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			// Get the yield and uncertainty values
			Double_t binVal = mapCosThetaPhi->GetBinContent(iCosTheta + 1, iPhi + 1);

			Double_t binUnc = mapCosThetaPhi->GetBinError(iCosTheta + 1, iPhi + 1);

			// Get the bin center coordinates
			double x = mapCosThetaPhi->GetXaxis()->GetBinCenter(iCosTheta + 1);

			double y = mapCosThetaPhi->GetYaxis()->GetBinCenter(iPhi + 1);

			// Create a TLatex object to write the signal extraction yield uncertainties on each bin
			TLatex latex;
			latex.SetTextSize(0.05); // Adjust text size as needed
			latex.SetTextAlign(22);  // Center alignment
			latex.SetTextColor(kWhite);

			if (displayError)
				latex.DrawLatex(x, y, Form("%.1f%%", binUnc / binVal * 100));

			else
				latex.DrawLatex(x, y, Form("%.2f", binVal));
		}
	}
}

double readYieldDiffMean(Int_t ptMin = 0, Int_t ptMax = 2, const char* refFrameName = "HX", Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, Int_t phiMin = 0, Int_t phiMax = 180, const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ExpTimesErr", int nPseudoExperiments = 10) {
	const char* fileLocation = "../SystematicUncertainties/YieldDifferencePlots_mass7to11p5_final/";

	// Open the ROOT file
	TFile* file = openFile(Form("%s%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%d.root", fileLocation, signalShapeName, bkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments));

	// Retrieve the histogram
	TH1* hist = (TH1*)file->Get("mergedHist"); // Replace "histogram_name" with your actual histogram's name
	if (!hist) {
		std::cerr << "Error: Could not find the histogram in the file!" << std::endl;
		file->Close();
		exit(1);
	}

	// Get the mean value
	double mean = hist->GetMean();
	std::cout << "The mean value of the histogram is: " << mean << std::endl;

	// Close the file
	file->Close();

	return mean;
}
