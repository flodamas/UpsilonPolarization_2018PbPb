// pre-difinced cosTheta bin edges for the polarization 1D fit
vector<Double_t> setCosThetaBinEdges(Int_t nCosThetaBins, Double_t cosThetaMin, Double_t cosThetaMax, Bool_t isUniform = kTRUE){

	vector<Double_t> cosThetaBinEdges = {};

	// define the bin edges along the cosTheta axis depending on the number of bins
	for (Int_t iCosTheta = 0; iCosTheta <= nCosThetaBins; iCosTheta++) {

		Double_t cosThetaBinWidth = (cosThetaMax - cosThetaMin) / nCosThetaBins;

		cosThetaBinEdges.push_back(cosThetaMin +  cosThetaBinWidth * iCosTheta); 

	}

	// // for the case that the bin width varies
	// if (!isUniform) {}

	// else if (!cosThetaBinEdges) {
	// 	cout << "Pre-defined binning not found. Check the input nCosThetaBins" << endl;
	// 	exit(1);
	// }

	return cosThetaBinEdges;
}

vector<Double_t> setPhiBinEdges(Int_t nPhiBins, Int_t phiMin, Int_t phiMax, Bool_t isUniform = kTRUE){

	// define the bin edges along the cosTheta axis depending on the number of bins
	vector<Double_t> phiBinEdges = {};

	// define the bin edges along the cosTheta axis depending on the number of bins
	for (Int_t iPhi = 0; iPhi <= nPhiBins; iPhi++) {

		Double_t phiBinWidth = (phiMax - phiMin) / nPhiBins;

		phiBinEdges.push_back(phiMin +  phiBinWidth * iPhi); 

	}

	return phiBinEdges;
}

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

TCanvas* drawUncertaintyPlot(const char* refFrameName, TH1D* uncPlot1, TH1D* uncPlot2, TH1D* uncPlot3, TH1D* uncPlot4, TH1D* uncPlot5, TH1D* uncPlot6, TH1D* uncPlot7){

	TCanvas *errCanvas = new TCanvas("errCanvas", "errCanvas", 650, 600);

	uncPlot1->GetYaxis()->SetRangeUser(0, 1);

	uncPlot1->SetXTitle(Form("cos #theta_{%s}", refFrameName));
	uncPlot1->SetYTitle("Relative Uncertainty");

	uncPlot1->GetXaxis()->CenterTitle();
	uncPlot1->GetYaxis()->CenterTitle();

	Int_t lineWidth = 6;

	uncPlot1->SetLineWidth(lineWidth);
	uncPlot1->SetLineColor(kRed+1);

	uncPlot1->Draw();

	uncPlot2->SetLineWidth(lineWidth);
	uncPlot2->SetLineColor(kRed-7);

	uncPlot2->Draw("SAME");

	uncPlot3->SetLineWidth(lineWidth);
	uncPlot3->SetLineColor(kOrange);

	uncPlot3->Draw("SAME");

	uncPlot4->SetLineWidth(lineWidth);
	uncPlot4->SetLineColor(kAzure-9);

	uncPlot4->Draw("SAME");	

	uncPlot5->SetLineWidth(lineWidth);
	uncPlot5->SetLineColor(kBlue-3);

	uncPlot5->Draw("SAME");

	uncPlot6->SetLineWidth(lineWidth);
	uncPlot6->SetLineColor(kViolet-8);

	uncPlot6->Draw("SAME");

	uncPlot7->SetLineWidth(lineWidth);
	uncPlot7->SetLineColor(kBlack);

	uncPlot7->Draw("SAME");

	return errCanvas;
}

// get maximum y value of the contour plot
double getMaxYValue(TGraph *graph) {
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
double getMinYValue(TGraph *graph) {
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

	TH2D *contourPlotFrame = new TH2D("contourPlotFrame", ";#lambda_{#theta};Normalization Factor", 20, -2, 2, 100, 0, 7000);

	// contour1->SetTitle(";#lambda#theta;Normalization Factor");

	contourPlotFrame->GetXaxis()->CenterTitle();
	contourPlotFrame->GetYaxis()->CenterTitle();

	contourPlotFrame->GetYaxis()->SetTitleOffset(1.4);
	
	contourPlotFrame->GetXaxis()->SetRangeUser(-2, 2);
	contourPlotFrame->GetYaxis()->SetRangeUser(0, 7000);

	contourPlotFrame->Draw();

	if (contour3){
		contour3->SetFillColorAlpha(kGreen-8, 0.2);
		contour3->SetLineColor(kGreen-8);
		contour3->SetLineWidth(2);
		contour3->Draw("FL SAME");
	}

	contour2->SetFillColorAlpha(kAzure-9, 0.7);
	contour2->SetLineColor(kAzure-9);
	contour2->SetLineWidth(2);
	contour2->Draw("FL SAME");

	contour1->SetFillColorAlpha(kOrange, 0.9);
	contour1->SetLineColor(kOrange+1);
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

// display the uncertainties signal extraction yield on each bin of 2D yield map
void displayYieldUncertainties(TH2D* yieldMap, Int_t nCosThetaBins = 10, Int_t nPhiBins = 6){

	if(!yieldMap) {
		cout << "no yield map found!!!" << endl;
		exit(1);
	}

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {

		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {

			// Get the yield and uncertainty values
			Double_t yield1SVal = yieldMap->GetBinContent(iCosTheta + 1, iPhi +1);

			Double_t yield1SUnc = yieldMap->GetBinError(iCosTheta + 1, iPhi +1);

            // Get the bin center coordinates
            double x = yieldMap->GetXaxis()->GetBinCenter(iCosTheta + 1);

            double y = yieldMap->GetYaxis()->GetBinCenter(iPhi + 1);

            // Create a TLatex object to write the signal extraction yield uncertainties on each bin
            TLatex latex;
            latex.SetTextSize(0.02);  // Adjust text size as needed
            latex.SetTextAlign(22);   // Center alignment
            latex.SetTextColor(kWhite);
            latex.DrawLatex(x, y, Form("%.2f%%", yield1SUnc / yield1SVal * 100));			
		}
	}
}