// pre-difinced cosTheta bin edges for the polarization 1D fit
vector<Double_t> setCosThetaBinEdges(Int_t nCosThetaBins){

	vector<Double_t> cosThetaBinEdges;

	// define the bin edges along the cosTheta axis depending on the number of bins
	if (nCosThetaBins == 5) cosThetaBinEdges = {-0.5, -0.3, -0.1, 0.1, 0.3, 0.5};
	else if (nCosThetaBins == 6) cosThetaBinEdges = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
	else if (nCosThetaBins == 7) cosThetaBinEdges = {-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7};
	else if (nCosThetaBins == 8) cosThetaBinEdges = {-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8};
	else if (nCosThetaBins == 9) cosThetaBinEdges = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9};
	else if (nCosThetaBins == 10) cosThetaBinEdges = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
	else {
		cout << "Pre-defined binning not found. Check the input nCosThetaBins" << endl;
		exit(1);
	}

	return cosThetaBinEdges;
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
TCanvas* drawContourPlots(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", TGraph* contour1 = nullptr, TGraph* contour2 = nullptr, TGraph* contour3 = nullptr) {
	TCanvas* contourCanvas = new TCanvas(contour1->GetName(), "", 650, 600);

	contourCanvas->SetLeftMargin(0.17);

	contour1->SetTitle(Form(";cos #theta_{%s};Normalization Factor", refFrameName));

	contour1->GetXaxis()->CenterTitle();
	contour1->GetYaxis()->CenterTitle();

	contour1->GetYaxis()->SetTitleOffset(1.4);

	contour1->GetYaxis()->SetRangeUser(getMinYValue(contour1) * 0.95, getMaxYValue(contour1) * 1.1);

	contour1->SetFillColorAlpha(kBlue-8, 0.2);
	contour1->SetLineColor(kBlue-8);
	contour1->SetLineWidth(2);
	contour1->Draw("ALF");

	contour2->SetFillColorAlpha(kGreen-8, 0.2);
	contour2->SetLineColor(kGreen-8);
	contour2->SetLineWidth(2);
	contour2->Draw("FL SAME");

	contour3->SetFillColorAlpha(kBlue+2, 0.2);
	contour3->SetLineColor(kBlue+2);
	contour3->SetLineWidth(2);
	contour3->Draw("FL SAME");

	TPaveText* text = new TPaveText(0.20, 0.80, 0.57, 0.90, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(CentralityRangeText(gCentralityBinMin, gCentralityBinMax));
	text->AddText(DimuonPtRangeText(ptMin, ptMax));

	text->SetAllWith("", "align", 12);
	text->Draw("SAME");

	TLegend contourLegend(.65, .75, .90, .91, NULL, "brNDC");
	contourLegend.SetTextSize(.05);
	
	contourLegend.AddEntry(contour3, "GMinuit, 1#sigma", "F");
   	contourLegend.AddEntry(contour2, "fitResult, 1#sigma", "F");
   	contourLegend.AddEntry(contour1, "GMinuit, 2#sigma", "F");

	contourLegend.DrawClone();

	gPad->Update();	

	return contourCanvas;
}