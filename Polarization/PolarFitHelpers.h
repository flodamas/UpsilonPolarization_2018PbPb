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