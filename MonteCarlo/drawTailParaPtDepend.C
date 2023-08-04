using namespace RooFit;

void drawTailParaPtDepend() {
	
	/// Define variables
	RooRealVar alphaInf("alphaInf", "", 0.1, 10);
	RooRealVar orderInf("orderInf", "", 0.1, 10);
	RooRealVar alphaSup("alphaSup", "", 0.1, 10);
	RooRealVar orderSup("orderSup", "", 0.1, 40);

	RooArgSet tailParams(alphaInf, orderInf, alphaSup, orderSup);


	Int_t centMin = 0;
	Int_t centMax = 90; 
	const Float_t ptCuts[8] = {0, 2, 4, 6, 8, 12, 16, 20};
	Int_t arrayLength = sizeof(ptCuts)/sizeof(Int_t);
	Int_t nBins = arrayLength -1 ;

	/// Define histogram
	TH1D* haH = new TH1D("haH", ";p_{T} (GeV); #alpha",nBins, ptCuts);
	TH1D* haL = new TH1D("haL", ";p_{T} (GeV); #alpha",nBins, ptCuts);
	TH1D* hnH = new TH1D("hnH", ";p_{T} (GeV); n",nBins, ptCuts);
	TH1D* hnL = new TH1D("hnL", ";p_{T} (GeV); n",nBins, ptCuts);

	// Read tail parameters from the file depending on pT
	for(int ievent = 0; ievent < nBins; ievent++){
	const char* mcFileName = Form("../MonteCarlo/SignalParameters/symCoreDSCB_cent%dto%d_pt%dto%d.txt", centMin, centMax, (int)ptCuts[ievent], (int)ptCuts[ievent+1]);

		if (fopen(mcFileName, "r")) {
			cout << endl
			     << "Found " << mcFileName << " file, will read the signal tail parameters from it" << endl;
			tailParams.readFromFile(mcFileName);
		} else {
			cout << endl
			     << mcFileName << " file does not seem to exist, you need to extract the signal tail paramaters from MC fit first!" << endl;
		}

		// fix the tail parameters
		alphaInf.setConstant();
		orderInf.setConstant();
		alphaSup.setConstant();
		orderSup.setConstant();

		// Fill the histogram with the tail parameters read from the file
		haH -> SetBinContent(ievent+1, alphaInf.getVal());
		haH -> SetBinError(ievent+1, alphaInf.getError());
		haL -> SetBinContent(ievent+1, alphaSup.getVal());
		haL -> SetBinError(ievent+1, alphaSup.getError());
		hnH -> SetBinContent(ievent+1, orderInf.getVal());
		hnH -> SetBinError(ievent+1, orderInf.getError());
		hnL -> SetBinContent(ievent+1, orderSup.getVal());
		hnL -> SetBinError(ievent+1, orderSup.getError());

		cout << endl
	     << "Tail parameters fixed to the following MC signal values:" << endl;

		tailParams.Print("v");

	}

	/// Draw graphs with legend

	/// alpha graph
	TCanvas* ca = new TCanvas("ca", "ca", 600, 600);
	// gStyle -> SetPadRightMargin(0.2);
	haH  -> GetYaxis() -> SetRangeUser(0, 4);
	haH -> SetLineColor(kRed-4);
	haH -> SetMarkerColor(kRed-4);
	haH -> SetMarkerSize(1);
	haH -> SetMarkerStyle(20);	
	haH -> Draw("P");
	haL -> SetLineColor(kAzure+2);
	haL -> SetMarkerColor(kAzure+2);
	haL -> SetMarkerSize(1);
	haL -> SetMarkerStyle(20);
	haL -> Draw("SAME P");

	/// Legend
	TLegend *lega = new TLegend(0.6, 0.6, 0.8, 0.8);
	lega -> SetBorderSize(0);
	lega -> SetTextSize(0.04);
	lega -> AddEntry(haH->GetName(), "#alpha_{H}", "lPE");
	lega -> AddEntry(haL->GetName(), "#alpha_{L}", "lPE");
	lega -> Draw();

	ca -> SaveAs("./SignalParameters/Alpha_pT_dependence.png");

	/// n graph
	TCanvas* cn = new TCanvas("cn", "cn", 600, 600);
	hnH -> GetYaxis() -> SetRangeUser(0, 25);
	hnH -> SetLineColor(kBlue-4);
	hnH -> SetMarkerColor(kBlue-4);
	hnH -> SetMarkerSize(1);
	hnH -> SetMarkerStyle(20);
	hnH -> Draw("P");
	hnL -> SetLineColor(kGreen+3);
	hnL -> SetMarkerColor(kGreen+3);
	hnL -> SetMarkerSize(1);
	hnL -> SetMarkerStyle(20);
	hnL -> Draw("SAME P");
	
	/// Legend
	TLegend *legn = new TLegend(0.6, 0.6, 0.8, 0.8);
	legn -> SetBorderSize(0);
	legn -> SetTextSize(0.04);
	legn -> AddEntry(hnH->GetName(), "n_{H}", "lPE");
	legn -> AddEntry(hnL->GetName(), "n_{L}", "lPE");
	legn -> Draw();

	cn -> SaveAs("./SignalParameters/n_pT_dependence.png");
}