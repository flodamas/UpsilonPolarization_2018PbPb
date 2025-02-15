#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

std::vector<std::vector<std::vector<double>>> readSystematicUncertainties() {
    std::ifstream inFile("../SystematicUncertainties/systematic_errors_total.txt");
    std::vector<std::vector<std::vector<double>>> sysError(2);
	std::string line;

    if (!inFile) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return sysError;
    }

    /// Skip the header lines
    std::getline(inFile, line); // Skip the first header line

    /// Read the remaining lines
    while (std::getline(inFile, line)) {

        std::istringstream iss(line);
        std::string refFrameName, ptRange;
        double lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde;

        // Parse the line
        iss >> refFrameName >> ptRange >> lambdaTheta >> lambdaPhi >> lambdaThetaPhi >> lambdaTilde;

		cout << "refFrameName: " << refFrameName << endl;
		cout << (refFrameName == "CS") << endl;
		cout << (refFrameName == "HX") << endl;

        // Save the data to the 2D vector
        if (refFrameName == "CS") sysError[0].push_back({lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde});
		else sysError[1].push_back({lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde});
    }

    inFile.close();
    return sysError;
}

void createSysErrorBox(TGraphErrors* gSys) {
	int nPoints = gSys->GetN();

	for (int i = 1; i < nPoints; ++i) {
		double x, y;
		gSys->GetPoint(i, x, y);
		double ex = 0.5;                         // x half-width, as set in SetPointError
		double eyHigh = gSys->GetErrorYhigh(i);    // upper y error
		double eyLow  = gSys->GetErrorYlow(i);     // lower y error

		TBox *box = new TBox(x - ex, y - eyLow, x + ex, y + eyHigh);
		box->SetLineColorAlpha(kAzure-9, 0.6);                // Outline color
		box->SetLineWidth(1);
		box->Draw("same");
	}
}

void finalResults(Bool_t isCSframe = true, const char* bkgShapeName = "ExpTimesErr", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", 0);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", 0);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", 0);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", 0);

	TH1D* lambdaThetaHist = new TH1D("lambdaThetaHist", "; p_{T} (GeV/c); #lambda_{#theta}", NPtBins, gPtBinning);
	TH1D* lambdaPhiHist = new TH1D("lambdaPhiHist", "; p_{T} (GeV/c); #lambda_{#varphi}", NPtBins, gPtBinning);
	TH1D* lambdaThetaPhiHist = new TH1D("lambdaThetaPhiHist", "; p_{T} (GeV/c); #lambda_{#theta#varphi}", NPtBins, gPtBinning);
	TH1D* lambdaTildeHist = new TH1D("lambdaTildeHist", "; p_{T} (GeV/c); #tilde{#lambda}", NPtBins, gPtBinning);

	const char* refFrameName = (isCSframe ? "CS" : "HX");
	const char* signalShapeName = "SymDSCB";
	const char* methodName = "rootFit";
	// const char* bkgShapeName = "ExpTimesErr";
	// const char* bkgShapeName = "Chebychev";

	int refFrameIndex = (isCSframe ? 0 : 1);

	std::vector<std::vector<std::vector<double>>> totalSysError = readSystematicUncertainties();

	TGraphErrors *gSysLambTheta = new TGraphErrors(NPtBins);
	TGraphErrors *gSysLambPhi = new TGraphErrors(NPtBins);
	TGraphErrors *gSysLambThetaPhi = new TGraphErrors(NPtBins);
	TGraphErrors *gSysLambTilde = new TGraphErrors(NPtBins);

	// for (Int_t ibin = 1; ibin <= NPtBins; ibin++) {
	for (Int_t ibin = 2; ibin <= NPtBins; ibin++) {
		// get polarization parameters
		const char* fitModelName = GetFitModelName(signalShapeName, gPtBinning[ibin - 1], gPtBinning[ibin], refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

		RooArgSet polarParams = GetPolarParams(lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde, methodName, fitModelName, bkgShapeName);

		double lambdaThetaVal = lambdaTheta->getVal();
		double lambdaThetaUnc = lambdaTheta->getError();

		double lambdaPhiVal = lambdaPhi->getVal();
		double lambdaPhiUnc = lambdaPhi->getError();

		double lambdaThetaPhiVal = lambdaThetaPhi->getVal();
		double lambdaThetaPhiUnc = lambdaThetaPhi->getError();

		double lambdaTildeVal = lambdaTilde->getVal();
		double lambdaTildeUnc = lambdaTilde->getError();

		cout << "lambdaThetaVal: " << lambdaThetaVal << endl;

		lambdaThetaHist->SetBinContent(ibin, lambdaThetaVal);
		lambdaThetaHist->SetBinError(ibin, lambdaThetaUnc);

		lambdaPhiHist->SetBinContent(ibin, lambdaPhiVal);
		lambdaPhiHist->SetBinError(ibin, lambdaPhiUnc);

		lambdaThetaPhiHist->SetBinContent(ibin, lambdaThetaPhiVal);
		lambdaThetaPhiHist->SetBinError(ibin, lambdaThetaPhiUnc);

		lambdaTildeHist->SetBinContent(ibin, lambdaTildeVal);
		lambdaTildeHist->SetBinError(ibin, lambdaTildeUnc);

		/// Fill the systematic uncertainties
		double x = (gPtBinning[ibin - 1] + gPtBinning[ibin]) / 2;
		double y = lambdaThetaVal;
		double xError = 0.5;

		gSysLambTheta->SetPoint(ibin - 1, x, lambdaThetaVal);
		gSysLambTheta->SetPointError(ibin - 1, xError, totalSysError[refFrameIndex][ibin - 1][0]);

		gSysLambPhi->SetPoint(ibin - 1, x, lambdaPhiVal);
		gSysLambPhi->SetPointError(ibin - 1, xError, totalSysError[refFrameIndex][ibin - 1][1]);

		gSysLambThetaPhi->SetPoint(ibin - 1, x, lambdaThetaPhiVal);
		gSysLambThetaPhi->SetPointError(ibin - 1, xError, totalSysError[refFrameIndex][ibin - 1][2]);

		gSysLambTilde->SetPoint(ibin - 1, x, lambdaTildeVal);
		gSysLambTilde->SetPointError(ibin - 1, xError, totalSysError[refFrameIndex][ibin - 1][3]);
	}

	TCanvas* polarParamsCanvas = new TCanvas("polarParamsCanvas", "polarParamsCanvas", 500, 900);

	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.66, 1, 1.0);

	pad1->SetTopMargin(0.12);
	pad1->SetBottomMargin(0.0);
	pad1->SetLeftMargin(0.17);
	pad1->Draw();
	pad1->cd();

	lambdaThetaHist->Draw("PL");
	lambdaThetaHist->SetMarkerStyle(20);
	lambdaThetaHist->SetMarkerSize(1.4);
	lambdaThetaHist->SetMarkerColor(kAzure+3);
	
	lambdaThetaHist->SetLineColor(kAzure+3);
	lambdaThetaHist->SetLineWidth(3);

	// cosmetics
	lambdaThetaHist->GetXaxis()->CenterTitle();
	lambdaThetaHist->GetXaxis()->SetLabelSize(0);
	lambdaThetaHist->GetXaxis()->SetTitleSize(0);

	lambdaThetaHist->GetYaxis()->SetRangeUser(-1.2, 1.2);
	lambdaThetaHist->GetYaxis()->CenterTitle();
	lambdaThetaHist->GetYaxis()->SetTitleSize(0.1);
	lambdaThetaHist->GetYaxis()->SetTitleOffset(0.8);

	lambdaThetaHist->GetYaxis()->SetLabelSize(0.075);

	gSysLambTheta->SetFillColorAlpha(kAzure-9, 0.3);

	gSysLambTheta->Draw("E2 SAME");
	
	createSysErrorBox(gSysLambTheta);

	TLatex text;
	text.SetTextAlign(22);
	text.SetTextSize(0.08);
	text.DrawLatexNDC(.31, .8, "#varUpsilon(1S) #rightarrow #mu^{+}#mu^{-}");

	TLatex refFrameText;
	refFrameText.SetTextAlign(22);
	refFrameText.SetTextSize(0.08);
	if (isCSframe)
		refFrameText.DrawLatexNDC(0.82, 0.8, "Collins-Soper");
	else
		refFrameText.DrawLatexNDC(0.86, 0.8, "Helicity");

	TLine* zeroLine = new TLine(gPtBinning[0], 0, gPtBinning[NPtBins], 0);
	zeroLine->Draw("SAME");
	zeroLine->SetLineStyle(kDashed);

	CMS_lumi(pad1, gCMSLumiText);

	polarParamsCanvas->cd();

	TPad* pad2 = new TPad("pad2", "pad2", 0, 0.35, 1, 0.66);

	pad2->SetTopMargin(0.0);
	pad2->SetBottomMargin(0.0);
	pad2->SetLeftMargin(0.17);
	pad2->Draw();
	pad2->cd();

	lambdaPhiHist->Draw();
	lambdaPhiHist->SetMarkerStyle(20);
	lambdaPhiHist->SetMarkerSize(1.4);
	lambdaPhiHist->SetMarkerColor(kAzure+3);
	
	lambdaPhiHist->SetLineColor(kAzure+3);
	lambdaPhiHist->SetLineWidth(3);

	// cosmetics
	lambdaPhiHist->GetXaxis()->CenterTitle();

	lambdaPhiHist->GetYaxis()->SetRangeUser(-1.2, 1.2);
	lambdaPhiHist->GetYaxis()->CenterTitle();
	lambdaPhiHist->GetYaxis()->SetTitleSize(0.11);
	lambdaPhiHist->GetYaxis()->SetTitleOffset(0.7);

	lambdaPhiHist->GetYaxis()->SetLabelSize(0.08);

	gSysLambPhi->SetFillColorAlpha(kAzure-9, 0.3);

	gSysLambPhi->Draw("E2 SAME");
	
	createSysErrorBox(gSysLambPhi);

	TLatex kinematicText;
	kinematicText.SetTextAlign(13);
	kinematicText.SetTextSize(0.085);
	kinematicText.DrawLatexNDC(.25, .85, Form("%s, %s", CentralityRangeText(), DimuonRapidityRangeText(0, 2.4)));

	zeroLine->Draw("SAME");

	polarParamsCanvas->cd();

	TPad* pad3 = new TPad("pad3", "pad3", 0, 0.00, 1, 0.35);

	pad3->SetTopMargin(0.0);
	pad3->SetBottomMargin(0.25);
	pad3->SetLeftMargin(0.17);
	pad3->Draw();
	pad3->cd();

	lambdaThetaPhiHist->Draw();
	lambdaThetaPhiHist->SetMarkerStyle(20);
	lambdaThetaPhiHist->SetMarkerSize(1.4);
	lambdaThetaPhiHist->SetMarkerColor(kAzure+3);
	
	lambdaThetaPhiHist->SetLineColor(kAzure+3);
	lambdaThetaPhiHist->SetLineWidth(3);

	// cosmetics
	lambdaThetaPhiHist->GetXaxis()->CenterTitle();
	lambdaThetaPhiHist->GetXaxis()->SetTitleSize(0.1);

	lambdaThetaPhiHist->GetXaxis()->SetLabelSize(0.075);

	lambdaThetaPhiHist->GetYaxis()->SetRangeUser(-1.2, 1.2);
	lambdaThetaPhiHist->GetYaxis()->CenterTitle();
	lambdaThetaPhiHist->GetYaxis()->SetTitleSize(0.1);
	lambdaThetaPhiHist->GetYaxis()->SetTitleOffset(0.75);

	lambdaThetaPhiHist->GetYaxis()->SetLabelSize(0.075);

	gSysLambThetaPhi->SetFillColorAlpha(kAzure-9, 0.3);

	gSysLambThetaPhi->Draw("E2 SAME");
	
	createSysErrorBox(gSysLambThetaPhi);

	zeroLine->Draw("SAME");

	const char* fitModelName2 = GetFitModelName(signalShapeName, gPtBinning[0], gPtBinning[NPtBins], refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	gSystem->mkdir("ParametersResults", kTRUE);
	polarParamsCanvas->SaveAs(Form("ParametersResults/polarParams_%s_%s_%s.png", methodName, fitModelName2, bkgShapeName), "RECREATE");

	TCanvas* invPolarParamsCanvas = new TCanvas("invPolarParamsCanvas", "invPolarParamsCanvas", 500, 500);

	lambdaTildeHist->Draw();
	
	lambdaTildeHist->SetMarkerStyle(20);
	lambdaTildeHist->SetMarkerSize(1.4);
	lambdaTildeHist->SetMarkerColor(kAzure+3);
	
	lambdaTildeHist->SetLineColor(kAzure+3);
	lambdaTildeHist->SetLineWidth(3);

	// cosmetics
	lambdaTildeHist->GetXaxis()->CenterTitle();

	lambdaTildeHist->GetYaxis()->SetRangeUser(-1.4, 1.4);
	lambdaTildeHist->GetYaxis()->CenterTitle();

	gSysLambTilde->SetFillColorAlpha(kAzure-9, 0.3);

	gSysLambTilde->Draw("E2 SAME");

	createSysErrorBox(gSysLambTilde);
	
	zeroLine->Draw("SAME");

	CMS_lumi(invPolarParamsCanvas, gCMSLumiText);

	text.SetTextSize(0.058);
	text.DrawLatexNDC(.32, .85, "#varUpsilon(1S) #rightarrow #mu^{+}#mu^{-}");

	refFrameText.SetTextSize(0.058);
	if (isCSframe)
		refFrameText.DrawLatexNDC(0.78, 0.85, "Collins-Soper");
	else
		refFrameText.DrawLatexNDC(0.85, 0.85, "Helicity");

	invPolarParamsCanvas->SaveAs(Form("ParametersResults/invariantParameter_%s_%s_%s.png", methodName, fitModelName2, bkgShapeName), "RECREATE");

	return;
}
