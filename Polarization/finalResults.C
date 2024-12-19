#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

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

	for (Int_t ibin = 1; ibin <= NPtBins; ibin++) {
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
	}

	TCanvas* polarParamsCanvas = new TCanvas("polarParamsCanvas", "polarParamsCanvas", 500, 900);

	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.66, 1, 1.0);

	pad1->SetTopMargin(0.12);
	pad1->SetBottomMargin(0.0);
	pad1->SetLeftMargin(0.17);
	pad1->Draw();
	pad1->cd();

	lambdaThetaHist->Draw();

	// cosmetics
	lambdaThetaHist->GetXaxis()->CenterTitle();
	lambdaThetaHist->GetXaxis()->SetLabelSize(0);
	lambdaThetaHist->GetXaxis()->SetTitleSize(0);

	lambdaThetaHist->GetYaxis()->SetRangeUser(-1.2, 1.2);
	lambdaThetaHist->GetYaxis()->CenterTitle();
	lambdaThetaHist->GetYaxis()->SetTitleSize(0.1);
	lambdaThetaHist->GetYaxis()->SetTitleOffset(0.8);

	lambdaThetaHist->GetYaxis()->SetLabelSize(0.075);

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

	// cosmetics
	lambdaPhiHist->GetXaxis()->CenterTitle();

	lambdaPhiHist->GetYaxis()->SetRangeUser(-1.2, 1.2);
	lambdaPhiHist->GetYaxis()->CenterTitle();
	lambdaPhiHist->GetYaxis()->SetTitleSize(0.11);
	lambdaPhiHist->GetYaxis()->SetTitleOffset(0.7);

	lambdaPhiHist->GetYaxis()->SetLabelSize(0.08);

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

	// cosmetics
	lambdaThetaPhiHist->GetXaxis()->CenterTitle();
	lambdaThetaPhiHist->GetXaxis()->SetTitleSize(0.1);

	lambdaThetaPhiHist->GetXaxis()->SetLabelSize(0.075);

	lambdaThetaPhiHist->GetYaxis()->SetRangeUser(-1.2, 1.2);
	lambdaThetaPhiHist->GetYaxis()->CenterTitle();
	lambdaThetaPhiHist->GetYaxis()->SetTitleSize(0.1);
	lambdaThetaPhiHist->GetYaxis()->SetTitleOffset(0.75);

	lambdaThetaPhiHist->GetYaxis()->SetLabelSize(0.075);

	zeroLine->Draw("SAME");

	const char* fitModelName2 = GetFitModelName(signalShapeName, gPtBinning[0], gPtBinning[NPtBins], refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	gSystem->mkdir("ParametersResults", kTRUE);
	polarParamsCanvas->SaveAs(Form("ParametersResults/polarParams_%s_%s_%s.png", methodName, fitModelName2, bkgShapeName), "RECREATE");

	TCanvas* invPolarParamsCanvas = new TCanvas("invPolarParamsCanvas", "invPolarParamsCanvas", 500, 500);

	lambdaTildeHist->Draw();

	// cosmetics
	lambdaTildeHist->GetXaxis()->CenterTitle();

	lambdaTildeHist->GetYaxis()->SetRangeUser(-1.4, 1.4);
	lambdaTildeHist->GetYaxis()->CenterTitle();

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
}
