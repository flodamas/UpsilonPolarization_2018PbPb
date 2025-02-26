// when the header is included several times, to avoid the redefinition error
#ifndef LEGENDS
#define LEGENDS

#include "TPaveText.h"

#include "RooRealVar.h"
#include "../../AnalysisParameters.h"

TPaveText* KinematicsText(Int_t centMin, Int_t centMax, Int_t ptMin, Int_t ptMax, float x1 = 0.57, float y1 = 0.9, float x2 = 0.95, float y2 = 0.6) {
	TPaveText* text = new TPaveText(x1, y1, x2, y2, "NDCNB");
	// TPaveText* text = new TPaveText(0.65, 0.90, 0.95, 0.60, "NDCNB");

	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(CentralityRangeText(centMin, centMax));
	text->AddText(gMuonPtCutText);
	text->AddText(DimuonRapidityRangeText(gRapidityMin, gRapidityMax));
	text->AddText(DimuonPtRangeText(ptMin, ptMax));

	text->SetAllWith("", "align", 32);
	return text;
}

TPaveText* RefFrameText(Bool_t isCSframe = true, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180, float x1 = 0.61, float y1 = 0.34, float x2 = 0.97, float y2 = 0.56) {
	TPaveText* text = new TPaveText(x1, y1, x2, y2, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	// text->AddText(Form("%d < p_{T}^{#mu#mu} < %d GeV/c", ptMin, ptMax));
	text->AddText(isCSframe ? "Collins-Soper frame" : "Helicity frame");
	text->AddText(CosThetaRangeText(isCSframe ? "CS" : "HX", cosThetaMin, cosThetaMax));
	text->AddText(PhiRangeText(isCSframe ? "CS" : "HX", phiMin, phiMax));

	text->SetAllWith("", "align", 32);
	return text;
}

TPaveText* RefFrameTextPhiFolded(Bool_t isCSframe = true, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180, float x1 = 0.14, float y1 = 0.65, float x2 = 0.50, float y2 = 0.90, int alignment = 12) {
	TPaveText* text = new TPaveText(x1, y1, x2, y2, "NDCNB"); // on the left side
	// TPaveText* text = new TPaveText(0.61, 0.34, 0.97, 0.56, "NDCNB"); // on the right side
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	// text->AddText(Form("%d < p_{T}^{#mu#mu} < %d GeV/c", ptMin, ptMax));
	text->AddText(isCSframe ? "Collins-Soper frame" : "Helicity frame");
	text->AddText(CosThetaRangeText(isCSframe ? "CS" : "HX", cosThetaMin, cosThetaMax));
	text->AddText(AbsPhiRangeText(isCSframe ? "CS" : "HX", phiMin, phiMax));

	text->SetAllWith("", "align", alignment); // on the left side
	// text->SetAllWith("", "align", 32); // on the right side
	return text;
}

TPaveText* FitResultText(RooRealVar n1S, Float_t signif1S, RooRealVar n2S, Float_t signif2S /*, RooRealVar nBkg*/, float x1 = 0.57, float y1 = 0.35, float x2 = 0.95, float y2 = 0.08, int alignment = 12) {
	// TPaveText* text = new TPaveText(0.6, 0.85, 0.95, 0.5, "NDCNB");
	TPaveText* text = new TPaveText(x1, y1, x2, y2, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	if (DoAsymptoticError) {
		text->AddText(Form("N(#varUpsilon(1S)) = %.0f #pm %.0f", n1S.getVal(), n1S.getError()));
		text->AddText(Form("S / #sqrt{S+B} (3#sigma) = %.1f", signif1S));
		text->AddText(Form("N(#varUpsilon(2S)) = %.0f #pm %.0f", n2S.getVal(), n2S.getError()));
		text->AddText(Form("S / #sqrt{S+B} (3#sigma) = %.1f", signif2S));
	} else { // assuming Minos is ON, print asymmetric errors
		text->AddText(Form("N(#varUpsilon(1S)) = %.0f^{ #plus%.0f}_{ %.0f}", n1S.getVal(), n1S.getErrorHi(), n1S.getErrorLo()));
		text->AddText(Form("S / #sqrt{S+B} (3#sigma) = %.1f", signif1S));
		text->AddText(Form("N(#varUpsilon(2S)) = %.0f^{ #plus%.0f}_{ %.0f}", n2S.getVal(), n2S.getErrorHi(), n2S.getErrorLo()));
		text->AddText(Form("S / #sqrt{S+B} (3#sigma) = %.1f", signif2S));
	}

	//	text->AddText(Form("N(bkg) = %.0f^{ #plus%.0f}_{ %.0f}", nBkg.getVal(), nBkg.getErrorHi(), nBkg.getErrorLo()));
	text->SetAllWith("", "align", alignment);
	return text;
}

TPaveText* AsymDoubleCBParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar sigmaR, RooRealVar alphaR, RooRealVar orderR) {
	auto* text = new TPaveText(0.15, 0.57, 0.50, 0.1, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Asymmetrical double-sided CB");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma_{L} = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha_{L} = %.3f #pm %.3f", alpha.getVal(), alpha.getError()));
	text->AddText(Form("n_{L} = %.3f #pm %.3f", order.getVal(), order.getError()));
	text->AddText(Form("#sigma_{H} = %.2f #pm %.2f MeV", 1000 * sigmaR.getVal(), 1000 * sigmaR.getError()));
	text->AddText(Form("#alpha_{H} = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n_{H} = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TPaveText* SymCoreDoubleCBParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar alphaR, RooRealVar orderR) {
	auto* text = new TPaveText(0.16, 0.8, 0.52, 0.4, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Symmetric double-sided CB");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha_{L} = %.3f #pm %.3f", alpha.getVal(), alpha.getError()));
	text->AddText(Form("n_{L}= %.3f #pm %.3f", order.getVal(), order.getError()));
	text->AddText(Form("#alpha_{H} = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n_{H} = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TPaveText* SymCoreDoubleCBDataParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar alphaR, RooRealVar orderR, RooRealVar turnon, RooRealVar err_sigma, RooRealVar exp_lambda) {
	auto* text = new TPaveText(0.15, 0.53, 0.46, 0.1, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("DSCB          ");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha_{L} = %.3f (fixed)", alpha.getVal()));
	text->AddText(Form("n_{L}= %.3f (fixed)", order.getVal()));
	text->AddText(Form("#alpha_{H} = %.3f (fixed)", alphaR.getVal()));
	text->AddText(Form("n_{H} = %.3f (fixed)", orderR.getVal()));
	text->AddText(Form("#mu_{err} = %.3f #pm %.3f", turnon.getVal(), turnon.getError()));
	text->AddText(Form("#sigma_{err} = %.3f #pm %.3f", err_sigma.getVal(), err_sigma.getError()));
	text->AddText(Form("#lambda_{exp} = %.3f #pm %.3f", exp_lambda.getVal(), exp_lambda.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TPaveText* SymCoreDoubleCBGaussParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar alphaR, RooRealVar orderR, RooRealVar sigmaG, RooRealVar fraction /*, RooRealVar ratio_sigma*/) {
	auto* text = new TPaveText(0.16, 0.8, 0.52, 0.4, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Symmetric DSCB + Gaussian");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha_{L} = %.3f #pm %.3f", alpha.getVal(), alpha.getError()));
	text->AddText(Form("n_{L}= %.3f #pm %.3f", order.getVal(), order.getError()));
	text->AddText(Form("#alpha_{H} = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n_{H} = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));
	text->AddText(Form("#sigma_{Gauss} = %.2f #pm %.2f MeV", 1000 * sigmaG.getVal(), 1000 * sigmaG.getError()));
	text->AddText(Form("Norm frac = %.2f #pm %.2f", fraction.getVal(), fraction.getError()));
	// text->AddText(Form("#frac{#sigma_{DSCB}}{#sigma_{g}} = %.2f ", ratio_sigma.getVal()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TPaveText* SymCoreDoubleCBGaussDataParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar alphaR, RooRealVar orderR, RooRealVar sigmaG, RooRealVar fraction, RooRealVar turnon, RooRealVar err_sigma, RooRealVar exp_lambda) {
	auto* text = new TPaveText(0.15, 0.56, 0.48, 0.07, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Symmetric DSCB + Gaussian");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	// text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000*sigma.getError()));
	text->AddText(Form("#sigma = %.2f MeV (fixed)", 1000 * sigma.getVal()));
	text->AddText(Form("#alpha_{L} = %.3f (fixed)", alpha.getVal()));
	text->AddText(Form("n_{L}= %.3f (fixed)", order.getVal()));
	text->AddText(Form("#alpha_{H} = %.3f (fixed)", alphaR.getVal()));
	text->AddText(Form("n_{H} = %.3f (fixed)", orderR.getVal()));
	// text->AddText(Form("#sigma_{g} = %.2f #pm %.2f MeV", 1000*sigmaG.getVal(), 1000*sigmaG.getError()));
	// text->AddText(Form("#sigma_{g} = %.2f MeV", 1000*sigmaG.getVal()));
	text->AddText(Form("#sigma_{g} = %.2f MeV (fixed)", 1000 * sigmaG.getVal()));
	text->AddText(Form("Frac_{DSCB, g} = %.2f #pm %.2f", fraction.getVal(), fraction.getError()));
	// text->AddText(Form("Frac_{DSCB, g} = %.2f (fixed)", fraction.getVal()));
	text->AddText(Form("#mu_{err} = %.3f #pm %.3f", turnon.getVal(), turnon.getError()));
	text->AddText(Form("#sigma_{err} = %.3f #pm %.3f", err_sigma.getVal(), err_sigma.getError()));
	text->AddText(Form("#lambda_{exp} = %.3f #pm %.3f", exp_lambda.getVal(), exp_lambda.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TPaveText* HypatiaParamsText(RooRealVar mean, RooRealVar lambda, RooRealVar zeta, RooRealVar beta, RooRealVar sigma, RooRealVar alphaL, RooRealVar orderL, RooRealVar alphaR, RooRealVar orderR) {
	auto* text = new TPaveText(0.16, 0.8, 0.5, 0.2, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Hypatia");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#lambda = %.3f #pm %.3f", lambda.getVal(), lambda.getError()));
	text->AddText(Form("#zeta = (%.2f #pm %.2f) x10^{-3}", 1000 * zeta.getVal(), 1000 * zeta.getError()));
	text->AddText(Form("#beta = %.3f #pm %.3f", beta.getVal(), beta.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha_{L} = %.3f #pm %.3f", alphaL.getVal(), alphaL.getError()));
	text->AddText(Form("n_{L}= %.3f #pm %.3f", orderL.getVal(), orderL.getError()));
	text->AddText(Form("#alpha_{H} = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n_{H} = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));

	text->SetAllWith(" ", "align", 12);

	return text;
}

TPaveText* JohnsonParamsText(RooRealVar mean, RooRealVar lambda, RooRealVar gamma, RooRealVar delta) {
	auto* text = new TPaveText(0.16, 0.8, 0.48, 0.48, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Johnson      ");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#lambda = %.3f #pm %.3f", lambda.getVal(), lambda.getError()));
	text->AddText(Form("#gamma = %.3f #pm %.3f", gamma.getVal(), gamma.getError()));
	text->AddText(Form("#delta = %.3f #pm %.3f", delta.getVal(), delta.getError()));

	text->SetAllWith(" ", "align", 12);

	return text;
}

TPaveText* TwoCBParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar sigmaFraction, RooRealVar normFraction) {
	auto* text = new TPaveText(0.35, 0.45, 0.75, 0.05, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Sum of two CB");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma_{1} = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha = %.3f #pm %.3f", alpha.getVal(), alpha.getError()));
	text->AddText(Form("n = %.3f #pm %.3f", order.getVal(), order.getError()));
	text->AddText(Form("#sigma_{2} / #sigma_{1} = %.3f #pm %.3f", sigmaFraction.getVal(), sigmaFraction.getError()));
	text->AddText(Form("norm factor = %.3f #pm %.3f", normFraction.getVal(), normFraction.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TPaveText* PolarParamsText(double lambThetaIn, double lambPhiIn, RooRealVar normCosThetaFit, RooRealVar lambThetaFit, RooRealVar normPhiFit, RooRealVar lambPhiFit, Bool_t isSecondOrder) {
	auto* text = new TPaveText(0.27, 0.09, 0.86, 0.36, "NDCNB");
	text->SetBorderSize(0);
	text->AddText(isSecondOrder ? Form("#lambda_{#varphi, 2D input} = %.3f", lambPhiIn) : Form("#lambda_{#theta, 2D input} = %.3f", lambThetaIn));
	text->AddText(isSecondOrder ? Form("#lambda_{#varphi, 1D fit} = %.3f #pm %.3f", lambPhiFit.getVal(), lambPhiFit.getError()) : Form("#lambda_{#theta, 1D fit} = %.3f #pm %.3f", lambThetaFit.getVal(), lambThetaFit.getError()));
	text->AddText(isSecondOrder ? Form("n_{1D fit} = %.2e #pm %.2e", normPhiFit.getVal(), normPhiFit.getError()) : Form("n_{1D fit} = %.2e #pm %.1e", normCosThetaFit.getVal(), normCosThetaFit.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TLegend* LegendBlock(Float_t x, Float_t y, int nEntries = 1) {
	TLegend* block = new TLegend(x, y, x + 0.25, y - nEntries * .05);

	block->SetTextSize(.05);

	return block;
}

void drawSingleLegend(TObject* graph, const char* text, Float_t x, Float_t y = .84, const char* opt = "p") {
	TLegend* legend = new TLegend(x, y + 0.03, x + 0.3, y - 0.03);
	// legend->SetTextAlign(12);
	legend->SetTextSize(.05);

	legend->AddEntry(graph, text, opt);

	legend->Draw();
}

#endif
