// when the header is included several times, to avoid the redefinition error
#ifndef LEGENDS
#define LEGENDS

#include "TPaveText.h"

#include "RooRealVar.h"
#include "../../AnalysisParameters.h"

TPaveText* KinematicsText(Int_t centMin, Int_t centMax, Int_t ptMin, Int_t ptMax) {
	TPaveText* text = new TPaveText(0.17, 0.9, 0.45, 0.6, "NDCNB");
	// TPaveText* text = new TPaveText(0.65, 0.90, 0.95, 0.60, "NDCNB");

	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(Form("Centrality %d-%d%%", centMin, centMax));
	text->AddText("p_{T}^{#mu} > 3.5 GeV/c");
	text->AddText(Form("%1.1f < |y^{#mu#mu}| < %1.1f", gRapidityMin, gRapidityMax));
	text->AddText(Form("%d < p_{T}^{#mu#mu} < %d GeV/c", ptMin, ptMax));

	text->SetAllWith("", "align", 12);
	return text;
}

TPaveText* RefFrameText(Bool_t isCSframe = true, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = 0, Int_t phiMax = 180) {
	TPaveText* text = new TPaveText(0.63, 0.9, 0.95, 0.7, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	// text->AddText(Form("%d < p_{T}^{#mu#mu} < %d GeV/c", ptMin, ptMax));
	text->AddText(isCSframe ? "Collins-Soper frame" : "Helicity frame");
	text->AddText(Form("%.2f < cos #theta < %.2f", cosThetaMin, cosThetaMax));
	text->AddText(Form("%d#circ < #varphi < %d#circ", phiMin, phiMax));

	text->SetAllWith("", "align", 32);
	return text;
}

TPaveText* FitResultText(RooRealVar n1S, Float_t signif1S, RooRealVar n2S, Float_t signif2S /*, RooRealVar nBkg*/) {
	// TPaveText* text = new TPaveText(0.6, 0.85, 0.95, 0.5, "NDCNB");
	TPaveText* text = new TPaveText(0.6, 0.67, 0.95, 0.42, "NDCNB");
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
	text->SetAllWith("", "align", 32);
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
	auto* text = new TPaveText(0.15, 0.53, 0.50, 0.1, "NDCNB");
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
	auto* text = new TPaveText(0.6, 0.85, 0.93, 0.25, "NDCNB");
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

	text->SetAllWith(" ", "align", 32);
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
	auto* text = new TPaveText(0.15, 0.58, 0.45, 0.07, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Hypatia");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#lambda = %.3f #pm %.3f GeV", lambda.getVal(), lambda.getError()));
	text->AddText(Form("#zeta = (%.2f #pm %.2f) x10^{-3}", 1000 * zeta.getVal(), 1000 * zeta.getError()));
	text->AddText(Form("#beta = %.3f #pm %.3f GeV", beta.getVal(), beta.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha_{L} = %.3f #pm %.3f", alphaL.getVal(), alphaL.getError()));
	text->AddText(Form("n_{L}= %.3f #pm %.3f", orderL.getVal(), orderL.getError()));
	text->AddText(Form("#alpha_{H} = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n_{H} = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));

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

TPaveText* PolarParamsText(double lambThetaIn, double lambPhiIn, RooRealVar lambThetaFit, RooRealVar lambPhiFit, Bool_t isSecondOrder) {
	auto* text = new TPaveText(0.35, 0.35, 0.75, 0.15, "NDCNB");
	text->SetBorderSize(0);
	text->AddText(isSecondOrder ? Form("#lambda_{#phi, 2D input} = %.3f", lambPhiIn) : Form("#lambda_{#theta, 2D input} = %.3f", lambThetaIn));
	text->AddText(isSecondOrder ? Form("#lambda_{#phi, 1D fit} = %.3f #pm %.3f", lambPhiFit.getVal(), lambPhiFit.getError()) : Form("#lambda_{#theta, 1D fit} = %.3f #pm %.3f", lambThetaFit.getVal(), lambThetaFit.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

#endif
