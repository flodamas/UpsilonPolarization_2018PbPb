TPaveText* KinematicsText(Int_t centMin, Int_t centMax, Int_t ptMin, Int_t ptMax) {
	TPaveText* text = new TPaveText(0.17, 0.9, 0.4, 0.65, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(Form("Centrality %d-%d%%", centMin, centMax));
	text->AddText("p_{T}^{#mu} > 3.5 GeV");
	text->AddText("|y^{#mu#mu}| < 2.4");
	text->AddText(Form("%d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));

	text->SetAllWith("", "align", 12);
	return text;
}

TPaveText* RefFrameText(Bool_t isCSframe, Float_t cosThetaMin, Float_t cosThetaMax, Int_t phiMin, Int_t phiMax) {
	TPaveText* text = new TPaveText(0.6, 0.9, 0.95, 0.65, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(isCSframe ? "Collins-Soper frame" : "Helicity frame");
	text->AddText(Form("%.2f < cos #theta < %.2f", cosThetaMin, cosThetaMax));
	text->AddText(Form("%d#circ < #varphi < %d#circ", phiMin, phiMax));

	text->SetAllWith("", "align", 32);
	return text;
}

TPaveText* FitResultText(RooRealVar n1S, Float_t signif1S, RooRealVar n2S, Float_t signif2S) {
	TPaveText* text = new TPaveText(0.6, 0.6, 0.95, 0.3, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(Form("N(#varUpsilon(1S)) = %.0f^{ #plus%.0f}_{ %.0f}", n1S.getVal(), n1S.getErrorHi(), n1S.getErrorLo()));
	text->AddText(Form("S / #sqrt{S+B} (3#sigma) = %.1f", signif1S));
	text->AddText(Form("N(#varUpsilon(2S)) = %.0f^{ #plus%.0f}_{ %.0f}", n2S.getVal(), n2S.getErrorHi(), n2S.getErrorLo()));
	text->AddText(Form("S / #sqrt{S+B} (3#sigma) = %.1f", signif2S));
	text->SetAllWith("", "align", 32);
	return text;
}

TPaveText* AsymDoubleCBParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar sigmaR, RooRealVar alphaR, RooRealVar orderR) {
	auto* text = new TPaveText(0.35, 0.45, 0.75, 0.05, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Asymmetrical double-sided CB");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha = %.3f #pm %.3f", alpha.getVal(), alpha.getError()));
	text->AddText(Form("n = %.3f #pm %.3f", order.getVal(), order.getError()));
	text->AddText(Form("#sigma' = %.2f #pm %.2f MeV", 1000 * sigmaR.getVal(), 1000 * sigmaR.getError()));
	text->AddText(Form("#alpha' = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n' = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TPaveText* SymCoreDoubleCBParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar alphaR, RooRealVar orderR) {
	auto* text = new TPaveText(0.35, 0.45, 0.75, 0.05, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Double-sided CB (symmetric Gaussian core)");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha = %.3f #pm %.3f", alpha.getVal(), alpha.getError()));
	text->AddText(Form("n = %.3f #pm %.3f", order.getVal(), order.getError()));
	text->AddText(Form("#alpha' = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n' = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));

	text->SetAllWith(" ", "align", 12);
	return text;
}

TPaveText* ExtendedCBParamsText(RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar order, RooRealVar alphaR, RooRealVar orderR) {
	auto* text = new TPaveText(0.35, 0.4, 0.7, 0.05, "NDCNB");
	text->SetBorderSize(0);

	text->AddText("Extended Crystal Ball");
	text->AddText(Form("m = %.3f #pm %.3f GeV", mean.getVal(), mean.getError()));
	text->AddText(Form("#sigma = %.2f #pm %.2f MeV", 1000 * sigma.getVal(), 1000 * sigma.getError()));
	text->AddText(Form("#alpha = %.3f #pm %.3f", alpha.getVal(), alpha.getError()));
	text->AddText(Form("n = %.3f #pm %.3f", order.getVal(), order.getError()));
	text->AddText(Form("#alpha' = %.3f #pm %.3f", alphaR.getVal(), alphaR.getError()));
	text->AddText(Form("n' = %.3f #pm %.3f", orderR.getVal(), orderR.getError()));

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
