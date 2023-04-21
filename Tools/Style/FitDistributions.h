
// create the pull distribution from the frame where the fit is performed
TPad* GetPadPullDistribution(RooPlot* frame, const int nFitPars) {
	TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .25);
	bottomPad->SetTopMargin(0.015);
	bottomPad->SetBottomMargin(0.34);
	bottomPad->SetTicks(1, 1);
	bottomPad->Draw();
	bottomPad->cd();

	RooPlot* pullFrame = (frame->getPlotVar())->frame(frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
	pullFrame->addPlotable(frame->pullHist(), "PZ");
	pullFrame->SetTitle(" ");
	pullFrame->GetYaxis()->SetTitleOffset(0.4);
	pullFrame->GetYaxis()->SetTitle("Pull");
	pullFrame->GetYaxis()->SetTitleSize(0.17);
	pullFrame->GetYaxis()->SetLabelSize(0.15);
	pullFrame->GetYaxis()->CenterTitle();

	//pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{#minus}} (GeV)");
	//pullFrame->GetXaxis()->SetTitleOffset(1.20);
	//pullFrame->GetXaxis()->SetLabelOffset(0.1);
	pullFrame->GetXaxis()->SetLabelSize(0.15);
	pullFrame->GetXaxis()->SetTitleSize(0.17);
	pullFrame->GetXaxis()->CenterTitle();
	// pullFrame->GetXaxis()->SetTitleFont(43);
	// pullFrame->GetYaxis()->SetTitleFont(43);

	pullFrame->GetYaxis()->SetTickSize(0.03);
	//pullFrame->GetYaxis()->SetNdivisions(505);
	pullFrame->GetXaxis()->SetTickSize(0.1);
	pullFrame->Draw();

	pullFrame->SetMaximum(3.6);
	pullFrame->SetMinimum(-3.6);

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.15);
	textChi2.DrawLatexNDC(0.7, 0.85, Form("#chi^{2} / n_{d.o.f.} = %.1f", frame->chiSquare(nFitPars)));

	return bottomPad;
}

void SaveMCSignalTailParameters(RooArgSet* params, const char* outputName) {
	params->writeToFile(Form("../MonteCarlo/SignalParameters/%s.txt", outputName));
}

RooArgSet* GetMCSignalParameters(const char* mcFileName) {
	RooRealVar alphaInf("alphaInf", "", 1);
	RooRealVar orderInf("orderInf", "", 1);
	RooRealVar alphaSup("alphaSup", "", 1);
	RooRealVar orderSup("orderSup", "", 1);

	// read the values from the corresponding txt files -> mind the ordering!!
	RooArgSet* tailParams = new RooArgSet(alphaInf, orderInf, alphaSup, orderSup);

	tailParams->readFromFile(Form("../MonteCarlo/SignalParameters/%s.txt", mcFileName));
	cout << endl
	     << "Tail parameters fixed to the following MC signal values:" << endl;
	tailParams->Print("v");

	// fix the tail parameters
	alphaInf.setConstant();
	orderInf.setConstant();
	alphaSup.setConstant();
	orderSup.setConstant();

	return tailParams;
}
