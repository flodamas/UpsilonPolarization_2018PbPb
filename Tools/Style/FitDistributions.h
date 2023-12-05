#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooAbsRealLValue.h"
#include "RooHist.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH2.h"

using namespace RooFit;

// make and draw the invariant mass distribution with fit results
RooPlot* InvariantMassRooPlot(RooWorkspace& wspace, RooDataSet* dataset) {
	RooPlot* frame = (*wspace.var("mass")).frame(Title(" "), Range(MassBinMin, MassBinMax));
	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad
	dataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"));

	auto* fitModel = wspace.pdf("fitModel");
	fitModel->plotOn(frame, Components(*wspace.pdf("bkgPDF")), LineColor(kGray + 2), LineStyle(kDashed));
	fitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_1S")), LineColor(kRed));
	fitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_2S")), LineColor(kRed));
	fitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_3S")), LineColor(kRed));
	fitModel->plotOn(frame, LineColor(kBlue));

	frame->GetYaxis()->SetMaxDigits(3);

	return frame;
}

// create the pull distribution from the frame where the fit is performed
TPad* GetPadPullDistribution(RooPlot* frame, const int nFitPars) {
	TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .25);
	bottomPad->SetTopMargin(0.015);
	bottomPad->SetBottomMargin(0.4);
	bottomPad->SetTicks(1, 1);
	bottomPad->Draw();
	bottomPad->cd();

	RooPlot* pullFrame = (frame->getPlotVar())->frame(frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
	pullFrame->addPlotable(frame->pullHist(), "PZ");
	pullFrame->SetTitle(" ");
	pullFrame->GetYaxis()->SetTitleOffset(0.35);
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

	pullFrame->SetMaximum(2.6);
	pullFrame->SetMinimum(-2.6);

	// pullFrame->SetMaximum(5.5);
	// pullFrame->SetMinimum(-5.5);

	TLine zeroLine(bottomPad->GetUxmin(), 0, bottomPad->GetUxmax(), 0);
	zeroLine.SetLineStyle(kDashed);
	zeroLine.Draw("SAME");

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.15);
	textChi2.DrawLatexNDC(0.75, 0.12, Form("#chi^{2} / n_{d.o.f.} = %.1f", frame->chiSquare(nFitPars)));

	return bottomPad;
}

// create a new 2D pull distribution
// TH2* GetPullDistribution2D(TH2* h2Ddata, TF2* h2Dfit, const int nFitPars) {
TH2* GetPullDistribution2D(TH2* h2Ddata, TH2* h2Dfit) {
	TH2* h2DPull = (TH2*)h2Dfit->Clone("h2DPull");
	//Fill in the histogram.
	for (int iXbin = 1; iXbin <= h2DPull->GetNbinsX(); iXbin++) {
		for (int iYbin = 1; iYbin <= h2DPull->GetNbinsY(); iYbin++) {
			float dataVal = h2Ddata->GetBinContent(iXbin, iYbin);
			float dataErr = h2Ddata->GetBinError(iXbin, iYbin);
			float fitVal = h2Dfit->GetBinContent(iXbin, iYbin);

			float pullVal = (dataVal - fitVal) / dataErr;
			h2DPull->SetBinContent(iXbin, iYbin, pullVal);
		}
	}

	// RooPlot* pullFrame = (frame->getPlotVar())->frame(frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
	// pullFrame->addPlotable(frame->pullHist(), "PZ");
	// pullFrame->SetTitle(" ");
	// pullFrame->GetYaxis()->SetTitleOffset(0.35);
	// pullFrame->GetYaxis()->SetTitle("Pull");
	// pullFrame->GetYaxis()->SetTitleSize(0.17);
	// pullFrame->GetYaxis()->SetLabelSize(0.15);
	// pullFrame->GetYaxis()->CenterTitle();

	// //pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{#minus}} (GeV)");
	// //pullFrame->GetXaxis()->SetTitleOffset(1.20);
	// //pullFrame->GetXaxis()->SetLabelOffset(0.1);
	// pullFrame->GetXaxis()->SetLabelSize(0.15);
	// pullFrame->GetXaxis()->SetTitleSize(0.17);
	// pullFrame->GetXaxis()->CenterTitle();
	// // pullFrame->GetXaxis()->SetTitleFont(43);
	// // pullFrame->GetYaxis()->SetTitleFont(43);

	// pullFrame->GetYaxis()->SetTickSize(0.03);
	// //pullFrame->GetYaxis()->SetNdivisions(505);
	// pullFrame->GetXaxis()->SetTickSize(0.1);
	// pullFrame->Draw();

	// pullFrame->SetMaximum(3.6);
	// pullFrame->SetMinimum(-3.6);

	// TLatex textChi2;
	// textChi2.SetTextAlign(12);
	// textChi2.SetTextSize(0.15);
	// textChi2.DrawLatexNDC(0.7, 0.85, Form("#chi^{2} / n_{d.o.f.} = %.1f", frame->chiSquare(nFitPars)));

	return h2DPull;
}

TCanvas* DrawMassFitDistributions(RooWorkspace& wspace, RooDataSet* dataset, Int_t nFloatParams = 1, Int_t ptMin = 0, Int_t ptMax = 30) {
	// one pad for the invariant mass data distribution with fit components, one for the pull distribution
	TCanvas* canvas = new TCanvas("massCanvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
	pad1->Draw();
	pad1->cd();

	RooPlot* frame = InvariantMassRooPlot(wspace, dataset);

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	//frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->addObject(FitResultText(*wspace.var("yield1S"), ComputeSignalSignificance(wspace, 1), *wspace.var("yield2S"), ComputeSignalSignificance(wspace, 2)));

	frame->Draw();

	gPad->RedrawAxis();

	// pull distribution
	canvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, nFloatParams);

	//canvas->Modified();
	//canvas->Update();
	canvas->cd();
	pad1->Draw();
	pad2->Draw();

	return canvas; // to be saved outside
}

void SaveMCSignalTailParameters(RooArgSet* params, const char* outputName) {
	params->writeToFile(Form("../MonteCarlo/SignalParameters/%s.txt", outputName));
}

void SaveMCSignalParameters(RooArgSet* params, const char* outputName) {
	params->writeToFile(Form("../MonteCarlo/SignalParameters/%s.txt", outputName));
}
