#include "../AnalysisParameters.h"

// plot the signal raw yield distribution for reco MC and data and fit the ratio to reweight the MC reco spectra

void nRecoSpectra() {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	// Raw data spectra
	Double_t nSignal1S_absy0to1p2[] = {797, 1680, 1606, 1348, 2506, 2265, 1831, 1344, 1321, 832};
	Double_t statSignal1S_absy0to1p2[] = {75, 88, 85, 96, 80, 113, 103, 72, 60, 38};

	TH1D* hRawData_absy0to1p2 = new TH1D("hRawData_absy0to1p2", ";p_{T} (GeV/c);dN(#varUpsilon(1S)) / dp_{T} (GeV/c)^{-1}", NPtBins, gPtBinning);

	Double_t nSignal1S_absy1p2to2p4[] = {889, 1980, 2909, 1562, 2232, 1964, 1261, 709, 746, 514};
	Double_t statSignal1S_absy1p2to2p4[] = {71, 123, 197, 140, 62, 95, 89, 65, 48, 35};

	TH1D* hRawData_absy1p2to2p4 = new TH1D("hRawData_absy1p2to2p4", ";p_{T} (GeV/c);dN(#varUpsilon(1S)) / dp_{T} (GeV/c)^{-1}", NPtBins, gPtBinning);

	for (int i = 0; i < NPtBins; i++) {
		hRawData_absy0to1p2->SetBinContent(i + 1, nSignal1S_absy0to1p2[i]);
		hRawData_absy0to1p2->SetBinError(i + 1, statSignal1S_absy0to1p2[i]);

		hRawData_absy1p2to2p4->SetBinContent(i + 1, nSignal1S_absy1p2to2p4[i]);
		hRawData_absy1p2to2p4->SetBinError(i + 1, statSignal1S_absy1p2to2p4[i]);
	}

	hRawData_absy0to1p2->Scale(1. / hRawData_absy0to1p2->Integral(), "width");

	hRawData_absy1p2to2p4->Scale(1. / hRawData_absy1p2to2p4->Integral(), "width");

	// MC reco distribution
	TFile* fileMC = TFile::Open("../Files/MCUpsilonSkimmedWeightedDataset.root", "READ");

	auto* hRecoMC_absy0to1p2 = (TH1D*)fileMC->Get("hReco1S_absy0to1p2");

	auto* hRecoMC_absy1p2to2p4 = (TH1D*)fileMC->Get("hReco1S_absy1p2to2p4");

	// divide the two
	TH1D* hRatio_absy0to1p2 = (TH1D*)hRawData_absy0to1p2->Clone("hRatio_absy0to1p2");

	hRatio_absy0to1p2->Divide(hRecoMC_absy0to1p2);

	TH1D* hRatio_absy1p2to2p4 = (TH1D*)hRawData_absy1p2to2p4->Clone("hRatio_absy1p2to2p4");

	hRatio_absy1p2to2p4->Divide(hRecoMC_absy1p2to2p4);

	gStyle->SetPadLeftMargin(.18);

	auto canvas = new TCanvas("canvas", "", 600, 650);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
	pad1->SetBottomMargin(0.04);
	pad1->SetTopMargin(0.1);
	pad1->Draw();
	pad1->cd();

	hRecoMC_absy1p2to2p4->SetLineColor(kRed + 1);
	hRecoMC_absy1p2to2p4->SetMarkerColor(kRed + 1);
	hRecoMC_absy1p2to2p4->SetLineWidth(2);

	hRawData_absy1p2to2p4->SetLineColor(kAzure + 1);
	hRawData_absy1p2to2p4->SetMarkerColor(kAzure + 1);
	hRawData_absy1p2to2p4->SetLineWidth(2);

	hRawData_absy1p2to2p4->GetXaxis()->SetLabelOffset(1);
	hRawData_absy1p2to2p4->GetYaxis()->SetTitleOffset(1);
	hRawData_absy1p2to2p4->GetYaxis()->SetTitleSize(0.085);
	hRawData_absy1p2to2p4->GetYaxis()->SetLabelSize(0.075);

	hRawData_absy1p2to2p4->Draw("PZ");

	hRecoMC_absy1p2to2p4->Draw("SAME PZ");

	// few cosmetics

	TPaveText* header = new TPaveText(.32, .85, .9, .75, "NDCNB");
	header->SetFillColor(4000);
	header->SetBorderSize(0);
	header->SetTextSize(.07);
	header->AddText("Centrality 0#minus90%, 1.2 < |y| < 2.4");
	header->SetAllWith("", "align", 12);
	header->Draw();

	TLegend* legend = new TLegend(.6, .7, .9, .55);
	legend->SetTextSize(.07);
	legend->AddEntry(hRawData_absy1p2to2p4, "Raw yield", "lep");
	legend->AddEntry(hRecoMC_absy1p2to2p4, "Reco MC", "lep");

	legend->Draw();

	// ratio plot
	canvas->cd();

	TPad* pad2 = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .4);
	pad2->SetTopMargin(0.025);
	pad2->SetBottomMargin(0.3);
	pad2->SetTicks(1, 1);
	pad2->Draw();
	pad2->cd();

	hRatio_absy1p2to2p4->SetTitle(" ");
	hRatio_absy1p2to2p4->GetYaxis()->SetTitleOffset(0.65);
	hRatio_absy1p2to2p4->GetYaxis()->SetTitle("Data / MC");
	hRatio_absy1p2to2p4->GetYaxis()->SetTitleSize(0.13);
	hRatio_absy1p2to2p4->GetYaxis()->SetLabelSize(0.11);
	hRatio_absy1p2to2p4->GetYaxis()->CenterTitle();

	hRatio_absy1p2to2p4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hRatio_absy1p2to2p4->GetXaxis()->SetLabelSize(0.11);
	hRatio_absy1p2to2p4->GetXaxis()->SetTitleSize(0.13);
	hRatio_absy1p2to2p4->GetXaxis()->SetTickSize(0.06);

	hRatio_absy1p2to2p4->GetYaxis()->SetNdivisions(505);

	hRatio_absy1p2to2p4->Draw("PZ");

	hRatio_absy1p2to2p4->SetMinimum(0.2);
	hRatio_absy1p2to2p4->SetMaximum(2.7);

	// fit the ratio

	TF1* fitFunc = new TF1("fitFunc", "([0]  + [1]*x*x) / ( x - [2])^3", gPtBinning[0], gPtBinning[NPtBins]);

	//TF1* fitFunc = new TF1("fitFunc", "([0]  + [1]*x +[2]*x*x) / ( x - [3])^2", gPtBinning[0], gPtBinning[NPtBins]);

	auto fitResult = hRatio_absy1p2to2p4->Fit(fitFunc, "QEMSR");
	fitResult->Print("v");

	// legend with fit result info

	TLegend* fitLegend = new TLegend(.32, .88, .62, .65);
	fitLegend->SetTextSize(.09);
	fitLegend->AddEntry(fitFunc, Form("#frac{A + B p_{T}^{2}}{(p_{T} - C)^{3}}  fit (#chi^{2} / d.o.f = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");

	fitLegend->Draw();

	TPaveText* fitHeader = new TPaveText(.65, .7, .9, .4, "NDCNB");
	fitHeader->SetFillColor(4000);
	fitHeader->SetBorderSize(0);
	fitHeader->SetTextSize(.09);
	fitHeader->AddText(Form("A = %.1f #pm %.1f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
	fitHeader->AddText(Form("B = %.1f #pm %.1f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
	fitHeader->AddText(Form("C = %.1f #pm %.1f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));

	fitHeader->SetAllWith("", "align", 12);
	fitHeader->Draw();

	canvas->Modified();
	canvas->Update();
	canvas->cd();

	pad1->Draw();

	pad2->Draw();

	CMS_lumi(canvas, lumi_PbPb);

	canvas->SaveAs("plots/recoPtSpectra_absy1p2to2p4.pdf", "RECREATE");
	/*
	// print the data / MC factors
	auto ratioGraph = ratioPlot->GetLowerRefGraph();

	auto yValues = ratioGraph->GetY();

	for (int i = 0; i < ratioGraph->GetN(); i++) cout << yValues[i] << ", ";
	*/
}
