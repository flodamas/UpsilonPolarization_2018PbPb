#include "../AnalysisParameters.h"

// plot the signal raw yield distribution for reco MC and data and fit the ratio to reweight the MC reco spectra

void nRecoSpectra() {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	// Raw data spectra
	Double_t nSignal1S_absy0to1p2[] = {794, 1961, 2009, 1898, 3442, 2597, 1946, 1500, 1349, 512, 346};
	Double_t statSignal1S_absy0to1p2[] = {76, 123, 147, 125, 176, 208, 128, 90, 64, 33, 26};

	const char* histoTitle = Form(";%s;dN(#varUpsilon(1S)) / d%s (%s)^{-1}", gPtAxisTitle, gPtVarName, gPtUnit);

	TH1D* hRawData_absy0to1p2 = new TH1D("hRawData_absy0to1p2", histoTitle, NPtFineBins, gPtFineBinning);

	Double_t nSignal1S_absy1p2to2p4[] = {824, 1988, 2549, 2436, 4206, 2513, 1244, 747, 808, 330, 192};
	Double_t statSignal1S_absy1p2to2p4[] = {70, 106, 127, 153, 190, 147, 146, 61, 53, 26, 22};

	TH1D* hRawData_absy1p2to2p4 = new TH1D("hRawData_absy1p2to2p4", histoTitle, NPtFineBins, gPtFineBinning);

	TFile* HEPDataFile = TFile::Open("./HEPData-ins1674529-v2-Table_1.root", "READ");
    if (!HEPDataFile || HEPDataFile->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }

    TDirectory* dir = (TDirectory*)HEPDataFile->Get("Table 1");
    if (!dir) {
        std::cerr << "Error: Could not find directory 'Table 1'!" << std::endl;
        HEPDataFile->Close();
        return;
    }
    
    // Retrieve the histogram from the directory
    TH1D* hppXSec = (TH1D*)dir->Get("Hist1D_y1");
    if (!hppXSec) {
        std::cerr << "Error: Could not find histogram 'Hist1D_y1'!" << std::endl;
        HEPDataFile->Close();
        return;
    }

	for (int i = 0; i < NPtFineBins; i++) {
		hRawData_absy0to1p2->SetBinContent(i + 1, nSignal1S_absy0to1p2[i]);
		hRawData_absy0to1p2->SetBinError(i + 1, statSignal1S_absy0to1p2[i]);

		hRawData_absy1p2to2p4->SetBinContent(i + 1, nSignal1S_absy1p2to2p4[i]);
		hRawData_absy1p2to2p4->SetBinError(i + 1, statSignal1S_absy1p2to2p4[i]);
	}

	hRawData_absy0to1p2->Scale(1. / hRawData_absy0to1p2->Integral(), "width");

	hRawData_absy1p2to2p4->Scale(1. / hRawData_absy1p2to2p4->Integral(), "width");

	// MC reco distribution
	TFile* fileMC = TFile::Open(Form("../Files/Y1SReconstructedMCDataset%s.root", gMuonAccName), "READ");

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
	// hRawData_absy0to1p2->GetYaxis()->SetRangeUser(0, 0.3);
	hRawData_absy1p2to2p4->Draw("PZ");

	hRecoMC_absy1p2to2p4->Draw("SAME PZ");

	// hppXSec->Draw("SAME P");
	// few cosmetics

	TPaveText* header = new TPaveText(.3, .85, .9, .75, "NDCNB");
	header->SetFillColor(4000);
	header->SetBorderSize(0);
	header->SetTextSize(.07);
	header->AddText(Form("%s, %s", CentralityRangeText(), DimuonRapidityRangeText(1.2, 2.4)));
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

	hRatio_absy1p2to2p4->GetXaxis()->SetTitle(gPtAxisTitle);
	hRatio_absy1p2to2p4->GetXaxis()->SetLabelSize(0.11);
	hRatio_absy1p2to2p4->GetXaxis()->SetTitleSize(0.13);
	hRatio_absy1p2to2p4->GetXaxis()->SetTickSize(0.06);

	hRatio_absy1p2to2p4->GetYaxis()->SetNdivisions(505);

	hRatio_absy1p2to2p4->Draw("PZ");

	hRatio_absy1p2to2p4->SetMinimum(0.1);
	hRatio_absy1p2to2p4->SetMaximum(2.3);

	// fit the ratio

	//TF1* fitFunc = new TF1("fitFunc", "([0]  + [1]*x*x) / ( x - [2])^3", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);

	// TF1* fitFunc = new TF1("fitFunc", "[0]/([1] + x)", gPtBinning[0], gPtBinning[NPtBins]);
	TF1* fitFunc = new TF1("fitFunc", "[0]/([1] + x)", gPtFineBinning[0], gPtFineBinning[NPtFineBins]);

	auto fitResult = hRatio_absy1p2to2p4->Fit(fitFunc, "QEMSR");
	fitResult->Print("v");

	// legend with fit result info

	TLegend* fitLegend = new TLegend(.4, .9, .7, .65);
	fitLegend->SetTextSize(.09);
	fitLegend->AddEntry(fitFunc, Form("#frac{A}{(p_{T} + B)}  fit (#chi^{2} / n_{dof} = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");
	//	fitLegend->AddEntry(fitFunc, Form("#frac{A + B p_{T}^{2}}{(p_{T} - C)^{3}}  fit (#chi^{2} / n_{dof} = %.1f / %d)", fitResult->Chi2(), fitResult->Ndf()), "l");

	fitLegend->Draw();

	TPaveText* fitHeader = new TPaveText(.65, .7, .9, .45, "NDCNB");
	fitHeader->SetFillColor(4000);
	fitHeader->SetBorderSize(0);
	fitHeader->SetTextSize(.09);
	fitHeader->AddText(Form("A = %.1f #pm %.1f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
	fitHeader->AddText(Form("B = %.1f #pm %.1f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
	//fitHeader->AddText(Form("C = %.1f #pm %.1f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));

	fitHeader->SetAllWith("", "align", 12);
	fitHeader->Draw();

	canvas->Modified();
	canvas->Update();
	canvas->cd();

	pad1->Draw();

	pad2->Draw();

	CMS_lumi(pad1, gCMSLumiText);

	canvas->SaveAs("plots/recoPtSpectra_absy1p2to2p4_bis.png", "RECREATE");
	/*
	// print the data / MC factors
	auto ratioGraph = ratioPlot->GetLowerRefGraph();

	auto yValues = ratioGraph->GetY();

	for (int i = 0; i < ratioGraph->GetN(); i++) cout << yValues[i] << ", ";
	*/
}
