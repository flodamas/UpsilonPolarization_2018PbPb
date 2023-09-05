#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "../AnalysisParameters.h"

TH1D* VertexZDistribution(const char* inputFileName = "OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root") {
	TFile* infile = TFile::Open(Form("../Files/%s", inputFileName), "READ");
	TTree* OniaTree = (TTree*)infile->Get("hionia/myTree");

	/// OniaTree variables
	Float_t zVtx;
	Int_t Centrality;
	ULong64_t HLTriggers;

	OniaTree->SetBranchAddress("zVtx", &zVtx);
	OniaTree->SetBranchAddress("Centrality", &Centrality);
	OniaTree->SetBranchAddress("HLTriggers", &HLTriggers);

	auto histo = new TH1D("histo", ";vertex z (cm);Self-normalized distribution", 2 * 40, -20, 20); // 0.5 cm wide bins
	histo->SetDirectory(0);

	Long64_t totEntries = OniaTree->GetEntries();

	for (Long64_t iEvent = 0; iEvent < (totEntries); iEvent++) {
		if (iEvent % 10000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, totEntries, 100. * iEvent / totEntries) << flush;
		}

		OniaTree->GetEntry(iEvent);

		if (!((HLTriggers & (ULong64_t)(1 << (gUpsilonHLTBit - 1))) == (ULong64_t)(1 << (gUpsilonHLTBit - 1)))) continue; // must fire the upsilon HLT path

		if (Centrality >= 2 * gCentralityBinMax) continue; // discard events with centrality >= 90% in 2018 data

		histo->Fill(zVtx);
	}

	infile->Close();

	histo->Scale(1. / histo->Integral());

	return histo;
}

void compareVertexZPosition() {
	writeExtraText = true; // if extra text
	//extraText = "      Internal";

	auto canvas = new TCanvas("canvas", "", 600, 650);

	auto distribMC = VertexZDistribution("OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root");
	distribMC->SetLineColor(kRed + 1);
	distribMC->SetMarkerColor(kRed + 1);

	auto distribData = VertexZDistribution("OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root");
	distribData->SetLineColor(kAzure + 1);

	auto ratioPlot = new TRatioPlot(distribData, distribMC);
	//ratioPlot->SetTicks(1, 1);
	ratioPlot->SetSeparationMargin(0.03);

	std::vector<double> lines = {-3, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3};
	ratioPlot->SetGridlines(lines);
	ratioPlot->SetGraphDrawOpt("APZ");
	ratioPlot->Draw("APZ");

	ratioPlot->GetLowYaxis()->SetNdivisions(505);
	ratioPlot->GetLowerRefYaxis()->SetTitle("Data / MC");
	ratioPlot->SetUpTopMargin(.08);
	ratioPlot->SetLeftMargin(.18);
	ratioPlot->SetRightMargin(.025);
	ratioPlot->SetLowBottomMargin(.4);
	ratioPlot->GetLowerRefGraph()->SetMinimum(0.2);
	ratioPlot->GetLowerRefGraph()->SetMaximum(1.8);

	ratioPlot->GetUpperPad()->cd();
	TLegend* legend = new TLegend(.72, .85, .92, .65);
	legend->SetTextSize(.07);

	legend->AddEntry(distribData, "Data", "l");
	legend->AddEntry(distribMC, "MC", "lep");

	legend->Draw();

	canvas->cd();
	CMS_lumi(canvas, lumi_PbPb);

	gPad->Update();

	canvas->SaveAs("plots/zVertexPosition.pdf", "RECREATE");

	// print the data / MC factors
	auto ratioGraph = ratioPlot->GetLowerRefGraph();

	auto yValues = ratioGraph->GetY();

	for (int i = 0; i < ratioGraph->GetN(); i++) cout << yValues[i] << ", ";
}
