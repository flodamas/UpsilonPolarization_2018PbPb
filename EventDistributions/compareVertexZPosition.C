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

	gStyle->SetLabelSize(0.045, "XYZ");
	gStyle->SetTitleSize(0.055, "XYZ");

	gStyle->SetTitleXOffset(1.);
	gStyle->SetTitleYOffset(1.5);

	auto canvas = new TCanvas("canvas", "", 700, 700);

	auto distribMC = VertexZDistribution("OniaTree_Y1S_pThat2_HydjetDrumMB_miniAOD.root");
	distribMC->SetLineColor(kRed + 1);
	distribMC->SetMarkerColor(kRed + 1);
	distribMC->SetLineWidth(2);

	auto distribData = VertexZDistribution("OniaTree_miniAOD_PbPbPrompt_112X_DATA_ep.root");
	distribData->SetMarkerColor(kAzure + 1);
	distribData->SetLineWidth(2);
	distribData->SetLineColor(kAzure + 1);

	auto ratioPlot = new TRatioPlot(distribData, distribMC);
	//ratioPlot->SetTicks(1, 1);
	ratioPlot->SetSeparationMargin(0.03);

	std::vector<double> lines = {-3, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3};
	ratioPlot->SetGridlines(lines);
	ratioPlot->SetGraphDrawOpt("PZ");
	ratioPlot->Draw("APZ");

	ratioPlot->GetLowYaxis()->SetNdivisions(505);
	ratioPlot->GetLowerRefYaxis()->SetTitle("Data / MC");
	ratioPlot->SetUpTopMargin(.08);
	ratioPlot->SetLeftMargin(.16);
	ratioPlot->SetRightMargin(.02);
	ratioPlot->SetLowBottomMargin(.4);
	ratioPlot->GetLowerRefGraph()->SetMinimum(0.2);
	ratioPlot->GetLowerRefGraph()->SetMaximum(1.8);
	ratioPlot->GetLowerRefGraph()->SetLineWidth(2);

	ratioPlot->GetUpperPad()->cd();
	((TH1D*)ratioPlot->GetUpperRefObject())->SetMaximum(0.054);

	TLegend* legendData = new TLegend(.2, .85, .45, .6);
	legendData->SetTextSize(.065);

	legendData->AddEntry(distribData, "Data", "l");
	legendData->AddEntry((TObject*)0, Form("#mu = %.2f cm", distribData->GetMean()), " ");
	legendData->AddEntry((TObject*)0, Form("#sigma = %.2f cm", distribData->GetStdDev()), " ");
	legendData->Draw();

	TLegend* legendMC = new TLegend(.65, .85, .9, .6);
	legendMC->SetTextSize(.065);
	legendMC->AddEntry(distribMC, "MC", "lep");
	legendMC->AddEntry((TObject*)0, Form("#mu = %.2f cm", distribMC->GetMean()), " ");
	legendMC->AddEntry((TObject*)0, Form("#sigma = %.2f cm", distribMC->GetStdDev()), " ");

	legendMC->Draw();

	canvas->cd();
	CMS_lumi(canvas, lumi_PbPb);

	gPad->Update();

	canvas->SaveAs("plots/zVertexPosition.pdf", "RECREATE");

	// print the data / MC factors
	auto ratioGraph = ratioPlot->GetLowerRefGraph();

	auto yValues = ratioGraph->GetY();

	for (int i = 0; i < ratioGraph->GetN(); i++) cout << yValues[i] << ", ";
}
