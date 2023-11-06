//#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"

void weightedInvMass(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root";
	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	writeExtraText = true; // if extra text
	extraText = "       Internal";

	Float_t binMin = 7, binMax = 14;
	Int_t nBins = 100;

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	RooDataSet* allDataset = (RooDataSet*)f->Get(isCSframe ? "datasetCS" : "datasetHX");

	RooWorkspace* wspace = new RooWorkspace("workspace");
	wspace->import(*allDataset);

	RooDataSet* massDataset = (isCSframe) ? ReducedWeightedMassDatasetCS(allDataset, wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax) : ReducedWeightedMassDatasetHX(allDataset, wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax);

	Long64_t nEntries = massDataset->sumEntries();

	auto* canvas = new TCanvas("canvas", "", 600, 600);

	RooPlot* frame = wspace->var("mass")->frame(Title(" "), Range(binMin, binMax));
	massDataset->plotOn(frame, Name("data"), Binning(nBins), DrawOption("P0Z"));

	frame->Draw();

	frame->GetYaxis()->SetMaxDigits(3);

	//frame->SetMaximum(nEntries / 150);
	//frame->SetMinimum(0.8);

	TPaveText* pt = new TPaveText(0.5, 0.9, 0.95, 0.7, "NDCNB");
	pt->SetFillColor(4000);
	pt->SetBorderSize(0);
	pt->AddText(Form("Centrality %d-%d%%", gCentralityBinMin, gCentralityBinMax));
	pt->AddText("|#eta^{#mu}| < 2.4, p_{T}^{#mu} > 3.5 GeV");
	pt->AddText(Form("|y^{#mu#mu}| < 2.4, %d < p_{T}^{#mu#mu} < %d GeV", ptMin, ptMax));

	pt->SetAllWith("", "align", 12);
	pt->Draw();

	//CMS_lumi(canvas, "2018 PbPb miniAOD, DoubleMuon PD");

	canvas->Modified();
	canvas->Update();

	canvas->SaveAs(Form("mass_distrib/weighted_cent%dto%d_pt%dto%dGeV.png", gCentralityBinMin, gCentralityBinMax, ptMin, ptMax), "RECREATE");
}
