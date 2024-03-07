// #include "../Tools/Style/tdrStyle.C"
// #include "../Tools/Style/CMS_lumi.C"

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

void plotWeightedInvMass(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = 0, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/WeightedUpsilonSkimmedDataset.root";

	TFile* f = TFile::Open(filename, "READ");
	if (!f) {
		cout << "File " << filename << " not found. Check the directory of the file." << endl;
		return;
	}

	cout << "File " << filename << " opened" << endl;

	const char* refFrameName = isCSframe ? "CS" : "HX";

	const char* datasetName = Form("dataset%s", refFrameName);
	RooDataSet* allDataset = (RooDataSet*)f->Get(datasetName);

	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooDataSet* reducedDataset = InvMassCosThetaPhiDataset(allDataset, wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	Long64_t nEntries = reducedDataset->sumEntries();

	auto* canvas = new TCanvas("canvas", "", 650, 600);

	RooPlot* frame = wspace.var("mass")->frame(Title(" "), Range(MassBinMin, MassBinMax));
	reducedDataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"), DataError(RooAbsData::SumW2));

	frame->GetYaxis()->SetMaxDigits(3);

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->Draw();

	//frame->SetMaximum(nEntries / 150);
	//frame->SetMinimum(0.8);

	//CMS_lumi(canvas, "2018 PbPb miniAOD, DoubleMuon PD");

	canvas->Modified();
	canvas->Update();

	gSystem->mkdir("mass_distrib", kTRUE);
	canvas->SaveAs(Form("mass_distrib/%s_cent%dto%d_pt%dto%dGeV_cosTheta%.1fto%.1f_phi%.1dto%.1d_test.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax), "RECREATE");
}

void scanPlotWeightedInvMass() {
	Int_t ptEdges[9] = {0, 2, 4, 6, 8, 12, 16, 20, 30};
	// Int_t ptEdges[4] = {2, 4, 6, 8};
	Float_t cosThetaEdges[11] = {-1., -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.};
	Int_t phiEdges[7] = {-180, -120, -60, 0, 60, 120, 180};
	// Int_t phiEdges[13] = {0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180};
	Int_t numPtEle = sizeof(ptEdges) / sizeof(Int_t);
	Int_t numCosThetaEle = sizeof(cosThetaEdges) / sizeof(Float_t);
	Int_t numPhiEle = sizeof(phiEdges) / sizeof(Int_t);
	for (Int_t ptIdx = 0; ptIdx < numPtEle - 1; ptIdx++) {
		for (Int_t cosThetaIdx = 0; cosThetaIdx < numCosThetaEle - 1; cosThetaIdx++) {
			for (Int_t phiIdx = 0; phiIdx < numPhiEle - 1; phiIdx++) {
				plotWeightedInvMass(ptEdges[ptIdx], ptEdges[ptIdx + 1], kFALSE, cosThetaEdges[cosThetaIdx], cosThetaEdges[cosThetaIdx + 1], phiEdges[phiIdx], phiEdges[phiIdx + 1]);
			}
		}
	}
}
