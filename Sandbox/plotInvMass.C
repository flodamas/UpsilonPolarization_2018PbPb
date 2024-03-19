// #include "../Tools/Style/tdrStyle.C"
// #include "../Tools/Style/CMS_lumi.C"

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"
#include "../Tools/Datasets/RooDataSetHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

void plotInvMass(Int_t ptMin = 0, Int_t ptMax = 30, Bool_t isCSframe = kTRUE, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/UpsilonSkimmedDataset.root";

	const char* refFrameName = isCSframe ? "CS" : "HX";

	RooWorkspace wspace = SetUpWorkspace(filename, refFrameName);

	RooDataSet* reducedDataset = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax, refFrameName, phiMin, phiMax);

	Long64_t nEntries = reducedDataset->sumEntries();

	auto* canvas = new TCanvas("canvas", "", 650, 600);

	RooPlot* frame = wspace.var("mass")->frame(Title(" "), Range(MassBinMin, MassBinMax));
	reducedDataset->plotOn(frame, Name("data"), Binning(NMassBins), DrawOption("P0Z"));

	frame->GetYaxis()->SetMaxDigits(3);

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	//frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->Draw();

	//frame->SetMaximum(nEntries / 150);
	//frame->SetMinimum(0.8);

	canvas->Modified();
	canvas->Update();

	CMS_lumi(canvas, gCMSLumiText);

	gSystem->mkdir("mass_distrib", kTRUE);
	canvas->SaveAs(Form("mass_distrib/RawData_%s_cent%dto%d_pt%dto%dGeV_cosTheta%.1fto%.1f_phi%dto%d.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax), "RECREATE");
}

void scanPlotInvMass() {
	Int_t ptEdges[9] = {0, 2, 4, 6, 8, 12, 16, 20, 30};
	// Int_t ptEdges[2] = {2, 4};
	Float_t cosThetaEdges[11] = {-1., -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.};
	Int_t phiEdges[7] = {-180, -120, -60, 0, 60, 120, 180};
	// Int_t phiEdges[2] = {0, 180};
	Int_t numPtEle = sizeof(ptEdges) / sizeof(Int_t);
	Int_t numCosThetaEle = sizeof(cosThetaEdges) / sizeof(Float_t);
	Int_t numPhiEle = sizeof(phiEdges) / sizeof(Int_t);
	for (Int_t ptIdx = 0; ptIdx < numPtEle - 1; ptIdx++) {
		for (Int_t cosThetaIdx = 0; cosThetaIdx < numCosThetaEle - 1; cosThetaIdx++) {
			for (Int_t phiIdx = 0; phiIdx < numPhiEle - 1; phiIdx++) {
				plotInvMass(ptEdges[ptIdx], ptEdges[ptIdx + 1], kFALSE, cosThetaEdges[cosThetaIdx], cosThetaEdges[cosThetaIdx + 1], phiEdges[phiIdx], phiEdges[phiIdx + 1]);
			}
		}
	}
}
