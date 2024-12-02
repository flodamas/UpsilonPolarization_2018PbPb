// #include "../Tools/Style/tdrStyle.C"
// #include "../Tools/Style/CMS_lumi.C"

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"
#include "../Tools/Datasets/RooDataSetHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

// reduce the whole raw dataset (N dimensions) to (invariant mass, cos theta, phi)
RooDataSet TempReducedDataset(RooWorkspace& wspace, const char* pairSign = "OS", Int_t ptMin = 0, Int_t ptMax = 30) {
	RooDataSet* allDataset = (RooDataSet*)wspace.data(Form("%sdataset", pairSign));

	if (allDataset == nullptr) {
		std::cerr << "Null RooDataSet provided to the reducer method!!" << std::endl;
		return RooDataSet();
	}

	const char* kinematicCut = Form("(centrality >= %d && centrality < %d) && (rapidity > %f && rapidity < %f) && (pt > %d && pt < %d)", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax);

	RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass")), *(wspace.var(CosThetaVarName("CS"))), *(wspace.var(PhiVarName("CS"))), *(wspace.var(CosThetaVarName("HX"))), *(wspace.var(PhiVarName("HX")))), kinematicCut);

	wspace.import(reducedDataset, RooFit::Rename(Form("reduced%sdataset", pairSign)));

	return reducedDataset;
}

void compareDimuonMassShapes(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* filename = "../Files/DimuonSkimmedDataset.root";

	TFile* f = TFile::Open(filename, "READ");

	RooDataSet* datasetOS = (RooDataSet*)f->Get("OSdataset");
	RooDataSet* datasetSS = (RooDataSet*)f->Get("SSdataset");

	RooWorkspace wspace("workspace");

	wspace.import(*datasetOS);
	wspace.import(*datasetSS);

	f->Close();

	RooDataSet reducedDatasetOS = TempReducedDataset(wspace, "OS", ptMin, ptMax);
	RooDataSet reducedDatasetSS = TempReducedDataset(wspace, "SS", ptMin, ptMax);

	//Long64_t nEntries = reducedDataset.sumEntries();

	auto* canvas = new TCanvas("canvas", "", 650, 600);

	RooPlot* frame = wspace.var("mass")->frame(Title(" "), Range(MassBinMin, MassBinMax));
	reducedDatasetOS.plotOn(frame, Name("dataOS"), Binning(NMassBins), DrawOption("P0Z"), MarkerColor(kRed + 1));
	reducedDatasetSS.plotOn(frame, Name("dataSS"), Binning(NMassBins), DrawOption("P0Z"), MarkerColor(kBlue + 1));

	frame->GetYaxis()->SetMaxDigits(3);

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	//frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->Draw();

	//frame->SetMaximum(nEntries / 150);
	//frame->SetMinimum(0.8);

	canvas->Modified();
	canvas->Update();

	CMS_lumi(canvas, gCMSLumiText);

	gSystem->mkdir("BkgShapes", kTRUE);
	canvas->SaveAs(Form("BkgShapes/pt%dto%dGeV_cosTheta%.2fto%.2f_phi%dto%d_%s.png", ptMin, ptMax, cosThetaMin, cosThetaMax, phiMin, phiMax, refFrameName), "RECREATE");
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
				compareDimuonMassShapes(ptEdges[ptIdx], ptEdges[ptIdx + 1], "CS", cosThetaEdges[cosThetaIdx], cosThetaEdges[cosThetaIdx + 1], phiEdges[phiIdx], phiEdges[phiIdx + 1]);
			}
		}
	}
}
