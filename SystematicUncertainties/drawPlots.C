#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"
#include "../Tools/Style/FitDistributions.h"

#include "../MonteCarlo/AccEffHelpers.h"

TPaveText* KinematicsText_YieldDiff(Int_t centMin, Int_t centMax, Int_t ptMin, Int_t ptMax) {
	TPaveText* text = new TPaveText(0.14, 0.9, 0.49, 0.65, "NDCNB");
	// TPaveText* text = new TPaveText(0.65, 0.90, 0.95, 0.60, "NDCNB");

	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(CentralityRangeText(centMin, centMax));
	text->AddText(gMuonPtCutText);
	text->AddText(DimuonRapidityRangeText(gRapidityMin, gRapidityMax));
	text->AddText(DimuonPtRangeText(ptMin, ptMax));

	text->SetAllWith("", "align", 12);
	return text;
}

TPaveText* RefFrameTextPhiFolded_YieldDiff(Bool_t isCSframe = true, Float_t cosThetaMin = -1, Float_t cosThetaMax = 1, Int_t phiMin = -180, Int_t phiMax = 180) {
	TPaveText* text = new TPaveText(0.14, 0.55, 0.45, 0.35, "NDCNB"); // on the left side
	// TPaveText* text = new TPaveText(0.61, 0.34, 0.97, 0.56, "NDCNB"); // on the right side
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	// text->AddText(Form("%d < p_{T}^{#mu#mu} < %d GeV/c", ptMin, ptMax));
	text->AddText(isCSframe ? "Collins-Soper frame" : "Helicity frame");
	text->AddText(CosThetaRangeText(isCSframe ? "CS" : "HX", cosThetaMin, cosThetaMax));
	text->AddText(AbsPhiRangeText(isCSframe ? "CS" : "HX", phiMin, phiMax));

	text->SetAllWith("", "align", 12); // on the left side
	// text->SetAllWith("", "align", 32); // on the right side
	return text;
}

void CustomizeStatBox(TH1D *hist, double x1 = 0.7, double x2 = 0.9, double y1 = 0.7, double y2 = 0.9) {
    if (!hist) return; // Ensure the histogram exists

    gPad->Update(); // Make sure the stat box is drawn
    TPaveStats *stats = (TPaveStats*)hist->FindObject("stats");
    if (stats) {
        stats->SetName("");        // Remove the title
        stats->SetX1NDC(x1);       // Left x
        stats->SetX2NDC(x2);       // Right x
        stats->SetY1NDC(y1);       // Bottom y
        stats->SetY2NDC(y2);       // Top y
        gPad->Modified();          // Mark the pad as modified
        gPad->Update();            // Update the pad
    }
}

void drawPlots(const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ExpTimesErr", 
                int ptMin = 0, int ptMax = 2, 
                Bool_t isCSframe = false,
                float cosThetaMin = -0.7, float cosThetaMax = -0.42, 
                int phiMin = 0, int phiMax = 60, 
                Long64_t nPseudoExperiments = 10) {
    
    writeExtraText = true; // if extra text
	extraText = "      Internal";

    const char* refFrameName = isCSframe? "CS":"HX";

    /// open files under the each directory (/YieldDifferencePlot_mass7to11p5/signalshape_bkgShape_ptMintoMax_refFrameName_cosThetaMintoMax_phiMintoMax_n/merged_output.root)
    TFile *file = TFile::Open(Form("YieldDifferencePlots_mass7to11p5/%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld/merged_output.root", signalShapeName, bkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments));
    
    if (file) cout << "Opening file: " << Form("YieldDifferencePlots_mass7to11p5/%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld/merged_output.root", signalShapeName, bkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments) << endl;
    else {
        cout << "Error: Failed to open file" << Form("YieldDifferencePlots_mass7to11p5/%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld/merged_output.root", signalShapeName, bkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments) << "!" << endl;
        return;
    }

    /// get histogram (mergedHist)
    TH1D *yield1SdiffHist = (TH1D*)file->Get("mergedHist");

    /// set up cosmetics
    auto* yield1SdiffCanvas = new TCanvas("yield1SdiffCanvas", "", 600, 400);
    yield1SdiffCanvas->SetRightMargin(0.035);

    gStyle->SetOptStat(1111);

    yield1SdiffHist->SetStats(1);

    yield1SdiffHist->SetMarkerStyle(20);  // Set marker style (e.g., 20 = solid circle)
    yield1SdiffHist->SetMarkerSize(1.2);  // Set marker size
    yield1SdiffHist->SetMarkerColor(kBlue);  // Set marker color (e.g., blue)
   
    yield1SdiffHist->SetXTitle("(Input - Fit) yield 1S");
    yield1SdiffHist->SetYTitle("Entries / 10");

    /// draw histogram
    yield1SdiffHist->Draw();

    TPaveText *KinematicsLegend = KinematicsText_YieldDiff(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);
    KinematicsLegend->Draw();

    TPaveText *RefFrameLegendPhiFolded = RefFrameTextPhiFolded_YieldDiff(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);
    RefFrameLegendPhiFolded->Draw();

    CustomizeStatBox(yield1SdiffHist, 0.7, 0.93, 0.73, 0.90);

    yield1SdiffCanvas->Update();  

    /// save under /YieldDifferencePlot_mass7to11p5/
    gSystem->mkdir("YieldDifferencePlots_mass7to11p5", kTRUE);

    yield1SdiffCanvas->SaveAs(Form("YieldDifferencePlots_mass7to11p5/yield1Sdiff_%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld.png", signalShapeName, bkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments));

    return;
}

void multipleDrawPlots(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ExpTimesErr") {

    Int_t nCosThetaBins = 5; 
    Double_t cosThetaMin = -0.7; 
    Double_t cosThetaMax = 0.7;

    Int_t nPhiBins = 3;
    Int_t phiMin = 0; 
    Int_t phiMax = 180; 

    std::vector<Double_t> cosThetaEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	std::vector<Double_t> phiEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

    for (Int_t i = 0; i < nCosThetaBins; ++i) {
        for (Int_t j = 0; j < nPhiBins; ++j) {
            drawPlots(signalShapeName, bkgShapeName, ptMin, ptMax, isCSframe, cosThetaEdges[i], cosThetaEdges[i+1], phiEdges[j], phiEdges[j+1], 10);
        }
    }
}