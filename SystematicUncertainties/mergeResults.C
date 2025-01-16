/// Merge the result histograms from several files into a single file

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"

#include "../MonteCarlo/AccEffHelpers.h"

// void mergeResults(const char* altSignalShapeName = "SymDSCB", const char* altBkgShapeName = "ChebychevOrder3", Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.70, Float_t cosThetaMax = -0.42, Int_t phiMin = 0, Int_t phiMax = 60, Long64_t nPseudoExperiments = 100) {
void mergeResults(const char* altSignalShapeName = "Johnson", const char* altBkgShapeName = "ExpTimesErr", Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.7, Float_t cosThetaMax = -0.42, Int_t phiMin = 0, Int_t phiMax = 60, Long64_t nPseudoExperiments = 100) { 

    int nFiles = 1000;

    // const char* altSignalShapeName = "SymDSCB";

    // // const char* altBkgShapeName = "ExpTimesErr";
    // const char* altBkgShapeName = "ChebychevOrder3";
    
    const char* refFrameName = isCSframe ? "CS" : "HX";

    TH1D *mergedHist = nullptr;

    char* fileLocation; 

    int startIdx = 0;

    for (int i = 0; i < nFiles; ++i) {

        fileLocation = Form("YieldDifferencePlots_mass7to11p5_final/%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld/", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments);
        // cout << fileLocation << endl;

        TString filename = Form("%soutput_%d.root", fileLocation, i);
        // TString filename = Form("%soutput_%d_cosTheta_%.2f_%.2f_phi_%.2f_%.2f.root", fileLocation, i, cosThetaMin, cosThetaMax, (double)phiMin, (double)phiMax);
        TFile *file = TFile::Open(filename, "READ");

        cout << filename << endl;

        if (!file || file->IsZombie()) {
            std::cerr << "Error: Failed to open file" << filename << "!" << std::endl;
            if (file) file->Close();
            startIdx++;
            continue;
        }

        TH1D *hist = (TH1D*)file->Get(Form("yield1Sdiff_%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments));

        if(i == startIdx){
            mergedHist = (TH1D*)hist->Clone("mergedHist");
        }
        else{
            mergedHist->Add(hist);
        }

        // file->Close(); // it causes the error 
    }

    // Save the merged histogram into a new file
    TFile *outFile = TFile::Open(Form("%s/merged_output.root", fileLocation), "RECREATE");
    
    if (mergedHist) {
        mergedHist->Write();
    } else {
        std::cerr << "Error: mergedHist is null. Cannot write to file." << std::endl;
    }
        
    outFile->Close();
    
    std::cout << "Histograms successfully merged into merged_output.root!" << std::endl;
}

void multipleMergeResults(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ChebychevOrder3") {

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
            mergeResults(signalShapeName, bkgShapeName, ptMin, ptMax, isCSframe, cosThetaEdges[i], cosThetaEdges[i+1], phiEdges[j], phiEdges[j+1], 10);
        }
    }
}