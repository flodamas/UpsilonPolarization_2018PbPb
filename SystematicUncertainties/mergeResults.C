/// Merge the result histograms from several files into a single file

void mergeResults(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.42, Float_t cosThetaMax = -0.14, Int_t phiMin = 60, Int_t phiMax = 120, Long64_t nPseudoExperiments = 100) {
    
    int nFiles = 1000;

    const char* altSignalShapeName = "SymDSCB";

    const char* altBkgShapeName = "ExpTimesErr";
    // const char* altBkgShapeName = "ChebychevOrder3";
    
    const char* refFrameName = isCSframe ? "CS" : "HX";

    TH1D *mergedHist = nullptr;

    char* fileLocation; 

    for (int i = 0; i < nFiles; ++i) {

        fileLocation = Form("YieldDifferencePlots/%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld/", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments);
        // cout << fileLocation << endl;

        TString filename = Form("%soutput_%d.root", fileLocation, i);
        TFile *file = TFile::Open(filename, "READ");

        cout << filename << endl;

        if (!file || file->IsZombie()) {
            std::cerr << "Error: Failed to open file" << filename << "!" << std::endl;
            if (file) file->Close();
            continue;
        }

        TH1D *hist = (TH1D*)file->Get(Form("yield1Sdiff_%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments));

        if(i == 0){
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