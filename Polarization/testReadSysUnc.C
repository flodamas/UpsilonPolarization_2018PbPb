#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
//#include "../sPlot/SPlotHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "PolarFitHelpers.h"

#include "../Tools/RooFitPDFs/InvariantMassModels.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../ReferenceFrameTransform/Transformations.h"

double testReadSysUnc(Int_t ptMin = 0, Int_t ptMax = 2, const char* refFrameName = "HX", Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, Int_t phiMin = 0, Int_t phiMax = 180, const char* signalShapeName = "SymDSCB", const char* bkgShapeName = "ExpTimesErr", int nPseudoExperiments = 10) {

    const char* fileLocation = "../SystematicUncertainties/YieldDifferencePlots_mass7to11p5_final/";

    // Open the ROOT file
    TFile *file = openFile(Form("%s%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%d.root", fileLocation, signalShapeName, bkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments));

    // Retrieve the histogram
    TH1 *hist = (TH1*)file->Get("mergedHist"); // Replace "histogram_name" with your actual histogram's name
    if (!hist) {
        std::cerr << "Error: Could not find the histogram in the file!" << std::endl;
        file->Close();
        return;
    }

    // Get the mean value
    double mean = hist->GetMean();
    std::cout << "The mean value of the histogram is: " << mean << std::endl;

    // Close the file
    file->Close();

    return mean;
}