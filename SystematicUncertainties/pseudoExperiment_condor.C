// This code performs pseudo-experiments to calculate the systematic uncertainties from the choice of the background model.
// 1. Generate a pseudo-sample using the nominal signal and background PDFs (DSCB and exp*erf) with the already obtained parameter values. 
// 2. Fit the pseudo-sample with the nominal signal PDF and an alternative background PDF (Chebyshev orJohnson's distribution) to obtain the signal yield.
// 3. Repeat the above steps for a large number of pseudo-samples and check how much it is biased.

// #include "../Tools/BasicHeaders.h"
// #include "BasicHeaders.h"

// #include "AnalysisParameters.h"

// #include "FitShortcuts.h"
// #include "Legends.h"
// #include "FitDistributions.h"

// #include "RooDataSetHelpers.h"
// #include "InvariantMassModels.h"

// #include "AccEffHelpers.h"
// #include "PolarFitHelpers.h"

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../Tools/RooFitPDFs/InvariantMassModels.h"

#include "../MonteCarlo/AccEffHelpers.h"
#include "../Polarization/PolarFitHelpers.h"

using namespace std;

RooDataSet* generatePseudoData(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.42, Float_t cosThetaMax = -0.14, Int_t phiMin = 60, Int_t phiMax = 120, Double_t* yield1SInput = nullptr, Int_t jobIndex = 1){

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* refFrameName = isCSframe ? "CS" : "HX";

    /// define nominal signal and background shape names
    const char* signalShapeName = "SymDSCB";

    const char* bkgShapeName = "ExpTimesErr";

    const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

    RooWorkspace wspace(Form("workspace%s", refFrameName));

    Float_t lowMassCut = 7, highMassCut = 13;
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);

    wspace.import(massVar);

    /// set the tentative limit of the yields variables for the nominal fit model
    Double_t nEntries = 1e6;

	BuildInvariantMassModel(wspace, signalShapeName, bkgShapeName, fitModelName, nEntries, true);

	RooAddPdf invMassModel = *((RooAddPdf*)wspace.pdf("invMassModel"));

    /// signal model parameters
    RooRealVar* mean_1S = (RooRealVar*)wspace.var("mean_1S");

    RooRealVar* yield1S = (RooRealVar*)wspace.var("yield1S");
    RooRealVar* yield2S = (RooRealVar*)wspace.var("yield2S");
    RooRealVar* yield3S = (RooRealVar*)wspace.var("yield3S");

    /// background model parameters
    RooRealVar* err_mu = (RooRealVar*)wspace.var("err_mu");
    RooRealVar* err_sigma = (RooRealVar*)wspace.var("err_sigma");
    RooRealVar* exp_lambda = (RooRealVar*)wspace.var("exp_lambda");

    RooRealVar* yieldBkg = (RooRealVar*)wspace.var("yieldBkg");

    /// set the parameter values using the fit results (now the example is from HX, pT 0 to 2, cosTheta -0.42 to -0.14, phi 60 to 120)
    /********************* enter the values using the already obtained fit results *******************************/
 
    /// HX frame, pT 0 to 2, cosTheta -0.42 to -0.14, phi 60 to 120
    // mean_1S->setVal(9.4572);
    // yield1S->setVal(5.0721e2);
    // yield2S->setVal(2.1923e1);
    // yield3S->setVal(3.9287e1);

    /// HX frame, pT 0 to 2, cosTheta -0.42 to -0.14, phi 60 to 120
    // mean_1S->setVal(9.4572);
    // yield1S->setVal(5.0721e2);
    // yield2S->setVal(2.1923e1);
    // yield3S->setVal(3.9287e1);

    // err_mu->setVal(7.2495e0);
    // err_sigma->setVal(9.1318e-1);
    // exp_lambda->setVal(1.5132e0);
    // yieldBkg->setVal(1.1140e4);

    /// HX frame, pT 2 to 6, cosTheta -0.42 to -0.14, phi 60 to 120
    mean_1S->setVal(9.4422);
    yield1S->setVal(1.8837e3);
    yield2S->setVal(2.1522e2);
    yield3S->setVal(7.0108e1);

    err_mu->setVal(7.1678e0);
    err_sigma->setVal(8.9041e-1);
    exp_lambda->setVal(1.5231e0);
    yieldBkg->setVal(3.8490e4);

    /************************************************************************************************************/

    /// record the input yield1S and pass it to the main function
    *yield1SInput = (Double_t)(yield1S->getVal());
    
    /// print out invMassModel and wspace for confirmation
    // invMassModel.Print();

    // wspace.Print();

    // ((RooRealVar*)wspace.var("mean_1S"))->Print();
    // invMassModel.Print("t");

    /// get the nominal fit model from the workspace
	auto* nominalFitModel = wspace.pdf("invMassModel");

    /// calculate the total yield
    double yieldTot = yield1S->getVal() + yield2S->getVal() + yield3S->getVal() + yieldBkg->getVal();

    /// radoimize the seed
    RooRandom::randomGenerator()->SetSeed(time(0) + jobIndex * 1000);

    /// generate the pseudo-data
    RooDataSet* pseudoData = nominalFitModel->generate(*wspace.var("mass"), yieldTot);

    /// draw the generated pseudo-data
    // auto* pseudoDataCanvas = new TCanvas("pseudoDataCanvas", "", 600, 600);

	// RooPlot* frame = (*wspace.var("mass")).frame(Title(" "), Range("MassFitRange"));

    // pseudoData->plotOn(frame, Name("data"), Binning(75), DrawOption("P0Z"));

	// nominalFitModel->plotOn(frame, Components(*wspace.pdf("bkgPDF")), LineColor(gColorBkg), LineStyle(kDashed), Range("MassFitRange"), NormRange("MassFitRange"));
	// nominalFitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_1S")), LineColor(gColor1S), Range("MassFitRange"));
	// nominalFitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_2S")), LineColor(gColor2S), Range("MassFitRange"));
	// nominalFitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_3S")), LineColor(gColor3S), Range("MassFitRange"));
	// nominalFitModel->plotOn(frame, LineColor(gColorTotalFit), Range("MassFitRange"), NormRange("MassFitRange"));

	// frame->GetYaxis()->SetMaxDigits(3);
	// gStyle->SetExponentOffset(-0.07, 0.005, "Y");
    
    // frame->Draw();

    return pseudoData;
}

void pseudoExperiment_condor(const char* outputfileName = "yield1Sdiff.root", Int_t jobIndex = 1, Int_t ptMin = 2, Int_t ptMax = 6, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.42, Float_t cosThetaMax = -0.14, Int_t phiMin = 60, Int_t phiMax = 120, Bool_t isPhiFolded = kTRUE, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax){

    /// Start measuring time
	clock_t start, end, cpu_time;
	start = clock();

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    /// define alternative signal and background shape names
    /********************* enter the alternavtive signal or background shape *******************************/

    const char* altSignalShapeName = "SymDSCB";

    const char* altBkgShapeName = "ExpTimesErr";
    // const char* altBkgShapeName = "ChebychevOrder3";

    /*******************************************************************************************************/

	const char* refFrameName = isCSframe ? "CS" : "HX";

    const char* altFitModelName = GetFitModelName(altSignalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

    RooWorkspace wspace(Form("workspace%s", refFrameName));

	Float_t lowMassCut = 6.5, highMassCut = 14.5;
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);

    massVar.setRange("MassFitRange", massMin, massMax); // in this way, the bin width is fixed to 0.08 GeV/c2

    wspace.import(massVar);

    // cout << "wspace print " << endl;
    // wspace.Print(); 

    Double_t nEntries = 1e5;
    
    /// Build the invariant mass model with the alternative signal and background shapes
	BuildInvariantMassModel(wspace, altSignalShapeName, altBkgShapeName, altFitModelName, nEntries, true);

    RooAddPdf invMassModel = *((RooAddPdf*)wspace.pdf("invMassModel"));

    Double_t yield1SInput = 0.;

    Long64_t nPseudoExperiments = 1e2;

    TH1D* yield1Sdiff = new TH1D(Form("yield1Sdiff_%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments), "", 60, -300, 300);

    RooDataSet* pseudoDataset = nullptr;

    RooFitResult* fitResult = nullptr;

    RooPlot* frame = nullptr;

    for (Long64_t iEvent = 0; iEvent < nPseudoExperiments; iEvent++){

        if (iEvent % 1000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, nPseudoExperiments, 100. * iEvent / nPseudoExperiments) << flush;
		}

        /// generate the pseudo-data
        pseudoDataset = generatePseudoData(ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax, &yield1SInput, jobIndex);
        wspace.import(*pseudoDataset);

        /// fit the pseudo-data with the alternative signal and background shapes
        fitResult = RawInvariantMassFit(wspace, *pseudoDataset);

        /// get the fitted yield1S
        RooRealVar* yield1Svar = dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find("yield1S"));

        frame = InvariantMassRooPlot(wspace, *pseudoDataset);

        double chi2 = frame->chiSquare(fitResult->floatParsFinal().getSize());
        // cout << "chi2/Ndf: " << chi2 << endl;

        if (chi2 > 5) {
            cout << "chi2/Ndf > 5, skipping this pseudo-experiment" << endl;
            continue;
        }

        /// fill the histogram with the difference between the input and fitted yield1S
        yield1Sdiff->Fill(yield1SInput - yield1Svar->getVal());
    }

    /// draw the yield difference histogram
    auto* yield1SdiffCanvas = new TCanvas("yield1SdiffCanvas", "", 600, 600);
    yield1SdiffCanvas->SetRightMargin(0.08);

    gStyle->SetOptStat(1111);

    yield1Sdiff->SetStats(1);

    yield1Sdiff->SetMarkerStyle(20);  // Set marker style (e.g., 20 = solid circle)
    yield1Sdiff->SetMarkerSize(1.2);  // Set marker size
    yield1Sdiff->SetMarkerColor(kBlue);  // Set marker color (e.g., blue)
   
    yield1Sdiff->SetXTitle("(Input - Fit) yield 1S");
    yield1Sdiff->SetYTitle("Entries / 10");

    yield1Sdiff->Draw();

    yield1SdiffCanvas->Update();  

    gSystem->mkdir(Form("YieldDifferencePlots/%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld/", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments), kTRUE);
    TFile* yield1SdiffFile = new TFile(Form("YieldDifferencePlots/%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld/%s", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments, outputfileName), "RECREATE");
    // TFile* yield1SdiffFile = new TFile(Form("%s", outputfileName), "RECREATE");
   
    yield1Sdiff->Write();

    yield1SdiffFile->Close();

    yield1SdiffCanvas->SaveAs(Form("YieldDifferencePlots/yield1Sdiff_%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld.png", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments));

    /// draw fitted invariant mass distribution
    auto* massCanvas = new TCanvas("massCanvas", "", 600, 600);
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);

	pad1->Draw();
	pad1->cd();

    /// draw the generated pseudo-data

	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	if (!isPhiFolded)
		frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	else
		frame->addObject(RefFrameTextPhiFolded(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	frame->addObject(FitResultText(*wspace.var("yield1S"), ComputeSignalSignificance(wspace, 1), *wspace.var("yield2S"), ComputeSignalSignificance(wspace, 2)));

	frame->Draw();

	gPad->RedrawAxis();

	// pull distribution
	massCanvas->cd();

	TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	//canvas->Modified();
	//canvas->Update();
	massCanvas->cd();
	pad1->Draw();
	pad2->Draw();

    /// End measuring time
	end = clock();
	cpu_time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << cpu_time / 60. << "minutes" << endl;
}