// This code performs pseudo-experiments to calculate the systematic uncertainties from the choice of the background model.
// 1. Generate a pseudo-sample using the nominal signal and background PDFs (DSCB and exp*erf) with the already obtained parameter values. 
// 2. Fit the pseudo-sample with the nominal signal PDF and an alternative background PDF (Chebyshev orJohnson's distribution) to obtain the signal yield.
// 3. Repeat the above steps for a large number of pseudo-samples and check how much it is biased.

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"
#include "../Tools/Style/FitDistributions.h"

#include "../Tools/Datasets/RooDataSetHelpers.h"
#include "../Tools/RooFitPDFs/InvariantMassModels.h"

#include "../MonteCarlo/AccEffHelpers.h"
#include "../Polarization/PolarFitHelpers.h"

RooDataSet* generatePseudoData(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.42, Float_t cosThetaMax = -0.14, Int_t phiMin = 60, Int_t phiMax = 120/*, Double_t& yield1SInput = */){

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* refFrameName = isCSframe ? "CS" : "HX";

    const char* signalShapeName = "SymDSCB";

    const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

    const char* bkgShapeName = "ExpTimesErr";

    RooWorkspace wspace(Form("workspace%s", refFrameName));

	Float_t lowMassCut = 7, highMassCut = 13;
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);

    wspace.import(massVar);

    Double_t nEntries = 1e6;

	BuildInvariantMassModel(wspace, signalShapeName, bkgShapeName, fitModelName, nEntries, true);

	RooAddPdf invMassModel = *((RooAddPdf*)wspace.pdf("invMassModel"));

    // signal model parameters
    RooRealVar* mean_1S = (RooRealVar*)wspace.var("mean_1S");

    RooRealVar* yield1S = (RooRealVar*)wspace.var("yield1S");
    RooRealVar* yield2S = (RooRealVar*)wspace.var("yield2S");
    RooRealVar* yield3S = (RooRealVar*)wspace.var("yield3S");

    // background model parameters
    RooRealVar* err_mu = (RooRealVar*)wspace.var("err_mu");
    RooRealVar* err_sigma = (RooRealVar*)wspace.var("err_sigma");
    RooRealVar* exp_lambda = (RooRealVar*)wspace.var("exp_lambda");

    RooRealVar* yieldBkg = (RooRealVar*)wspace.var("yieldBkg");

    // set the parameter values using the fit results (now the example is from HX, pT 0 to 2, cosTheta -0.42 to -0.14, phi 60 to 120)
    mean_1S->setVal(9.4572);
    yield1S->setVal(5.0721e2);
    yield2S->setVal(2.1923e1);
    yield3S->setVal(3.9287e1);

    err_mu->setVal(7.2495e0);
    err_sigma->setVal(9.1318e-1);
    exp_lambda->setVal(1.5132e0);
    yieldBkg->setVal(1.1140e4);

    // yield1SInput = yield1S->getVal();
    
    // invMassModel.Print();

    // wspace.Print();

    // ((RooRealVar*)wspace.var("mean_1S"))->Print();
    // invMassModel.Print("t");

    auto* pseudoDataCanvas = new TCanvas("pseudoDataCanvas", "", 600, 600);

	RooPlot* frame = (*wspace.var("mass")).frame(Title(" "), Range("MassFitRange"));

	auto* nominalFitModel = wspace.pdf("invMassModel");

    double yieldTot = yield1S->getVal() + yield2S->getVal() + yield3S->getVal() + yieldBkg->getVal();

    RooRandom::randomGenerator()->SetSeed(time(0));

    RooDataSet* pseudoData = nominalFitModel->generate(*wspace.var("mass"), yieldTot);

    pseudoData->plotOn(frame, Name("data"), Binning(75), DrawOption("P0Z"));

	// nominalFitModel->plotOn(frame, Components(*wspace.pdf("bkgPDF")), LineColor(gColorBkg), LineStyle(kDashed), Range("MassFitRange"), NormRange("MassFitRange"));
	// nominalFitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_1S")), LineColor(gColor1S), Range("MassFitRange"));
	// nominalFitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_2S")), LineColor(gColor2S), Range("MassFitRange"));
	// nominalFitModel->plotOn(frame, Components(*wspace.pdf("signalPDF_3S")), LineColor(gColor3S), Range("MassFitRange"));
	// nominalFitModel->plotOn(frame, LineColor(gColorTotalFit), Range("MassFitRange"), NormRange("MassFitRange"));

	frame->GetYaxis()->SetMaxDigits(3);
	// gStyle->SetExponentOffset(-0.07, 0.005, "Y");
    frame->Draw();

    return pseudoData;
}

void pseudoExperiment(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.42, Float_t cosThetaMax = -0.14, Int_t phiMin = 60, Int_t phiMax = 120, Bool_t isPhiFolded = kTRUE){

	writeExtraText = true; // if extra text
	extraText = "      Internal";

	/// Set up the data
	using namespace RooFit;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	const char* refFrameName = isCSframe ? "CS" : "HX";

    const char* signalShapeName = "SymDSCB";

    const char* altBkgShapeName = "ExpTimesErr";
    // const char* altBkgShapeName = "ChebychevOrder2";

    const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

    RooWorkspace wspace(Form("workspace%s", refFrameName));

	Float_t lowMassCut = 7, highMassCut = 13;
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);

    wspace.import(massVar);

    // Double_t nEntries = pseudoDataset->sumEntries();
    Double_t nEntries = 1e5;

	BuildInvariantMassModel(wspace, signalShapeName, altBkgShapeName, fitModelName, nEntries, true);

    RooAddPdf invMassModel = *((RooAddPdf*)wspace.pdf("invMassModel"));

    TH1D* yield1Sdiff = new TH1D("yield1Sdiff", "", 20, -100, 100);

    Double_t yield1SInput = 5.0721e2;

    for (int i = 0; i < 1000; i++){
        RooDataSet* pseudoDataset = generatePseudoData(ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax/*, yield1SInput*/);
        wspace.import(*pseudoDataset);

        auto* fitResult = RawInvariantMassFit(wspace, *pseudoDataset);

        RooRealVar* yield1Svar = dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find("yield1S"));

        yield1Sdiff->Fill(yield1SInput - yield1Svar->getVal());
    }

    auto* yield1SdiffCanvas = new TCanvas("yield1SdiffCanvas", "", 600, 600);

    gStyle->SetOptStat(1111);

    yield1Sdiff->SetStats(1);

    yield1Sdiff->SetMarkerStyle(20);  // Set marker style (e.g., 20 = solid circle)
    yield1Sdiff->SetMarkerSize(1.2);  // Set marker size
    yield1Sdiff->SetMarkerColor(kBlue);  // Set marker color (e.g., blue)

    yield1Sdiff->Draw();

    yield1SdiffCanvas->Update();  

    // auto* massCanvas = new TCanvas("massCanvas", "", 600, 600);
	// TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	// pad1->SetBottomMargin(0.03);

	// pad1->Draw();
	// pad1->cd();

	// RooPlot* frame = InvariantMassRooPlot(wspace, *pseudoDataset);
	// frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad

	// frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	// if (!isPhiFolded)
	// 	frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	// else
	// 	frame->addObject(RefFrameTextPhiFolded(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	// frame->addObject(FitResultText(*wspace.var("yield1S"), ComputeSignalSignificance(wspace, 1), *wspace.var("yield2S"), ComputeSignalSignificance(wspace, 2)));

	// frame->Draw();

	// gPad->RedrawAxis();

	// // pull distribution
	// massCanvas->cd();

	// TPad* pad2 = GetPadPullDistribution(frame, fitResult->floatParsFinal().getSize());

	// //canvas->Modified();
	// //canvas->Update();
	// massCanvas->cd();
	// pad1->Draw();
	// pad2->Draw();

	// // const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

	// RooArgSet* signalYields = new RooArgSet(*wspace.var("yield1S"), *wspace.var("yield2S"), *wspace.var("yield3S"));

	// const char* totalFitModelName = GetTotalFitModelName(altBkgShapeName, signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);


}