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

using namespace std;

RooDataSet InvMassDataset(RooWorkspace& wspace, Int_t ptMin = 0, Int_t ptMax = 30, Float_t cosThetaMin = -0.1, Float_t cosThetaMax = 0.1, const char* refFrameName = "CS", Int_t phiMin = -180, Int_t phiMax = 180, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax) {
	RooDataSet* allDataset = (RooDataSet*)wspace.data(RawDatasetName(""));

	if (allDataset == nullptr) {
		std::cerr << "Null RooDataSet provided to the reducer method!!" << std::endl;
		return RooDataSet();
	}

	// const char* kinematicCut = Form("centrality >= %d && centrality < %d && rapidity > %f && rapidity < %f && pt > %d && pt < %d && cosTheta%s > %f && cosTheta%s < %f && %s > %d && %s < %d", 2 * gCentralityBinMin, 2 * gCentralityBinMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);
	const char* kinematicCut = Form("centrality >= %d && centrality < %d && mass >= %f && mass < %f && rapidity > %f && rapidity < %f && pt > %d && pt < %d && cosTheta%s > %f && cosTheta%s < %f && fabs(%s) > %d && fabs(%s) < %d", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);

	RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(RooArgSet(*(wspace.var("mass"))), kinematicCut);

	wspace.import(reducedDataset, Rename(Form("dataset_cosTheta_%.2fto%.2f_absphi_%dto%d", cosThetaMin, cosThetaMax, phiMin, phiMax)));

	return reducedDataset;
}

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

RooDataSet* generatePseudoData(Int_t ptMin = 2, Int_t ptMax = 6, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.42, Float_t cosThetaMax = -0.14, Int_t phiMin = 60, Int_t phiMax = 120, Double_t* yield1SInput = nullptr, Float_t massMin = MassBinMin, Float_t massMax = MassBinMax){

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

    // Float_t lowMassCut = 6.5, highMassCut = 14.5;
    // Float_t lowMassCut = 7, highMassCut = 14;
    Float_t lowMassCut = massMin, highMassCut = massMax;
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);

    massVar.setRange("MassFitRange", massMin, massMax);

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
	
    /// get the fit results from the nominal fit model
    const char* totalFitModelName = GetTotalFitModelName(bkgShapeName, signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

    RooFitResult* fitResults = GetFitResults(totalFitModelName, gMuonAccName);

    RooArgSet* fitParams = new RooArgSet(fitResults->floatParsFinal());

    mean_1S->setVal(((RooRealVar*)fitParams->find("mean_1S"))->getVal());
    yield1S->setVal(((RooRealVar*)fitParams->find("yield1S"))->getVal());
    yield2S->setVal(((RooRealVar*)fitParams->find("yield2S"))->getVal());
    yield3S->setVal(((RooRealVar*)fitParams->find("yield3S"))->getVal());

    err_mu->setVal(((RooRealVar*)fitParams->find("err_mu"))->getVal());
    err_sigma->setVal(((RooRealVar*)fitParams->find("err_sigma"))->getVal());
    exp_lambda->setVal(((RooRealVar*)fitParams->find("exp_lambda"))->getVal());
    yieldBkg->setVal(((RooRealVar*)fitParams->find("yieldBkg"))->getVal());

    // cout << "yield1S: " << yield1S->getVal() << endl;

    /// record the input yield1S and pass it to the main function
    *yield1SInput = (Double_t)(yield1S->getVal());
    
    /// print out invMassModel and wspace for confirmation
    // invMassModel.Print();

    // wspace.Print();

    // ((RooRealVar*)wspace.var("mean_1S"))->Print();
    // invMassModel.Print("t");

    /// get the nominal fit model from the workspace
	auto* nominalFitModel = wspace.pdf("invMassModel");
    nominalFitModel->setNormRange("MassFitRange");
 
    /// calculate the total yield
    double yieldTot = yield1S->getVal() + yield2S->getVal() + yield3S->getVal() + yieldBkg->getVal();

    /// radoimize the seed
    RooRandom::randomGenerator()->SetSeed(time(0));

    /// generate the pseudo-data
    // RooDataSet* pseudoData = nominalFitModel->generate(*wspace.var("mass"), NumEvents(yieldTot));
    // RooDataSet* pseudoData = nominalFitModel->generate(*wspace.var("mass"), yieldTot);
    RooDataSet* pseudoData = nominalFitModel->generate(*wspace.var("mass"));

    double generatedYield = pseudoData->sumEntries();
    std::cout << "Input Yield: " << yieldTot << ", Generated Yield: " << generatedYield << std::endl;

    /// draw the generated pseudo-data
    auto* pseudoDataCanvas = new TCanvas("pseudoDataCanvas", "", 600, 600);

    RooPlot* frame = InvariantMassRooPlot(wspace, *pseudoData);

    std::cout << "NominalFitModel Total Yield: " 
              << nominalFitModel->expectedEvents(*wspace.var("mass")) << std::endl;

    // Check ranges
    std::cout << "NominalFitModel Total Yield: " 
              << nominalFitModel->expectedEvents(*wspace.var("mass")) << std::endl;

    std::cout << "Generated Events: " << pseudoData->numEntries() << std::endl;

    // Integrals over default range and fit range
    double integralDefault = nominalFitModel->createIntegral(RooArgSet(*wspace.var("mass")))->getVal();
    double integralFitRange = nominalFitModel->createIntegral(RooArgSet(*wspace.var("mass")), RooFit::Range("MassFitRange"))->getVal();

    std::cout << "Integral over default range: " << integralDefault << std::endl;
    std::cout << "Integral over fit range: " << integralFitRange << std::endl;

    // Generated events in range
    std::cout << "Generated Events in Range: " 
              << pseudoData->sumEntries("mass > 6.5 && mass < 14.5") << std::endl;

    // Draw the frame
    frame->Draw();
    // nominalFitModel->Print("t");

    return pseudoData;
}

void pseudoExperiment(Int_t ptMin = 0, Int_t ptMax = 2, Bool_t isCSframe = kFALSE, Float_t cosThetaMin = -0.7, Float_t cosThetaMax = -0.42, Int_t phiMin = 0, Int_t phiMax = 60, Bool_t isPhiFolded = kTRUE, Float_t massMin = MassBinMin, Float_t massMax = 11.5){

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
    // const char* altSignalShapeName = "Johnson";

    const char* altBkgShapeName = "ExpTimesErr";
    // const char* altBkgShapeName = "ChebychevOrder3";

    /*******************************************************************************************************/

	const char* refFrameName = isCSframe ? "CS" : "HX";

    const char* altFitModelName = GetFitModelName(altSignalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

    RooWorkspace wspace(Form("workspace%s", refFrameName));

	// Float_t lowMassCut = 6.5, highMassCut = 14.5;
    Float_t lowMassCut = MassBinMin, highMassCut = MassBinMax;
	RooRealVar massVar("mass", gMassVarTitle, lowMassCut, highMassCut, gMassUnit);

    massVar.setRange("MassFitRange", massMin, massMax); // in this way, the bin width is fixed to 0.08 GeV/c2

    wspace.import(massVar);

    // cout << "wspace print " << endl;
    // wspace.Print(); 

    /// test with the real dataset********************************************************
 	// const char* filename = Form("../Files/UpsilonSkimmedDataset%s.root", gMuonAccName);
    
    // RooWorkspace wspace = SetUpWorkspace(filename, "");

    // RooRealVar massVar = *wspace.var("mass");
    // // wspace.Print("t");
	// massVar.setRange("MassFitRange", massMin, massMax);
    
    // auto allDataset = InvMassCosThetaPhiDataset(wspace, ptMin, ptMax, "");
  	
    // RooDataSet reducedDataset = InvMassDataset(wspace, ptMin, ptMax, cosThetaMin, cosThetaMax, refFrameName, phiMin, phiMax);
	
    // Long64_t nEntries = reducedDataset.sumEntries();
    /// **********************************************************************************

    Double_t nEntries = 1e5;

    Double_t yield1SInput = 0.;

    Long64_t nPseudoExperiments = 1e0;

    // TH1D* yield1Sdiff = new TH1D(Form("yield1Sdiff_%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments), "", 60, -300, 300);
    TH1D* yield1Sdiff = new TH1D(Form("%s+%s", altSignalShapeName, altBkgShapeName), "", 60, -300, 300);
    // TH1D* yield1Sdiff = new TH1D(Form("yield1Sdiff_%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments), "", 400, -2000, 2000);

    RooDataSet* pseudoDataset = nullptr;

    RooFitResult* fitResult = nullptr;

    RooPlot* frame = nullptr;

    for (Long64_t iEvent = 0; iEvent < nPseudoExperiments; iEvent++){

        if (iEvent % 1000 == 0) {
			cout << Form("\rProcessing event %lld / %lld (%.0f%%)", iEvent, nPseudoExperiments, 100. * iEvent / nPseudoExperiments) << flush;
		}

        /// generate the pseudo-data
        pseudoDataset = generatePseudoData(ptMin, ptMax, isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax, &yield1SInput, massMin, massMax);
        wspace.import(*pseudoDataset);

        altFitModelName = GetFitModelName(altSignalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax); // lost the information of the name so define again

        /// Build the invariant mass model with the alternative signal and background shapes
        BuildInvariantMassModel(wspace, altSignalShapeName, altBkgShapeName, altFitModelName, nEntries, true);
        // BuildInvariantMassModel(wspace, altSignalShapeName, altBkgShapeName, altFitModelName, nEntries, false);

        RooAddPdf invMassModel = *((RooAddPdf*)wspace.pdf("invMassModel"));
        invMassModel.setNormRange("MassFitRange");

        /// fit the pseudo-data with the alternative signal and background shapes
        fitResult = RawInvariantMassFit(wspace, *pseudoDataset);

        /// test with the real dataset
        // fitResult = RawInvariantMassFit(wspace, reducedDataset);

        /// get the fitted yield1S
        RooRealVar* yield1Svar = dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find("yield1S"));

        frame = InvariantMassRooPlot(wspace, *pseudoDataset);
        // frame = InvariantMassRooPlot(wspace, reducedDataset);

        int status = fitResult->status();
        // cout << "fit result status: " << status << endl;

        int covQual = fitResult->covQual();
        // cout << "covQual: " << covQual << endl;

        double chi2 = frame->chiSquare(fitResult->floatParsFinal().getSize());
        // cout << "chi2/Ndf: " << chi2 << endl;

        cout << "correlation matrix: " << endl;
        fitResult->correlationMatrix().Print("v");

        if (chi2 > 5 || status != 0 || covQual < 2) {
            if (chi2 > 5) cout << "chi2/Ndf > 5, skipping this pseudo-experiment" << endl;
            if (status != 0) cout << "status != 0, skipping this pseudo-experiment" << endl;
            if (covQual < 2) cout << "covQual < 2, skipping this pseudo-experiment" << endl;
            continue;
        }

        /// fill the histogram with the difference between the input and fitted yield1S
        yield1Sdiff->Fill(yield1SInput - yield1Svar->getVal());
    }

    /// draw the yield difference histogram
    auto* yield1SdiffCanvas = new TCanvas("yield1SdiffCanvas", "", 600, 400);
    yield1SdiffCanvas->SetRightMargin(0.035);

    gStyle->SetOptStat(1111);

    yield1Sdiff->SetStats(1);

    yield1Sdiff->SetMarkerStyle(20);  // Set marker style (e.g., 20 = solid circle)
    yield1Sdiff->SetMarkerSize(1.2);  // Set marker size
    yield1Sdiff->SetMarkerColor(kBlue);  // Set marker color (e.g., blue)
   
    yield1Sdiff->SetXTitle("(Input - Fit) yield 1S");
    yield1Sdiff->SetYTitle("Entries / 10");

    yield1Sdiff->Draw();

    TPaveText *KinematicsLegend = KinematicsText_YieldDiff(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax);
    KinematicsLegend->Draw();

    TPaveText *RefFrameLegendPhiFolded = RefFrameTextPhiFolded_YieldDiff(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax);
    RefFrameLegendPhiFolded->Draw();

    CustomizeStatBox(yield1Sdiff, 0.7, 0.93, 0.73, 0.90);

    yield1SdiffCanvas->Update();  

    gSystem->mkdir("YieldDifferencePlots", kTRUE);
    TFile* yield1SdiffFile = new TFile("YieldDifferencePlots/yield1Sdiff.root", "UPDATE");
   
    yield1Sdiff->Write();

    yield1SdiffFile->Close();

    yield1SdiffCanvas->SaveAs(Form("YieldDifferencePlots/yield1Sdiff_%s_%s_pt%dto%d_%s_cosTheta%.2fto%.2f_phi%dto%d_n%lld.png", altSignalShapeName, altBkgShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax, nPseudoExperiments));

    /// draw fitted invariant mass distribution
    auto* massCanvas = new TCanvas("massCanvas", "", 600, 600);
	
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.03);
    pad1->SetRightMargin(0.035);

	pad1->Draw();
	pad1->cd();

    /// draw the generated pseudo-data

	frame->GetXaxis()->SetLabelOffset(1); // to make it disappear under the pull distribution pad

	frame->addObject(KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin, ptMax));

	if (!isPhiFolded)
		frame->addObject(RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));
	else
		frame->addObject(RefFrameTextPhiFolded(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax));

	if (strcmp(altSignalShapeName, "SymDSCB") == 0)
        frame->addObject(FitResultText(*wspace.var("yield1S"), ComputeSignalSignificance(wspace, 1), *wspace.var("yield2S"), ComputeSignalSignificance(wspace, 2)));
    else
        frame->addObject(FitResultText(*wspace.var("yield1S"), 0, *wspace.var("yield2S"), 0));

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