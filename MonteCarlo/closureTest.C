#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/Parameters/PhaseSpace.h"

#include "AccEffHelpers.h"

#include "../Polarization/PolarFitHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/Parameters/CentralityValues.h"
#include "../Tools/Parameters/EfficiencyWeights.h"
#include "../Tools/Parameters/MuonScaleFactors.h"

#include "../ReferenceFrameTransform/Transformations.h"

#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.cxx"

#include "../Tools/RooFitPDFs/PhiPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PhiPolarizationPDF.cxx"

#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.h"
#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.cxx"

using namespace RooFit;

/// This code extracts the polarization parameters from the 2D angular distribution function generated using the official reco MC simulation
/// Input polarization parameters are already applied to the MC sample in File/skimReconstructedMCWeighted.C

TEfficiency* getAcceptance3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE) {

	/// get acceptance and efficiency in 1D
	const char* nominalMapName = NominalTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

    TString fileName;
    if (isPhiFolded) {
        fileName = Form("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults%s.root", gMuonAccName);
    } else {
        fileName = Form("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults%s_fullPhi.root", gMuonAccName);
    }

	// get acceptance maps
	TFile* acceptanceFile = openFile(fileName.Data());

    if (!acceptanceFile) {
        std::cerr << "Error: acceptanceFile is null." << std::endl;
        return nullptr;
    }

	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);

    if (!accMap) {
        std::cerr << "Error: accMap is null." << std::endl;
        return nullptr;
    }

	return accMap;
}

TEfficiency* getEfficiency3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE) {

	/// get acceptance and efficiency in 1D
	const char* nominalMapName = NominalTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

    TString fileName;
    if (isPhiFolded) {
        fileName = Form("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults%s.root", gMuonAccName);
    } else {
        fileName = Form("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults%s_fullPhi.root", gMuonAccName);
    }

	// get efficiency maps
	TFile* efficiencyFile = openFile(fileName.Data());

    if (!efficiencyFile) {
        std::cerr << "Error: efficiencyFile is null." << std::endl;
        return nullptr;
    }

	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);

    if (!effMap) {
        std::cerr << "Error: effMap is null." << std::endl;
        return nullptr;
    }

	return effMap;
}

TH3D* getSysEff3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE) {

	/// get acceptance and efficiency in 1D
	const char* nominalMapName = SystTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

    TString fileName;
    if (isPhiFolded) {
        fileName = Form("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults%s.root", gMuonAccName);
    } else {
        fileName = Form("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults%s_fullPhi.root", gMuonAccName);
    }
	// get relative systematic uncertainty of efficiency
	TFile* efficiencyFile = openFile(fileName.Data());

    if (!efficiencyFile) {
        std::cerr << "Error: efficiencyFile is null." << std::endl;
        return nullptr;
    }

	auto* systEff = (TH3D*)efficiencyFile->Get(nominalMapName);

    if (!systEff) {
        std::cerr << "Error: systEff is null." << std::endl;
        return nullptr;
    }

	return systEff;
}

void extractPolarParam(TH2D* correctedHist, TString refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {}, Bool_t isPhiFolded = kFALSE) {
    writeExtraText = true; // if extra text
	// extraText = "       Internal";
    extraText = "       Simulation Preliminary";
    
    TH2D* fitHist = (TH2D*)correctedHist->Clone("fitHist");

    Float_t maxYield = fitHist->GetMaximum();

	TF2* polarFunc2D = getGeneralPolarFunc(maxYield);

    /// Polarization fit
	/// with Root Fit function (E: minos on, S: save results, V: Verbose, I: integral, M: imporve algorithm, R: specificed range)
    TFitResultPtr fitResults = fitHist->Fit("generalPolarFunc", "ESVIMR"); // chi2 fit to the integrated bin 
   
    /// Fit results
    double chi2 = fitResults->Chi2();
    double nDOF = nCosThetaBins * nPhiBins - polarFunc2D->GetNpar();
    
    double pValue = TMath::Prob(chi2, nDOF); // https://root.cern.ch/doc/master/namespaceTMath.html#a3aba6abaf2bc605680c7be8fc5fd98aa

    double normVal = fitResults->Parameter(0);
    double normErr = fitResults->ParError(0);

    double lambdaThetaVal = fitResults->Parameter(1);
    double lambdaThetaErr = fitResults->ParError(1);

    double lambdaPhiVal = fitResults->Parameter(2);
    double lambdaPhiErr = fitResults->ParError(2);

    double lambdaThetaPhiVal = fitResults->Parameter(3);
    double lambdaThetaPhiErr = fitResults->ParError(3);

    double lambdaTildeVal = (lambdaThetaVal + 3. * lambdaPhiVal) / (1. - lambdaPhiVal);
    double lambdaTildeErr = TMath::Hypot(1. / (1. - lambdaPhiVal) * lambdaPhiErr, (3. - lambdaThetaVal - 6. * lambdaPhiVal) / TMath::Power((1. - lambdaPhiVal), 2) * lambdaPhiErr);

    RooRealVar lambdaTheta("lambdaTheta", "lambdaTheta", -3., 3.);
    RooRealVar lambdaPhi("lambdaPhi", "lambdaPhi", -3., 3.);
    RooRealVar lambdaThetaPhi("lambdaThetaPhi", "lambdaThetaPhi", -3., 3.);
    RooRealVar lambdaTilde("lambdaTilde", "lambdaTilde", -3., 3.);

    lambdaTheta.setVal(lambdaThetaVal);
    lambdaTheta.setError(lambdaThetaErr);

    lambdaPhi.setVal(lambdaPhiVal);
    lambdaPhi.setError(lambdaPhiErr);

    lambdaThetaPhi.setVal(lambdaThetaPhiVal);
    lambdaThetaPhi.setError(lambdaThetaPhiErr);

    lambdaTilde.setVal(lambdaTildeVal);
    lambdaTilde.setError(lambdaTildeErr);

    TCanvas* mc2DCanvas = draw2DMap(fitHist, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kTRUE, kFALSE, 1, isPhiFolded);

    fitHist->GetZaxis()->SetRangeUser(0, maxYield * 2);
    cout << "maxYield: " << maxYield << endl;

	fitHist->GetZaxis()->SetTitle("Corrected #varUpsilon(1S) Yields");
	fitHist->GetZaxis()->SetTitleSize(0.055);
	fitHist->GetZaxis()->SetTitleOffset(1.6);
	fitHist->GetZaxis()->SetLabelSize(0.044);
	
	fitHist->GetYaxis()->SetTitleSize(0.055);
	fitHist->GetYaxis()->SetTitleOffset(1.6);
	fitHist->GetYaxis()->SetLabelSize(0.044);
	fitHist->GetYaxis()->SetLabelOffset(0);

	fitHist->GetXaxis()->SetTitleSize(0.055);	
	fitHist->GetXaxis()->SetTitleOffset(1.2);
	fitHist->GetXaxis()->SetLabelSize(0.044);

	// polarFunc2D->Draw("SURFACE SAME");

    TPaveText* kinematicsText = new TPaveText(0.15, 0.86, 0.77, 0.95, "NDCNB");
	kinematicsText->SetFillColor(4000);
	kinematicsText->SetBorderSize(0);
	kinematicsText->AddText(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	kinematicsText->SetAllWith("", "align", 12);
	kinematicsText->Draw("SAME");

    TLegend legend2(.17, .76, .23, .87);
    legend2.SetTextSize(.045);
    legend2.SetFillColor(0);
    legend2.SetFillStyle(1001);

    legend2.AddEntry(fitHist, "#varUpsilon(1S) corrected yield", "lp");
    // legend2.AddEntry(polarFunc2D, Form("fit: #lambda_{#theta}  = %.2f #pm %.2f    #lambda_{#varphi} = %.2f #pm %.2f", lambdaThetaVal, lambdaThetaErr, lambdaPhiVal, lambdaPhiErr), "l");
    // legend2.AddEntry((TObject*)0, Form("     #lambda_{#theta#varphi} = %.2f #pm %.2f  #tilde{#lambda}  = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr, lambdaTildeVal, lambdaTildeErr), "");
    legend2.AddEntry(polarFunc2D, "fit: ", "l");

    legend2.DrawClone();

    TPaveText *resultTextRight = new TPaveText(0.23, 0.70, 0.60, 0.81, "NDC"); // Adjust coordinates
    resultTextRight->SetFillColor(0);   // White background
    resultTextRight->SetFillStyle(1001); // Solid fill
    resultTextRight->SetBorderSize(0);  // Optional: Thin border
    resultTextRight->SetTextSize(0.045);
    resultTextRight->SetTextAlign(12);  // Align text left
    resultTextRight->SetMargin(0.03);
    resultTextRight->AddText(Form("#lambda_{#theta}  = %.2f #pm %.2f ", lambdaThetaVal, lambdaThetaErr));
    resultTextRight->AddText(Form("#lambda_{#theta#varphi} = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr));
    resultTextRight->Draw();

    TPaveText *resultTextLeft = new TPaveText(0.52, 0.72, 0.78, 0.805, "NDC"); // Adjust coordinates
    resultTextLeft->SetFillColor(0);   // White background
    resultTextLeft->SetFillStyle(1001); // Solid fill
    resultTextLeft->SetBorderSize(0);  // Optional: Thin border
    resultTextLeft->SetTextSize(0.045);
    resultTextLeft->SetTextAlign(12);  // Align text left
    resultTextLeft->SetMargin(0.03);
    resultTextLeft->AddText(Form("#lambda_{#varphi} = %.2f #pm %.2f", lambdaPhiVal, lambdaPhiErr));
    resultTextLeft->AddText(Form("#tilde{#lambda}  = %.2f #pm %.2f", lambdaTildeVal, lambdaTildeErr));
    resultTextLeft->Draw();

    TLatex textChi2;
    textChi2.SetTextAlign(22);
    textChi2.SetTextSize(0.04);
    textChi2.DrawLatexNDC(0.72, 0.044, Form("#chi^{2} / n_{dof} = %.2f, p-value = %.2f", chi2 / nDOF, pValue));

    gPad->Update();

    gSystem->mkdir("closureTest", kTRUE);
    if (isPhiFolded) mc2DCanvas->SaveAs(Form("closureTest/2Dfit_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaThetaVal, lambdaPhiVal, lambdaThetaPhiVal), "RECREATE");
    else mc2DCanvas->SaveAs(Form("closureTest/2Dfit_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaThetaVal, lambdaPhiVal, lambdaThetaPhiVal), "RECREATE");
}

void correctMC2DHist(TH2D* polarizedHist, TH2D* correctedHist, TString refFrameName = "CS", Float_t lambdaTheta0 = 0, Float_t lambdaPhi0 = 0, Float_t lambdaThetaPhi0 = 0, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 20, const vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 18, const vector<Double_t>& phiBinEdges = {}, Bool_t drawPlot = kFALSE, Bool_t isPhiFolded = kTRUE, Bool_t applyAcc = kTRUE, Bool_t applyEff = kFALSE) {
    writeExtraText = true; // if extra text
	// extraText = "       Internal";
    extraText = "       Simulation Preliminary";

	/// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

    // set zero polarization parameters for acceptance and efficiency MC samples
    float lambdaTheta = 0;
    float lambdaPhi = 0;
    float lambdaThetaPhi = 0;

    // float lambdaTheta = lambdaTheta0;
    // float lambdaPhi = lambdaPhi0;
    // float lambdaThetaPhi = lambdaThetaPhi0;

	/// get acceptance and efficiency 3D map
	auto* accMap = getAcceptance3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded); 
	auto* effMap = getEfficiency3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded); 
	auto* systEff = getSysEff3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);	

	/// rebin acceptance and efficiency, efficiency systematic Uncertainty
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

    TH2D* hTotalCosThetaPhi = (TH2D*)accMapCosThetaPhi->GetTotalHistogram();

	/// draw acceptance and efficiency graph for check
	TCanvas* accCanvas = DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kTRUE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);
	TCanvas* effCanvas = DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kFALSE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	/// apply acc x eff correction weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			Double_t weight = 1;

			// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotalCosThetaPhi->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotalCosThetaPhi->GetYaxis()->GetBinCenter(iPhi + 1);

			Int_t iGlobalBin = hTotalCosThetaPhi->FindFixBin(binCenterCosTheta, binCenterPhi);

			// get the corresponding weights
			double acceptance = accMapCosThetaPhi->GetEfficiency(iGlobalBin);
			double efficiency = effMapCosThetaPhi->GetEfficiency(iGlobalBin);
            // double efficiency = 1;

            cout << "acceptance: " << acceptance << endl;
            cout << "efficiency: " << efficiency << endl;

			// calculate weight
            if (applyEff) {
                if (acceptance == 0 || efficiency == 0) weight = 0.;
                else weight = 1. / (acceptance * efficiency);
            }
            
            else if (applyAcc) {
                if (acceptance == 0) weight = 0.;
                else weight = 1. / (acceptance);
            }
            
            else weight = 1.;

			// weightMap->SetBinContent(iCosTheta + 1, iPhi + 1, weight);

			// // propagate both scale factor uncertainties and efficiency stat errors to the weight
			// double relSystUnc = systEff->GetBinContent(iCosTheta + 1, iPhi + 1) / efficiency;

			// double relEffUncHigh = effMap->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
			// double relEffUncLow = effMap->GetEfficiencyErrorLow(iGlobalBin) / efficiency;

			// double relAccUncHigh = accMap->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
			// double relAccUncLow = accMap->GetEfficiencyErrorLow(iGlobalBin) / acceptance;

			// totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

			// totalRelUncHigh = TMath::Hypot(relEffUncHigh, relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(relEffUncLow, relAccUncLow);

			// totalUncHigh = totalRelUncHigh * efficiency * acceptance;
			// totalUncLow = totalRelUncLow * efficiency * acceptance;

			double recoMCVal = polarizedHist->GetBinContent(iCosTheta + 1, iPhi + 1);

			double recoMCUnc = polarizedHist->GetBinError(iCosTheta + 1, iPhi + 1);

			cout << "recoMCVal: " << recoMCVal << endl;
			cout << "recoMCUnc: " << recoMCUnc << endl;

			// set the bin contents reflecting weights
			// yield with acceptance x efficiency correction
			correctedHist->SetBinContent(iCosTheta + 1, iPhi + 1, recoMCVal * weight);
					
			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);
		
			correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, recoMCUnc * weight);
			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(recoMCUnc / recoMCVal, totalRelUncHigh));

            cout << "correctedHist->GetBinContent(" << iCosTheta + 1 << ", " << iPhi + 1 << ") = " << correctedHist->GetBinContent(iCosTheta + 1, iPhi + 1) << endl;
            cout << "correctedHist->GetBinError(" << iCosTheta + 1 << ", " << iPhi + 1 << ") = " << correctedHist->GetBinError(iCosTheta + 1, iPhi + 1) << endl;
			// // fill uncertainty histograms
			// relSystEffCosThetaPhi->SetBinContent(iCosTheta +1, iPhi + 1, relSystUnc);
			// statHighEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relEffUncHigh);
			// statLowEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relEffUncLow);
			// statHighAccCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relAccUncHigh);
			// statLowAccCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relAccUncLow);
			// yield1SUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SUnc / yield1SVal);
			// totalRelUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

			// if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
		}
	}
    
    if (drawPlot) {
        TCanvas *correctedCanvas = draw2DMap(correctedHist, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, kFALSE);

        // Styles of the texts in the plot
        TLatex* legend1 = new TLatex();
        legend1->SetTextAlign(22);
        legend1->SetTextSize(0.05);

        // Put texts inside the plot
        legend1->DrawLatexNDC(.50, .88, "Corrected MC");
        legend1->DrawLatexNDC(.48, .80, Form("Input: #lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta0, lambdaPhi0, lambdaThetaPhi0));

        gPad->Update();

        gSystem->mkdir("closureTest", kTRUE);
        correctedCanvas->SaveAs(Form("closureTest/2Dcorrected_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta0, lambdaPhi0, lambdaThetaPhi0));
    }

    gSystem->mkdir("closureTest", kTRUE);
    if (isPhiFolded) {
        accCanvas->SaveAs(Form("closureTest/2Dacc_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
        effCanvas->SaveAs(Form("closureTest/2Deff_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
    }
    else {
        accCanvas->SaveAs(Form("closureTest/2Dacc_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
        effCanvas->SaveAs(Form("closureTest/2Deff_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
    }

    return;
}

/// read the input polarized MC file and fill the 2D angular distribution histogram within the specified binning
void getPolarizedMCHist(TH2D* angDistHist2D, TString refFrameName = "CS", Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {}, Bool_t isPhiFolded = kFALSE) {
    
    writeExtraText = true; // if extra text
	// extraText = "       Internal";
    extraText = "       Simulation Preliminary";

    /// read the input MC file
    const char* inputFileName = Form("Y1SReconstructedMCWeightedDataset_TriggerAcc_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", lambdaTheta, lambdaPhi, lambdaThetaPhi);
    TFile* infile = TFile::Open(Form("../Files/%s", inputFileName), "READ");

    if (!infile) {
        cout << "File " << inputFileName << " not found. Check the directory of the file." << endl;
        return;
    }

    cout << "File " << inputFileName << " opened" << endl;

    /// get the dataset
    TString datasetName = Form("MCdataset%s", refFrameName.Data());
    RooDataSet* allDataset = (RooDataSet*)infile->Get(datasetName.Data());

    if (!allDataset) {
        cout << "Dataset " << datasetName << " not found. Check the dataset name." << endl;
        return;
    }

    allDataset->Print("V");

    const char* kinematicCut = Form("(pt > %d) && (pt < %d) && (cosTheta%s > %f) && (cosTheta%s < %f) && (%s > %d) && (%s < %d)", ptMin, ptMax, refFrameName.Data(), cosThetaBinEdges[0], refFrameName.Data(), cosThetaBinEdges[nCosThetaBins], PhiVarName(refFrameName.Data()), (int)phiBinEdges[0], PhiVarName(refFrameName.Data()), (int)phiBinEdges[nPhiBins]);
    RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(kinematicCut);

    /// import the dataset to a workspace
    RooWorkspace wspace("workspace");
    wspace.import(reducedDataset);

    /// define the varibles
    RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName.Data()));
    RooRealVar phi = *wspace.var(PhiVarName(refFrameName.Data()));
    RooRealVar phiTilde = *wspace.var(PhiTildeVarName(refFrameName.Data()));

    wspace.Print("V");

    /// dummy histogram to adjust the plot range
    TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ); Number of generated #varUpsilons", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins + 5, phiBinEdges[0] - 20, phiBinEdges[nPhiBins] + 100);
  
    /// histogram to store errors on the weights
    TH2D* errorHist = new TH2D("hdummy", ";cos #theta; #varphi (#circ); Number of generated #varUpsilons", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins]); 

    /// temporary arrays to store errors on the weights
	Double_t tempErrorUp[100][100] = {0};
	Double_t tempErrorDown[100][100] = {0};
    
    /// loop over the dataset and fill 2D angular distribution histogram
    ULong64_t totEvents = reducedDataset.numEntries();

    for (Int_t iEvent = 0; iEvent < totEvents; iEvent++) {

   		// if (iEvent > 100) break;
        
        /// get the event
        const RooArgSet* iRooArgSet = reducedDataset.get(iEvent);
        
        /// get the angular variables
        Double_t cosThetaVal = ((RooRealVar*) iRooArgSet->find(CosThetaVarName(refFrameName.Data())))->getVal();
        Double_t phiVal = ((RooRealVar*) iRooArgSet->find(PhiVarName(refFrameName.Data())))->getVal();
        // Double_t phiTildeVal = ((RooRealVar*) iRooArgSet->find(PhiTildeVarName(refFrameName.Data())))->getVal();
        
        /// get the weight
        Double_t weight = reducedDataset.weight(); // weight = nColl * Gen_weight * dimuonPtWeight * dimuWeight_nominal * polarWeight
        // Double_t weightErrorUp = dynamic_cast<RooRealVar*> (iRooArgSet->find(Form("eventWeight%s", refFrameName.Data())))->getErrorHi();
        // Double_t weightErrorDown = ((RooRealVar*) iRooArgSet->find("errorWeightDown"))->getVal();

        /// get the bin numbers
        Int_t iCosTheta = angDistHist2D->GetXaxis()->FindBin(cosThetaVal);
        Int_t iPhi = angDistHist2D->GetYaxis()->FindBin(phiVal);

        /// calculate the total errors on the weights
        // tempErrorUp[iCosTheta - 1][iPhi - 1] = TMath::Hypot(tempErrorUp[iCosTheta - 1][iPhi - 1], weightErrorUp);
        // tempErrorDown[iCosTheta - 1][iPhi - 1] = TMath::Hypot(tempErrorDown[iCosTheta - 1][iPhi - 1], weightErrorDown);

        /// print out the values
 		// cout << "weight: " << weight << endl;

		// cout << "weight Error up: " << weightErrorUp << endl;
		// cout << "weight Error down: " << weightErrorDown << endl;

		// cout << "costTheta: " << cosThetaVal << endl;
		// cout << "phi: " << phiVal << endl;

		// cout << "iCosTheta: " << iCosTheta << endl;
		// cout << "iPhi: " << iPhi << endl;

		// cout << "tempErrorUp: " << tempErrorUp[iCosTheta - 1][iPhi - 1] << endl;
		// cout << "tempErrorDown: " << tempErrorDown[iCosTheta - 1][iPhi - 1] << endl;

		// cout << " " << endl;
        
        /// fill the histogram
        angDistHist2D->Fill(cosThetaVal, phiVal, weight);
        // generalPolarTildeHist->Fill(cosThetaVal, phiTildeVal, weight);
    }
    
    // /// fill the error histogram
    // for (int iCosTheta = 1; iCosTheta <= nCosThetaBins; iCosTheta++) {
	// 	for (int iPhi = 1; iPhi <= nPhiBins; iPhi++) { 
			// errorHist->SetBinContent(iCosTheta, iPhi, tempErrorUp[iCosTheta - 1][iPhi - 1]);

			// angDistHist2D->SetBinError(iCosTheta, iPhi, tempErrorUp[iCosTheta - 1][iPhi - 1]);
			// cout << "BinError: " << angDistHist2D->GetBinError(iCosTheta, iPhi) << endl;
			// cout << "BinValue: " << angDistHist2D->GetBinContent(iCosTheta, iPhi) << endl;
	// 	}
	// }

	/// draw polarized MC histogram
	draw2DMap(angDistHist2D, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
    // display2DMapContents(angDistHist2D, nCosThetaBins, nPhiBins, kFALSE);

	//draw error histogram
	// draw2DMap(errorHist, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
	// display2DMapContents(errorHist, nCosThetaBins, nPhiBins, kFALSE);

	// Styles of the texts in the plot
	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(22);
	legend1->SetTextSize(0.05);

	// Put texts inside the plot
	legend1->DrawLatexNDC(.50, .88, "Reconstructed MC");
	legend1->DrawLatexNDC(.50, .80, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	gPad->Update();

	// // save the plot
	// gSystem->mkdir("DistributionFitsMC", kTRUE);
	// mc2DCanvas->SaveAs(Form("DistributionFitsMC/RecoMC_Theta%.2f_Phi%.2f_ThetaPhi%.2f.png", lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");

    return;
}

void getGENMCHist(TH2D* angDistHist2D, TString refFrameName = "CS", Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {}, Bool_t applyAcc = kTRUE, Bool_t isPhiFolded = kFALSE) {
    
    writeExtraText = true; // if extra text
    // extraText = "       Internal";
    extraText = "       Simulation Preliminary";

    /// read the input MC file
    const char* inputFileName;
    if (applyAcc) inputFileName = Form("Y1SGenNoFilterMCDataset%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", gMuonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi); // Gen + Acc
    else inputFileName = Form("Y1SGenNoFilterMCDataset_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", lambdaTheta, lambdaPhi, lambdaThetaPhi); // Pure Gen
    
    TFile* infile = TFile::Open(Form("../Files/%s", inputFileName), "READ");

    if (!infile) {
        cout << "File " << inputFileName << " not found. Check the directory of the file." << endl;
        return;
    }

    cout << "File " << inputFileName << " opened" << endl;

    /// get the dataset
    TString datasetName = Form("MCdataset%s", refFrameName.Data());
    RooDataSet* allDataset = (RooDataSet*)infile->Get(datasetName.Data());

    if (!allDataset) {
        cout << "Dataset " << datasetName << " not found. Check the dataset name." << endl;
        return;
    }

    allDataset->Print("V");

    // const char* kinematicCut = Form("centrality >= %d && centrality < %d && mass >= %f && mass < %f && rapidity > %f && rapidity < %f && pt > %d && pt < %d && cosTheta%s > %f && cosTheta%s < %f && fabs(%s) > %d && fabs(%s) < %d", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);
    // const char* kinematicCut = Form("(pt > %d) && (pt < %d) && (cosTheta%s > %f) && (cosTheta%s < %f) && (fabs(%s) > %d) && (fabs(%s) < %d)", ptMin, ptMax, refFrameName.Data(), cosThetaBinEdges[0], refFrameName.Data(), cosThetaBinEdges[5], PhiVarName(refFrameName.Data()), (int)phiBinEdges[0], PhiVarName(refFrameName.Data()), (int)phiBinEdges[nPhiBins]);
    const char* kinematicCut = Form("(pt > %d) && (pt < %d) && (cosTheta%s > %f) && (cosTheta%s < %f) && ((%s) > %d) && ((%s) < %d)", ptMin, ptMax, refFrameName.Data(), cosThetaBinEdges[0], refFrameName.Data(), cosThetaBinEdges[nCosThetaBins], PhiVarName(refFrameName.Data()), (int)phiBinEdges[0], PhiVarName(refFrameName.Data()), (int)phiBinEdges[nPhiBins]);
    RooDataSet reducedDataset = *(RooDataSet*)allDataset->reduce(kinematicCut);

    /// import the dataset to a workspace
    RooWorkspace wspace("workspace");
    wspace.import(reducedDataset);

    /// define the varibles
    RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName.Data()));
    RooRealVar phi = *wspace.var(PhiVarName(refFrameName.Data()));
    // RooRealVar phiTilde = *wspace.var(PhiTildeVarName(refFrameName.Data()));

    wspace.Print("V");

    /// dummy histogram to adjust the plot range
    TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ); Number of generated #varUpsilons", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins + 5, phiBinEdges[0] - 20, phiBinEdges[nPhiBins] + 100);
  
    /// loop over the dataset and fill 2D angular distribution histogram
    ULong64_t totEvents = reducedDataset.numEntries();

    for (Int_t iEvent = 0; iEvent < totEvents; iEvent++) {
   		// if (iEvent > 100) break;
        
        /// get the event
        const RooArgSet* iRooArgSet = reducedDataset.get(iEvent);
        
        /// get the angular variables
        Double_t cosThetaVal = ((RooRealVar*) iRooArgSet->find(CosThetaVarName(refFrameName.Data())))->getVal();
        Double_t phiVal = ((RooRealVar*) iRooArgSet->find(PhiVarName(refFrameName.Data())))->getVal();
        // Double_t phiTildeVal = ((RooRealVar*) iRooArgSet->find(PhiTildeVarName(refFrameName.Data())))->getVal();
        
        // /// get the weight
        Double_t weight = reducedDataset.weight(); // weight = nColl * Gen_weight * dimuonPtWeight * dimuWeight_nominal * polarWeight
        // Double_t weightErrorUp = dynamic_cast<RooRealVar*> (iRooArgSet->find(Form("eventWeight%s", refFrameName.Data())))->getErrorHi();
        // Double_t weightErrorDown = ((RooRealVar*) iRooArgSet->find("errorWeightDown"))->getVal();

        /// get the bin numbers
        Int_t iCosTheta = angDistHist2D->GetXaxis()->FindBin(cosThetaVal);
        Int_t iPhi = angDistHist2D->GetYaxis()->FindBin(phiVal);

        /// calculate the total errors on the weights
        // tempErrorUp[iCosTheta - 1][iPhi - 1] = TMath::Hypot(tempErrorUp[iCosTheta - 1][iPhi - 1], weightErrorUp);
        // tempErrorDown[iCosTheta - 1][iPhi - 1] = TMath::Hypot(tempErrorDown[iCosTheta - 1][iPhi - 1], weightErrorDown);

        /// print out the values
 		// cout << "weight: " << weight << endl;

		// cout << "weight Error up: " << weightErrorUp << endl;
		// cout << "weight Error down: " << weightErrorDown << endl;

		// cout << "costTheta: " << cosThetaVal << endl;
		// cout << "phi: " << phiVal << endl;

		// cout << "iCosTheta: " << iCosTheta << endl;
		// cout << "iPhi: " << iPhi << endl;

		// cout << "tempErrorUp: " << tempErrorUp[iCosTheta - 1][iPhi - 1] << endl;
		// cout << "tempErrorDown: " << tempErrorDown[iCosTheta - 1][iPhi - 1] << endl;

		// cout << " " << endl;
        
        /// fill the histogram
        angDistHist2D->Fill(cosThetaVal, phiVal, weight);
        
        // generalPolarTildeHist->Fill(cosThetaVal, phiTildeVal, weight);
    }
    
    // /// fill the error histogram
    // for (int iCosTheta = 1; iCosTheta <= nCosThetaBins; iCosTheta++) {
	// 	for (int iPhi = 1; iPhi <= nPhiBins; iPhi++) { 
			// errorHist->SetBinContent(iCosTheta, iPhi, tempErrorUp[iCosTheta - 1][iPhi - 1]);

			// angDistHist2D->SetBinError(iCosTheta, iPhi, tempErrorUp[iCosTheta - 1][iPhi - 1]);
			// cout << "BinError: " << angDistHist2D->GetBinError(iCosTheta, iPhi) << endl;
			// cout << "BinValue: " << angDistHist2D->GetBinContent(iCosTheta, iPhi) << endl;
	// 	}
	// }

	/// draw polarized MC histogram
	TCanvas *mc2DCanvas = draw2DMap(angDistHist2D, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
    // display2DMapContents(angDistHist2D, nCosThetaBins, nPhiBins, kFALSE);

	//draw error histogram
	// draw2DMap(errorHist, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
	// display2DMapContents(errorHist, nCosThetaBins, nPhiBins, kFALSE);

	// Styles of the texts in the plot
	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(22);
	legend1->SetTextSize(0.05);

	// Put texts inside the plot
	if (applyAcc) legend1->DrawLatexNDC(.50, .88, "GEN MC + Acc");
    else legend1->DrawLatexNDC(.50, .88, "No Filter GEN MC");
	legend1->DrawLatexNDC(.48, .80, Form("Input: #lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	gPad->Update();

	// save the plot
	gSystem->mkdir("ClosureTest", kTRUE);
	if (isPhiFolded) mc2DCanvas->SaveAs(Form("ClosureTest/GenMC_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_Theta%.2f_Phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
    else mc2DCanvas->SaveAs(Form("ClosureTest/GenMC_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_Theta%.2f_Phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");

    return;
}

void closureTest(TString refFrameName = "CS", 
                Float_t lambdaTheta0 = 1, Float_t lambdaPhi0 = 0, Float_t lambdaThetaPhi0 = 0, 
                Int_t ptMin = 2, Int_t ptMax = 6, 
                const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, 
                const Int_t nPhiBins = 6, Int_t phiMin = -180, Int_t phiMax = 180, 
                Bool_t isPhiFolded = kFALSE, Bool_t applyAcc = kTRUE, Bool_t applyEff = kFALSE) {

    /// set bin edges and width 
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);
	
	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	/// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;
            
	/// create 2D angular distribution histogram
	TH2D* angDistHist2D = new TH2D("angDistHist2D", "; cos #theta; #varphi (#circ); Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 
 
    getGENMCHist(angDistHist2D, refFrameName, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, applyAcc, isPhiFolded);
	// getPolarizedMCHist(angDistHist2D, refFrameName, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
    
    /// apply acceptance and efficiency correction
    TH2D *correctedHist = new TH2D("correctedHist", "; cos #theta; #varphi (#circ); Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
    correctMC2DHist(angDistHist2D, correctedHist, refFrameName, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kTRUE, isPhiFolded, applyAcc, applyEff);

    TH2D *fittedHist = new TH2D("fittedHist", "; cos #theta; #varphi (#circ); Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
    correctMC2DHist(angDistHist2D, fittedHist, refFrameName, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, isPhiFolded, applyAcc, applyEff);   

    extractPolarParam(fittedHist, refFrameName, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, isPhiFolded);

    return;
}
