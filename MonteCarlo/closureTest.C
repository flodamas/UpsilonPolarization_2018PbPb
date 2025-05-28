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

#include "../MonteCarlo/acceptanceMap_noGenFilter.C"
#include "../MonteCarlo/weightedEfficiencyMaps.C"

using namespace RooFit;

/// This code extracts the polarization parameters from the 2D angular distribution function generated using the official reco MC simulation
/// The input polarization parameters are already applied to the MC sample in File/skimReconstructedMCWeighted.C or in Files/skimGenUpsilonMC.C
/// The corresponding acceptance and efficiency should be skimmed at least once in advance as well by using MonteCarlo/acceptanceMap_noGenFilter.C and MonteCarlo/weightedEfficiencyMaps.C (this can be improved in the future)

/// 1. get the polarizaed gen MC or reco MC histograms
/// 2. correct the angular distribution histogram with the acceptance or/and efficiency maps
/// 3. extract the polarization parameters from the corrected angular distribution histogram
/// 4. iterate the steps 2 and 3 until the polarization parameters converge (lambda diff < 0.01) (or the number of iterations is reached)

TEfficiency* getAcceptance3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE) {
	/// get acceptance and efficiency in 1D
	TString nominalMapName = NominalTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TString fileName = Form("%s/AcceptanceResults%s.root", AcceptanceResultsPath(gMuonAccName), isPhiFolded ? "" : "_fullPhi");

	// get acceptance maps
	// TFile* acceptanceFile = openFile(fileName.Data());
	TFile* acceptanceFile = TFile::Open(fileName, "READ");

	if (!acceptanceFile) {
		std::cerr << "Error: acceptanceFile is null." << std::endl;
		acceptanceFile->Close();
		delete acceptanceFile;

		// acceptanceMap_noGenFilter(0, 30, gUpsilonState, lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded, "MuonUpsilonTriggerAcc");

		return nullptr;
	}

	acceptanceFile = openFile(fileName.Data());

	/// get acceptance maps
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName.Data());

	/// check if the acceptance map is null and create it if it doesn't exist
	if (!accMap) {
		std::cerr << "Error: accMap is null." << std::endl;
		acceptanceFile->Close();
		delete acceptanceFile;
		delete accMap;

		// create the acceptance map
		cout << Form("Acceptance map not found. Creating a new one with (lambdaTheta, phi, thetaPhi) = (%.2f, %.2f, %.2f)....", lambdaTheta, lambdaPhi, lambdaThetaPhi) << endl;

		acceptanceMap_noGenFilter(0, 30, isPhiFolded, "UpsilonTriggerThresholds", lambdaTheta, lambdaPhi, lambdaThetaPhi, gUpsilonState);

        // force ROOT to forget cached file
        gROOT->GetListOfFiles()->Remove(gROOT->GetFile(fileName.Data()));
        delete TFile::Open(fileName.Data());  // flush disk buffers once

		TFile* newAcceptanceFile = openFile(fileName.Data());
        // newAcceptanceFile->cd();
        // newAcceptanceFile->ls();  // now you should see your object

        cout << "nominal map name: " << nominalMapName << endl;
		auto* newAccMap = (TEfficiency*)newAcceptanceFile->Get(nominalMapName.Data());
		cout << "Acceptance map is loaded." << endl;

		/// check if the acceptance map is still null after creating it
		if (!newAccMap) {
			std::cerr << "Error: accMap is still null after creating it." << std::endl;
			// accMap->ls();
			return nullptr;
		}

		return newAccMap;
	}

	return accMap;
}

TEfficiency* getEfficiency3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE) {
	/// get acceptance and efficiency in 1D
	TString nominalMapName = NominalTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TString fileName = Form("%s/EfficiencyResults%s.root", EfficiencyResultsPath(gMuonAccName), isPhiFolded ? "" : "_fullPhi");

	// get efficiency maps
	TFile* efficiencyFile = openFile(fileName.Data());

	if (!efficiencyFile) {
		std::cerr << "Error: efficiencyFile is null." << std::endl;

		return nullptr;
	}

	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName.Data());

	if (!effMap) {
		std::cerr << "Error: effMap is null." << std::endl;
		weightedEfficiencyMaps(0, 30, "UpsilonTriggerThresholds", lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded, gUpsilonState);

        // force ROOT to forget cached file
        gROOT->GetListOfFiles()->Remove(gROOT->GetFile(fileName.Data()));
        delete TFile::Open(fileName.Data());  // flush disk buffers once

		efficiencyFile = openFile(fileName.Data());
        efficiencyFile->cd();
        efficiencyFile->ls();  // now you should see your object

        cout << "nominal map name: " << nominalMapName.Data() << endl;
		effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName.Data());
		cout << "Efficiency map is loaded." << endl;

		if (!effMap) {
			std::cerr << "Error: effMap is still null after creating it." << std::endl;
			// effMap->ls();
			return nullptr;
		}
	}

	return effMap;
}

TH3D* getSysEff3DMap(const char* refFrameName, Double_t lambdaTheta = 0, Double_t lambdaPhi = 0, Double_t lambdaThetaPhi = 0, Bool_t isPhiFolded = kFALSE) {
	/// get acceptance and efficiency in 1D
	const char* nominalMapName = SystTEfficiency3DName(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TString fileName = Form("%s/EfficiencyResults%s.root", EfficiencyResultsPath(gMuonAccName), isPhiFolded ? "" : "_fullPhi");

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

TLine* drawLine(double x1, double y1, double x2, double y2) {
	TLine* line = new TLine(x1, y1, x2, y2);
	line->SetLineStyle(2);
	line->SetLineColor(kBlack);

	line->Draw("SAME");
	return line;
}
/// fit the corrected histo and extract polarization parameters
RooArgSet extractPolarParam(TH2D* correctedHist, TString refFrameName = "CS",
                            Int_t ptMin = 0, Int_t ptMax = 30,
                            Int_t nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {},
                            Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {},
                            Bool_t isPhiFolded = kFALSE, Bool_t applyEff = kTRUE) {
	
	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";

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

	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", -3., 3.);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", -3., 3.);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", -3., 3.);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", -3., 3.);

	lambdaTheta->setVal(lambdaThetaVal);
	lambdaTheta->setError(lambdaThetaErr);

	lambdaPhi->setVal(lambdaPhiVal);
	lambdaPhi->setError(lambdaPhiErr);

	lambdaThetaPhi->setVal(lambdaThetaPhiVal);
	lambdaThetaPhi->setError(lambdaThetaPhiErr);

	lambdaTilde->setVal(lambdaTildeVal);
	lambdaTilde->setError(lambdaTildeErr);

	/// store fit results
	RooArgSet fittedParams(*lambdaTheta, *lambdaPhi, *lambdaThetaPhi, *lambdaTilde);

	TCanvas* mc2DCanvas = draw2DMap(fitHist, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kTRUE, kFALSE, 1, isPhiFolded);

	fitHist->GetZaxis()->SetRangeUser(0, maxYield * 2);
	cout << "maxYield: " << maxYield << endl;

	fitHist->GetZaxis()->SetTitle("#varUpsilon(1S) weighted events");
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
	TPaveText* kinematicsText = new TPaveText(0.28, 0.89, 0.85, 0.93, "NDCNB"); //0.375,0.8521739,0.9459877,0.9043478
	kinematicsText->SetFillColor(4000);
	kinematicsText->SetBorderSize(0);
	kinematicsText->AddText(Form("%s, %s", DimuonRapidityRangeText(gRapidityMin, gRapidityMax), DimuonPtRangeText(ptMin, ptMax)));
	kinematicsText->SetAllWith("", "align", 12);
	kinematicsText->Draw("SAME");

	/// add results to the plot
	TLegend legend2(.16, .78, .46, .89); // 0.159292,0.7791304,0.460177,0.8886957
	legend2.SetTextSize(.045);
	legend2.SetFillColor(0);
	legend2.SetFillStyle(1001);

	legend2.AddEntry(fitHist, "#varUpsilon(1S) corrected yield", "lp");
	// legend2.AddEntry(polarFunc2D, Form("fit: #lambda_{#theta}  = %.2f #pm %.2f    #lambda_{#varphi} = %.2f #pm %.2f", lambdaThetaVal, lambdaThetaErr, lambdaPhiVal, lambdaPhiErr), "l");
	// legend2.AddEntry((TObject*)0, Form("     #lambda_{#theta#varphi} = %.2f #pm %.2f  #tilde{#lambda}  = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr, lambdaTildeVal, lambdaTildeErr), "");
	legend2.AddEntry(polarFunc2D, "fit: ", "l");

	legend2.DrawClone();

	TPaveText* resultTextRight = new TPaveText(0.28, 0.72, 0.65, 0.83, "NDC"); // Adjust coordinates
	resultTextRight->SetFillColor(0);                                          // White background
	resultTextRight->SetFillStyle(1001);                                       // Solid fill
	resultTextRight->SetBorderSize(0);                                         // Optional: Thin border
	resultTextRight->SetTextSize(0.045);
	resultTextRight->SetTextAlign(12); // Align text left
	resultTextRight->SetMargin(0.03);
	resultTextRight->AddText(Form("#lambda_{#theta}  = %.2f #pm %.2f ", lambdaThetaVal, lambdaThetaErr));
	resultTextRight->AddText(Form("#lambda_{#theta#varphi} = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr));
	resultTextRight->Draw();

	TPaveText* resultTextLeft = new TPaveText(0.57, 0.74, 0.83, 0.825, "NDC"); // Adjust coordinates
	resultTextLeft->SetFillColor(0);                                           // White background
	resultTextLeft->SetFillStyle(1001);                                        // Solid fill
	resultTextLeft->SetBorderSize(0);                                          // Optional: Thin border
	resultTextLeft->SetTextSize(0.045);
	resultTextLeft->SetTextAlign(12); // Align text left
	resultTextLeft->SetMargin(0.03);
	resultTextLeft->AddText(Form("#lambda_{#varphi} = %.2f #pm %.2f", lambdaPhiVal, lambdaPhiErr));
	resultTextLeft->AddText(Form("#tilde{#lambda}  = %.2f #pm %.2f", lambdaTildeVal, lambdaTildeErr));
	resultTextLeft->Draw();

	TLatex textChi2;
	textChi2.SetTextAlign(22);
	textChi2.SetTextSize(0.04);
	textChi2.DrawLatexNDC(0.72, 0.044, Form("#chi^{2} / n_{dof} = %.2f, p-value = %.2f", chi2 / nDOF, pValue));

	// gPad->RedrawAxis();

	gPad->Update();

	gSystem->mkdir(Form("closureTest/%s", applyEff ? "HydjetAccEff": "PythiaAcc"), kTRUE);
	if (isPhiFolded)
		mc2DCanvas->SaveAs(Form("closureTest/%s/2Dfit_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaThetaVal, lambdaPhiVal, lambdaThetaPhiVal));
	else
		mc2DCanvas->SaveAs(Form("closureTest/%s/2Dfit_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaThetaVal, lambdaPhiVal, lambdaThetaPhiVal));

	return fittedParams;
}

/// correct the polarized MC with weights (1 / (acceptance and/or efficiency))
void correctMC2DHist(TH2D* polarizedHist, TH2D* correctedHist, TString refFrameName = "CS",
                     RooArgSet polarParams = RooArgSet(),
                     Int_t ptMin = 0, Int_t ptMax = 30,
                     Int_t nCosThetaBins = 20, const vector<Double_t>& cosThetaBinEdges = {},
                     Int_t nPhiBins = 18, const vector<Double_t>& phiBinEdges = {},
                     Bool_t drawPlot = kFALSE, Bool_t isPhiFolded = kTRUE, Bool_t applyAcc = kTRUE, Bool_t applyEff = kFALSE,
					 float lambdaTheta0 = 0, float lambdaPhi0 = 0, float lambdaThetaPhi0 = 0) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";

	/// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	RooRealVar* lambdaThetaVar = (RooRealVar*)polarParams.find("lambdaTheta");
	RooRealVar* lambdaPhiVar = (RooRealVar*)polarParams.find("lambdaPhi");
	RooRealVar* lambdaThetaPhiVar = (RooRealVar*)polarParams.find("lambdaThetaPhi");

	/// set the polarization parameters for acceptance and efficiency MC samples to the input polarization parameters
	float lambdaTheta = lambdaThetaVar->getVal();
	float lambdaPhi = lambdaPhiVar->getVal();
	float lambdaThetaPhi = lambdaThetaPhiVar->getVal();

	/// get acceptance and efficiency 3D map
	TEfficiency* accMap;
	TEfficiency* effMap;
	TH3D* systEff;
	
	if (!applyEff) {
		/// apply underlying polarization parameters to only acceptance
		accMap = getAcceptance3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);
    	effMap = getEfficiency3DMap(refFrameName.Data(), 0, 0, 0, isPhiFolded);
		systEff = getSysEff3DMap(refFrameName.Data(), 0, 0, 0, isPhiFolded);
	}
	
	else {
		/// apply underlying polarization parameters to only efficiency
		// accMap = getAcceptance3DMap(refFrameName.Data(), 0, 0, 0, isPhiFolded);
		// effMap = getEfficiency3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);
    	// systEff = getSysEff3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);

		/// apply underlying polarization parameters to both acceptance and efficiency maps
		accMap = getAcceptance3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);
    	effMap = getEfficiency3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);
		systEff = getSysEff3DMap(refFrameName.Data(), lambdaTheta, lambdaPhi, lambdaThetaPhi, isPhiFolded);
	}

	/// rebin acceptance and efficiency, efficiency systematic Uncertainty
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	TH2D* hTotalCosThetaPhiAcc = (TH2D*)accMapCosThetaPhi->GetTotalHistogram();

	TH2D* hPassedCosThetaPhiAcc = (TH2D*)accMapCosThetaPhi->GetPassedHistogram();

	TH2D* hTotalCosThetaPhiEff = (TH2D*)effMapCosThetaPhi->GetTotalHistogram();

	TH2D* hPassedCosThetaPhiEff = (TH2D*)effMapCosThetaPhi->GetPassedHistogram();

	TH2D* hRatioCosThetaPhi = (TH2D*)hTotalCosThetaPhiAcc->Clone("hRatioCosThetaPhi");

	hRatioCosThetaPhi->Divide(hPassedCosThetaPhiAcc, hTotalCosThetaPhiEff, 1, 1);
	
	/// draw acceptance and efficiency graph for check
	TCanvas* accCanvas = nullptr;
	TCanvas* effCanvas = nullptr;
	
	TCanvas* accTotalCanvas = nullptr;
	TCanvas* accPassedCanvas = nullptr;
	TCanvas* effTotalCanvas = nullptr;
	TCanvas* effPassedCanvas = nullptr;

	TCanvas* ratioCanvas = nullptr;

	if (!applyEff) {
		accCanvas = DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kTRUE, kFALSE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);
		effCanvas = DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kFALSE, kFALSE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, 0, 0, 0);
	}

	else {
		accCanvas = DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kTRUE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);
		effCanvas = DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kFALSE, kTRUE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);
        
		ratioCanvas = draw2DMap(hRatioCosThetaPhi, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
		display2DMapContents(hRatioCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE, 0.04, kBlack, 4);

		// accTotalCanvas = draw2DMap(hTotalCosThetaPhiAcc, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
        // display2DMapContents(hTotalCosThetaPhiAcc, nCosThetaBins, nPhiBins, kFALSE);
       
        accPassedCanvas = draw2DMap(hPassedCosThetaPhiAcc, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
        display2DMapContents(hPassedCosThetaPhiAcc, nCosThetaBins, nPhiBins, kFALSE);
        
        effTotalCanvas = draw2DMap(hTotalCosThetaPhiEff, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
        display2DMapContents(hTotalCosThetaPhiEff, nCosThetaBins, nPhiBins, kFALSE);

        effPassedCanvas = draw2DMap(hPassedCosThetaPhiEff, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
        display2DMapContents(hPassedCosThetaPhiEff, nCosThetaBins, nPhiBins, kFALSE);

        // accCanvas = DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kTRUE, kFALSE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, lambdaTheta, lambdaPhi, lambdaThetaPhi);
        // effCanvas = DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, gUpsilonState, kFALSE, kFALSE, kFALSE, "_TriggerAcc", isPhiFolded, kTRUE, 0, 0, 0);
	}
	
	TCanvas* dummyCanvas = new TCanvas("dummyCanvas", "dummyCanvas", 600, 600);

	/// apply acc x eff correction weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			Double_t weight = 1;

			// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotalCosThetaPhiAcc->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotalCosThetaPhiAcc->GetYaxis()->GetBinCenter(iPhi + 1);

			Int_t iGlobalBin = hTotalCosThetaPhiAcc->FindFixBin(binCenterCosTheta, binCenterPhi);

			// get the corresponding weights
			double acceptance = accMapCosThetaPhi->GetEfficiency(iGlobalBin);
			double efficiency = effMapCosThetaPhi->GetEfficiency(iGlobalBin);
			double residual = hRatioCosThetaPhi->GetBinContent(iCosTheta + 1, iPhi + 1);
			// double efficiency = 1;

			cout << "acceptance: " << acceptance << endl;
			cout << "efficiency: " << efficiency << endl;

			double relEffUncHigh = 0;
			double relEffUncLow = 0;

			double relAccUncHigh = 0;
			double relAccUncLow = 0;

			/// calculate weight
			/// apply both eff and acc
			if (applyAcc && applyEff) {
				if (acceptance == 0 || efficiency == 0) {
					weight = 0.;
					relEffUncHigh = 0.;
					relAccUncHigh = 0.;
				}

				else {
					// weight = 1. / (acceptance * efficiency);
					weight = 1. / (acceptance * efficiency) * residual;

					relEffUncHigh = effMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
					relAccUncHigh = accMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
				}
			}

			/// apply only eff
			else if (!applyAcc && applyEff) {
				if (efficiency == 0) {
					weight = 0.;
					relEffUncHigh = 0.;
				}
				else {
					weight = 1. / (efficiency);
					relEffUncHigh = effMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
				}
			}

			/// apply only acc
			else if (applyAcc && !applyEff) {
				if (acceptance == 0) {
					weight = 0.;
					relAccUncHigh = 0.;
				}
				else {
					weight = 1. / (acceptance);
					relAccUncHigh = accMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
				}
			}

			/// no detector effects
			else
				weight = 1.;

			// weightMap->SetBinContent(iCosTheta + 1, iPhi + 1, weight);

			// propagate both scale factor uncertainties and efficiency stat errors to the weight
			// double relSystUnc = systEff->GetBinContent(iCosTheta + 1, iPhi + 1) / efficiency;

			// totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

			// totalRelUncHigh = TMath::Hypot(relEffUncHigh, relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(relEffUncLow, relAccUncLow);

			// totalUncHigh = totalRelUncHigh * efficiency * acceptance;
			// totalUncLow = totalRelUncLow * efficiency * acceptance;

			// full relative uncertainty on weight
			double relWeightUnc = TMath::Hypot(relAccUncHigh, relEffUncHigh);
			double deltaW = weight * relWeightUnc;

			double recoMCVal = polarizedHist->GetBinContent(iCosTheta + 1, iPhi + 1);

			double recoMCUnc = polarizedHist->GetBinError(iCosTheta + 1, iPhi + 1);

			double residualVal = hRatioCosThetaPhi->GetBinContent(iCosTheta + 1, iPhi + 1);
			double residualUnc = hRatioCosThetaPhi->GetBinError(iCosTheta + 1, iPhi + 1);

			cout << "recoMCVal: " << recoMCVal << endl;
			cout << "recoMCUnc: " << recoMCUnc << endl;

			// now propagate both sources of uncertainty:
			double totalError = TMath::Sqrt(pow(weight * recoMCUnc, 2) + pow(recoMCVal * deltaW, 2)); // error from acc × eff)

			cout << "weight: " << weight << endl;
			cout << "relEffUncHigh: " << relEffUncHigh << endl;
			cout << "relEffUncLow: " << relEffUncLow << endl;
			cout << "relAccUncHigh: " << relAccUncHigh << endl;
			cout << "relAccUncLow: " << relAccUncLow << endl;
			cout << "relWeightUnc: " << relWeightUnc << endl;
			cout << "deltaW: " << deltaW << endl;
			cout << "residualVal: " << residualVal << endl;
			cout << "residualUnc: " << residualUnc << endl;
			cout << "recoMCUnc: " << recoMCUnc << endl;
			cout << "recoMCVal: " << recoMCVal << endl;

			cout << "totalError: " << totalError << endl;
			/// set the bin contents reflecting weights
			/// yield with acceptance x efficiency correction
			correctedHist->SetBinContent(iCosTheta + 1, iPhi + 1, recoMCVal * weight);
			// correctedHist->SetBinContent(iCosTheta + 1, iPhi + 1, residualVal);

			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, recoMCUnc * weight);
			correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, totalError);
			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(recoMCUnc / recoMCVal, totalRelUncHigh));
			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, residualUnc);

			cout << "correctedHist->GetBinContent(" << iCosTheta + 1 << ", " << iPhi + 1 << ") = " << correctedHist->GetBinContent(iCosTheta + 1, iPhi + 1) << endl;
			cout << "correctedHist->GetBinError(" << iCosTheta + 1 << ", " << iPhi + 1 << ") = " << correctedHist->GetBinError(iCosTheta + 1, iPhi + 1) << endl;
			// /// fill uncertainty histograms
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
		TCanvas* correctedCanvas = draw2DMap(correctedHist, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, kFALSE);
		display2DMapContents(correctedHist, nCosThetaBins, nPhiBins, kFALSE, 0.035, kWhite, 0);

		/// Styles of the texts in the plot
		TLatex* legend1 = new TLatex();
		legend1->SetTextAlign(22);
		legend1->SetTextSize(0.05);

		/// Put texts inside the plot
		legend1->DrawLatexNDC(.50, .88, "Corrected MC");
		legend1->DrawLatexNDC(.48, .80, Form("Input: #lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta0, lambdaPhi0, lambdaThetaPhi0));

		gPad->RedrawAxis();

		gPad->Update();

		gSystem->mkdir(Form("closureTest/%s", applyEff ? "HydjetAccEff": "PythiaAcc"), kTRUE);
		correctedCanvas->SaveAs(Form("closureTest/%s/2Dcorrected_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.pdf", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		correctedCanvas->SaveAs(Form("closureTest/%s/2Dcorrected_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
	}

	gPad->RedrawAxis();

	/// save plots
	gSystem->mkdir(Form("closureTest/%s", applyEff ? "HydjetAccEff": "PythiaAcc"), kTRUE);
	if (isPhiFolded) {
		/// png
		accCanvas->SaveAs(Form("closureTest/%s/2Dacc_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		effCanvas->SaveAs(Form("closureTest/%s/2Deff_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
	} else {
		/// pdf
		accCanvas->SaveAs(Form("closureTest/%s/2Dacc_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.pdf", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		effCanvas->SaveAs(Form("closureTest/%s/2Deff_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.pdf", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		///png
		accCanvas->SaveAs(Form("closureTest/%s/2Dacc_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
		effCanvas->SaveAs(Form("closureTest/%s/2Deff_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_lambdaTheta%.2f_phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi));
	}

	return;
}

/// read the input polarized MC file (recoMC) and fill the 2D angular distribution histogram within the specified binning
void getPolarizedMCHist(TH2D* angDistHist2D, TString refFrameName = "CS",
                        Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0,
                        Int_t ptMin = 0, Int_t ptMax = 30,
                        Int_t nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {},
                        Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {},
                        Bool_t isPhiFolded = kFALSE) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";

	/// read the input MC file

	// const char* inputFileName = Form("../Files/Y1SReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_noMassCut.root", gMuonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi);
	const char* inputFileName = Form("../Files/Y1SReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_gen.root", gMuonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi);
	// const char* inputFileName = Form("../Files/Y1SReconstructedMCWeightedDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f_HydjetWeight.root", gMuonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi);

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

	/// apply kinematic cuts on the dataset
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
	TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins + 5, phiBinEdges[0] - 20, phiBinEdges[nPhiBins] + 100);

	/// histogram to store errors on the weights
	TH2D* errorHist = new TH2D("hdummy", ";cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins]);

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
		Double_t cosThetaVal = ((RooRealVar*)iRooArgSet->find(CosThetaVarName(refFrameName.Data())))->getVal();
		Double_t phiVal = ((RooRealVar*)iRooArgSet->find(PhiVarName(refFrameName.Data())))->getVal();
		// Double_t phiTildeVal = ((RooRealVar*) iRooArgSet->find(PhiTildeVarName(refFrameName.Data())))->getVal();

		/// get the weight and its errors
		Double_t weight = reducedDataset.weight(); // weight = nColl * Gen_weight * dimuonPtWeight * dimuWeight_nominal * polarWeight

		double weightErrorUp = ((RooRealVar*)iRooArgSet->find(Form("errorWeightUp%s", refFrameName.Data())))->getVal();
		double weightErrorDown = ((RooRealVar*)iRooArgSet->find(Form("errorWeightUp%s", refFrameName.Data())))->getVal();

		/// get the bin numbers
		Int_t iCosTheta = angDistHist2D->GetXaxis()->FindBin(cosThetaVal);
		Int_t iPhi = angDistHist2D->GetYaxis()->FindBin(phiVal);

		/// calculate the total errors on the weights
		tempErrorUp[iCosTheta - 1][iPhi - 1] = TMath::Hypot(tempErrorUp[iCosTheta - 1][iPhi - 1], weightErrorUp);
		tempErrorDown[iCosTheta - 1][iPhi - 1] = TMath::Hypot(tempErrorDown[iCosTheta - 1][iPhi - 1], weightErrorDown);

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
		// angDistHist2D->SetBinError(iCosTheta, iPhi, tempErrorUp[iCosTheta - 1][iPhi - 1]);
		// generalPolarTildeHist->Fill(cosThetaVal, phiTildeVal, weight);
	}

	/// fill the error histogram
	for (int iCosTheta = 1; iCosTheta <= nCosThetaBins; iCosTheta++) {
		for (int iPhi = 1; iPhi <= nPhiBins; iPhi++) {
			angDistHist2D->SetBinError(iCosTheta, iPhi, tempErrorUp[iCosTheta - 1][iPhi - 1]);
			// cout << "Bin: (" << iCosTheta << ", " << iPhi << ")" << endl;
			// cout << "BinError: " << angDistHist2D->GetBinError(iCosTheta, iPhi) << endl;
			// cout << "BinValue: " << angDistHist2D->GetBinContent(iCosTheta, iPhi) << endl;
		}
	}

	/// draw polarized MC histogram
	draw2DMap(angDistHist2D, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
	display2DMapContents(angDistHist2D, nCosThetaBins, nPhiBins, kFALSE);

	/// draw error histogram
	// draw2DMap(errorHist, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
	// display2DMapContents(errorHist, nCosThetaBins, nPhiBins, kFALSE);

	/// Styles of the texts in the plot
	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(22);
	legend1->SetTextSize(0.05);

	/// Put texts inside the plot
	legend1->DrawLatexNDC(.50, .88, "Reconstructed MC");
	legend1->DrawLatexNDC(.50, .80, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	gPad->Update();

	// // save the plot
	// gSystem->mkdir("DistributionFitsMC", kTRUE);
	// mc2DCanvas->SaveAs(Form("DistributionFitsMC/RecoMC_Theta%.2f_Phi%.2f_ThetaPhi%.2f.png", lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");

	return;
}

/// get polarized GEN sample using noGENFilter sample (can switch on and off applying acceptance)
void getGENMCHist(TH2D* angDistHist2D, TString refFrameName = "CS",
                  Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0,
                  Int_t ptMin = 0, Int_t ptMax = 30,
                  Int_t nCosThetaBins = 5, const vector<Double_t>& cosThetaBinEdges = {},
                  Int_t nPhiBins = 6, const vector<Double_t>& phiBinEdges = {},
                  Bool_t applyAcc = kTRUE, Bool_t isPhiFolded = kFALSE) {
	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";

	/// read the input MC file (choose applying acceptance or not)
	const char* inputFileName;
	if (applyAcc)
		inputFileName = Form("Y1SGenNoFilterMCDataset_%s_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", gMuonAccName, lambdaTheta, lambdaPhi, lambdaThetaPhi); // Gen + Acc
	else
		inputFileName = Form("Y1SGenNoFilterMCDataset_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", lambdaTheta, lambdaPhi, lambdaThetaPhi); // Pure Gen

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

	/// apply kinematic cuts
	char* kinematicCut;

	// const char* kinematicCut = Form("centrality >= %d && centrality < %d && mass >= %f && mass < %f && rapidity > %f && rapidity < %f && pt > %d && pt < %d && cosTheta%s > %f && cosTheta%s < %f && fabs(%s) > %d && fabs(%s) < %d", 2 * gCentralityBinMin, 2 * gCentralityBinMax, massMin, massMax, gRapidityMin, gRapidityMax, ptMin, ptMax, refFrameName, cosThetaMin, refFrameName, cosThetaMax, PhiVarName(refFrameName), phiMin, PhiVarName(refFrameName), phiMax);
	if (isPhiFolded)
		kinematicCut = Form("(pt > %d) && (pt < %d) && (cosTheta%s > %f) && (cosTheta%s < %f) && (fabs(%s) > %d) && (fabs(%s) < %d)", ptMin, ptMax, refFrameName.Data(), cosThetaBinEdges[0], refFrameName.Data(), cosThetaBinEdges[5], PhiVarName(refFrameName.Data()), (int)phiBinEdges[0], PhiVarName(refFrameName.Data()), (int)phiBinEdges[nPhiBins]);
	else
		kinematicCut = Form("(pt > %d) && (pt < %d) && (cosTheta%s > %f) && (cosTheta%s < %f) && ((%s) > %d) && ((%s) < %d)", ptMin, ptMax, refFrameName.Data(), cosThetaBinEdges[0], refFrameName.Data(), cosThetaBinEdges[nCosThetaBins], PhiVarName(refFrameName.Data()), (int)phiBinEdges[0], PhiVarName(refFrameName.Data()), (int)phiBinEdges[nPhiBins]);

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
	TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins + 5, phiBinEdges[0] - 20, phiBinEdges[nPhiBins] + 100);

	/// loop over the dataset and fill 2D angular distribution histogram
	ULong64_t totEvents = reducedDataset.numEntries();

	for (Int_t iEvent = 0; iEvent < totEvents; iEvent++) {
		// if (iEvent > 100) break;

		/// get the event
		const RooArgSet* iRooArgSet = reducedDataset.get(iEvent);

		/// get the angular variables
		Double_t cosThetaVal = ((RooRealVar*)iRooArgSet->find(CosThetaVarName(refFrameName.Data())))->getVal();
		Double_t phiVal = ((RooRealVar*)iRooArgSet->find(PhiVarName(refFrameName.Data())))->getVal();
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
	TCanvas* mc2DCanvas = draw2DMap(angDistHist2D, refFrameName.Data(), nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
	display2DMapContents(angDistHist2D, nCosThetaBins, nPhiBins, kFALSE);

	/// draw error histogram
	// draw2DMap(errorHist, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1, isPhiFolded);
	// display2DMapContents(errorHist, nCosThetaBins, nPhiBins, kFALSE);

	/// Styles of the texts in the plot
	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(22);
	legend1->SetTextSize(0.05);

	/// Put texts inside the plot
	if (applyAcc)
		legend1->DrawLatexNDC(.50, .88, gMuonPtCutText);
	else
		legend1->DrawLatexNDC(.50, .88, "No muon kinematic selection");

	legend1->DrawLatexNDC(.48, .80, Form("Input: #lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	gPad->RedrawAxis();

	gPad->Update();

	/// save the plot
	gSystem->mkdir("closureTest/PythiaAcc", kTRUE);
	if (isPhiFolded)
		mc2DCanvas->SaveAs(Form("ClosureTest/PythiaAcc/GenMC_%s_pt%dto%d_cosTheta%.2fto%.2f_absphi%dto%d_Theta%.2f_Phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
	else {
		/// pdf
		mc2DCanvas->SaveAs(Form("ClosureTest/PythiaAcc/GenMC_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_Theta%.2f_Phi%.2f_ThetaPhi%.2f.pdf", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
		/// png
		mc2DCanvas->SaveAs(Form("ClosureTest/PythiaAcc/GenMC_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_Theta%.2f_Phi%.2f_ThetaPhi%.2f.png", refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
	}

	return;
}

void closureTest(TString refFrameName = "CS",
                 Float_t lambdaTheta0 = 1, Float_t lambdaPhi0 = 0, Float_t lambdaThetaPhi0 = 0,
                 Int_t ptMin = 2, Int_t ptMax = 6,
                 const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7,
                 const Int_t nPhiBins = 6, Int_t phiMin = -180, Int_t phiMax = 180,
                 Bool_t isPhiFolded = kFALSE, Bool_t applyAcc = kTRUE, Bool_t applyEff = kFALSE, int totNItrs = 0) {
	
	writeExtraText = true; // if extra text
	extraText = "       Internal";
	// extraText = "       Simulation Preliminary";
	
	/// set bin edges and width
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	/// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	/// create 2D angular distribution histogram
	TH2D* angDistHist2D = new TH2D("angDistHist2D", "; cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	/// define polarization parameters that will be used for correction (acceptance and efficiency)
	RooRealVar* lambdaTheta = new RooRealVar("lambdaTheta", "lambdaTheta", -3., 3.);
	RooRealVar* lambdaPhi = new RooRealVar("lambdaPhi", "lambdaPhi", -3., 3.);
	RooRealVar* lambdaThetaPhi = new RooRealVar("lambdaThetaPhi", "lambdaThetaPhi", -3., 3.);
	RooRealVar* lambdaTilde = new RooRealVar("lambdaTilde", "lambdaTilde", -3., 3.);

	/// always start with null polarization assumption in acceptance and efficiency
	lambdaTheta->setVal(0);
	lambdaPhi->setVal(0);
	lambdaThetaPhi->setVal(0);
	lambdaTilde->setVal(0);

	/// place holder of the polarization parameters before and after iterations
	vector<double> lambdaThetaArr, lambdaPhiArr, lambdaThetaPhiArr, lambdaTildeArr, numItrArr;
	vector<double> lambdaThetaErrArr, lambdaPhiErrArr, lambdaThetaPhiErrArr, lambdaTildeErrArr, xErrArr;

	/// set the initial values of the polarization parameters for acceptance and efficiency correction
	lambdaThetaArr.push_back(lambdaTheta->getVal());
	lambdaPhiArr.push_back(lambdaPhi->getVal());
	lambdaThetaPhiArr.push_back(lambdaThetaPhi->getVal());
	lambdaTildeArr.push_back(lambdaTilde->getVal());
	numItrArr.push_back(-1);

	lambdaThetaErrArr.push_back(lambdaTheta->getError());
	lambdaPhiErrArr.push_back(lambdaPhi->getError());
	lambdaThetaPhiErrArr.push_back(lambdaThetaPhi->getError());
	lambdaTildeErrArr.push_back(lambdaTilde->getError());
	
	xErrArr.push_back(0.);

	// cout << "lambdaTheta arr size: " << lambdaThetaArr.size() << endl;
	// cout << lambdaThetaArr[0] << ", " << lambdaThetaArr[1] << endl;
	// cout << "lambdaPhi arr size: " << lambdaPhiArr.size() << endl;
	// cout << lambdaPhiArr[0] << ", " << lambdaPhiArr[1] << endl;
	// cout << "lambdaThetaPhi arr size: " << lambdaThetaPhiArr.size() << endl;
	// cout << lambdaThetaPhiArr[0] << ", " << lambdaThetaPhiArr[1] << endl;
	// cout << "lambdaTilde arr size: " << lambdaTildeArr.size() << endl;
	// cout << lambdaTildeArr[0] << ", " << lambdaTildeArr[1] << endl;
	// cout << "numItr arr size: " << numItrArr.size() << endl;
	// cout << numItrArr[0] << ", " << numItrArr[1] << endl;

	/// store lambda parameters
	RooArgSet polarParams(*lambdaTheta, *lambdaPhi, *lambdaThetaPhi, *lambdaTilde);

	/// get MC histo with polarization (genMC without any filters or recoMC)
	if (!applyEff)
		getGENMCHist(angDistHist2D, refFrameName, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, applyAcc, isPhiFolded); /// apply polarization to NoGENFilter sample
	else
		getPolarizedMCHist(angDistHist2D, refFrameName, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges); /// apply polarization to recoMC sample

	/// define histograms
	TH2D* correctedHist = new TH2D("correctedHist", "; cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);
	TH2D* fittedHist = new TH2D("fittedHist", "; cos #theta; #varphi (#circ); Number of generated #varUpsilon(1S) events", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax);

	/// loop over the number of iterations
	int nItrs = 0; /// number of iterations

	// while (!((lambdaTheta0 == (round(lambdaTheta->getVal() * 100.) / 100.)) && (lambdaPhi0 == (round(lambdaPhi->getVal() * 100.) / 100.)) && (lambdaThetaPhi0 == round(lambdaThetaPhi->getVal() * 100. / 100.)))) {
	for (int iItr = 0; iItr <= totNItrs; iItr++) {
		/// correct the MC
		correctMC2DHist(angDistHist2D, correctedHist, refFrameName, polarParams, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kTRUE, isPhiFolded, applyAcc, applyEff, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0);

		/// due to drawing issue, separate correctedHist and fittedHist
		correctMC2DHist(angDistHist2D, fittedHist, refFrameName, polarParams, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, isPhiFolded, applyAcc, applyEff, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0);

		/// extract parameters from the corrected hist
		polarParams = extractPolarParam(fittedHist, refFrameName, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, isPhiFolded, applyEff);
		cout << "extracted polarization parameters: " << endl;
		polarParams.Print("v");

		lambdaThetaArr.push_back(lambdaTheta->getVal());
		lambdaPhiArr.push_back(lambdaPhi->getVal());
		lambdaThetaPhiArr.push_back(lambdaThetaPhi->getVal());
		lambdaTildeArr.push_back(lambdaTilde->getVal());

		numItrArr.push_back(iItr);

		lambdaThetaErrArr.push_back(lambdaTheta->getError());
		lambdaPhiErrArr.push_back(lambdaPhi->getError());
		lambdaThetaPhiErrArr.push_back(lambdaThetaPhi->getError());
		lambdaTildeErrArr.push_back(lambdaTilde->getError());
	
		xErrArr.push_back(0.);

		// cout << "lambdaTheta arr size: " << lambdaThetaArr.size() << endl;
		// cout << lambdaThetaArr[0] << ", " << lambdaThetaArr[1] << endl;
		// cout << "lambdaPhi arr size: " << lambdaPhiArr.size() << endl;
		// cout << lambdaPhiArr[0] << ", " << lambdaPhiArr[1] << endl;
		// cout << "lambdaThetaPhi arr size: " << lambdaThetaPhiArr.size() << endl;
		// cout << lambdaThetaPhiArr[0] << ", " << lambdaThetaPhiArr[1] << endl;
		// cout << "lambdaTilde arr size: " << lambdaTildeArr.size() << endl;
		// cout << lambdaTildeArr[0] << ", " << lambdaTildeArr[1] << endl;
		// cout << "numItr arr size: " << numItrArr.size() << endl;
		// cout << numItrArr[0] << ", " << numItrArr[1] << endl;
	
		
		/// if the difference between the new paramter and the previous one is less than 0.01, stop the iterative procedure
		if ((fabs(lambdaThetaArr[iItr] - lambdaThetaArr[iItr + 1]) < 0.01) && (fabs(lambdaPhiArr[iItr] - lambdaPhiArr[iItr + 1]) < 0.01) && (fabs(lambdaThetaPhiArr[iItr] - lambdaThetaPhiArr[iItr + 1]) < 0.01)) {
			// cout << "converged after " << iItr << " iterations" << endl;
			// cout << "lambdaTheta diff: " << lambdaThetaArr[iItr] << ", " << lambdaThetaArr[iItr + 1] << ", " << fabs(lambdaThetaArr[iItr + 1] - lambdaThetaArr[iItr]) << endl;
			// cout << "lambdaPhi diff: " << lambdaPhiArr[iItr] << ", " << fabs(lambdaPhiArr[iItr + 1] - lambdaPhiArr[iItr]) << endl;
			// cout << "lambdaThetaPhi diff: " << lambdaThetaPhiArr[iItr] << ", " << fabs(lambdaThetaPhiArr[iItr + 1] - lambdaThetaPhiArr[iItr]) << endl;

			// break;
		}

		nItrs++;

	}

	/// draw plots for polarization parameters vs # of Iterations
	TCanvas* iterCanvas = new TCanvas("iterCanvas", "iterCanvas", 650, 600);
	iterCanvas->SetRightMargin(0.05);
	iterCanvas->cd();

	TGraphErrors* lambdaThetaGraph = new TGraphErrors(lambdaThetaArr.size(), numItrArr.data(), lambdaThetaArr.data(), xErrArr.data(), lambdaThetaErrArr.data());
	lambdaThetaGraph->SetTitle("; Number of Iterations; Polarization Parameter Values");
	lambdaThetaGraph->GetXaxis()->CenterTitle();
	lambdaThetaGraph->GetYaxis()->CenterTitle();
	lambdaThetaGraph->SetMarkerStyle(20);
	lambdaThetaGraph->SetMarkerColor((TColor::GetColor("#FF1F5B")));
	lambdaThetaGraph->SetLineColor((TColor::GetColor("#FF1F5B")));
	lambdaThetaGraph->SetLineWidth(2);
	lambdaThetaGraph->SetMarkerSize(1);
	lambdaThetaGraph->SetMinimum(-1.2);
	lambdaThetaGraph->SetMaximum(1.5);

	gPad->Update(); // ensures the histogram behind the graph is created

	lambdaThetaGraph->GetXaxis()->SetNdivisions(nItrs + 2);
	lambdaThetaGraph->Draw("APL");

	TGraphErrors* lambdaPhiGraph = new TGraphErrors(lambdaPhiArr.size(), numItrArr.data(), lambdaPhiArr.data(), xErrArr.data(), lambdaPhiErrArr.data());
	lambdaPhiGraph->SetMarkerStyle(20);
	lambdaPhiGraph->SetMarkerColor(TColor::GetColor("#009ADE"));
	lambdaPhiGraph->SetLineColor(TColor::GetColor("#009ADE"));
	lambdaPhiGraph->SetLineWidth(2);
	lambdaPhiGraph->SetMarkerSize(1);
	lambdaPhiGraph->Draw("PL same");

	TGraphErrors* lambdaThetaPhiGraph = new TGraphErrors(lambdaThetaPhiArr.size(), numItrArr.data(), lambdaThetaPhiArr.data(), xErrArr.data(), lambdaThetaPhiErrArr.data());
	lambdaThetaPhiGraph->SetMarkerStyle(20);
	lambdaThetaPhiGraph->SetMarkerColor(TColor::GetColor("#bababa"));
	lambdaThetaPhiGraph->SetLineColor(TColor::GetColor("#bababa"));
	lambdaThetaPhiGraph->SetLineWidth(2);
	lambdaThetaPhiGraph->SetMarkerSize(1);
	lambdaThetaPhiGraph->Draw("PL same");

	if (nItrs == 0) {
		drawLine(-1.1, 0, 0, 0);
		drawLine(-1.1, lambdaTheta0, 0, lambdaTheta0);
		drawLine(-1.1, lambdaPhi0, 0, lambdaPhi0);
		drawLine(-1.1, lambdaThetaPhi0, 0, lambdaThetaPhi0);
	}

	else {
		cout << "nItrs: " << nItrs << endl;
		drawLine(-1.57, 0, (nItrs - 1) * 1.1, 0);
		drawLine(-1.57, lambdaTheta0, (nItrs - 1) * 1.1, lambdaTheta0);
		drawLine(-1.57, lambdaPhi0, (nItrs - 1) * 1.1, lambdaPhi0);
		drawLine(-1.57, lambdaThetaPhi0, (nItrs - 1) * 1.1, lambdaThetaPhi0);
	}

	TLegend* legend = new TLegend(0.16, 0.78, 0.46, 0.91);
	legend->SetBorderSize(0);
	legend->SetFillColor(0);
	legend->SetTextSize(0.05);

	legend->AddEntry(lambdaThetaGraph, "#lambda_{#theta}", "p");
	legend->AddEntry(lambdaPhiGraph, "#lambda_{#varphi}", "p");
	legend->AddEntry(lambdaThetaPhiGraph, "#lambda_{#theta#varphi}", "p");
	legend->Draw();

	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(23);
	legend1->SetTextSize(0.04);
	legend1->DrawLatexNDC(.65, .84, Form("Input: (#lambda_{#theta},  #lambda_{#varphi}, #lambda_{#theta#varphi}) = (%.2f, %.2f, %.2f)", lambdaTheta0, lambdaPhi0, lambdaThetaPhi0));

	Bool_t isCSframe;

	if (refFrameName == (TString) "CS")
		isCSframe = true;
	else
		isCSframe = false;

	TPaveText* text = new TPaveText(0.38, 0.85, 0.95, 0.90, "NDCNB");
	text->SetFillColor(4000);
	text->SetBorderSize(0);
	text->AddText(Form("%s, %s", isCSframe ? "Collins-Soper frame" : "Helicity frame", DimuonPtRangeText(ptMin, ptMax)));
	// text->AddText(DimuonPtRangeText(ptMin, ptMax));
	text->SetAllWith("", "align", 32);
	text->Draw("same");

	gPad->Update();   // important to get pad and axis ranges

	// Get y-axis range (used for placing the box slightly below the axis)
	double yMin = gPad->GetUymin();
	double yMax = gPad->GetUymax();
	double yRange = yMax - yMin;
	
	// Define box size and position around x = -1
	double box_x1 = -1.2;
	double box_x2 = -0.8;
	double box_y1 = yMin - 0.06 * yRange;  // slightly below visible axis
	double box_y2 = yMin + 0.02 * yRange;  // near the bottom of plot
	
	// Draw box to cover the default label
	TBox* cover = new TBox(box_x1, box_y1, box_x2, box_y2);
	cover->SetFillColor(kWhite);
	cover->SetLineColor(kWhite);
	cover->SetFillStyle(1001);   // solid fill
	cover->Draw("same");

	// Add custom label at x = -1
	TLatex* t = new TLatex(-1, -1.37, "Initial");
	t->SetTextAlign(21);  // center align
	t->SetTextSize(0.047);
	t->Draw("same");

	// TPaveText* refFrameText = RefFrameText(isCSframe, cosThetaMin, cosThetaMax, phiMin, phiMax, 0.61, 0.34, 0.97, 0.56);
	// TPaveText* kinematicsText = KinematicsText(gCentralityBinMin, gCentralityBinMax, ptMin,  ptMax, 0.57, 0.95, 0.95, 0.65);

	// TPaveText* refFrameText = new TPaveText( 0.61, 0.34, 0.97, 0.56, "NDCNB");
	// refFrameText->SetFillColor(4000);
	// refFrameText->SetBorderSize(0);
	// refFrameText->AddText(isCSframe ? "Collins-Soper frame" : "Helicity frame");

	// refFrameText->Draw("same");

	// TPaveText* kinematicsText = new TPaveText(0.57, 0.95, 0.95, 0.65, "NDCNB");
	// kinematicsText->SetFillColor(4000);
	// kinematicsText->SetBorderSize(0);
	// // kinematicsText->AddText(CentralityRangeText(gCentralityBinMin, gCentralityBinMax));
	// // kinematicsText->AddText(gMuonPtCutText);
	// // kinematicsText->AddText(DimuonRapidityRangeText(gRapidityMin, gRapidityMax));
	// kinematicsText->AddText(DimuonPtRangeText(ptMin, ptMax));

	// kinematicsText->SetAllWith("", "align", 32);
	// kinematicsText->Draw("same");

	gPad->RedrawAxis();

	iterCanvas->SaveAs(Form("ClosureTest/%s/Iterations_totItr%d_%s_pt%dto%d_cosTheta%.2fto%.2f_phi%dto%d_LambdaTheta%.2f_Phi%.2f_ThetaPhi%.2f.png", applyEff ? "HydjetAccEff": "PythiaAcc", totNItrs, refFrameName.Data(), ptMin, ptMax, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], (int)phiBinEdges[0], (int)phiBinEdges[nPhiBins], lambdaTheta0, lambdaPhi0, lambdaThetaPhi0), "RECREATE");

	return;
}
