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

void rawYield_2D_customizedFits_sysUnc(Int_t ptMin = 0, Int_t ptMax = 2, const char* refFrameName = "CS", const Int_t nCosThetaBins = 5, Double_t cosThetaMin = -0.7, Double_t cosThetaMax = 0.7, const Int_t nPhiBins = 3, Int_t phiMin = 0, Int_t phiMax = 180, Int_t iState = gUpsilonState, Bool_t LEGOplot = kTRUE, const char* defaultBkgShapeName = "ExpTimesErr", int chooseYield = 0, double errPercent = 1, int chooseEff = 0) {
	writeExtraText = true; // if extra text
	extraText = "      Internal";

	using namespace RooFit;
	//using namespace RooStats;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	/// Bin edges and width
	// Set the bin edges along the cosTheta/phi axis depending on the number of bins, min and max values
	// (If want to use non-uniform bin width, the bin edges should be pre-defined in PolarFitHelpers.h)
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);

	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

    /// Read the yields from the signal extraction / pseudo-experiment file
	RooRealVar* yield1S = new RooRealVar("yield1S", "", 1000);
	RooRealVar* yield2S = new RooRealVar("yield2S", "", 100);
	RooRealVar* yield3S = new RooRealVar("yield3S", "", 10);

    /// assign signal and background shape name to read the file for the pseudo-experiment yield difference results
	const char* nominalSignalShapeName = "SymDSCB";
    const char* nominalBkgShapeName = "ExpTimesErr";

    const char* altSignalShapeName = "Johnson";
    const char* altBkgShapeName = "ChebychevOrder3";

	/// assign signal and background shape name to read the file for the yield extraction results
    const char* signalShapeName = "SymDSCB";

	/// background shape array: ChebychevOrderN or ExpTimesErr (this part is for exception backround shape handling)
	const Int_t nCosThetaBinsMax = 20;
	const Int_t nPhiBinsMax = 10;

	std::string bkgShapeName[nCosThetaBinsMax][nPhiBinsMax];

	/// fill the background shape array with ChebychevOrder2
	if (strcmp(defaultBkgShapeName, "Chebychev") == 0) {
		std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ChebychevOrder2");
	}

    /// fill the background shape array with ExpTimesErr
	if (strcmp(defaultBkgShapeName, "ExpTimesErr") == 0) {
		// fill the background shape array with ExpTimesErr
		std::fill(&bkgShapeName[0][0], &bkgShapeName[0][0] + nCosThetaBinsMax * nPhiBinsMax, "ExpTimesErr");
	}

	/// define a TH2D for the yield map before any correction
	TH2D* yieldMap = new TH2D("yieldMap", " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
    
    /// define a TH2D for the yield map after correction
    /// "Standard" procedure means extract the yields per bin
	TH2D* standardCorrectedMap = new TH2D("standardCorrectedMap", " ", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

    /// Read accepatance and efficiecy files for correction
	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

    /// polarization weight of the acceptance and efficiency maps (for now all 0)
	Double_t lambdaTheta = 0., lambdaPhi = 0., lambdaThetaPhi = 0.;

	/// get acceptance maps
	TFile* acceptanceFile = openFile(Form("../MonteCarlo/AcceptanceMaps/%dS/AcceptanceResults%s.root", iState, gMuonAccName));
	auto* accMap = (TEfficiency*)acceptanceFile->Get(nominalMapName);

	/// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	/// get efficiency maps
	TFile* efficiencyFile = openFile(Form("../MonteCarlo/EfficiencyMaps/%dS/EfficiencyResults%s.root", iState, gMuonAccName));
	auto* effMap = (TEfficiency*)efficiencyFile->Get(nominalMapName);

    auto* trkSysUpMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_trk_systUp_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
    auto* trkSysDownMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_trk_systDown_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	auto* muIDSysUpMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_muId_systUp_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
    auto* muIDSysDownMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_muId_systDown_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
    auto* trigSysUpMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_trig_systUp_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
    auto* trigSysDownMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_trig_systDown_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	
	auto* trkStatUpMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_trk_statUp_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	auto* trkStatDownMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_trk_statDown_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	auto* muIDStatUpMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_muId_statUp_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	auto* muIDStatDownMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_muId_statDown_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	auto* trigStatUpMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_trig_statUp_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	auto* trigStatDownMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_trig_statDown_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));

    auto* totalSysUpMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_total_systUp_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));   
    auto* totalSysDownMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_total_systDown_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	auto* totalStatUpMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_total_statUp_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));
	auto* totalStatDownMap = (TEfficiency*)efficiencyFile->Get(Form("h%s_total_statDown_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi));

	/// rebin efficiency maps based on costheta, phi, and pT selection
	cout << "rebinning nominal efficiency map" << endl;
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	cout << "rebinning trkSysUp efficiency map" << endl;
    TEfficiency* trkSysUpCosThetaPhi = rebinTEff3DMap(trkSysUpMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
    
	cout << "rebinning trkSysDown efficiency map" << endl;
	TEfficiency* trkSysDownCosThetaPhi = rebinTEff3DMap(trkSysDownMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
    TEfficiency* muIDSysUpCosThetaPhi = rebinTEff3DMap(muIDSysUpMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
    TEfficiency* muIDSysDownCosThetaPhi = rebinTEff3DMap(muIDSysDownMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
    TEfficiency* trigSysUpCosThetaPhi = rebinTEff3DMap(trigSysUpMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
    TEfficiency* trigSysDownCosThetaPhi = rebinTEff3DMap(trigSysDownMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	TEfficiency* trkStatUpCosThetaPhi = rebinTEff3DMap(trkStatUpMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* trkStatDownCosThetaPhi = rebinTEff3DMap(trkStatDownMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* muIDStatUpCosThetaPhi = rebinTEff3DMap(muIDStatUpMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* muIDStatDownCosThetaPhi = rebinTEff3DMap(muIDStatDownMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* trigStatUpCosThetaPhi = rebinTEff3DMap(trigStatUpMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* trigStatDownCosThetaPhi = rebinTEff3DMap(trigStatDownMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

    TEfficiency* totalSysUpCosThetaPhi = rebinTEff3DMap(totalSysUpMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
    TEfficiency* totalSysDownCosThetaPhi = rebinTEff3DMap(totalSysDownMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	TEfficiency* totalStatUpCosThetaPhi = rebinTEff3DMap(totalStatUpMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* totalStatDownCosThetaPhi = rebinTEff3DMap(totalStatDownMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

    cout << "nCosThetaBins: " << nCosThetaBins << endl;
    cout << "nPhiBins: " << nPhiBins << endl;
    cout << "cosThetaBinEdges: " << cosThetaBinEdges[0] << endl;
    cout << "cosThetaBinEdges: " << cosThetaBinEdges[1] << endl;
    cout << "cosThetaBinEdges: " << cosThetaBinEdges[2] << endl;
    // cout << "cosThetaBinEdges: " << cosThetaBinEdges[4] << endl;
    cout << "phiBinEdges: " << phiBinEdges[0] << endl;
    cout << "phiBinEdges: " << phiBinEdges[1] << endl;

	// get relative systematic uncertainty of efficiency
	auto* systEff = (TH3D*)efficiencyFile->Get(SystTEfficiency3DName(refFrameName));

	// rebin uncertainty map based on costheta, phi, and pT selection
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// draw acc and eff histograms to check if the rebinning works well
	DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kTRUE, kTRUE);

	DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
   
    DrawEfficiency2DHist(trkSysUpCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
    DrawEfficiency2DHist(trkSysDownCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
    DrawEfficiency2DHist(muIDSysUpCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
    DrawEfficiency2DHist(muIDSysDownCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
    DrawEfficiency2DHist(trigSysUpCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
    DrawEfficiency2DHist(trigSysDownCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);

	DrawEfficiency2DHist(trkStatUpCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
	DrawEfficiency2DHist(trkStatDownCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
	DrawEfficiency2DHist(muIDStatUpCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
	DrawEfficiency2DHist(muIDStatDownCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
	DrawEfficiency2DHist(trigStatUpCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
	DrawEfficiency2DHist(trigStatDownCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);

    DrawEfficiency2DHist(totalSysUpCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
    DrawEfficiency2DHist(totalSysDownCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);

	DrawEfficiency2DHist(totalStatUpCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);
	DrawEfficiency2DHist(totalStatDownCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kTRUE);

    /// get the total histogram of the acceptance map for the binning in the loop below
	TH2D* hTotalCosThetaPhi = (TH2D*)accMapCosThetaPhi->GetTotalHistogram();

	/// define a histogram to draw weight map
	TH2D* weightMap = new TH2D("weightMap", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	/// define histograms to draw uncertainty plots
	TH2D* relSystEffCosThetaPhi = new TH2D("relSystEffCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* statHighEffCosThetaPhi = new TH2D("statHighEffCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* statLowEffCosThetaPhi = new TH2D("statLowEffCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* statHighAccCosThetaPhi = new TH2D("statHighAccCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* statLowAccCosThetaPhi = new TH2D("statLowAccCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* yield1SUncCosThetaPhi = new TH2D("yield1SUncCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());
	TH2D* totalRelUncCosThetaPhi = new TH2D("totalRelUncCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	TH2D* yield1SSysUncCosThetaPhi = new TH2D("yield1SSysUncCosThetaPhi", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	const char* efficiencyName = "";

    const char* pseudoExpModelName = "";

	/// Apply weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			Double_t weight = 1;

			/// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotalCosThetaPhi->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotalCosThetaPhi->GetYaxis()->GetBinCenter(iPhi + 1);

			Int_t iGlobalBin = hTotalCosThetaPhi->FindFixBin(binCenterCosTheta, binCenterPhi);

			/// get the acceptance values
			double acceptance = accMapCosThetaPhi->GetEfficiency(iGlobalBin);
    
            /// get the efficiency values
			double efficiency = 0;

            if (chooseEff == 0) {efficiency = effMapCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "nominal";}
            
            else if (chooseEff == 1) {efficiency = trkSysUpCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "trkSysUp";} 
            else if (chooseEff == 2) {efficiency = trkSysDownCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "trkSysDown";}
            else if (chooseEff == 3) {efficiency = muIDSysUpCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "muIDSysUp";}
            else if (chooseEff == 4) {efficiency = muIDSysDownCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "muIDSysDown";}
            else if (chooseEff == 5) {efficiency = trigSysUpCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "trigSysUp";}
            else if (chooseEff == 6) {efficiency = trigSysDownCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "trigSysDown";}

            else if (chooseEff == 7) {efficiency = totalSysUpCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "totalSysUp";}
            else if (chooseEff == 8) {efficiency = totalSysDownCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "totalSysDown";}

            else if (chooseEff == 9) {efficiency = trkStatUpCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "trkStatUp";} 
            else if (chooseEff == 10) {efficiency = trkStatDownCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "trkStatDown";}
            else if (chooseEff == 11) {efficiency = muIDStatUpCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "muIDStatUp";}
            else if (chooseEff == 12) {efficiency = muIDStatDownCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "muIDStatDown";}
            else if (chooseEff == 13) {efficiency = trigStatUpCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "trigStatUp";}
            else if (chooseEff == 14) {efficiency = trigStatDownCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "trigStatDown";}

            else if (chooseEff == 15) {efficiency = totalStatUpCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "totalStatUp";}
            else if (chooseEff == 16) {efficiency = totalStatDownCosThetaPhi->GetEfficiency(iGlobalBin); efficiencyName = "totalStatDown";}

			/// calculate weight
			if (acceptance == 0 || efficiency == 0)
                /// if acceptance or efficiency is 0, set the weight to 0
				weight = 0.;
			else
				weight = 1. / (acceptance * efficiency);
			
			cout << "weight: " << weight << endl;

            /// fill the weight map
			weightMap->SetBinContent(iCosTheta + 1, iPhi + 1, weight);

            /// Error propagation
			/// propagate both scale factor uncertainties and efficiency stat errors to the weight
			double relSystUnc = systEffCosThetaPhi->GetBinContent(iCosTheta + 1, iPhi + 1) / efficiency;

			double relEffUncHigh = effMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
			double relEffUncLow = effMapCosThetaPhi->GetEfficiencyErrorLow(iGlobalBin) / efficiency;

			double relAccUncHigh = accMapCosThetaPhi->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
			double relAccUncLow = accMapCosThetaPhi->GetEfficiencyErrorLow(iGlobalBin) / acceptance;

			// totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

            /// apply only the statistical uncertainties of the efficiency and acceptance
			totalRelUncHigh = TMath::Hypot(relEffUncHigh, relAccUncHigh);
			totalRelUncLow = TMath::Hypot(relEffUncLow, relAccUncLow);

			totalUncHigh = totalRelUncHigh * efficiency * acceptance;
			totalUncLow = totalRelUncLow * efficiency * acceptance;

            /// Read yield difference uncertainties (mean values) from pseudo-experiments
            double yieldDiffMean_nominal = readYieldDiffMean(ptMin, ptMax, refFrameName, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[iPhi], (Int_t)phiBinEdges[iPhi + 1], nominalSignalShapeName, nominalBkgShapeName);
            double yieldDiffMean_altBkg = readYieldDiffMean(ptMin, ptMax, refFrameName, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[iPhi], (Int_t)phiBinEdges[iPhi + 1], nominalSignalShapeName, altBkgShapeName);
            double yieldDiffMean_altSig = readYieldDiffMean(ptMin, ptMax, refFrameName, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[iPhi], (Int_t)phiBinEdges[iPhi + 1], altSignalShapeName, nominalBkgShapeName);

            /// calculate the systematic uncertainty
            // double systUnc_fitModels = sqrt(pow(yieldDiffMean_nominal, 2) + pow(yieldDiffMean_altBkg, 2) + pow(yieldDiffMean_altSig, 2));

            // cout << "bin: " << refFrameName << ", pT " << ptMin << " to " << ptMax << ", cosTheta "<< cosThetaBinEdges[iCosTheta] << " to " << cosThetaBinEdges[iCosTheta + 1] << ", phi " << phiBinEdges[iPhi] << " to " << phiBinEdges[iPhi + 1] << endl;
            // cout << "systUnc from fit model: " << systUnc_fitModels << endl;

			/// Get yields and their uncertainties
			Int_t absiPhi = yieldMap->GetYaxis()->FindBin(fabs(binCenterPhi)) - 1;
            // cout << "absiPhi: " << absiPhi << endl;

			const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaBinEdges[iCosTheta], cosThetaBinEdges[iCosTheta + 1], (Int_t)phiBinEdges[absiPhi], (Int_t)phiBinEdges[absiPhi + 1]);

			RooArgSet signalYields = GetSignalYields(yield1S, yield2S, yield3S, Form("%s", bkgShapeName[iCosTheta][absiPhi].c_str()), fitModelName, Form("%s", gMuonAccName));
			// cout << "signalYields: " << signalYields << endl;

            /// get the yield value
			double yield1SVal = 0;

            /// shift the yield value based on the pseudo-experiment model
			if (chooseYield == 0) yield1SVal = (yield1S->getVal());
            else if (chooseYield == 1) {yield1SVal = (yield1S->getVal()) - yieldDiffMean_nominal; pseudoExpModelName = "nominal";}
            else if (chooseYield == 2) {yield1SVal = (yield1S->getVal()) - yieldDiffMean_altBkg; pseudoExpModelName = "altBkg";}
            else if (chooseYield == 3) {yield1SVal = (yield1S->getVal()) - yieldDiffMean_altSig; pseudoExpModelName = "altSig";}
            // cout << "yield1SVal: " << yield1SVal << endl;

            /// get the yield uncertainty
			double yield1SUnc = (yield1S->getError()) * errPercent * 0.01; // apply relartive error of statistical uncertainty to the yield

			/// Set the bin contents reflecting weights
			/// only raw yield itself before correction
			yieldMap->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal);

            /// Set the bin error of the raw yield
  			// yieldMap->SetBinError(iCosTheta + 1, yield1SErr * weight);
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh));
			// yieldMap->SetBinError(iCosTheta + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));

			yieldMap->SetBinError(iCosTheta + 1, iPhi + 1, yield1SUnc);

			// yieldMap->SetBinError(iCosTheta + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

			// yieldMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			/// yield with acceptance x efficiency correction
			standardCorrectedMap->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SVal * weight);

            /// set the bin error of the corrected yield
			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);
			
            standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);

            // standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, systUnc_fitModels * weight);

			// fill uncertainty histograms
			relSystEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relSystUnc);
			// statHighEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relEffUncHigh);
			// statLowEffCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relEffUncLow);
			// statHighAccCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relAccUncHigh);
			// statLowAccCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, relAccUncLow);
			// yield1SUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, yield1SUnc / yield1SVal);
			// totalRelUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

            // yield1SSysUncCosThetaPhi->SetBinContent(iCosTheta + 1, iPhi + 1, systUnc_fitModels / yield1SVal);

			// if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight / 2.;
			if ((yield1S->getVal()) * weight > maxYield) maxYield = (yield1S->getVal()) * weight;
		}
	}

	/// Polarization fit
	/// with Root Fit function (E: minos on, S: save results, V: Verbose, I: integral, M: imporve algorithm, R: specificed range)

	// TVirtualFitter::SetDefaultFitter("Minuit");

	/// draw 1/(acceptance x efficiency) map

	TCanvas* weightCanvas = draw2DMap(weightMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE);

	weightMap->GetZaxis()->SetTitle("1 / (acc x eff)");
	weightMap->GetZaxis()->SetTitleOffset(1.);

	display2DMapContents(weightMap, nCosThetaBins, nPhiBins, kFALSE);

	TPaveText* kinematicsText = new TPaveText(0.15, 0.86, 0.77, 0.95, "NDCNB");
	kinematicsText->SetFillColor(4000);
	kinematicsText->SetBorderSize(0);
	kinematicsText->AddText(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	kinematicsText->SetAllWith("", "align", 12);
	kinematicsText->Draw("SAME");

	weightCanvas->Modified();
	weightCanvas->Update();

	/// draw yield map before applying corrections

	TCanvas* yieldCanvas = draw2DMap(yieldMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, LEGOplot);

	yieldMap->GetZaxis()->SetTitle(Form("#varUpsilon(%dS) yields", iState));
	yieldMap->GetZaxis()->SetTitleOffset(1.);

	yieldMap->GetZaxis()->SetRangeUser(0, maxYield * 2);
	// yieldMap->SetMinimum(0);

	// yieldCanvas->SetLogz(); // useful when the one of the bins has an exceptionally high value :')

	TPaveText* kinematicsText_2D = new TPaveText(0.17, 0.78, 0.79, 0.88, "NDCNB");
	kinematicsText_2D ->SetFillColor(4000);
	kinematicsText_2D ->SetBorderSize(0);
	kinematicsText_2D ->AddText(Form("%s, %s", CentralityRangeText(gCentralityBinMin, gCentralityBinMax), DimuonPtRangeText(ptMin, ptMax)));
	kinematicsText_2D ->SetAllWith("", "align", 12);

	if (LEGOplot) kinematicsText->Draw("SAME");

	else if (!LEGOplot) {
		display2DMapContents(yieldMap, nCosThetaBins, nPhiBins, kTRUE);
		kinematicsText_2D ->Draw("SAME");
	}

	yieldCanvas->Modified();
	yieldCanvas->Update();

	/// draw yield map corrected by acceptance and efficiency

	TCanvas* correctedMapCanvas = draw2DMap(standardCorrectedMap, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, LEGOplot);

	standardCorrectedMap->GetZaxis()->SetTitle("Corrected #varUpsilon(1S) Yields");
	standardCorrectedMap->GetZaxis()->SetTitleSize(0.055);
	standardCorrectedMap->GetZaxis()->SetTitleOffset(1.4);
	standardCorrectedMap->GetZaxis()->SetLabelSize(0.044);
	
	standardCorrectedMap->GetYaxis()->SetTitleSize(0.055);
	standardCorrectedMap->GetYaxis()->SetTitleOffset(1.6);
	standardCorrectedMap->GetYaxis()->SetLabelSize(0.044);
	standardCorrectedMap->GetYaxis()->SetLabelOffset(0);

	standardCorrectedMap->GetXaxis()->SetTitleSize(0.055);	
	standardCorrectedMap->GetXaxis()->SetTitleOffset(1.2);
	standardCorrectedMap->GetXaxis()->SetLabelSize(0.044);

	if (!LEGOplot) standardCorrectedMap->GetZaxis()->SetTitleOffset(1.2);

	standardCorrectedMap->GetZaxis()->SetRangeUser(0, maxYield * 2);

	// standardCorrectedMap->SetMinimum(0);

	if (LEGOplot) kinematicsText->Draw("SAME");

	else if (!LEGOplot) {
		display2DMapContents(standardCorrectedMap, nCosThetaBins, nPhiBins, kTRUE);
		kinematicsText_2D->Draw("SAME");
	}
	/// Polarization fit!!!

	TF2* polarFunc2D = getGeneralPolarFunc(maxYield);

	if (LEGOplot) {
		TFitResultPtr fitResults = standardCorrectedMap->Fit("generalPolarFunc", "ESVIMR");

		// Fit results

		double chi2 = fitResults->Chi2();
		double nDOF = nCosThetaBins * nPhiBins - polarFunc2D->GetNpar();

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

		RooArgSet* savedParams = new RooArgSet(lambdaTheta, lambdaPhi, lambdaThetaPhi, lambdaTilde);

		const char* fitModelName = GetFitModelName(signalShapeName, ptMin, ptMax, refFrameName, cosThetaMin, cosThetaMax, phiMin, phiMax);

		SavePolarizationFitParameters(savedParams, Form("rootFit_SysFitModels_%s_%s", pseudoExpModelName, efficiencyName), fitModelName, defaultBkgShapeName);

		// TLegend legend2(.17, .55, .28, .84);
		// legend2.SetTextSize(.05);
		// legend2.AddEntry(standardCorrectedMap, "#varUpsilon(1S) corrected yield", "lp");
		// legend2.AddEntry(polarFunc2D, Form("distribution fit: #lambda_{#theta} = %.2f #pm %.2f", lambdaThetaVal, lambdaThetaErr), "l");
		// legend2.AddEntry((TObject*)0, Form("                       #lambda_{#varphi} = %.2f #pm %.2f", lambdaPhiVal, lambdaPhiErr), "");
		// legend2.AddEntry((TObject*)0, Form("                       #lambda_{#theta#varphi} = %.2f #pm %.2f", lambdaThetaPhiVal, lambdaThetaPhiErr), "");
		// legend2.AddEntry((TObject*)0, Form("                       #tilde{#lambda} = %.2f #pm %.2f", lambdaTildeVal, lambdaTildeErr), "");
		// legend2.AddEntry((TObject*)0, Form("                       n  = %.2f #pm %.2f", normVal, normErr), "");

		// legend2.DrawClone();

		TLegend legend2(.17, .76, .23, .87);
		legend2.SetTextSize(.045);
		legend2.SetFillColor(0);
		legend2.SetFillStyle(1001);

		legend2.AddEntry(standardCorrectedMap, "#varUpsilon(1S) corrected yield", "lp");
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
		textChi2.SetTextAlign(12);
		textChi2.SetTextSize(0.045);
		textChi2.DrawLatexNDC(0.70, 0.043, Form("#chi^{2} / n_{dof} = %.2f", chi2 / nDOF));
	}

	gPad->Update();

	// // /// draw uncertainty 2D plots

	// // // statistical uncertainty of acceptance up
	// // TCanvas* statHighAccCanvas = draw2DMap(statHighAccCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	// // statHighAccCosThetaPhi->GetZaxis()->SetTitle("stat uncer high of acceptance");
	// // statHighAccCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	// // kinematicsText->Draw("SAME");

	// // display2DMapContents(statHighAccCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	// // statHighAccCanvas->Modified();
	// // statHighAccCanvas->Update();

	// // // statistical uncertainty of acceptance down
	// // TCanvas* statLowAccCanvas = draw2DMap(statLowAccCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	// // statLowAccCosThetaPhi->GetZaxis()->SetTitle("stat uncer low of acceptance");
	// // statLowAccCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	// // kinematicsText->Draw("SAME");

	// // display2DMapContents(statLowAccCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	// // statLowAccCanvas->Modified();
	// // statLowAccCanvas->Update();

	// // // statistical uncertainty of efficiency up
	// // TCanvas* statHighEffCanvas = draw2DMap(statHighEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	// // statHighEffCosThetaPhi->GetZaxis()->SetTitle("stat uncer high of efficiency");
	// // statHighEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	// // kinematicsText->Draw("SAME");

	// // display2DMapContents(statHighEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	// // statHighEffCanvas->Modified();
	// // statHighEffCanvas->Update();

	// // /// statistical uncertainty of efficiency down
	// // TCanvas* statLowEffCanvas = draw2DMap(statLowEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	// // statLowEffCosThetaPhi->GetZaxis()->SetTitle("stat uncer low of efficiency");
	// // statLowEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	// // kinematicsText->Draw("SAME");

	// // display2DMapContents(statLowEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	// // statLowEffCanvas->Modified();
	// // statLowEffCanvas->Update();

	// /// systematic uncertainty of efficiency (muon scale factor)
	// TCanvas* systEffCanvas = draw2DMap(relSystEffCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	// relSystEffCosThetaPhi->GetZaxis()->SetTitle("syst uncer of efficiency (muon SF)");
	// relSystEffCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	// kinematicsText->Draw("SAME");

	// display2DMapContents(relSystEffCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	// systEffCanvas->Modified();
	// systEffCanvas->Update();

	// // /// systematic uncertainty of efficiency (muon scale factor)
	// // TCanvas* systEffFitModelCanvas = draw2DMap(yield1SSysUncCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

    // // yield1SSysUncCosThetaPhi->GetZaxis()->SetTitle("syst uncertainty (fit models)");
    // // yield1SSysUncCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	// // kinematicsText->Draw("SAME");

	// // display2DMapContents(yield1SSysUncCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	// // systEffFitModelCanvas->Modified();
	// // systEffFitModelCanvas->Update();

	// // /// yield extraction uncertainty
	// // TCanvas* yieldUncCanvas = draw2DMap(yield1SUncCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	// // yield1SUncCosThetaPhi->GetZaxis()->SetTitle("raw yield statistical uncertainty");
	// // //yield1SUncCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	// // kinematicsText->Draw("SAME");

	// // display2DMapContents(yield1SUncCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	// // yieldUncCanvas->Modified();
	// // yieldUncCanvas->Update();

	// // /// total uncertatinty

	// // TCanvas* totalUncCanvas = draw2DMap(totalRelUncCosThetaPhi, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kTRUE);

	// // totalRelUncCosThetaPhi->GetZaxis()->SetTitle("total uncertainty");
	// // totalRelUncCosThetaPhi->GetZaxis()->SetRangeUser(0, 1);

	// // kinematicsText->Draw("SAME");

	// // display2DMapContents(totalRelUncCosThetaPhi, nCosThetaBins, nPhiBins, kFALSE);

	// // totalUncCanvas->Modified();
	// // totalUncCanvas->Update();

	/// save canvas
	gSystem->mkdir(Form("EfficiencyMaps/%dS", iState), kTRUE);
	weightCanvas->SaveAs(Form("EfficiencyMaps/%dS/WeightsMapCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");

	gSystem->mkdir("YieldMap/2D_sysUnc", kTRUE);
	if (!LEGOplot) yieldCanvas->SaveAs(Form("YieldMap/2D_sysUnc/YieldMapCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s_%s_err%.2f_%sEff.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName, pseudoExpModelName, errPercent, efficiencyName), "RECREATE");
    else yieldCanvas->SaveAs(Form("YieldMap/2D_sysUnc/YieldMapCosTheta%s_fit_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s_%s_err%.2f_%sEff.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName, pseudoExpModelName, errPercent, efficiencyName), "RECREATE");
 
	if (!LEGOplot) correctedMapCanvas->SaveAs(Form("YieldMap/2D_sysUnc/CorrectedMapCosTheta%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s_%s_err%.2f_%sEff.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName, pseudoExpModelName, errPercent, efficiencyName), "RECREATE");
    else correctedMapCanvas->SaveAs(Form("YieldMap/2D_sysUnc/CorrectedMapCosTheta%s_fit_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s_%s_err%.2f_%sEff.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName, pseudoExpModelName, errPercent, efficiencyName), "RECREATE");

	// gSystem->mkdir("UncertaintyPlots/2D_sysUnc", kTRUE);
	// // statHighAccCanvas->SaveAs(Form("UncertaintyPlots/2D/statHighAcc%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// // statLowAccCanvas->SaveAs(Form("UncertaintyPlots/2D/statLowAcc%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// // statHighEffCanvas->SaveAs(Form("UncertaintyPlots/2D/statHighEff%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// // statLowEffCanvas->SaveAs(Form("UncertaintyPlots/2D/statLowEff%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// systEffCanvas->SaveAs(Form("UncertaintyPlots/2D/sysEff%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// // yieldUncCanvas->SaveAs(Form("UncertaintyPlots/2D/yieldUnc%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// // totalUncCanvas->SaveAs(Form("UncertaintyPlots/2D/totalUnc%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
    // // systEffFitModelCanvas->SaveAs(Form("UncertaintyPlots/2D_sysUnc/systUncFitModels%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.png", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");

	// // // save the histograms to the root files to see them with the root viewer
	// // TFile* EfficiencyOutFile = new TFile(Form("EfficiencyMaps/%dS/efficiencyHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.root", iState, refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// // accMapCosThetaPhi->Write();
	// // effMapCosThetaPhi->Write();
	// // systEffCosThetaPhi->Write();
	// // weightMap->Write();
	// // EfficiencyOutFile->Close();

	// TFile* ResultOutFile = new TFile(Form("YieldMap/2D_sysUnc/resultsHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.root", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// yieldMap->Write();
	// standardCorrectedMap->Write();
	// ResultOutFile->Close();

	// // TFile* uncOutFile = new TFile(Form("UncertaintyPlots/2D/uncertaintyHistos%s_cent%dto%d_pt%dto%dGeV_phi%dto%d_costheta%.1fto%.1f%s.root", refFrameName, gCentralityBinMin, gCentralityBinMax, ptMin, ptMax, (Int_t)phiBinEdges[0], (Int_t)phiBinEdges[nPhiBins], cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], gMuonAccName), "RECREATE");
	// // statLowAccCosThetaPhi->Write();
	// // statLowEffCosThetaPhi->Write();
	// // statHighAccCosThetaPhi->Write();
	// // statHighEffCosThetaPhi->Write();
	// // relSystEffCosThetaPhi->Write();
	// // yield1SUncCosThetaPhi->Write();
	// // totalRelUncCosThetaPhi->Write();
	// // uncOutFile->Close();
}

void scanRawYield_2D_customizedFits_sysUnc(const char* refFrameName = "CS", double errPercent = 1) {
	
	// /// loop over all pseudo-experiments
 	// Int_t chooseEff = 0; // fix the efficiency to the nominal value (no shift)
    
	// for (Int_t ptIdx = 0; ptIdx < NPtBins; ptIdx++) {
	// // for (Int_t ptIdx = 1; ptIdx < 2; ptIdx++) {
    //     for (Int_t chooseYield = 0; chooseYield <= 0; chooseYield++) {
    //         for (Int_t idx = 0; idx < 2; idx++) {
    //             if (idx == 0) rawYield_2D_customizedFits_sysUnc(gPtBinning[ptIdx], gPtBinning[ptIdx + 1], refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kFALSE, "ExpTimesErr", chooseYield, errPercent);
    //             else rawYield_2D_customizedFits_sysUnc(gPtBinning[ptIdx], gPtBinning[ptIdx + 1], refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kTRUE, "ExpTimesErr", chooseYield, errPercent);
    //             // if (idx == 0) rawYield_2D_customizedFits_sysUnc(0, 2, refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kFALSE, "ExpTimesErr", chooseYield);
    //             // else rawYield_2D_customizedFits_sysUnc(0, 2, refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kTRUE, "ExpTimesErr", chooseYield);
    //         }
    //     }
    // }
	
	/// loop over different efficiency scenarios
	Int_t chooseYield = 0; // fix the yield to the nominal value (no shift)

	for (Int_t ptIdx = 0; ptIdx < NPtBins; ptIdx++) {
		for (Int_t chooseEff = 1; chooseEff <= 16; chooseEff++) {
            for (Int_t idx = 0; idx < 2; idx++) {
                if (idx == 0) rawYield_2D_customizedFits_sysUnc(gPtBinning[ptIdx], gPtBinning[ptIdx + 1], refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kFALSE, "ExpTimesErr", chooseYield, errPercent, chooseEff);
                else rawYield_2D_customizedFits_sysUnc(gPtBinning[ptIdx], gPtBinning[ptIdx + 1], refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kTRUE, "ExpTimesErr", chooseYield, errPercent, chooseEff);
                // if (idx == 0) rawYield_2D_customizedFits_sysUnc(0, 2, refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kFALSE, "ExpTimesErr", chooseYield);
                // else rawYield_2D_customizedFits_sysUnc(0, 2, refFrameName, 5, -0.7, 0.7, 3, 0, 180, 1, kTRUE, "ExpTimesErr", chooseYield);
            }
        }
    }
}

