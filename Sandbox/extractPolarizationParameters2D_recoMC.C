// This code extracts the polarization parameters 
// from the 2D angular distribution function (ideal MC distribution for the test) in the costheta phi space 
// using 1D fit

#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../MonteCarlo/AccEffHelpers.h"

#include "../Polarization/PolarFitHelpers.h"

#include "../Tools/FitShortcuts.h"
#include "../Tools/Style/Legends.h"

#include "../Tools/Style/FitDistributions.h"

#include "../Tools/RooFitPDFs/PolarFunc.h"

#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.h"
#include "../Tools/RooFitPDFs/CosThetaPolarizationPDF.cxx"

#include "../Tools/RooFitPDFs/PhiPolarizationPDF.h"
#include "../Tools/RooFitPDFs/PhiPolarizationPDF.cxx"

#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.h"
#include "../Tools/RooFitPDFs/GeneralPolarizationPDF.cxx"

using namespace RooFit;

TEfficiency* getAcceptance3DMap(const char* refFrameName, Double_t lambdaTheta = 0.88, Double_t lambdaPhi = -0.8, Double_t lambdaThetaPhi = 0) {

	/// get acceptance and efficiency in 1D
	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get acceptance maps
	TFile* acceptanceFile = openFile("../MonteCarlo/AcceptanceMaps/1S/AcceptanceResults.root");
	auto* accMap = (TEfficiency*)acceptanceFile->Get(Form("%s_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", nominalMapName, lambdaTheta, lambdaPhi, lambdaThetaPhi));

	return accMap;
}

TEfficiency* getEfficiency3DMap(const char* refFrameName, Double_t lambdaTheta = 0.88, Double_t lambdaPhi = -0.8, Double_t lambdaThetaPhi = 0) {

	/// get acceptance and efficiency in 1D
	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get efficiency maps
	TFile* efficiencyFile = openFile("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root");
	auto* effMap = (TEfficiency*)efficiencyFile->Get(Form("%s_LambdaTheta%.2fPhi%.2fThetaPhi%.2f", nominalMapName, lambdaTheta, lambdaPhi, lambdaThetaPhi));

	return effMap;
}

TH3D* getSysEff3DMap(const char* refFrameName) {

	/// get acceptance and efficiency in 1D
	const char* nominalMapName = NominalTEfficiency3DName(refFrameName);

	// get relative systematic uncertainty of efficiency
	TFile* efficiencyFile = openFile("../MonteCarlo/EfficiencyMaps/1S/EfficiencyResults.root");
	auto* systEff = (TH3D*)efficiencyFile->Get(RelativeSystTEfficiency3DName(refFrameName));

	return systEff;
}

void correctMC2DHist(TH2D* polarizationHist, TH2D* correctedHist, TEfficiency* accMap, TEfficiency* effMap, TH2D* systEff, Int_t nCosThetaBins, Int_t nPhiBins) {

	TH2D* hTotalCosThetaPhi = (TH2D*)accMap->GetTotalHistogram();

	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	Float_t maxYield = 0;

	TCanvas* massCanvas = 0;

	/// apply weights and errors to each costheta bin
	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			Double_t weight = 1;

			// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotalCosThetaPhi->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotalCosThetaPhi->GetYaxis()->GetBinCenter(iPhi + 1);

			Int_t iGlobalBin = hTotalCosThetaPhi->FindFixBin(binCenterCosTheta, binCenterPhi);

			// get the corresponding weights
			double acceptance = accMap->GetEfficiency(iGlobalBin);
			double efficiency = effMap->GetEfficiency(iGlobalBin);

			// calculate weight
			if (acceptance == 0 || efficiency == 0) weight = 1.;
			else weight = 1. / (acceptance * efficiency);

			// weightMap->SetBinContent(iCosTheta + 1, iPhi + 1, weight);

			// propagate both scale factor uncertainties and efficiency stat errors to the weight
			double relSystUnc = systEff->GetBinContent(iCosTheta + 1, iPhi + 1) / efficiency;

			double relEffUncHigh = effMap->GetEfficiencyErrorUp(iGlobalBin) / efficiency;
			double relEffUncLow = effMap->GetEfficiencyErrorLow(iGlobalBin) / efficiency;

			double relAccUncHigh = accMap->GetEfficiencyErrorUp(iGlobalBin) / acceptance;
			double relAccUncLow = accMap->GetEfficiencyErrorLow(iGlobalBin) / acceptance;

			// totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
			// totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

			totalRelUncHigh = TMath::Hypot(relEffUncHigh, relAccUncHigh);
			totalRelUncLow = TMath::Hypot(relEffUncLow, relAccUncLow);

			totalUncHigh = totalRelUncHigh * efficiency * acceptance;
			totalUncLow = totalRelUncLow * efficiency * acceptance;

			double recoMCVal = polarizationHist->GetBinContent(iCosTheta + 1, iPhi + 1);

			double recoMCUnc = polarizationHist->GetBinError(iCosTheta + 1, iPhi + 1);

			cout << "recoMCVal: " << recoMCVal << endl;
			cout << "recoMCUnc: " << recoMCUnc << endl;

			// set the bin contents reflecting weights
			// yield with acceptance x efficiency correction
			correctedHist->SetBinContent(iCosTheta + 1, iPhi + 1, recoMCVal * weight);
					
			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);
		
			correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(recoMCUnc / recoMCVal, totalRelUncHigh) * recoMCVal * weight);
			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(recoMCUnc / recoMCVal, totalRelUncHigh));

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

}

void getPolarizedMCHist(TH2D* generalPolarHist, TH2D* generalPolarTildeHist, const char* refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 20, const vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 18, const vector<Double_t>& phiBinEdges = {}, Int_t iState = 1, Double_t lambdaTheta = 0.88, Double_t lambdaPhi = -0.8, Double_t lambdaThetaPhi = 0){

	// reconstruced MC
	const char* mcFileName = Form("../Files/Y1SReconstructedMCWeightedDataset_Lambda_Theta%.2f_Phi%.2f_ThetaPhi%.2f.root", lambdaTheta, lambdaPhi, lambdaThetaPhi);

	TFile* mcFile = openFile(mcFileName);

	RooDataSet* allDataset = (RooDataSet*)mcFile->Get(Form("MCdataset%s", refFrameName));

	// import the dataset to a workspace
	RooWorkspace wspace("workspace");
	wspace.import(*allDataset);

	RooRealVar cosTheta = *wspace.var(CosThetaVarName(refFrameName));
	RooRealVar phi = *wspace.var(PhiVarName(refFrameName));
	RooRealVar phiTilde = *wspace.var(PhiTildeVarName(refFrameName));
	
	allDataset->Print("V");
	wspace.Print("V");

	/// dummy histogram to adjust the plot range
	TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins + 5, phiBinEdges[0] - 20, phiBinEdges[nPhiBins] + 100); 

	TH2D* errorHist = new TH2D("hdummy", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins]); 

	Double_t tempErrorUp[100][100] = {0};
	Double_t tempErrorDown[100][100] = {0};

	for (Int_t iEvent = 0; iEvent < allDataset->numEntries(); iEvent++) {

		// if (iEvent > 100) break;

		const RooArgSet* iRooArgSet = allDataset->get(iEvent);
		
		Double_t cosThetaVal = ((RooRealVar*) iRooArgSet->find(CosThetaVarName(refFrameName)))->getVal();
		Double_t phiVal = ((RooRealVar*) iRooArgSet->find(PhiVarName(refFrameName)))->getVal();
		Double_t phiTildeVal = ((RooRealVar*) iRooArgSet->find(PhiTildeVarName(refFrameName)))->getVal();
		
		Double_t weight = allDataset->weight();
		Double_t weightErrorUp = ((RooRealVar*) iRooArgSet->find("errorWeightUp"))->getVal();
		Double_t weightErrorDown = ((RooRealVar*) iRooArgSet->find("errorWeightDown"))->getVal();

		Int_t iCosTheta = generalPolarHist->GetXaxis()->FindBin(cosThetaVal);
		Int_t iPhi = generalPolarHist->GetYaxis()->FindBin(phiVal);

		tempErrorUp[iCosTheta - 1][iPhi - 1] = TMath::Hypot(tempErrorUp[iCosTheta - 1][iPhi - 1], weightErrorUp);
		tempErrorDown[iCosTheta - 1][iPhi - 1] = TMath::Hypot(tempErrorDown[iCosTheta - 1][iPhi - 1], weightErrorDown);

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

		generalPolarHist->Fill(cosThetaVal, phiVal, weight);
		generalPolarTildeHist->Fill(cosThetaVal, phiTildeVal, weight);
	}	

	for (int iCosTheta = 1; iCosTheta <= nCosThetaBins; iCosTheta++) {
		for (int iPhi = 1; iPhi <= nPhiBins; iPhi++) { 
			errorHist->SetBinContent(iCosTheta, iPhi, tempErrorUp[iCosTheta - 1][iPhi - 1]);

			generalPolarHist->SetBinError(iCosTheta, iPhi, tempErrorUp[iCosTheta - 1][iPhi - 1]);
			// cout << "BinError: " << generalPolarHist->GetBinError(iCosTheta, iPhi) << endl;
			// cout << "BinValue: " << generalPolarHist->GetBinContent(iCosTheta, iPhi) << endl;
		}
	}
	// get acceptance and efficiency 3D map
	auto* accMap = getAcceptance3DMap(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi); 
	auto* effMap = getEfficiency3DMap(refFrameName, lambdaTheta, lambdaPhi, lambdaThetaPhi); 
	auto* systEff = getSysEff3DMap(refFrameName);	

	// rebin acceptance and efficiency, efficiency systematic Uncertainty
	TEfficiency* accMapCosThetaPhi = rebinTEff3DMap(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TEfficiency* effMapCosThetaPhi = rebinTEff3DMap(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);
	TH2D* systEffCosThetaPhi = rebinRel3DUncMap(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges);

	// draw acceptance and efficiency graph for check
	DrawEfficiency2DHist(accMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kTRUE, kFALSE);
	DrawEfficiency2DHist(effMapCosThetaPhi, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, kFALSE, kFALSE);

	//draw error histogram
	draw2DMap(errorHist, refFrameName, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, kFALSE, kFALSE, 1);
	// display2DMapContents(errorHist, nCosThetaBins, nPhiBins, kFALSE);

	// apply acceptance and efficiency correction
	TH2D* correctedPolarHist = new TH2D("correctedPolarHist", "correctedPolarHist", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	correctMC2DHist(generalPolarHist, correctedPolarHist, accMapCosThetaPhi, effMapCosThetaPhi, systEffCosThetaPhi, nCosThetaBins, nPhiBins);

	Float_t maxYield = correctedPolarHist->GetEntries();

	TF2* generalPolarFuncFit = getGeneralPolarFunc(maxYield);

	TFitResultPtr fitResults = correctedPolarHist->Fit("generalPolarFunc", "ESVIM0"); // chi2 fit to the integrated bin 

	// draw plots

	TCanvas* mc2DCanvas = new TCanvas("mc2DCanvas", "mc2DCanvas", 1250, 600);

	mc2DCanvas->Divide(2);

	mc2DCanvas->cd(1);

    gPad->SetRightMargin(0.18);

	hdummy->Draw("COLZ");

	correctedPolarHist->Draw("COLZ SAME");

	// Styles of the texts in the plot
	TLatex* legend1 = new TLatex();
	legend1->SetTextAlign(22);
	legend1->SetTextSize(0.05);

	// Put texts inside the plot
	legend1->DrawLatexNDC(.50, .88, "Reconstructed MC");
	legend1->DrawLatexNDC(.50, .80, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	mc2DCanvas->cd(2);

	gPad->SetTopMargin(0.05);

	correctedPolarHist->Draw("LEGO");

	generalPolarFuncFit->Draw("SURFACE SAME");

	/// cosmetics	
	mc2DCanvas->SetTopMargin(.15);
	mc2DCanvas->SetLeftMargin(.1);

	hdummy->GetZaxis()->SetRangeUser(0, correctedPolarHist->GetMaximum());

	correctedPolarHist->SetStats(0);

	correctedPolarHist->GetZaxis()->SetMaxDigits(3);

	generalPolarFuncFit->SetRange(-1, -180, 0, 1, 180, generalPolarFuncFit->GetMaximum() * 1.1);
	generalPolarFuncFit->SetMaximum(generalPolarFuncFit->GetMaximum());
	generalPolarFuncFit->SetMinimum(0);

	// Styles of the texts in the plot
	TLatex* legend2 = new TLatex();

	// legend2->SetTextAlign(22);
	legend2->SetTextSize(0.05);

	// Put texts inside the plot
	legend2->DrawLatexNDC(.5, .90, Form("#lambda_{#theta, fit} = %.4f #pm %.4f", fitResults->Parameter(1), fitResults->ParError(1)));
	legend2->DrawLatexNDC(.5, .84, Form("#lambda_{#varphi, fit} = %.4f #pm %.4f", fitResults->Parameter(2), fitResults->ParError(2)));
	legend2->DrawLatexNDC(.5, .78, Form("#lambda_{#theta#varphi, fit} = %.4f #pm %.4f", fitResults->Parameter(3), fitResults->ParError(3)));

	// Set the plot styles
	hdummy->GetZaxis()->SetTitleOffset(1.);
	hdummy->GetZaxis()->SetMaxDigits(3);

	hdummy->GetXaxis()->SetTitleOffset(1.);
	hdummy->GetXaxis()->CenterTitle();

	hdummy->GetYaxis()->SetTitleOffset(1.5);
	hdummy->GetYaxis()->CenterTitle();

	SetColorPalette(gPreferredColorPaletteName);

	gPad->Update();

	mc2DCanvas->Update();

	// save the plot
	gSystem->mkdir("DistributionFitsMC", kTRUE);
	mc2DCanvas->SaveAs(Form("DistributionFitsMC/RecoMC_Theta%.2f_Phi%.2f_ThetaPhi%.2f.png", lambdaTheta, lambdaPhi, lambdaThetaPhi), "RECREATE");
}

void extractPolarizationParameters2D_recoMC(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", 
											const Int_t nCosThetaBins = 10, Double_t cosThetaMin = -1, Double_t cosThetaMax = 1, 
											const Int_t nPhiBins = 6, Int_t phiMin = -180, Int_t phiMax = 180, 
											Int_t iState = gUpsilonState,
											Double_t lambdaTheta0 = 0.88, Double_t lambdaPhi0 = -0.8, Double_t lambdaThetaPhi0 = 0) {  

	/// Bin edges and width 
	// Set the bin edges along the cosTheta/phi axis depending on the number of bins 
	vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaBins, cosThetaMin, cosThetaMax);
	
	vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiBins, phiMin, phiMax);

	// bin width
	Double_t cosThetaStep = (cosThetaBinEdges[nCosThetaBins] - cosThetaBinEdges[0]) / nCosThetaBins;

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	// get acceptance 3D map
	auto* accMap = getAcceptance3DMap(refFrameName); 

	// rebin acceptance maps based on costheta, phi, and pT selection
	TEfficiency* accHistCosTheta = rebinTEff3DMapCosTheta(accMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);
	TEfficiency* accHistPhi = rebinTEff3DMapPhi(accMap, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);
	
	// get efficiency 3D map
	auto* effMap = getEfficiency3DMap(refFrameName); 	

	// rebin efficiency maps based on costheta, phi, and pT selection
	TEfficiency* effHistCosTheta = rebinTEff3DMapCosTheta(effMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);
	TEfficiency* effHistPhi = rebinTEff3DMapPhi(effMap, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);
	
	auto* systEff = getSysEff3DMap(refFrameName);

	// rebin uncertainty map based on costheta, phi, and pT selection
	TH1D* systEffCosTheta = rebinRel3DUncCosTheta(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiMin, phiMax);
	TH1D* systEffPhi = rebinRel3DUncPhi(effMap, systEff, ptMin, ptMax, cosThetaMin, cosThetaMax, nPhiBins, phiBinEdges);

	// draw acc and eff histograms to check if the rebinning works well
	DrawEfficiency1DHist(accHistCosTheta, ptMin, ptMax, iState, kTRUE, kTRUE);
	DrawEfficiency1DHist(effHistCosTheta, ptMin, ptMax, iState, kFALSE, kTRUE);

	DrawEfficiency1DHist(accHistPhi, ptMin, ptMax, iState, kTRUE, kFALSE);
	DrawEfficiency1DHist(effHistPhi, ptMin, ptMax, iState, kFALSE, kFALSE);

	/// Generate a Toy Data (This part can be replaced by data)

	// (the sample from reconstructed MC)
	TH2D* generalPolarHist = new TH2D("generalPolarHist", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 
	TH2D* generalPolarTildeHist = new TH2D("generalPolarTildeHist", ";cos #theta; #tilde{#varphi} (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaMin, cosThetaMax, nPhiBins, phiMin, phiMax); 

	getPolarizedMCHist(generalPolarHist, generalPolarTildeHist, refFrameName, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, nPhiBins, phiBinEdges, iState, lambdaTheta0, lambdaPhi0, lambdaThetaPhi0);

	Double_t nEntries = generalPolarHist->GetEntries();

	cout << "--------------------------------------" << endl;
	cout << "number of entries: " <<  nEntries << endl;
	cout << "--------------------------------------" << endl;

	return;
}