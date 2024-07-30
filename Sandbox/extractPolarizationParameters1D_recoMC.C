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

TPaveStats* positionStatBox(TH1* hist) {

	// Access the stat box
    auto *statBox = dynamic_cast<TPaveStats*>(hist->FindObject("stats"));

    if (statBox) {
        statBox->SetX1NDC(0.26); // New X1 coordinate in NDC
        statBox->SetX2NDC(0.87); // New X2 coordinate in NDC
        statBox->SetY1NDC(0.09); // New Y1 coordinate in NDC
        statBox->SetY2NDC(0.32); // New Y2 coordinate in NDC
        statBox->SetLineColor(0);
    }
    else cout << "not found stat box!!" << endl;

    return statBox;
}

TPad* getPadPullDistribution(TH1* dataHist, TF1* fitFunction, TFitResultPtr fitResult) {

    Int_t nBins = dataHist->GetNbinsX();
    
    TH1D* pullHist = new TH1D("pullHist", "Pull Distribution", nBins, dataHist->GetXaxis()->GetXmin(), dataHist->GetXaxis()->GetXmax());

    for (int iBin=1; iBin<=nBins; iBin++) {
    	double xValue = dataHist->GetXaxis()->GetBinCenter(iBin);
        double observed = dataHist->GetBinContent(iBin);
        double fitted = fitFunction->Eval(xValue);
        double error = dataHist->GetBinError(iBin);

        if (error != 0) {
            double pull = (observed - fitted) / error;
            pullHist->SetBinContent(iBin, pull);
        } else {
            pullHist->SetBinContent(iBin, 0);
        }
    }

    TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0, 0.0, 1, .25);
    bottomPad->SetTopMargin(0.015);
	bottomPad->SetBottomMargin(0.4);
	bottomPad->SetRightMargin(0.01);
	bottomPad->SetLeftMargin(0.13);
	bottomPad->SetTicks(1, 1);
	bottomPad->Draw();
	bottomPad->cd();

	pullHist->SetTitle(" ");
	
	pullHist->GetYaxis()->SetTitleOffset(0.3);
	pullHist->GetYaxis()->SetTitle("Pull");
	
	pullHist->GetYaxis()->SetTitleSize(0.17);
	pullHist->GetYaxis()->SetLabelSize(0.15);
	pullHist->GetYaxis()->CenterTitle();

	pullHist->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle());
	pullHist->GetXaxis()->SetTitleSize(0.17);	
	pullHist->GetXaxis()->SetLabelSize(0.15);
	pullHist->GetXaxis()->CenterTitle();
	
	pullHist->GetYaxis()->SetTickSize(0.03);
	pullHist->GetXaxis()->SetTickSize(0.1);

	pullHist->SetMaximum(9.5);
	pullHist->SetMinimum(-9.5);
    
	pullHist->SetMarkerStyle(8);
	pullHist->SetMarkerSize(1);
	pullHist->SetMarkerColor(kBlack);

    pullHist->Draw("PE");

    TLine zeroLine(bottomPad->GetUxmin(), 0, bottomPad->GetUxmax(), 0);
	zeroLine.SetLineStyle(kDashed);
	zeroLine.Draw("SAME");

	bottomPad->Update();

	TLatex textChi2;
	textChi2.SetTextAlign(12);
	textChi2.SetTextSize(0.15);
	textChi2.DrawLatexNDC(0.75, 0.12, Form("#chi^{2} / n_{dof} = %.1f", fitResult->Chi2() / fitResult->Ndf()));

    return bottomPad;
}

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

/// apply weights and errors to each costheta bin
Float_t correctMC1DHist(TH1D* projectionHist, TH1D* correctedHist, TEfficiency* accMap, TEfficiency* effMap, TH1D* systEff, Int_t nBins = 10) {
	// initialize variables used in the loop below
	Double_t totalRelUncHigh = 0, totalRelUncLow = 0;
	Double_t totalUncHigh = 0, totalUncLow = 0;

	TCanvas* massCanvas = 0;

	Float_t maxMCVal = 0;

	for (Int_t iBin = 0; iBin < nBins; iBin++) {
		Double_t weight = 0;

		// get the corresponding weights
		double acceptance = accMap->GetEfficiency(iBin + 1);
		double efficiency = effMap->GetEfficiency(iBin + 1);

		// calculate weight
		weight = 1. / (acceptance * efficiency);
		// weight = 1. ;

		// propagate both scale factor uncertainties and efficiency stat errors to the weight
		double relSystUnc = systEff->GetBinContent(iBin + 1);

		double relEffUncHigh = effMap->GetEfficiencyErrorUp(iBin + 1) / efficiency;
		double relEffUncLow = effMap->GetEfficiencyErrorLow(iBin + 1) / efficiency;

		double relAccUncHigh = accMap->GetEfficiencyErrorUp(iBin + 1) / acceptance;
		double relAccUncLow = accMap->GetEfficiencyErrorLow(iBin + 1) / acceptance;

		totalRelUncHigh = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncHigh), relAccUncHigh);
		totalRelUncLow = TMath::Hypot(TMath::Hypot(relSystUnc, relEffUncLow), relAccUncLow);

		totalUncHigh = totalRelUncHigh * efficiency * acceptance;
		totalUncLow = totalRelUncLow * efficiency * acceptance;

		// get yields and their uncertainties
		double recoMCVal = projectionHist->GetBinContent(iBin + 1);

		cout << "iBin " << iBin << endl;
		cout << "recoMCval: " << recoMCVal << endl;	
		cout << "acceptance: " << acceptance << endl;
		cout << "efficiency: " << efficiency << endl;
		cout << "weight: " << weight << endl;

		double recoMCUnc = projectionHist->GetBinError(iBin + 1);

		// set the bin contents reflecting weights
		correctedHist->SetBinContent(iBin + 1, recoMCVal * weight);

		// standardCorrectedHist->SetBinError(iBin + 1, yield1SErr * weight);
		// standardCorrectedHist->SetBinError(iBin + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relAccUncHigh));
		// standardCorrectedHist->SetBinError(iBin + 1, yield1SVal * weight * TMath::Hypot(yield1SErr / yield1SVal, relEffUncHigh));

		// standardCorrectedHist->SetBinError(iBin + 1, yield1SUnc * weight);

		// standardCorrectedHist->SetBinError(iBin + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

		// correctedHist->SetBinError(iBin + 1, TMath::Hypot(recoMCUnc / recoMCVal, totalRelUncHigh) * recoMCVal * weight);

		// // fill uncertainty histograms
		// statHighEff->SetBinContent(iBin + 1, relEffUncHigh);
		// statLowEff->SetBinContent(iBin + 1, relEffUncLow);
		// statHighAcc->SetBinContent(iBin + 1, relAccUncHigh);
		// statLowAcc->SetBinContent(iBin + 1, relAccUncLow);
		// yield1SUnc->SetBinContent(iBin + 1, yield1SUnc / yield1SVal);
		// totalRelUnc->SetBinContent(iBin + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh));

		if (recoMCVal * weight > maxMCVal) maxMCVal = recoMCVal * weight;

	}

	return maxMCVal;
}

void correctMC2DHist(TH2D* projectionHist, TH2D* correctedHist, TEfficiency* accMap, TEfficiency* effMap, TH2D* systEff, Int_t nCosThetaBins, Int_t nPhiBins) {

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
			if (acceptance == 0 || efficiency == 0) weight = 0;
			else weight = 1. / (acceptance * efficiency);

			// else weight = 1;

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

			double recoMCVal = projectionHist->GetBinContent(iCosTheta + 1, iPhi + 1);

			double recoMCUnc = projectionHist->GetBinError(iCosTheta + 1, iPhi + 1);

			// set the bin contents reflecting weights
			// yield with acceptance x efficiency correction
			correctedHist->SetBinContent(iCosTheta + 1, iPhi + 1, recoMCVal * weight);
					
			// standardCorrectedMap->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(yield1SUnc / yield1SVal, totalRelUncHigh) * yield1SVal * weight);
		
			// correctedHist->SetBinError(iCosTheta + 1, iPhi + 1, TMath::Hypot(recoMCUnc / recoMCVal, totalRelUncHigh) * recoMCVal * weight);

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

	/// dummy histogram to adjust the plot range
	TH2F* hdummy = new TH2F("hdummy", ";cos #theta; #varphi (#circ);Number of generated #varUpsilons", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins + 5, phiBinEdges[0] - 20, phiBinEdges[nPhiBins] + 100); 

	for (Int_t iEvent = 0; iEvent < allDataset->numEntries(); iEvent++) {
		const RooArgSet* iRooArgSet = allDataset->get(iEvent);
		
		Double_t cosThetaVal = ((RooRealVar*)iRooArgSet->find(CosThetaVarName(refFrameName)))->getVal();
		Double_t phiVal = ((RooRealVar*)iRooArgSet->find(PhiVarName(refFrameName)))->getVal();
		Double_t phiTildeVal = ((RooRealVar*)iRooArgSet->find(PhiTildeVarName(refFrameName)))->getVal();
		Double_t weight = allDataset->weight();

		generalPolarHist->Fill(cosThetaVal, phiVal, weight);
		generalPolarTildeHist->Fill(cosThetaVal, phiTildeVal, weight);
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

	// hdummy->Draw("LEGO");

	correctedPolarHist->Draw("LEGO");

	generalPolarFuncFit->Draw("SURFACE SAME");

	/// cosmetics	
	mc2DCanvas->SetTopMargin(.15);
	mc2DCanvas->SetLeftMargin(.1);

	hdummy->GetZaxis()->SetRangeUser(0, correctedPolarHist->GetMaximum());

	correctedPolarHist->GetZaxis()->SetMaxDigits(3);

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

void extractPolarizationParameters1D_recoMC(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", 
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

	/// Make 2D histograms to 1D 
	// (integrate over phi, that is, costheta graph)
	TH1D* polarHistCosTheta = generalPolarHist->ProjectionX("cos #theta", 0, nPhiBins); // arguments: (name, firstybin, lastybin)

	// (integrate over cosTheta, that is, phi graph)
	TH1D* polarHistPhi = generalPolarHist->ProjectionY("#varphi", 0, nCosThetaBins);

	TH1D* polarHistPhiTilde = generalPolarTildeHist->ProjectionY("#tilde{#varphi}", 0, nCosThetaBins);

	/// apply weights and errors to each costheta bin

	TH1D* correctedHistCosTheta = new TH1D("correctedHistCosTheta", " ", nCosThetaBins, cosThetaBinEdges.data());

	TH1D* correctedHistPhi = new TH1D("correctedHistPhi", " ", nPhiBins, phiBinEdges.data());

	Float_t maxYieldCosTheta = correctMC1DHist(polarHistCosTheta, correctedHistCosTheta, accHistCosTheta, effHistCosTheta, systEffCosTheta, nCosThetaBins);
	
	Float_t maxYieldPhi = correctMC1DHist(polarHistPhi, correctedHistPhi, accHistPhi, effHistPhi, systEffPhi, nPhiBins);

	/// using ROOT::Fit()

	// perform fit of cosTheta graph 
	Float_t maxYield = polarHistCosTheta->GetEntries();

	TF1* cosThetaPolarFuncFit = getCosThetaPolarFunc(maxYield);

	TFitResultPtr cosThetafitResults = correctedHistCosTheta->Fit("cosThetaPolarFunc", "ESVIMR0"); // chi2 fit to the integrated bin 

	cout << "--------------------------------------" << endl;
	cout << "Done costheta fit!" << endl;
	cout << "normalization: " << cosThetafitResults->Parameter(0) << ", lambdaTheta: " << cosThetafitResults->Parameter(1) << endl;
	cout << "--------------------------------------" << endl;

	cosThetafitResults->Print("v");

	TF1* phiPolarFuncFit = getPhiPolarFunc(maxYield);

	// fix the extracted lambdaTheta value
	// phiPolarFuncFit->SetParameter(1, cosThetafitResults->Parameter(1));
	phiPolarFuncFit->FixParameter(1, cosThetafitResults->Parameter(1));

	// perform fit of phi graph 
	TFitResultPtr phifitResults = correctedHistPhi->Fit("phiPolarFunc", "ESVIMR0"); // chi2 fit to the integrated bin 

	cout << "--------------------------------------" << endl;
	cout << "Done phi fit!" << endl;
	cout << "normalization: " << phifitResults->Parameter(0) << ", lambdaPhi: " << phifitResults->Parameter(2) << endl;
	cout << "--------------------------------------" << endl;
	
	phifitResults->Print("v");

	TF1* phiTildePolarFuncFit = getPhiTildePolarFunc(maxYield);

	// fix the extracted lambdaTheta value
	phiTildePolarFuncFit->FixParameter(1, cosThetafitResults->Parameter(1));

	// perform fit of phi graph 
	TFitResultPtr phiTildefitResults = polarHistPhiTilde->Fit("phiTildePolarFunc", "ESVIMR0"); // chi2 fit to the integrated bin 

	cout << "--------------------------------------" << endl;
	cout << "Done phiTilde fit!" << endl;
	cout << "normalization: " << phiTildefitResults->Parameter(0) << ", lambdaPhiTilde: " << phiTildefitResults->Parameter(2) << endl;
	cout << "--------------------------------------" << endl;
	
	phiTildefitResults->Print("v");

	// draw the histogram and the fit
	TCanvas* polarCanvas1D = new TCanvas("polarCanvas1D", "", 1350, 500);

	polarCanvas1D->Divide(3, 1);

	polarCanvas1D->cd(1);

	TPad* padCosTheta = new TPad("padCosTheta", "padCosTheta", 0, 0.25, 1, 1.0);
	padCosTheta->SetBottomMargin(0.03);
	padCosTheta->SetRightMargin(0.01);
	padCosTheta->SetLeftMargin(0.13);
	padCosTheta->SetTicks(1, 1);
	padCosTheta->Draw();
	padCosTheta->cd();

	correctedHistCosTheta->SetMinimum(0);

	correctedHistCosTheta->SetMarkerStyle(8);
	correctedHistCosTheta->SetMarkerSize(1.5);
	correctedHistCosTheta->SetMarkerColor(kBlack);
	
	correctedHistCosTheta->GetXaxis()->SetLabelSize(0);

	correctedHistCosTheta->GetYaxis()->SetTitle(Form("Events / (%.1f)", (cosThetaMax - cosThetaMin) / nCosThetaBins));
	correctedHistCosTheta->GetYaxis()->SetTitleOffset(1.1);

	correctedHistCosTheta->Draw("PE");

	cosThetaPolarFuncFit->Draw("SAME");

  	padCosTheta->Draw();

  	TPaveStats* cosThetaStatbox = positionStatBox(correctedHistCosTheta);
  	cosThetaStatbox->Draw();

	polarCanvas1D->cd(1);

	// pull Distribution
	TPad* padCosThetaPull = getPadPullDistribution(correctedHistCosTheta, cosThetaPolarFuncFit, cosThetafitResults);

	// phi graph
	polarCanvas1D->cd(2);

	TPad* padPhi = new TPad("padPhi", "padPhi", 0, 0.25, 1, 1.0);
	padPhi->SetBottomMargin(0.03);
	padPhi->SetRightMargin(0.01);
	padPhi->SetLeftMargin(0.13);
	padPhi->SetTicks(1, 1);
	padPhi->Draw();
	padPhi->cd();

	correctedHistPhi->SetMinimum(0);

	correctedHistPhi->SetMarkerStyle(8);
	correctedHistPhi->SetMarkerSize(1.5);
	correctedHistPhi->SetMarkerColor(kBlack);
	
	correctedHistPhi->GetXaxis()->SetLabelSize(0);

	correctedHistPhi->GetYaxis()->SetTitle(Form("Events / (%.1f)", (Double_t)((phiMax - phiMin) / nPhiBins)));
	correctedHistPhi->GetYaxis()->SetTitleOffset(1.1);

	correctedHistPhi->Draw("PE");

	phiPolarFuncFit->Draw("SAME");
  	
  	padPhi->Draw();

  	TPaveStats* phiStatbox = positionStatBox(correctedHistPhi);
  	phiStatbox->Draw();

	polarCanvas1D->cd(2);

	// pull Distribution
	TPad* padPhiPull = getPadPullDistribution(correctedHistPhi, phiPolarFuncFit, phifitResults);

	// phi tilde graph
	polarCanvas1D->cd(3);

	TPad* padPhiTilde = new TPad("padPhiTilde", "padPhiTilde", 0, 0.25, 1, 1.0);
	padPhiTilde->SetBottomMargin(0.03);
	padPhiTilde->SetRightMargin(0.01);
	padPhiTilde->SetLeftMargin(0.13);
	padPhiTilde->SetTicks(1, 1);
	padPhiTilde->Draw();
	padPhiTilde->cd();

	polarHistPhiTilde->SetMinimum(0);

	polarHistPhiTilde->SetMarkerStyle(8);
	polarHistPhiTilde->SetMarkerSize(1.5);
	polarHistPhiTilde->SetMarkerColor(kBlack);
	
	polarHistPhiTilde->GetXaxis()->SetLabelSize(0);

	polarHistPhiTilde->GetYaxis()->SetTitle(Form("Events / (%.1f)", (Double_t)((phiMax - phiMin) / nPhiBins)));
	polarHistPhiTilde->GetYaxis()->SetTitleOffset(1.1);

	polarHistPhiTilde->Draw("PE");

	phiTildePolarFuncFit->Draw("SAME");
  	
  	padPhiTilde->Draw();

  	TPaveStats* phiTildeStatbox = positionStatBox(polarHistPhiTilde);
  	phiTildeStatbox->Draw();

	polarCanvas1D->cd(3);

	// pull Distribution
	TPad* padPhiTildePull = getPadPullDistribution(polarHistPhiTilde, phiTildePolarFuncFit, phiTildefitResults);

	gSystem->mkdir("DistributionFitsMC", kTRUE);
  	polarCanvas1D->SaveAs(Form("DistributionFitsMC/RecoMC_fit1D_lambdaCosTheta%.2fPhi%.2fThetaPhi%.2f.png", lambdaTheta0, lambdaPhi0, lambdaThetaPhi0), "RECREATE");

	return;
}