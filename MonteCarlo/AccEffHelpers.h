#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "../Tools/FitShortcuts.h"

#include "TEfficiency.h"

/// Definition of common helping functions for acceptance and efficiency distribution making codes

// Option to estimate the confidence intervals in TEfficiency
TEfficiency::EStatOption gTEffStatOption = TEfficiency::EStatOption::kFCP; /*default = Clopper-Pearson as recommended by PDG*/

// Define cosTheta bin edges for acceptance and efficiency 3D map grid
std::vector<Double_t> setCosThetaBinEdges(Int_t nCosThetaBins, Double_t cosThetaMin, Double_t cosThetaMax) {
	std::vector<Double_t> cosThetaBinEdges = {};

	// define the bin edges along the cosTheta axis depending on the number of bins
	for (Int_t iCosTheta = 0; iCosTheta <= nCosThetaBins; iCosTheta++) {
		Double_t cosThetaBinWidth = (cosThetaMax - cosThetaMin) / nCosThetaBins;

		cosThetaBinEdges.push_back(cosThetaMin + cosThetaBinWidth * iCosTheta);
	}

	return cosThetaBinEdges;
}

std::vector<Double_t> setPhiBinEdges(Int_t nPhiBins, Int_t phiMin, Int_t phiMax) {
	std::vector<Double_t> phiBinEdges = {};

	// define the bin edges along the phi axis depending on the number of bins
	for (Int_t iPhi = 0; iPhi <= nPhiBins; iPhi++) {
		Double_t phiBinWidth = (phiMax - phiMin) / nPhiBins;

		phiBinEdges.push_back(phiMin + phiBinWidth * iPhi);
	}

	return phiBinEdges;
}

/// Object names (because used here and there) same for acceptance AND efficiency (saved in different files anyway)

// common naming convention
const char* TEfficiencyEndName(const char* refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30, Float_t lambdaTheta = 0., Float_t lambdaPhi = 0., Float_t lambdaThetaPhi = 0.) {
	return Form("%s%s_pt%dto%d_%s", refFrameName, gMuonAccName, ptMin, ptMax, PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi));
}

const char* TEfficiencyMainTitle(int iState = gUpsilonState, const char* title = "total efficiency") {
	return Form("#varUpsilon(%dS) %s", iState, title);
}

const char* TEfficiencyAccMainTitle(int iState = gUpsilonState) {
	return Form("#varUpsilon(%dS) acceptance", iState);
}

/// TEfficiency constructors

// cos theta distribution for a given pT range
TEfficiency* CosThetaTEfficiency1D(const char* refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30, int iState = gUpsilonState, Bool_t isAcc = kFALSE, Float_t lambdaTheta = 0., Float_t lambdaPhi = 0., Float_t lambdaThetaPhi = 0.) {
	const char* name = Form("CosTheta%s", TEfficiencyEndName(refFrameName, ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi));

	char* title = nullptr;
	if (isAcc)
		title = Form(";%s;%s", CosThetaVarTitle(refFrameName), TEfficiencyAccMainTitle(iState));
	else
		title = Form(";%s;%s", CosThetaVarTitle(refFrameName), TEfficiencyMainTitle(iState));

	TEfficiency* tEff = nullptr;
	if (std::string(refFrameName) == "CS")
		tEff = new TEfficiency(name, title, NCosThetaBinsCS, CosThetaBinningCS);
	else
		tEff = new TEfficiency(name, title, NCosThetaBinsHX, CosThetaBinningHX);

	tEff->SetStatisticOption(gTEffStatOption);

	return tEff;
}

// (cos theta, phi) distribution for a given pT range
const char* SystTEfficiency2DName(const char* refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30, Float_t lambdaTheta = 0., Float_t lambdaPhi = 0., Float_t lambdaThetaPhi = 0.) {
	return Form("SystEff2D%s", TEfficiencyEndName(refFrameName, ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi));
}

const char* CosThetaPhiTEfficiency2DName(const char* refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30, Float_t lambdaTheta = 0., Float_t lambdaPhi = 0., Float_t lambdaThetaPhi = 0.) {
	return Form("CosThetaAbsPhi%s", TEfficiencyEndName(refFrameName, ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi));
}

TEfficiency* CosThetaPhiTEfficiency2D(const char* refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30, int iState = gUpsilonState, Float_t lambdaTheta = 0., Float_t lambdaPhi = 0., Float_t lambdaThetaPhi = 0.) {
	const char* name = CosThetaPhiTEfficiency2DName(refFrameName, ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	const char* title = Form(";%s;%s;%s", CosThetaVarTitle(refFrameName), PhiAxisTitle(refFrameName), TEfficiencyMainTitle(iState));

	TEfficiency* tEff = new TEfficiency(name, title, NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);

	tEff->SetStatisticOption(gTEffStatOption);

	return tEff;
}

TEfficiency* CosThetaPhiAcceptance2D(const char* refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30, Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0) {
	const char* name = CosThetaPhiTEfficiency2DName(refFrameName, ptMin, ptMax, lambdaTheta, lambdaPhi, lambdaThetaPhi);

	const char* title = Form(";%s;%s;acceptance", CosThetaVarTitle(refFrameName), PhiAxisTitle(refFrameName));

	TEfficiency* tEff = new TEfficiency(name, title, NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);

	tEff->SetStatisticOption(gTEffStatOption);

	return tEff;
}

// (cos theta, phi, pT) 3D maps

// naming function so that we can more easily get the TEfficiency from other codes
const char* NominalTEfficiency3DName(const char* refFrameName = "CS", Float_t lambdaTheta = 0., Float_t lambdaPhi = 0., Float_t lambdaThetaPhi = 0.) {
	return Form("Nominal3D%s%s_%s", refFrameName, gMuonAccName, PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi));
}

const char* SystTEfficiency3DName(const char* refFrameName = "CS", Float_t lambdaTheta = 0., Float_t lambdaPhi = 0., Float_t lambdaThetaPhi = 0.) {
	return Form("SystEff3D%s%s_%s", refFrameName, gMuonAccName, PolaWeightName(lambdaTheta, lambdaPhi, lambdaThetaPhi));
}

const char* TEfficiency3DAxisTitle(const char* refFrameName = "CS") {
	return Form("%s;%s;%s", CosThetaVarTitle(refFrameName), PhiAxisTitle(refFrameName), gPtAxisTitle);
}

const char* TEfficiency3DTitle(const char* refFrameName = "CS", int iState = gUpsilonState) {
	return Form("%s;%s", TEfficiencyMainTitle(iState), TEfficiency3DAxisTitle(refFrameName));
}

TEfficiency* TEfficiency3D(const char* name, const char* refFrameName = "CS", int iState = gUpsilonState, Bool_t isPhiFolded = kFALSE) {
	Int_t nCosThetaUltraFineBinning = 200; // cosTheta bin width: 0.01
	Int_t nPhiUltraFineBinning;            // phi bin width: 1 degree

	// Int_t nCosThetaUltraFineBinning = 4; // cosTheta bin width: 0.01
	// Int_t nPhiUltraFineBinning = 2;      // phi bin width: 1 degree
	float phiMin, phiMax;

	if (isPhiFolded) {
		phiMax = gPhiMax;
		phiMin = gPhiMin;
		nPhiUltraFineBinning = 180; // phi bin width: 2 degree
	} else {
		phiMax = gPhiMax;
		phiMin = -180;
		nPhiUltraFineBinning = 360; // phi bin width: 1 degree
	}

	std::vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(nCosThetaUltraFineBinning, gCosThetaMin, gCosThetaMax);
	std::vector<Double_t> phiBinEdges = setPhiBinEdges(nPhiUltraFineBinning, phiMin, phiMax);

	// std::vector<Double_t> cosThetaBinEdges = setCosThetaBinEdges(5, -0.7, 0.7);
	// std::vector<Double_t> phiBinEdges = setPhiBinEdges(3, 0, 180);

	TEfficiency* tEff = new TEfficiency(name, TEfficiency3DTitle(refFrameName, iState), nCosThetaUltraFineBinning, cosThetaBinEdges.data(), nPhiUltraFineBinning, phiBinEdges.data(), NPtFineBins, gPtFineBinning); // variable size binning for the stats
	// TEfficiency* tEff = new TEfficiency(name, TEfficiency3DTitle(refFrameName, iState), 5, cosThetaBinEdges.data(), 3, phiBinEdges.data(), NPtBins, gPtBinning); // variable size binning for the stats

	tEff->SetStatisticOption(gTEffStatOption);

	return tEff;
}

// rebin TEfficiency 3D maps to TEfficiency 1D cosTheta based on costheta, phi, and pT selection

TEfficiency* rebinTEff3DMapCosTheta(TEfficiency* TEff3DMap, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 10, const std::vector<Double_t>& cosThetaBinEdges = {}, Int_t phiMin = -180, Int_t phiMax = 180) {
	// extract the numerator and the denominator from the 3D TEfficiency Map
	TH3D* hPassed = (TH3D*)TEff3DMap->GetPassedHistogram();
	TH3D* hTotal = (TH3D*)TEff3DMap->GetTotalHistogram();

	// obtain the bin numbers of the boundaries on phi and pt
	Int_t iPhiMin = hPassed->GetYaxis()->FindBin(phiMin);
	Int_t iPhiMax = hPassed->GetYaxis()->FindBin(phiMax);

	Int_t iPtMin = hPassed->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = hPassed->GetZaxis()->FindBin(ptMax);

	// std::cout << iPhiMin << " " << iPhiMax << " " << iPtMin << " " << iPtMax << std::endl;
	// cout << "hPassed->GetXaxis()->GetNbins() = " << hPassed->GetXaxis()->GetNbins() << endl;
	// cout << "hPassed->GetYaxis()->GetNbins() = " << hPassed->GetYaxis()->GetNbins() << endl;
	// cout << "hPassed->GetZaxis()->GetNbins() = " << hPassed->GetZaxis()->GetNbins() << endl;

	// cout << "hPassed->GetBinContent(1, 1, 2) = " << hPassed->GetBinContent(1, 1, 2) << endl;
	// cout << "hPassed->GetBinContent(1, 2, 2) = " << hPassed->GetBinContent(1, 2, 2) << endl;
	// cout << "hPassed->GetBinContent(2, 1, 2) = " << hPassed->GetBinContent(2, 1, 2) << endl;
	// cout << "hPassed->GetBinContent(2, 2, 2) = " << hPassed->GetBinContent(2, 2, 2) << endl;

	// cout << "hTotal->GetBinContent(1, 1, 2) = " << hTotal->GetBinContent(1, 1, 2) << endl;
	// cout << "hTotal->GetBinContent(1, 2, 2) = " << hTotal->GetBinContent(1, 2, 2) << endl;
	// cout << "hTotal->GetBinContent(2, 1, 2) = " << hTotal->GetBinContent(2, 1, 2) << endl;
	// cout << "hTotal->GetBinContent(2, 2, 2) = " << hTotal->GetBinContent(2, 2, 2) << endl;

	// obtain the projection histogram along the costheta axis within boundaries of phi and pt
	// (option e: calculate errors, o: only bins inside the selected range will be filled)
	TH1D* hPassedCosTheta = (TH1D*)hPassed->ProjectionX("hPassedCosTheta", iPhiMin, iPhiMax - 1, iPtMin, iPtMax - 1, "eo");
	TH1D* hTotalCosTheta = (TH1D*)hTotal->ProjectionX("hTotalCosTheta", iPhiMin, iPhiMax - 1, iPtMin, iPtMax - 1, "eo");

	if (iPhiMin >= iPhiMax - 1) {
		hPassedCosTheta->Reset();
		hTotalCosTheta->Reset();
	}

	// rebin the projection histogram (default: 20 bins from -1 to 1)
	TH1D* hPassedCosTheta_Rebin = (TH1D*)hPassedCosTheta->Rebin(nCosThetaBins, "hPassedCosTheta_Rebin", cosThetaBinEdges.data());
	TH1D* hTotalCosTheta_Rebin = (TH1D*)hTotalCosTheta->Rebin(nCosThetaBins, "hTotalCosTheta_Rebin", cosThetaBinEdges.data());

	// define TEfficiency using the final numerator and denominator
	TEfficiency* TEffCosTheta = new TEfficiency("TEffCosTheta", "cos #theta; efficiency", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);

	TEffCosTheta->SetPassedHistogram(*hPassedCosTheta_Rebin, "f");
	TEffCosTheta->SetTotalHistogram(*hTotalCosTheta_Rebin, "f");

	return TEffCosTheta;
}

// rebin TEfficiency 3D maps to TEfficiency 1D phi based on costheta, phi, and pT selection

TEfficiency* rebinTEff3DMapPhi(TEfficiency* TEff3DMap, Int_t ptMin = 0, Int_t ptMax = 30, Double_t cosThetaMin = -1, Double_t cosThetaMax = 1, Int_t nPhiBins = 10, const std::vector<Double_t>& phiBinEdges = {}) {
	// extract the numerator and the denominator from the 3D TEfficiency Map
	TH3D* hPassed = (TH3D*)TEff3DMap->GetPassedHistogram();
	TH3D* hTotal = (TH3D*)TEff3DMap->GetTotalHistogram();

	// obtain the bin numbers of the boundaries on phi and pt
	Int_t iCosThetaMin = hPassed->GetXaxis()->FindBin(cosThetaMin);
	Int_t iCosThetaMax = hPassed->GetXaxis()->FindBin(cosThetaMax);

	Int_t iPtMin = hPassed->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = hPassed->GetZaxis()->FindBin(ptMax);

	// obtain the projection histogram along the costheta axis within boundaries of phi and pt
	// (option e: calculate errors, o: only bins inside the selected range will be filled)
	TH1D* hPassedPhi = (TH1D*)hPassed->ProjectionY("hPassedPhi", iCosThetaMin, iCosThetaMax - 1, iPtMin, iPtMax - 1, "eo");
	TH1D* hTotalPhi = (TH1D*)hTotal->ProjectionY("hTotalPhi", iCosThetaMin, iCosThetaMax - 1, iPtMin, iPtMax - 1, "eo");

	if (iCosThetaMin >= iCosThetaMax - 1) {
		hPassedPhi->Reset();
		hTotalPhi->Reset();
	}

	// rebin the projection histogram (default: 18 bins from -180 to 180)
	TH1D* hPassedPhi_Rebin = (TH1D*)hPassedPhi->Rebin(nPhiBins, "hPassedPhi_Rebin", phiBinEdges.data());
	TH1D* hTotalPhi_Rebin = (TH1D*)hTotalPhi->Rebin(nPhiBins, "hTotalPhi_Rebin", phiBinEdges.data());

	// define TEfficiency using the final numerator and denominator
	TEfficiency* TEffPhi = new TEfficiency("TEffPhi", "#varphi; efficiency", nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins]);

	TEffPhi->SetPassedHistogram(*hPassedPhi_Rebin, "f");
	TEffPhi->SetTotalHistogram(*hTotalPhi_Rebin, "f");

	return TEffPhi;
}

// rebin TEfficiency 3D maps to TEfficiency 2D map (cosTheta, phi) based on costheta, phi, and pT selection

TEfficiency* rebinTEff3DMap(TEfficiency* TEff3DMap, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 10, const std::vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 10, const std::vector<Double_t>& phiBinEdges = {}) {
	// define TEfficiency that will contain cosTheta-phi 2D TEfficiency map
	TEfficiency* TEffCosThetaPhi = new TEfficiency(TEff3DMap->GetName(), "cos #theta;#varphi; efficiency", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges.data());

	TH2D* hTotalCosThetaPhi = (TH2D*)TEffCosThetaPhi->GetTotalHistogram();

	// loop over the rebinning code along CosTheta and set the contents of the 2D map row by row
	for (Int_t iRow = 0; iRow < nPhiBins; iRow++) {
		TEfficiency* TEffCosTheta = rebinTEff3DMapCosTheta(TEff3DMap, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, (Int_t)phiBinEdges[iRow], (Int_t)phiBinEdges[iRow + 1]);
		// std::cout << "Phi: " << phiBinEdges[iRow] << " - " << phiBinEdges[iRow + 1] << std::endl;

		TH1D* hPassedCosTheta = (TH1D*)TEffCosTheta->GetPassedHistogram();
		TH1D* hTotalCosTheta = (TH1D*)TEffCosTheta->GetTotalHistogram();

		for (Int_t iColumn = 0; iColumn < nCosThetaBins; iColumn++) {
			Double_t passedEvents = hPassedCosTheta->GetBinContent(iColumn + 1);
			Double_t totalEvents = hTotalCosTheta->GetBinContent(iColumn + 1);

			Double_t binCenterCosTheta = hTotalCosThetaPhi->GetXaxis()->GetBinCenter(iColumn + 1);
			Double_t binCenterPhi = hTotalCosThetaPhi->GetYaxis()->GetBinCenter(iRow + 1);

			Int_t iGlobalBin = TEffCosThetaPhi->FindFixBin(binCenterCosTheta, binCenterPhi);

			TEffCosThetaPhi->SetTotalEvents(iGlobalBin, totalEvents);
			TEffCosThetaPhi->SetPassedEvents(iGlobalBin, passedEvents);

			// std::cout << "CosTheta: " << binCenterCosTheta << " Phi: " << binCenterPhi << " Passed: " << passedEvents << " Total: " << totalEvents << std::endl;
		}
	}

	return TEffCosThetaPhi;
}

// rebin TEfficiency 3D maps of efficiency systematic uncertainty to TEfficiency 1D cosTheta based on costheta, phi, and pT selection

TH1D* rebinRel3DUncCosTheta(TEfficiency* effMap, TH3D* systEff, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 10, const std::vector<Double_t>& cosThetaBinEdges = {}, Int_t phiMin = -180, Int_t phiMax = 180) {
	/// rebin efficiency maps based on costheta, phi, and pT selection
	// uncertainty addition is sqrt(pow(unc1, 2) + pow(unc2, 2)), so fold it manually

	TH1D* h1DSystEffCosTheta = new TH1D("h1DSystEffCosTheta", "", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);

	// obtain the bin numbers of the boundaries on phi and pt
	Int_t iPhiMin = systEff->GetYaxis()->FindBin(phiMin);
	Int_t iPhiMax = systEff->GetYaxis()->FindBin(phiMax) - 1;

	Int_t iPtMin = systEff->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = systEff->GetZaxis()->FindBin(ptMax) - 1;

	// calculate the systematic uncertainties merging phi and pt bins
	// and set the bin content of 1D costheta hist using the calculated value
	for (int ibin = 1; ibin <= nCosThetaBins; ibin++) {
		Double_t cosThetaSumSystEff = 0;

		Int_t iCosThetaMin = systEff->GetXaxis()->FindBin(cosThetaBinEdges[ibin - 1]);
		Int_t iCosThetaMax = systEff->GetXaxis()->FindBin(cosThetaBinEdges[ibin]) - 1;

		// merge bins along the cosTheta axis within the bin width
		for (int iCosTheta = iCosThetaMin; iCosTheta <= iCosThetaMax; iCosTheta++) {
			Double_t phiSumSystEff = 0;

			// sum uncertainties along the phi axis
			for (int iPhi = iPhiMin; iPhi <= iPhiMax; iPhi++) {
				Double_t ptSumSystEff = 0;

				// sum uncertainties along the pt axis
				for (int iPt = iPtMin; iPt <= iPtMax; iPt++) {
					Int_t globalBin = effMap->GetGlobalBin(iCosTheta, iPhi, iPt);

					ptSumSystEff = TMath::Hypot(ptSumSystEff, systEff->GetBinContent(iCosTheta, iPhi, iPt) /* effMap->GetEfficiency(globalBin)*/);
					// cout << "systEff->GetBinContent(" << iCosTheta << ", " << iPhi << ", " << iPt << ") = " << systEff->GetBinContent(iCosTheta, iPhi, iPt) * effMap->GetEfficiency(globalBin) << endl;
					// cout << "ptSumSystEff = " << ptSumSystEff << endl;
				}

				phiSumSystEff = TMath::Hypot(phiSumSystEff, ptSumSystEff);
				// cout << "phiSumSystEff = " << phiSumSystEff << endl;
			}

			cosThetaSumSystEff = TMath::Hypot(cosThetaSumSystEff, phiSumSystEff);
			// cout << "cosThetaSumSystEff = " << cosThetaSumSystEff << endl;
		}

		h1DSystEffCosTheta->SetBinContent(ibin, cosThetaSumSystEff);
		// cout << "ibin = " << ibin << " cosThetaSumSystEff = " << cosThetaSumSystEff << endl;
	}

	return h1DSystEffCosTheta;
}

// rebin TEfficiency 3D maps of efficiency systematic uncertainty to TEfficiency 1D phi based on costheta, phi, and pT selection

TH1D* rebinRel3DUncPhi(TEfficiency* effMap, TH3D* systEff, Int_t ptMin = 0, Int_t ptMax = 30, Double_t cosThetaMin = -1, Double_t cosThetaMax = 1, Int_t nPhiBins = 6, const std::vector<Double_t>& phiBinEdges = {}) {
	/// rebin efficiency maps based on costheta, phi, and pT selection
	// uncertainty addition is sqrt(pow(unc1, 2) + pow(unc2, 2)), so fold it manually

	TH1D* h1DSystEffPhi = new TH1D("h1DSystEffPhi", "", nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins]);

	// obtain the bin numbers of the boundaries on phi and pt
	Int_t iCosThetaMin = systEff->GetXaxis()->FindBin(cosThetaMin);
	Int_t iCosThetaMax = systEff->GetXaxis()->FindBin(cosThetaMax) - 1;

	Int_t iPtMin = systEff->GetZaxis()->FindBin(ptMin);
	Int_t iPtMax = systEff->GetZaxis()->FindBin(ptMax) - 1;

	// calculate the systematic uncertainties merging cosTheta and pt bins
	// and set the bin content of 1D phi hist using the calculated value
	for (int ibin = 1; ibin <= nPhiBins; ibin++) {
		Double_t phiSumSystEff = 0;

		Int_t iPhiMin = systEff->GetYaxis()->FindBin(phiBinEdges[ibin - 1]);
		Int_t iPhiMax = systEff->GetYaxis()->FindBin(phiBinEdges[ibin]) - 1;

		// merge bins along the phi axis within the bin width
		for (int iPhi = iPhiMin; iPhi <= iPhiMax; iPhi++) {
			Double_t cosThetaSumSystEff = 0;

			// sum uncertainties along the cosTheta axis
			for (int iCosTheta = iCosThetaMin; iCosTheta <= iCosThetaMax; iCosTheta++) {
				Double_t ptSumSystEff = 0;

				// sum uncertainties along the pt axis
				for (int iPt = iPtMin; iPt <= iPtMax; iPt++) {
					Int_t globalBin = effMap->GetGlobalBin(iCosTheta, iPhi, iPt);

					ptSumSystEff = TMath::Hypot(ptSumSystEff, systEff->GetBinContent(iCosTheta, iPhi, iPt) /* effMap->GetEfficiency(globalBin)*/);
				}

				cosThetaSumSystEff = TMath::Hypot(cosThetaSumSystEff, ptSumSystEff);
			}

			phiSumSystEff = TMath::Hypot(phiSumSystEff, cosThetaSumSystEff);
		}

		h1DSystEffPhi->SetBinContent(ibin, phiSumSystEff);
	}

	return h1DSystEffPhi;
}

// rebin TEfficiency 3D maps of efficiency systematic uncertainty to TEfficiency 2D map based on costheta, phi, and pT selection

TH2D* rebinRel3DUncMap(TEfficiency* effMap, TH3D* systEff, Int_t ptMin = 0, Int_t ptMax = 30, Int_t nCosThetaBins = 5, const std::vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 6, const std::vector<Double_t>& phiBinEdges = {}) {
	/// rebin efficiency maps based on costheta, phi, and pT selection
	// uncertainty addition is sqrt(pow(unc1, 2) + pow(unc2, 2)), so fold it manually

	TH2D* h2DSystEffMap = new TH2D("h2DSystEffMap", "", nCosThetaBins, cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins], nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins]);

	// calculate the systematic uncertainties merging cosTheta and pt bins
	// and set the bin content of 2D cosTheta-phi uncertainty hist using the calculated value

	for (int iPhi = 1; iPhi <= nPhiBins; iPhi++) {
		TH1D* h1DSystEffCosTheta = rebinRel3DUncCosTheta(effMap, systEff, ptMin, ptMax, nCosThetaBins, cosThetaBinEdges, phiBinEdges[iPhi - 1], phiBinEdges[iPhi]);

		for (int iCosTheta = 1; iCosTheta <= nCosThetaBins; iCosTheta++) {
			h2DSystEffMap->SetBinContent(iCosTheta, iPhi, h1DSystEffCosTheta->GetBinContent(iCosTheta));
			// cout << "iCosTheta = " << iCosTheta << " iPhi = " << iPhi << " h1DSystEffCosTheta->GetBinContent(iCosTheta) = " << h1DSystEffCosTheta->GetBinContent(iCosTheta) << endl;
		}
	}

	return h2DSystEffMap;
}

// display the uncertainties signal extraction yield on each bin of 2D yield map

void displayEfficiencies(TEfficiency* effMap, Int_t nCosThetaBins = 10, Int_t nPhiBins = 6, Bool_t displayErrors = kFALSE) {
	if (!effMap) {
		std::cout << "no efficiency map found!!!" << std::endl;
		exit(1);
	}

	TH2D* hTotal2D = (TH2D*)effMap->GetTotalHistogram();

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotal2D->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotal2D->GetYaxis()->GetBinCenter(iPhi + 1);

			Int_t iGlobalBin = hTotal2D->FindFixBin(binCenterCosTheta, binCenterPhi);

			// Get the yield and uncertainty values
			Double_t effVal = effMap->GetEfficiency(iGlobalBin);

			Double_t effUncUp = effMap->GetEfficiencyErrorUp(iGlobalBin);
			Double_t effUncDown = effMap->GetEfficiencyErrorLow(iGlobalBin);

			// Get the bin center coordinates
			double x = hTotal2D->GetXaxis()->GetBinCenter(iCosTheta + 1);

			double y = hTotal2D->GetYaxis()->GetBinCenter(iPhi + 1);

			// Create a TLatex object to write the signal extraction yield uncertainties on each bin
			TLatex latex;
			latex.SetTextAlign(22); // Center alignment
			latex.SetTextColor(kWhite);
			if (displayErrors == kTRUE) {
				latex.SetTextSize(0.05);
				latex.DrawLatex(x, y, Form(" ^{+%.4f}_{-%.4f}", effUncUp, effUncDown));
			} else {
				latex.SetTextSize(0.05);
				latex.DrawLatex(x, y, Form("%.2f", effVal));
			}
		}
	}
}

void displayYields(TEfficiency* effMap, Int_t nCosThetaBins = 10, Int_t nPhiBins = 6, std::string prefix = "", std::string postfix = "") {
	if (!effMap) {
		std::cout << "no efficiency map found!!!" << std::endl;
		exit(1);
	}

	TH2D* hTotal2D = (TH2D*)effMap->GetTotalHistogram();

	for (Int_t iCosTheta = 0; iCosTheta < nCosThetaBins; iCosTheta++) {
		for (Int_t iPhi = 0; iPhi < nPhiBins; iPhi++) {
			// get the global bin number of Efficiency
			Double_t binCenterCosTheta = hTotal2D->GetXaxis()->GetBinCenter(iCosTheta + 1);
			Double_t binCenterPhi = hTotal2D->GetYaxis()->GetBinCenter(iPhi + 1);

			Double_t binEdgeCosThetaMin = hTotal2D->GetXaxis()->GetBinLowEdge(iCosTheta + 1);
			Double_t binEdgeCosThetaMax = hTotal2D->GetXaxis()->GetBinUpEdge(iCosTheta + 1);

			Int_t binEdgePhiMin = (int)(hTotal2D->GetYaxis()->GetBinLowEdge(iPhi + 1));
			Int_t binEdgePhiMax = (int)(hTotal2D->GetYaxis()->GetBinUpEdge(iPhi + 1));

			// Get the yield and uncertainty values

			RooRealVar yield1S("yield1S", "", 1000);
			RooRealVar yield2S("yield2S", "", 100);
			RooRealVar yield3S("yield3S", "", 10);

			RooRealVar yield1S_freeSigma("yield1S", "", 1000);
			RooRealVar yield2S_freeSigma("yield2S", "", 100);
			RooRealVar yield3S_freeSigma("yield3S", "", 10);

			// RooRealVar* yield1S_freeSigma = new RooRealVar("yield1S", "", 1000);
			// RooRealVar* yield2S_freeSigma = new RooRealVar("yield2S", "", 100);
			// RooRealVar* yield3S_freeSigma = new RooRealVar("yield3S", "", 10);
			// cout << "prefix = " << prefix << endl;
			// cout << "postfix = " << postfix << endl;
			// cout << "effMap name = " << effMap->GetName() << endl;
			// cout << Form("%s_cosTheta%.2fto%.2f_phi%dto%d_%s", prefix, binEdgeCosThetaMin, binEdgeCosThetaMax, binEdgePhiMin, binEdgePhiMax, postfix) << endl;
			// const char* fileName = Form("%s_cosTheta%.2fto%.2f_phi%dto%d_%s.txt", prefix, binEdgeCosThetaMin, binEdgeCosThetaMax, binEdgePhiMin, binEdgePhiMax, postfix);
			// const char* fileName_freeSigma = Form("%s_cosTheta%.2fto%.2f_phi%dto%d_%s_%s.txt", prefix, binEdgeCosThetaMin, binEdgeCosThetaMax, binEdgePhiMin, binEdgePhiMax, postfix, "woSigmaConstraints");

			std::string fileName = Form("%s_cosTheta%.2fto%.2f_phi%dto%d_%s.txt", prefix.c_str(), binEdgeCosThetaMin, binEdgeCosThetaMax, binEdgePhiMin, binEdgePhiMax, postfix.c_str());
			std::string fileName_freeSigma = Form("%s_cosTheta%.2fto%.2f_phi%dto%d_%s_%s.txt", prefix.c_str(), binEdgeCosThetaMin, binEdgeCosThetaMax, binEdgePhiMin, binEdgePhiMax, postfix.c_str(), "woSigmaConstraints");

			// get yields and their uncertainties
			RooArgSet signalYields = GetSignalYields(&yield1S, &yield2S, &yield3S, fileName.c_str());
			cout << "fileName_freeSigma = " << fileName_freeSigma.c_str() << endl;
			RooArgSet signalYields_freeSigma = GetSignalYields(&yield1S_freeSigma, &yield2S_freeSigma, &yield3S_freeSigma, fileName_freeSigma.c_str());

			double yield1Sdiff = (yield1S.getVal()) - (yield1S_freeSigma.getVal());

			// Get the bin center coordinates
			double x = hTotal2D->GetXaxis()->GetBinCenter(iCosTheta + 1);

			double y = hTotal2D->GetYaxis()->GetBinCenter(iPhi + 1);

			// Create a TLatex object to write the signal extraction yield uncertainties on each bin
			TLatex latex;
			latex.SetTextSize(0.05); // Adjust text size as needed
			latex.SetTextAlign(22);  // Center alignment
			latex.SetTextColor(kWhite);
			latex.DrawLatex(x, y, Form("%.0f", yield1Sdiff));
		}
	}
}

// Draw 1D efficincy plot

void DrawEfficiency1DHist(TEfficiency* effHist, Int_t ptMin, Int_t ptMax, Int_t iState = gUpsilonState, Bool_t isAcc = kTRUE, Bool_t isCosTheta = kTRUE, const char* refFrameName = "CS", Float_t lambdaTheta = 0, Float_t lambdaPhi = 0, Float_t lambdaThetaPhi = 0, const char* extraString = "", Bool_t isPhiFolded = kTRUE) {
	TCanvas* canvas = new TCanvas(effHist->GetName(), "", 600, 600);
	canvas->SetRightMargin(0.05);

	// empty frame for the axes
	TH1D* frameHist;

	if (isCosTheta) {
		frameHist = new TH1D("frameHist", "", NCosThetaBinsHX, CosThetaBinningHX);
	} else {
		frameHist = new TH1D("frameHist", "", NPhiBinsHX, PhiBinningHX);
	}

	frameHist->Draw();

	// draw efficiency plot on top of the histogram frame
	effHist->SetLineWidth(3);
	effHist->Draw("PL E0 SAME");

	// cosmetics of the histogram
	CMS_lumi(canvas, Form("#varUpsilon(%dS) Pythia 8 (5.02 TeV)", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.55, .87, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));

	if (strcmp(extraString, "_TriggerAcc") == 0)
		legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	else if (strcmp(extraString, "_SimpleAcc") == 0)
		legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", iState));
	else if (strcmp(extraString, "_2018PbPbAcc") == 0)
		legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", iState));

	// legend.DrawLatexNDC(.55, .72, Form("#lambda_{#theta} = %.2f, #lambda_{#varphi} = %.2f, #lambda_{#theta#varphi} = %.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi));

	if (isCosTheta) {
		frameHist->SetXTitle(CosThetaVarTitle(refFrameName));

	}

	else {
		if (isPhiFolded == kTRUE)
			frameHist->SetXTitle(AbsPhiAxisTitle(refFrameName));
		else
			frameHist->SetXTitle(PhiVarTitle(refFrameName));
	}

	if (isAcc)
		frameHist->SetYTitle(TEfficiencyMainTitle(iState, "Acceptance"));
	else
		frameHist->SetYTitle(TEfficiencyMainTitle(iState));

	frameHist->GetXaxis()->CenterTitle();
	frameHist->GetYaxis()->CenterTitle();

	// frameHist->GetXaxis()->SetRangeUser(XBinning[0], XBinning[NXBins]);
	frameHist->GetYaxis()->SetRangeUser(0, 1);

	frameHist->GetXaxis()->SetNdivisions(510, kTRUE);

	// // save the plot
	// if (isAcc) {
	// 	gSystem->mkdir(Form("AcceptanceMaps/%dS", iState), kTRUE);
	// 	canvas->SaveAs(Form("AcceptanceMaps/%dS/1Dacc_%s%s_pt%dto%d%s.png", iState, effHist->GetName(), refFrameName, ptMin, ptMax, extraString), "RECREATE");
	// } else {
	// 	gSystem->mkdir(Form("EfficiencyMaps/%dS", iState), kTRUE);
	// 	canvas->SaveAs(Form("EfficiencyMaps/%dS/1Deff_%s%s_pt%dto%d%s.png", iState, effHist->GetName(), refFrameName, ptMin, ptMax, extraString), "RECREATE");
	// }

	// save the plot
	if (isAcc) {
		gSystem->mkdir(Form("AcceptanceMaps/%dS/analysisBin", iState), kTRUE);
		canvas->SaveAs(Form("AcceptanceMaps/%dS/analysisBin/1Dacc_%s%s_pt%dto%d%s.png", iState, effHist->GetName(), refFrameName, ptMin, ptMax, extraString), "RECREATE");
	} else {
		gSystem->mkdir(Form("EfficiencyMaps/%dS/analysisBin", iState), kTRUE);
		canvas->SaveAs(Form("EfficiencyMaps/%dS/analysisBin/1Deff_%s%s_pt%dto%d%s.png", iState, effHist->GetName(), refFrameName, ptMin, ptMax, extraString), "RECREATE");
	}
}

void DrawEfficiency2DHist(TEfficiency* effHist, Int_t ptMin, Int_t ptMax, Int_t nCosThetaBins = 5, const std::vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 5, const std::vector<Double_t>& phiBinEdges = {}, Int_t iState = gUpsilonState, Bool_t isAcc = kTRUE, Bool_t displayValues = kFALSE, Bool_t displayErrors = kFALSE, const char* extraString = "", Bool_t isPhiFolded = kTRUE) {
	TCanvas* canvas;
	if (isAcc)
		canvas = new TCanvas("accCosThetaPhi", "", 680, 600);
	else
		canvas = new TCanvas("effCosThetaPhi", "", 680, 600);
	canvas->SetRightMargin(0.18);

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	// empty frame for the axes
	TH2D* frameHist2D = new TH2D("frameHist2D", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins] + phiStep);

	frameHist2D->Draw("COLZ");

	// draw efficiency plot on top of the histogram frame
	effHist->SetLineWidth(3);
	effHist->Draw("COLZ SAME");

	// cosmetics of the histogram
	CMS_lumi(canvas, Form("#varUpsilon(%dS) Pythia 8 (5.02 TeV)", iState));

	if (displayValues) displayEfficiencies(effHist, nCosThetaBins, nPhiBins, displayErrors);

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.04);
	legend.DrawLatexNDC(.48, .86, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));

	if (strcmp(extraString, "_TriggerAcc") == 0)
		legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	else if (strcmp(extraString, "_SimpleAcc") == 0)
		legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", iState));
	else if (strcmp(extraString, "_2018PbPbAcc") == 0)
		legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", iState));

	if (strstr(effHist->GetName(), "CS")) {
		frameHist2D->SetXTitle(CosThetaVarTitle("CS"));
		if (isPhiFolded == kTRUE)
			frameHist2D->SetYTitle(AbsPhiAxisTitle("CS"));
		else
			frameHist2D->SetYTitle(PhiAxisTitle("CS"));
	} else {
		frameHist2D->SetXTitle(CosThetaVarTitle("HX"));
		if (isPhiFolded == kTRUE)
			frameHist2D->SetYTitle(AbsPhiAxisTitle("HX"));
		else
			frameHist2D->SetYTitle(PhiAxisTitle("HX"));
	}

	if (isAcc) {
		effHist->SetTitle(Form(";;;%s", TEfficiencyMainTitle(iState, "Acceptance")));
		SetColorPalette(gAcceptanceColorPaletteName);
	} else {
		effHist->SetTitle(Form(";;;%s", TEfficiencyMainTitle(iState)));
		SetColorPalette(gEfficiencyColorPaletteName);
	}

	frameHist2D->GetXaxis()->CenterTitle();
	frameHist2D->GetYaxis()->CenterTitle();

	frameHist2D->GetXaxis()->SetRangeUser(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	frameHist2D->GetYaxis()->SetRangeUser(phiBinEdges[0], phiBinEdges[nPhiBins] + phiStep);
	frameHist2D->GetZaxis()->SetRangeUser(0, 1);

	// frameHist2D->GetXaxis()->SetNdivisions(-500 - (nCosThetaBins / 2.));
	// frameHist2D->GetYaxis()->SetNdivisions(-500 - (nPhiBins / 2.));

	// frameHist2D->GetXaxis()->SetNdivisions(-500 - (5));
	// frameHist2D->GetYaxis()->SetNdivisions(-500 - (5 + 1));

	/// Ndivision setting for the analysis binning (costheta: 5 bins from -0.7 to 0.7, phi: 7 bins from -180 to 240)
	frameHist2D->GetXaxis()->SetNdivisions(-500 - (nCosThetaBins));
	frameHist2D->GetYaxis()->SetNdivisions(-500 - (nPhiBins)-1);

	// /// Ndivision setting for the finer binning (costheta: 20 bins from -1 to 1, phi: 23 bins from -180 to 280)
	// frameHist2D->GetXaxis()->SetNdivisions(-210); //cosTheta
	// frameHist2D->GetYaxis()->SetNdivisions(-308); //phi

	frameHist2D->GetXaxis()->SetLabelSize(0.045);
	frameHist2D->GetYaxis()->SetLabelSize(0.045);

	// // save the plot
	// gSystem->mkdir(Form("EfficiencyMaps/%dS", iState), kTRUE);
	// if (isAcc) {
	// 	if (isPhiFolded == kTRUE) canvas->SaveAs(Form("EfficiencyMaps/%dS/2Dacc_%s_pt%dto%d%s.png", iState, effHist->GetName(), ptMin, ptMax, extraString), "RECREATE");
	// 	else canvas->SaveAs(Form("EfficiencyMaps/%dS/2Dacc_%s_pt%dto%d%s_fullPhi.png", iState, effHist->GetName(), ptMin, ptMax, extraString), "RECREATE");
	// }
	// else {
	// 	if (isPhiFolded == kTRUE)  canvas->SaveAs(Form("EfficiencyMaps/%dS/2Deff_%s_pt%dto%d%s.png", iState, effHist->GetName(), ptMin, ptMax, extraString), "RECREATE");
	//     else canvas->SaveAs(Form("EfficiencyMaps/%dS/2Deff_%s_pt%dto%d%s_fullPhi.png", iState, effHist->GetName(), ptMin, ptMax, extraString), "RECREATE");

	// }

	// save the plot
	gSystem->mkdir(Form("EfficiencyMaps/%dS/analysisBin", iState), kTRUE);
	if (isAcc) {
		if (isPhiFolded == kTRUE)
			canvas->SaveAs(Form("EfficiencyMaps/%dS/analysisBin/2Dacc_%s_pt%dto%d%s.png", iState, effHist->GetName(), ptMin, ptMax, extraString), "RECREATE");
		else
			canvas->SaveAs(Form("EfficiencyMaps/%dS/analysisBin/2Dacc_%s_pt%dto%d%s_fullPhi.png", iState, effHist->GetName(), ptMin, ptMax, extraString), "RECREATE");
	} else {
		if (isPhiFolded == kTRUE)
			canvas->SaveAs(Form("EfficiencyMaps/%dS/analysisBin/2Deff_%s_pt%dto%d%s.png", iState, effHist->GetName(), ptMin, ptMax, extraString), "RECREATE");
		else
			canvas->SaveAs(Form("EfficiencyMaps/%dS/analysisBin/2Deff_%s_pt%dto%d%s_fullPhi.png", iState, effHist->GetName(), ptMin, ptMax, extraString), "RECREATE");
	}
}

TEfficiency* MultiplyTEfficiencies1D(const TEfficiency* eff1, const TEfficiency* eff2) {
	if (!eff1 || !eff2) return nullptr;

	auto* h1 = eff1->GetCopyTotalHisto();
	auto* h2 = eff2->GetCopyTotalHisto();

	if (h1->GetNbinsX() != h2->GetNbinsX()) {
		std::cerr << "Error: The efficiency histograms have different binning!" << std::endl;
		return nullptr;
	}

	auto* result = new TEfficiency(*eff1); // Clone structure

	for (int i = 1; i <= h1->GetNbinsX(); i++) {
		double e1 = eff1->GetEfficiency(h1->GetBin(i));
		// cout << "e1 = " << e1 << endl;
		// cout << "h1->GetBinContent(" << i << ", " << j << ") = " << h1->GetBinContent(i, j) << endl;
		double e2 = eff2->GetEfficiency(h2->GetBin(i));
		// cout << "e2 = " << e2 << endl;
		// cout << "h2->GetBinContent(" << i << ", " << j << ") = " << h2->GetBinContent(i, j) << endl;

		double err1 = eff1->GetEfficiencyErrorUp(h1->GetBin(i));
		double err2 = eff2->GetEfficiencyErrorUp(h2->GetBin(i));

		if (e1 == 0 || e2 == 0) continue; // Avoid division by zero

		double eff_product = e1 * e2;
		double err_product = eff_product * sqrt(pow(err1 / e1, 2) + pow(err2 / e2, 2)); // Error propagation

		int bin = h1->GetBin(i);
		result->SetTotalEvents(bin, h1->GetBinContent(bin));
		result->SetPassedEvents(bin, eff_product * h1->GetBinContent(bin));
		// cout << "accxeff value" << result->GetEfficiency(bin) << endl;
	}

	return result;
}

TEfficiency* MultiplyTEfficiencies2D(const TEfficiency* eff1, const TEfficiency* eff2) {
	if (!eff1 || !eff2) return nullptr;

	auto* h1 = eff1->GetCopyTotalHisto();
	auto* h2 = eff2->GetCopyTotalHisto();

	if (h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsY() != h2->GetNbinsY()) {
		std::cerr << "Error: The efficiency histograms have different binning!" << std::endl;
		return nullptr;
	}

	auto* result = new TEfficiency(*eff1); // Clone structure

	for (int i = 1; i <= h1->GetNbinsX(); i++) {
		for (int j = 1; j <= h1->GetNbinsY(); j++) {
			double e1 = eff1->GetEfficiency(h1->GetBin(i, j));
			// cout << "e1 = " << e1 << endl;
			// cout << "h1->GetBinContent(" << i << ", " << j << ") = " << h1->GetBinContent(i, j) << endl;
			double e2 = eff2->GetEfficiency(h2->GetBin(i, j));
			// cout << "e2 = " << e2 << endl;
			// cout << "h2->GetBinContent(" << i << ", " << j << ") = " << h2->GetBinContent(i, j) << endl;

			double err1 = eff1->GetEfficiencyErrorUp(h1->GetBin(i, j));
			double err2 = eff2->GetEfficiencyErrorUp(h2->GetBin(i, j));

			if (e1 == 0 || e2 == 0) continue; // Avoid division by zero

			double eff_product = e1 * e2;
			double err_product = eff_product * sqrt(pow(err1 / e1, 2) + pow(err2 / e2, 2)); // Error propagation

			int bin = h1->GetBin(i, j);
			result->SetTotalEvents(bin, h1->GetBinContent(bin));
			result->SetPassedEvents(bin, eff_product * h1->GetBinContent(bin));
		}
	}

	return result;
}

double GetMaxEfficiencyValue(TEfficiency& eff) {
	// Get the number of bins for both X and Y axes
	int nBinsX = eff.GetTotalHistogram()->GetNbinsX();
	int nBinsY = eff.GetTotalHistogram()->GetNbinsY();

	double maxEfficiency = 0.0;

	for (int i = 1; i <= nBinsX; ++i) {
		for (int j = 1; j <= nBinsY; ++j) {
			int globalbin = eff.GetGlobalBin(i, j);

			double efficiency = eff.GetEfficiency(globalbin); // Get efficiency for bin (i, j)
			if (efficiency > maxEfficiency) {
				maxEfficiency = efficiency;
			}
		}
	}
	std::cout << "Maximum efficiency: " << maxEfficiency << std::endl;

	return maxEfficiency;
}

TEfficiency* DrawEffxAcc1DHist(TEfficiency* accHist, TEfficiency* effHist, Int_t ptMin, Int_t ptMax, Int_t nBins = 5, const std::vector<Double_t>& BinEdges = {}, Double_t intRangeMin = -1, Double_t intRangeMax = 1, Bool_t isCosTheta = kTRUE, const char* refFrameName = "CS", Int_t iState = gUpsilonState, const char* extraString = "", Bool_t isPhiFolded = kFALSE) {
	TCanvas* canvas = new TCanvas("EffxAccCosThetaPhi", "", 600, 600);
	canvas->SetRightMargin(0.05);

	TEfficiency* EffxAccHist = MultiplyTEfficiencies1D(accHist, effHist);

	// empty frame for the axes
	TH1D* frameHist1D = new TH1D("frameHist1D", "", nBins, BinEdges.data());

	frameHist1D->Draw();

	// draw efficiency plot on top of the histogram frame
	EffxAccHist->SetLineWidth(3);
	EffxAccHist->Draw("PE SAME");
	EffxAccHist->SetMarkerStyle(20);
	EffxAccHist->SetMarkerSize(1.5);
	EffxAccHist->SetMarkerColor(TColor::GetColor("#A559AA"));

	EffxAccHist->SetLineColor(TColor::GetColor("#A559AA"));
	EffxAccHist->SetLineWidth(3);

	// cosmetics of the histogram
	CMS_lumi(canvas, Form("#varUpsilon(%dS) Pythia 8 (5.02 TeV)", iState));

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.05);
	legend.DrawLatexNDC(.55, .87, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));

	if (strcmp(extraString, "_TriggerAcc") == 0)
		legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	else if (strcmp(extraString, "_SimpleAcc") == 0)
		legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", iState));
	else if (strcmp(extraString, "_2018PbPbAcc") == 0)
		legend.DrawLatexNDC(.55, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", iState));

	if (isCosTheta == kTRUE)
		legend.DrawLatexNDC(.55, .72, Form("%.0f < #varphi < %.0f #circ", intRangeMin, intRangeMax));
	else
		legend.DrawLatexNDC(.55, .72, Form("%.2f < cos#theta < %.2f", intRangeMin, intRangeMax));

	if (isCosTheta) {
		frameHist1D->SetXTitle(CosThetaVarTitle(refFrameName));

	}

	else {
		if (isPhiFolded == kTRUE)
			frameHist1D->SetXTitle(AbsPhiAxisTitle(refFrameName));
		else
			frameHist1D->SetXTitle(PhiAxisTitle(refFrameName));
	}

	frameHist1D->SetYTitle(TEfficiencyMainTitle(iState, "acceptance x efficiency"));

	frameHist1D->GetXaxis()->CenterTitle();
	frameHist1D->GetYaxis()->CenterTitle();
	// frameHist1D->GetYaxis()->SetTitleOffset(1.2);

	frameHist1D->GetXaxis()->SetRangeUser(BinEdges[0], BinEdges[nBins]);
	frameHist1D->GetYaxis()->SetRangeUser(0, 1);

	// // save the plot
	// gSystem->mkdir(Form("AccxEffMaps/%dS", iState), kTRUE);

	// if (isPhiFolded == kTRUE) canvas->SaveAs(Form("AccxEffMaps/%dS/1DAccxEff_%s%s_pt%dto%d.png", iState, EffxAccHist->GetName(), refFrameName, ptMin, ptMax), "RECREATE");
	// else canvas->SaveAs(Form("AccxEffMaps/%dS/1DAccxEff_%s%s_pt%dto%d%s%s.png", iState, EffxAccHist->GetName(), refFrameName, ptMin, ptMax, extraString, "_fullPhi"), "RECREATE");

	// save the plot
	gSystem->mkdir(Form("AccxEffMaps/%dS/analysisBin", iState), kTRUE);

	if (isPhiFolded == kTRUE) {
		if (isCosTheta == kTRUE)
			canvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/1DAccxEff_%s%s_pt%dto%d_cosTheta%.2fto%.2f.png", iState, EffxAccHist->GetName(), refFrameName, ptMin, ptMax, intRangeMin, intRangeMax), "RECREATE");
		else
			canvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/1DAccxEff_%s%s_pt%dto%d_phi%.0fto%.0f.png", iState, EffxAccHist->GetName(), refFrameName, ptMin, ptMax, intRangeMin, intRangeMax), "RECREATE");
	} else {
		if (isCosTheta == kTRUE)
			canvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/1DAccxEff_%s%s_pt%dto%d_cosTheta%.2fto%.2f%s%s.png", iState, EffxAccHist->GetName(), refFrameName, ptMin, ptMax, intRangeMin, intRangeMax, extraString, "_fullPhi"), "RECREATE");
		else
			canvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/1DAccxEff_%s%s_pt%dto%d_phi%.0fto%.0f%s%s.png", iState, EffxAccHist->GetName(), refFrameName, ptMin, ptMax, intRangeMin, intRangeMax, extraString, "_fullPhi"), "RECREATE");
	}

	return EffxAccHist;
}

void DrawEffxAcc2DHist(TEfficiency* accHist, TEfficiency* effHist, Int_t ptMin, Int_t ptMax, Int_t nCosThetaBins = 5, const std::vector<Double_t>& cosThetaBinEdges = {}, Int_t nPhiBins = 5, const std::vector<Double_t>& phiBinEdges = {}, Int_t iState = gUpsilonState, Bool_t displayValues = kFALSE, Bool_t displayErrors = kFALSE, Bool_t displayYieldValues = kTRUE, const char* extraString = "", Bool_t isPhiFolded = kFALSE) {
	TCanvas* canvas = new TCanvas("EffxAccCosThetaPhi", "", 680, 600);

	canvas->SetRightMargin(0.18);
	canvas->SetLeftMargin(0.15);
	canvas->SetBottomMargin(0.15);

	Double_t phiStep = (phiBinEdges[nPhiBins] - phiBinEdges[0]) / nPhiBins;

	TEfficiency* EffxAccHist = MultiplyTEfficiencies2D(accHist, effHist);

	// empty frame for the axes
	TH2D* frameHist2D = new TH2D("frameHist2D", "", nCosThetaBins, cosThetaBinEdges.data(), nPhiBins, phiBinEdges[0], phiBinEdges[nPhiBins] + phiStep);

	frameHist2D->Draw("COLZ");

	// draw efficiency plot on top of the histogram frame
	EffxAccHist->SetLineWidth(3);
	EffxAccHist->Draw("COLZ SAME");

	// cosmetics of the histogram
	CMS_lumi(canvas, Form("#varUpsilon(%dS) Pythia 8 (5.02 TeV)", iState));

	if (displayValues) displayEfficiencies(EffxAccHist, nCosThetaBins, nPhiBins, displayErrors);

	TLatex legend;
	legend.SetTextAlign(22);
	legend.SetTextSize(0.04);
	legend.DrawLatexNDC(.48, .86, Form("%s < 2.4, %s", gDimuonRapidityVarTitle, DimuonPtRangeText(ptMin, ptMax)));

	if (strcmp(extraString, "_TriggerAcc") == 0)
		legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, %s", iState, gMuonPtCutText));
	else if (strcmp(extraString, "_SimpleAcc") == 0)
		legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}", iState));
	else if (strcmp(extraString, "_2018PbPbAcc") == 0)
		legend.DrawLatexNDC(.48, .80, Form("#varUpsilon(%dS) acc. for |#eta^{#mu}| < 2.4, #it{p}_{T}^{ #mu} > 2018PbPbAcc", iState));

	const char* refFrameName = "";

	if (strstr(EffxAccHist->GetName(), "CS")) {
		refFrameName = "CS";

		frameHist2D->SetXTitle(CosThetaVarTitle(refFrameName));
		if (isPhiFolded == kTRUE)
			frameHist2D->SetYTitle(AbsPhiAxisTitle(refFrameName));
		else
			frameHist2D->SetYTitle(PhiAxisTitle(refFrameName));
	} else {
		refFrameName = "HX" frameHist2D->SetXTitle(CosThetaVarTitle(refFrameName));
		if (isPhiFolded == kTRUE)
			frameHist2D->SetYTitle(AbsPhiAxisTitle(refFrameName));
		else
			frameHist2D->SetYTitle(PhiAxisTitle(refFrameName));
	}

	EffxAccHist->SetTitle(Form(";;;%s", TEfficiencyMainTitle(iState, "acceptance x efficiency")));
	// SetColorPalette(gAcceptanceColorPaletteName);
	SetColorPalette("sunset");

	frameHist2D->GetXaxis()->CenterTitle();
	frameHist2D->GetYaxis()->CenterTitle();
	frameHist2D->GetYaxis()->SetTitleOffset(1.2);

	frameHist2D->GetXaxis()->SetRangeUser(cosThetaBinEdges[0], cosThetaBinEdges[nCosThetaBins]);
	frameHist2D->GetYaxis()->SetRangeUser(phiBinEdges[0], phiBinEdges[nPhiBins] + phiStep);

	double maxEffxAcc = GetMaxEfficiencyValue(*EffxAccHist);
	frameHist2D->GetZaxis()->SetRangeUser(0, maxEffxAcc);
	// frameHist2D->GetZaxis()->SetRangeUser(0, 1);

	/// Ndivision setting for the analysis binning (costheta: 5 bins from -0.7 to 0.7, phi: 10 bins from -180 to 280)
	frameHist2D->GetXaxis()->SetNdivisions(-500 - (nCosThetaBins));
	frameHist2D->GetYaxis()->SetNdivisions(-500 - (nPhiBins)-1);

	// /// Ndivision setting for the finer binning (costheta: 20 bins from -1 to 1, phi: 23 bins from -180 to 280)
	// frameHist2D->GetXaxis()->SetNdivisions(-210); //cosTheta
	// frameHist2D->GetYaxis()->SetNdivisions(-308); //phi

	frameHist2D->GetXaxis()->SetLabelSize(0.045);
	frameHist2D->GetYaxis()->SetLabelSize(0.045);

	// if (displayYieldsValues) displayYields(EffxAccHist, nCosThetaBins, nPhiBins);
	std::string prefix = Form("RawYields/ExpTimesErr_SymDSCB_cent0to90_absy0p0to2p4_pt%dto%d", ptMin, ptMax);
	std::string postfix = Form("%s%s%s", refFrameName, gMuonAccName, extraString);

	if (displayYieldValues) displayYields(EffxAccHist, nCosThetaBins, nPhiBins, prefix, postfix);

	// // save the plot
	// gSystem->mkdir(Form("AccxEffMaps/%dS", iState), kTRUE);

	// if (isPhiFolded == kTRUE) canvas->SaveAs(Form("AccxEffMaps/%dS/2DAccxEff_%s_pt%dto%d.png", iState, EffxAccHist->GetName(), ptMin, ptMax), "RECREATE");
	// else canvas->SaveAs(Form("AccxEffMaps/%dS/2DAccxEff_%s_pt%dto%d%s%s.png", iState, EffxAccHist->GetName(), ptMin, ptMax, extraString, "_fullPhi"), "RECREATE");
	// // canvas->SaveAs(Form("AccxEffMaps/%dS/2DAccxEff_%s_pt%dto%d%s.png", iState, EffxAccHist->GetName(), ptMin, ptMax, "_diff"), "RECREATE");

	// save the plot
	gSystem->mkdir(Form("AccxEffMaps/%dS/analysisBin", iState), kTRUE);

	if (isPhiFolded == kTRUE)
		canvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/2DAccxEff_%s_pt%dto%d%s.png", iState, EffxAccHist->GetName(), ptMin, ptMax, displayErrors ? "_DisplayUnc" : "_DisplayEff"), "RECREATE");
	else
		canvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/2DAccxEff_%s_pt%dto%d%s%s%s.png", iState, EffxAccHist->GetName(), ptMin, ptMax, extraString, displayErrors ? "_DisplayUnc" : "_DisplayEff", "_fullPhi"), "RECREATE");
	// canvas->SaveAs(Form("AccxEffMaps/%dS/analysisBin/2DAccxEff_%s_pt%dto%d%s.png", iState, EffxAccHist->GetName(), ptMin, ptMax, "_diff"), "RECREATE");
}
