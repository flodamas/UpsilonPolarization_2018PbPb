#include "../Tools/BasicHeaders.h"

#include "../AnalysisParameters.h"

#include "TEfficiency.h"

/// Definition of common helping functions for acceptance and efficiency distribution making codes

// Option to estimate the confidence intervals in TEfficiency
TEfficiency::EStatOption gTEffStatOption = TEfficiency::EStatOption::kFCP; /*default = Clopper-Pearson as recommended by PDG*/

/// Object names (because used here and there) same for acceptance AND efficiency (saved in different files anyway)

// common naming convention
const char* TEfficiencyEndName(const char* refFrameName = "CS", Int_t ptMin = 0, Int_t ptMax = 30) {
	return Form("%s_pt%dto%d", refFrameName, ptMin, ptMax);
}

const char* TEfficiencyMainTitle(int iState = gUpsilonState) {
	return Form("#varUpsilon(%dS) total efficiency", iState);
}

/// TEfficiency constructors

// cos theta distribution for a given pT range
TEfficiency* CosThetaTEfficiency1D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", int iState = gUpsilonState) {
	const char* name = Form("CosTheta%s", TEfficiencyEndName(refFrameName, ptMin, ptMax));

	const char* title = Form(";%s;%s", CosThetaVarTitle(refFrameName), TEfficiencyMainTitle(iState));

	TEfficiency* tEff = new TEfficiency(name, title, NCosThetaBins, gCosThetaMin, gCosThetaMax);

	tEff->SetStatisticOption(gTEffStatOption);

	return tEff;
}

// (cos theta, phi) distribution for a given pT range
const char* RelativeSystTEfficiency2DName(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS") {
	return Form("RelativeSyst2D%s", TEfficiencyEndName(refFrameName, ptMin, ptMax));
}

const char* CosThetaPhiTEfficiency2DName(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS") {
	return Form("CosThetaPhi%s", TEfficiencyEndName(refFrameName, ptMin, ptMax));
}

TEfficiency* CosThetaPhiTEfficiency2D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS", int iState = gUpsilonState) {
	const char* name = CosThetaPhiTEfficiency2DName(ptMin, ptMax, refFrameName);

	const char* title = Form(";%s;%s;%s", CosThetaVarTitle(refFrameName), PhiAxisTitle(refFrameName), TEfficiencyMainTitle(iState));

	TEfficiency* tEff = new TEfficiency(name, title, NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);

	tEff->SetStatisticOption(gTEffStatOption);

	return tEff;
}

TEfficiency* CosThetaPhiAcceptance2D(Int_t ptMin = 0, Int_t ptMax = 30, const char* refFrameName = "CS") {
	const char* name = CosThetaPhiTEfficiency2DName(ptMin, ptMax, refFrameName);

	const char* title = Form(";%s;%s;acceptance", CosThetaVarTitle(refFrameName), PhiAxisTitle(refFrameName));

	TEfficiency* tEff = new TEfficiency(name, title, NCosThetaBins, gCosThetaMin, gCosThetaMax, NPhiBins, gPhiMin, gPhiMax);

	tEff->SetStatisticOption(gTEffStatOption);

	return tEff;
}

// (cos theta, phi, pT) 3D maps

// naming function so that we can more easily get the TEfficiency from other codes
const char* NominalTEfficiency3DName(const char* refFrameName = "CS") {
	return Form("Nominal3D%s", refFrameName);
}

const char* RelativeSystTEfficiency3DName(const char* refFrameName = "CS") {
	return Form("RelativeSyst3D%s", refFrameName);
}

const char* TEfficiency3DAxisTitle(const char* refFrameName = "CS") {
	return Form("%s;%s;%s", CosThetaVarTitle(refFrameName), PhiAxisTitle(refFrameName), gPtAxisTitle);
}

const char* TEfficiency3DTitle(const char* refFrameName = "CS", int iState = gUpsilonState) {
	return Form("%s;%s", TEfficiencyMainTitle(iState), TEfficiency3DAxisTitle(refFrameName));
}

TEfficiency* TEfficiency3D(const char* name, const char* refFrameName = "CS", int iState = gUpsilonState) {
	TEfficiency* tEff = new TEfficiency(name, TEfficiency3DTitle(refFrameName, iState), NCosThetaFineBins, gCosThetaFineBinning, NPhiFineBins, gPhiFineBinning, NPtFineBins, gPtFineBinning); // variable size binning for the stats

	tEff->SetStatisticOption(gTEffStatOption);

	return tEff;
}
