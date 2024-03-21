// when the header is included several times, to avoid the redefinition error
#ifndef Analysis_Parameters_GUARD
#define Analysis_Parameters_GUARD

#include "Tools/Style/tdrStyle.C"
#include "Tools/Style/CMS_lumi.C"

/// Definition of all global parameters to be used within this analysis framework

const char* gCMSLumiText = "PbPb 1.61 nb^{-1} (5.02 TeV)";

// CANNOT BE CHANGED!
const int gUpsilonHLTBit = 14;
const int gL2FilterBit = 19;
const int gL3FilterBit = 20;

const float gMuonPtCut = 3.5;

// One state at a time, patience is the key!

const int gUpsilonState = 1;

/// Kinematic intervals
const int gCentralityBinMin = 0;
const int gCentralityBinMax = 90;

const float gRapidityMin = 0.0;
const float gRapidityMax = 2.4;

// pT binning for the final results!
const int NPtBins = 7;
const double gPtBinning[NPtBins + 1] = {0, 2, 4, 6, 8, 12, 16, 30};

const int NPtFineBins = 10;
const double gPtFineBinning[NPtFineBins + 1] = {0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 30}; // can afford to split bins into two for the MC pt spectrum reweighting

// fine binning for (cos theta, phi) correction maps, same binning regardless of the reference frame
const int NCosThetaFineBins = 20;
const double gCosThetaFineBinning[NCosThetaFineBins + 1] = {-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};

const int NPhiFineBins = 18; // i.e. 20 degrees wide bins
const double gPhiFineBinning[NCosThetaFineBins + 1] = {-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180};

/// cos theta and phi binning, specific to each pt interval

const int NCosThetaBins = 20;
const float gCosThetaMin = -1;
const float gCosThetaMax = 1;

const int NPhiBins = 18;
const float gPhiMin = -180;
const float gPhiMax = 180;

// select one and comment out the others!!
/*
const int gPtMin = 0;
const int gPtMax = 2;

// 0 < pt < 2 GeV, Collins-Soper
const int NCosThetaBinsCS = 6;
//const double CosThetaBinningCS[NCosThetaBinsCS + 1] = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
const double CosThetaBinningCS[NCosThetaBinsCS + 1] = {-0.5, -0.3, -0.15, 0, 0.15, 0.3, 0.5};

const int NPhiBinsCS = 8;
const double PhiBinningCS[NPhiBinsCS + 1] = {-180, -120, -80, -40, 0, 40, 80, 120, 180};

// 0 < pt < 2 GeV, Helicity
const int NCosThetaBinsHX = 8;
const double CosThetaBinningHX[NCosThetaBinsHX + 1] = {-1, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 1};

const int NPhiBinsHX = 6;
const double PhiBinningHX[NPhiBinsHX + 1] = {-180, -120, -60, 0, 60, 120, 180};
*/

const int gPtMin = 0;
const int gPtMax = 30;

// 16 < pt < 30 GeV, Lab
const int NCosThetaBinsLab = 10;
//const double CosThetaBinningLab[NCosThetaBinsLab + 1] = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
const double CosThetaBinningLab[NCosThetaBinsLab + 1] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};

const int NPhiBinsLab = 1; // 1D analysis
const double PhiBinningLab[NPhiBinsLab + 1] = {-180, 180};

// 16 < pt < 30 GeV, Collins-Soper
const int NCosThetaBinsCS = 10;
//const double CosThetaBinningCS[NCosThetaBinsCS + 1] = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
const double CosThetaBinningCS[NCosThetaBinsCS + 1] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};

const int NPhiBinsCS = 1; // 1D analysis
const double PhiBinningCS[NPhiBinsCS + 1] = {-180, 180};

// 16 < pt < 30 GeV, Helicity
const int NCosThetaBinsHX = 10;
//const double CosThetaBinningHX[NCosThetaBinsHX + 1] = {-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6};
const double CosThetaBinningHX[NCosThetaBinsHX + 1] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};

const int NPhiBinsHX = 1; // 1D analysis
const double PhiBinningHX[NPhiBinsHX + 1] = {-180, 180};

/// Settings for invariant mass fits
const int NCPUs = 3;

const float MassBinMin = 7;
const float MassBinMax = 13;
const int NMassBins = 100;

// error estimate
bool DoMCWeightedError = true;

bool DoAsymptoticError = true;

/// Cosmetics

// colors (to be agreed on)
const Color_t Color1S = kRed + 1;
const Color_t Color2S = kGreen + 2;
const Color_t Color3S = kAzure + 1;
const Color_t ColorBkg = kGray + 1;

// consistent naming
const char* CentralityRangeText(int centMin = gCentralityBinMin, int centMax = gCentralityBinMax) {
	return Form("Centrality %d-%d%%", centMin, centMax);
}

const char* gPtVarName = "#it{p}_{T}";
const char* gPtUnit = "GeV/#it{c}";
const char* gPtAxisTitle = Form("%s (%s)", gPtVarName, gPtUnit);

const char* gMuonPtCutText = "#it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}"; // weird text other wise, don't know why...

const char* gDimuonPtVarTitle = "#it{p}_{T}^{ #mu#mu}";
const char* DimuonPtRangeText(int ptMin, int ptMax) {
	return Form("%d < %s < %d %s", ptMin, gDimuonPtVarTitle, ptMax, gPtUnit);
}

const char* gDimuonRapidityVarTitle = "|#it{y}^{ #mu#mu}|";
const char* DimuonRapidityRangeText(float rapidityMin, float rapidityMax) {
	return Form("%1.1f < %s < %1.1f", rapidityMin, gDimuonRapidityVarTitle, rapidityMax);
}

const char* gMassVarTitle = "#it{m}_{ #mu^{#plus}#mu^{#font[122]{\55}}}";
const char* gMassUnit = "GeV/#it{c}^{ 2}";
const char* gMassAxisTitle = Form("%s (%s)", gMassVarTitle, gMassUnit);

const char* CosThetaVarName(const char* refFrameName = "CS") {
	return Form("cosTheta%s", refFrameName);
}
const char* CosThetaVarTitle(const char* refFrameName = "CS") {
	return Form("cos #theta_{%s}", refFrameName);
}
const char* CosThetaRangeText(const char* refFrameName = "CS", float cosThetaMin = -1, float cosThetaMax = 1) {
	return Form("%.2f < %s < %.2f", cosThetaMin, CosThetaVarTitle(refFrameName), cosThetaMax);
}

const char* gPhiSymbol = "#varphi";
const char* gPhiUnit = "#circ";

const char* PhiVarName(const char* refFrameName = "CS") {
	return Form("phi%s", refFrameName);
}
const char* PhiVarTitle(const char* refFrameName = "CS") {
	return Form("%s_{%s}", gPhiSymbol, refFrameName);
}
const char* PhiAxisTitle(const char* refFrameName) {
	return Form("%s (%s)", PhiVarTitle(refFrameName), gPhiUnit);
}
const char* PhiRangeText(const char* refFrameName = "CS", int phiMin = -180, int phiMax = 180) {
	return Form("%d < %s < %d %s", phiMin, PhiVarTitle(refFrameName), phiMax, gPhiUnit);
}

const char* RawDatasetName(const char* refFrameName = "CS") {
	return Form("dataset%s", refFrameName);
}

#endif
