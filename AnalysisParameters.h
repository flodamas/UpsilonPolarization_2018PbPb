// when the header is included several times, to avoid the redefinition error
#ifndef Analysis_Parameters_GUARD
#define Analysis_Parameters_GUARD

#include "Tools/Style/tdrStyle.C"
#include "Tools/Style/CMS_lumi.C"
#include "Tools/Style/ColorSchemes.h"

bool BeVerbose = true;

/// Definition of all global parameters to be used within this analysis framework

const char* gCMSLumiText = "PbPb 1.61 nb^{#minus1} (5.02 TeV)";

const char* gMuonAccName = "UpsilonTriggerThresholds"; // possible names listed in PhaseSpace.h
// const char* gMuonAccName = "_2018Acc";
// const char* gMuonAccName = "";

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

const float gRapidityMin = 0.;
const float gRapidityMax = 2.4;

// pT binning for the final results!
const int NPtBins = 4;
const double gPtBinning[NPtBins + 1] = {0, 2, 6, 12, 20};

const int NPtFineBins = 11;
const double gPtFineBinning[NPtFineBins + 1] = {0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 30}; // can afford to split bins into two for the MC pt spectrum reweighting

// const int NPtFineBins = 2;
// const double gPtFineBinning[NPtFineBins + 1] = {0, 2, 30}; // can afford to split bins into two for the MC pt spectrum reweighting

// fine binning for (cos theta, phi) correction maps, same binning regardless of the reference frame
const int NCosThetaFineBins = 20;
const double gCosThetaFineBinning[NCosThetaFineBins + 1] = {-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};

const int NPhiFineBins = 19;
const double gPhiFineBinning[NPhiFineBins + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 240}; // last extra bin for legends

/// cos theta and phi binning, specific to each pt interval

const int NCosThetaBins = 5;
const float gCosThetaMin = -0.7;
const float gCosThetaMax = 0.7;
const double gCosThetaBinning[NCosThetaBins + 1] = {-0.7, -0.42, -0.14, 0.14, 0.42, 0.7};

const int NPhiBins = 3;
//const float gPhiMin = -180;
const float gPhiMin = 0;
const float gPhiMax = 180;
const double gPhiBinning[NPhiBins + 1] = {0, 60, 120, 180};

const int NFullPhiBins = 6;
const float gFullPhiMin = -180;
const float gFullPhiMax = 180;
const double gFullPhiBinning[NFullPhiBins + 1] = {-180, -120, -60, 0, 60, 120, 180};

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
const int gPtMax = 40;

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
// const float MassBinMax = 13;
const float MassBinMax = 12;
const int NMassBins = 100;

// error estimate
bool DoMCWeightedError = true;

bool DoAsymptoticError = true;

/// Cosmetics

// definition of the actual colors globally used in the code

std::string gPreferredColorPaletteName = "sunset";
std::string gAcceptanceColorPaletteName = "cividis";
std::string gEfficiencyColorPaletteName = "cividis";

const Color_t gColor1S = TamDragonRed;
const Color_t gColor2S = TamDragonOrange;
const Color_t gColor3S = TamDragonOrangeYellow;
const Color_t gColorBkg = kGray + 1;
const Color_t gColorTotalFit = TamDragonDarkPurple;

// consistent naming
const char* CentralityRangeText(int centMin = gCentralityBinMin, int centMax = gCentralityBinMax) {
	return Form("Centrality %d-%d%%", centMin, centMax);
}

const char* gPtVarName = "#it{p}_{T}";
const char* gPtUnit = "GeV/#it{c}";
const char* gPtAxisTitle = Form("%s (%s)", gPtVarName, gPtUnit);

const char* gMuonPtCutText = "#it{p}_{T}^{ #mu} > 4 and 3 GeV/#it{c}";
//const char* gMuonPtCutText = "#it{p}_{T}^{ #mu} > thresholds(#eta^{ #mu})"; // "#it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}"; // weird text other wise, don't know why...
// const char* gMuonPtCutText = "#it{p}_{T}^{ #mu} > 3.5 GeV/#it{c}"; // weird text other wise, don't know why...

const char* gDimuonPtVarTitle = "#it{p}_{T}";
const char* DimuonPtRangeText(int ptMin, int ptMax) {
	return Form("%d < %s < %d %s", ptMin, gDimuonPtVarTitle, ptMax, gPtUnit);
}

const char* gDimuonRapidityVarTitle = "|#it{y}|";
const char* DimuonRapidityRangeText(float rapidityMin = gRapidityMin, float rapidityMax = gRapidityMax) {
	return Form("%1.1f < %s < %1.1f", rapidityMin, gDimuonRapidityVarTitle, rapidityMax);
}

const char* gMassVarTitle = "#it{m}_{#mu^{#plus}#mu^{#font[122]{\55}}}";
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

const char* AbsPhiVarName(const char* refFrameName = "CS") {
	return Form("absPhi%s", refFrameName);
}
const char* AbsPhiVarTitle(const char* refFrameName = "CS") {
	return Form("|%s_{%s}|", gPhiSymbol, refFrameName);
}
const char* AbsPhiAxisTitle(const char* refFrameName) {
	return Form("%s (%s)", AbsPhiVarTitle(refFrameName), gPhiUnit);
}
const char* AbsPhiRangeText(const char* refFrameName = "CS", int phiMin = -180, int phiMax = 180) {
	return Form("%d < %s < %d %s", phiMin, AbsPhiVarTitle(refFrameName), phiMax, gPhiUnit);
}

const char* gPhiTildeSymbol = "#tilde{#varphi}";

const char* PhiTildeVarName(const char* refFrameName = "CS") {
	return Form("phiTilde%s", refFrameName);
}
const char* PhiTildeVarTitle(const char* refFrameName = "CS") {
	return Form("%s_{%s}", gPhiTildeSymbol, refFrameName);
}

const char* RawDatasetName(const char* refFrameName = "CS") {
	return Form("dataset%s", refFrameName);
}

const char* PolaWeightName(Float_t lambdaTheta = 0., Float_t lambdaPhi = 0., Float_t lambdaThetaPhi = 0.) {
	return Form("LambdaTheta%.2fPhi%.2fThetaPhi%.2f", lambdaTheta, lambdaPhi, lambdaThetaPhi);
}

#endif
